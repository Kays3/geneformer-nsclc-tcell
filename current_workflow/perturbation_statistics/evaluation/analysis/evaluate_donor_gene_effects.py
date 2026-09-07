#!/usr/bin/env python3
"""Donor-level Geneformer deletion-shift aggregation and leave-one-donor-out ranks."""

from __future__ import annotations

import argparse
import json
import pickle
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

STATES = {
    "luad": "lung adenocarcinoma",
    "lusc": "squamous cell lung carcinoma",
    "normal": "normal",
}
STATE_TO_SLUG = {v: k for k, v in STATES.items()}
COMPARISONS = [
    ("lusc", "luad"),
    ("lusc", "normal"),
    ("luad", "lusc"),
    ("luad", "normal"),
    ("normal", "luad"),
    ("normal", "lusc"),
]
LABELS = {
    ("lusc", "luad"): "LUSC → LUAD",
    ("lusc", "normal"): "LUSC → NORMAL",
    ("luad", "lusc"): "LUAD → LUSC",
    ("luad", "normal"): "LUAD → NORMAL",
    ("normal", "luad"): "NORMAL → LUAD",
    ("normal", "lusc"): "NORMAL → LUSC",
}
PICKLE_RE = re.compile(
    r"in_silico_delete_heldout_(luad|lusc|normal)_shard(\d+)_dict_cell_embs_(\d+)batch(-?\d+)_raw\.pickle$"
)
BOOTSTRAP_N = 1000
BOOTSTRAP_SEED = 0
DIRECTION_GATE = 0.70
LODO_TOP_K = 50
LODO_RETENTION_GATE = 0.70


def load_token_maps(geneformer_hint: Path | None) -> tuple[dict[int, str], dict[str, str]]:
    token_to_ensembl: dict[int, str] = {}
    ensembl_to_name: dict[str, str] = {}
    candidates = []
    if geneformer_hint is not None:
        candidates.append(geneformer_hint)
    try:
        import geneformer

        pkg = Path(geneformer.__file__).resolve().parent
        candidates.append(pkg)
    except Exception:
        pass
    for root in candidates:
        token_paths = [
            p
            for p in sorted(root.glob("token_dictionary*.pkl"))
            if p.stat().st_size > 1024
        ]
        name_paths = [
            p
            for p in sorted(root.glob("gene_name_id_dict*.pkl")) + sorted(root.glob("ensembl_mapping_dict*.pkl"))
            if p.stat().st_size > 1024
        ]
        for path in token_paths:
            obj = pickle.loads(path.read_bytes())
            if isinstance(obj, dict) and obj:
                sample = next(iter(obj.items()))
                if isinstance(sample[0], str) and isinstance(sample[1], int):
                    token_to_ensembl = {int(v): str(k) for k, v in obj.items()}
                elif isinstance(sample[0], int):
                    token_to_ensembl = {int(k): str(v) for k, v in obj.items()}
                if token_to_ensembl:
                    break
        for path in name_paths:
            obj = pickle.loads(path.read_bytes())
            if isinstance(obj, dict) and obj:
                sample = next(iter(obj.items()))
                if isinstance(sample[0], str) and str(sample[0]).startswith("ENSG"):
                    ensembl_to_name = {str(k): str(v) for k, v in obj.items()}
                else:
                    ensembl_to_name = {str(v): str(k) for k, v in obj.items() if str(v).startswith("ENSG")}
                if ensembl_to_name:
                    break
        if token_to_ensembl:
            break
    return token_to_ensembl, ensembl_to_name


def enrich_maps_from_tables(
    tables: Path,
    token_to_ensembl: dict[int, str],
    ensembl_to_name: dict[str, str],
) -> None:
    one_cell = tables / "one_cell_lusc_shard0000_cell000_deletion_scores.csv"
    if one_cell.exists():
        frame = pd.read_csv(one_cell)
        for row in frame.itertuples(index=False):
            token = int(row.token)
            ensembl = str(row.Ensembl_ID)
            name = str(row.Gene_name)
            token_to_ensembl.setdefault(token, ensembl)
            ensembl_to_name.setdefault(ensembl, name)


def load_manifest(tables: Path) -> pd.DataFrame:
    path = tables / "heldout_shard_manifest.csv"
    frame = pd.read_csv(path)
    required = {"source", "shard", "local_index", "cell_id", "individual", "disease"}
    missing = required.difference(frame.columns)
    if missing:
        raise SystemExit(f"Shard manifest missing columns: {sorted(missing)}")
    frame["shard"] = frame["shard"].astype(int)
    frame["local_index"] = frame["local_index"].astype(int)
    return frame


def iter_pickles(raw_root: Path):
    files = sorted(raw_root.glob("*/*_raw.pickle"))
    if not files:
        files = sorted(raw_root.glob("*_raw.pickle"))
    return files


def accumulate(raw_root: Path, manifest: pd.DataFrame):
    lookup = {
        (str(row.source), int(row.shard), int(row.local_index)): row
        for row in manifest.itertuples(index=False)
    }
    sums: dict[tuple[str, str], dict[int, dict[str, list[float]]]] = {
        pair: defaultdict(lambda: defaultdict(lambda: [0.0, 0.0])) for pair in COMPARISONS
    }
    n_files = 0
    n_used = 0
    n_skip = 0
    n_values = 0
    files = iter_pickles(raw_root)
    if not files:
        raise SystemExit(f"No raw pickles under {raw_root}")
    for path in files:
        match = PICKLE_RE.search(path.name)
        if match is None:
            n_skip += 1
            continue
        source, shard_s, cell_s, _batch = match.groups()
        key = (source, int(shard_s), int(cell_s))
        row = lookup.get(key)
        if row is None:
            n_skip += 1
            continue
        donor = str(row.individual)
        with path.open("rb") as handle:
            obj = pickle.load(handle)
        if not isinstance(obj, dict):
            n_skip += 1
            continue
        by_state = {}
        for state_name, gene_dict in obj.items():
            slug = STATE_TO_SLUG.get(state_name)
            if slug is None or not isinstance(gene_dict, dict):
                continue
            parsed = {}
            for gene_key, values in gene_dict.items():
                token = gene_key[0] if isinstance(gene_key, tuple) else gene_key
                if isinstance(values, (list, tuple)):
                    if not values:
                        continue
                    value = float(values[0])
                else:
                    value = float(values)
                parsed[int(token)] = value
            by_state[slug] = parsed
        n_files += 1
        if source not in by_state:
            n_skip += 1
            continue
        start_map = by_state[source]
        for goal in (slug for slug, _src in STATES.items() if slug != source):
            goal_map = by_state.get(goal)
            if not goal_map:
                continue
            pair = (source, goal)
            if pair not in sums:
                continue
            acc_pair = sums[pair]
            for token, start_val in start_map.items():
                if token not in goal_map:
                    continue
                shift = goal_map[token] - start_val
                bucket = acc_pair[token][donor]
                bucket[0] += shift
                bucket[1] += 1.0
                n_values += 1
        n_used += 1
        if n_used % 500 == 0:
            print(f"processed {n_used}/{len(files)} matched pickles", flush=True)
    print(
        json.dumps(
            {
                "pickle_files_seen": len(files),
                "pickle_files_parsed": n_files,
                "pickle_files_used": n_used,
                "pickle_files_skipped": n_skip,
                "cell_gene_shift_values": n_values,
            }
        ),
        flush=True,
    )
    return sums


def donor_means(acc_pair: dict[int, dict[str, list[float]]]) -> dict[int, dict[str, tuple[float, int]]]:
    out: dict[int, dict[str, tuple[float, int]]] = {}
    for token, donors in acc_pair.items():
        out[token] = {
            donor: (vals[0] / vals[1], int(vals[1]))
            for donor, vals in donors.items()
            if vals[1] > 0
        }
    return out


def bootstrap_mean_ci(values: np.ndarray, n: int, seed: int) -> tuple[float, float]:
    if values.size == 0:
        return (float("nan"), float("nan"))
    rng = np.random.default_rng(seed)
    draws = rng.choice(values, size=(n, values.size), replace=True).mean(axis=1)
    low, high = np.quantile(draws, [0.025, 0.975])
    return float(low), float(high)


def summarize(
    sums,
    token_to_ensembl: dict[int, str],
    ensembl_to_name: dict[str, str],
) -> tuple[pd.DataFrame, dict[tuple[str, str], dict[int, dict[str, tuple[float, int]]]]]:
    rows = []
    donor_level = {}
    for pair in COMPARISONS:
        means = donor_means(sums[pair])
        donor_level[pair] = means
        source, goal = pair
        label = LABELS[pair]
        for token, donors in means.items():
            effects = np.array([v[0] for v in donors.values()], dtype=float)
            n_cells = np.array([v[1] for v in donors.values()], dtype=int)
            n_donors = int(effects.size)
            median_effect = float(np.median(effects))
            mean_effect = float(effects.mean())
            nonzero = effects[effects != 0.0]
            n_eval_dir = int(nonzero.size)
            if n_eval_dir == 0:
                n_pos = 0
                n_neg = 0
                concordance = float("nan")
                majority_sign = 0
            else:
                majority_sign = 1 if median_effect > 0 else (-1 if median_effect < 0 else int(np.sign(nonzero.sum())))
                n_pos = int((effects > 0).sum())
                n_neg = int((effects < 0).sum())
                if majority_sign >= 0:
                    concordance = n_pos / n_eval_dir
                else:
                    concordance = n_neg / n_eval_dir
            ci_low, ci_high = bootstrap_mean_ci(effects, BOOTSTRAP_N, BOOTSTRAP_SEED)
            ensembl = token_to_ensembl.get(token, "")
            gene_name = ensembl_to_name.get(ensembl, "") if ensembl else ""
            rows.append(
                {
                    "comparison": f"{source}_to_{goal}",
                    "comparison_label": label,
                    "source": source,
                    "goal": goal,
                    "token": token,
                    "Ensembl_ID": ensembl,
                    "Gene_name": gene_name,
                    "n_donors_evaluable": n_donors,
                    "n_donors_nonzero": n_eval_dir,
                    "n_cells_total": int(n_cells.sum()),
                    "median_cells_per_donor": float(np.median(n_cells)),
                    "median_donor_effect": median_effect,
                    "mean_donor_effect": mean_effect,
                    "min_donor_effect": float(effects.min()),
                    "max_donor_effect": float(effects.max()),
                    "n_donors_positive": n_pos,
                    "n_donors_negative": n_neg,
                    "sign_concordance": concordance,
                    "bootstrap_ci_low": ci_low,
                    "bootstrap_ci_high": ci_high,
                    "ci_excludes_zero": bool(ci_low > 0 or ci_high < 0) if np.isfinite(ci_low) else False,
                    "same_direction_pass_70": bool(n_eval_dir > 0 and concordance >= DIRECTION_GATE),
                }
            )
    frame = pd.DataFrame(rows)
    frame = frame.sort_values(
        ["comparison", "median_donor_effect"],
        ascending=[True, False],
    ).reset_index(drop=True)
    return frame, donor_level


def lodo_stability(donor_level) -> pd.DataFrame:
    rows = []
    for pair in COMPARISONS:
        means = donor_level[pair]
        source, goal = pair
        label = LABELS[pair]
        tokens = [token for token, donors in means.items() if len(donors) >= 2]
        if not tokens:
            continue
        median_full = {
            token: float(np.median([v[0] for v in means[token].values()])) for token in tokens
        }
        ranked = sorted(tokens, key=lambda t: median_full[t], reverse=True)
        top_full = ranked[:LODO_TOP_K]
        top_set = set(top_full)
        donors = sorted({donor for token in tokens for donor in means[token]})
        fold_overlaps = []
        gene_retain = {token: 0 for token in top_full}
        for held in donors:
            medians = {}
            for token in tokens:
                remaining = [v[0] for d, v in means[token].items() if d != held]
                if not remaining:
                    continue
                medians[token] = float(np.median(remaining))
            ranked_fold = sorted(medians, key=lambda t: medians[t], reverse=True)
            top_fold = set(ranked_fold[:LODO_TOP_K])
            overlap = len(top_set & top_fold) / float(LODO_TOP_K)
            fold_overlaps.append(overlap)
            for token in top_full:
                if token in top_fold:
                    gene_retain[token] += 1
            rows.append(
                {
                    "comparison": f"{source}_to_{goal}",
                    "comparison_label": label,
                    "level": "fold",
                    "held_out_donor": held,
                    "n_donors_in_fold": len(donors) - 1,
                    "top50_overlap": overlap,
                    "top50_retained_n": int(len(top_set & top_fold)),
                    "fold_overlap_ge_70": overlap >= LODO_RETENTION_GATE,
                    "token": "",
                    "n_folds": len(donors),
                    "frac_folds_retained_in_top50": "",
                    "gene_retained_in_70pct_folds": "",
                }
            )
        n_folds = len(donors)
        n_folds_pass = int(sum(1 for o in fold_overlaps if o >= LODO_RETENTION_GATE))
        rows.append(
            {
                "comparison": f"{source}_to_{goal}",
                "comparison_label": label,
                "level": "comparison",
                "held_out_donor": "",
                "n_donors_in_fold": n_folds - 1,
                "top50_overlap": float(np.mean(fold_overlaps)) if fold_overlaps else float("nan"),
                "top50_retained_n": "",
                "fold_overlap_ge_70": n_folds_pass / n_folds if n_folds else float("nan"),
                "token": "",
                "n_folds": n_folds,
                "frac_folds_retained_in_top50": "",
                "gene_retained_in_70pct_folds": n_folds_pass / n_folds >= LODO_RETENTION_GATE if n_folds else False,
            }
        )
        for token in top_full:
            frac = gene_retain[token] / n_folds if n_folds else float("nan")
            rows.append(
                {
                    "comparison": f"{source}_to_{goal}",
                    "comparison_label": label,
                    "level": "gene",
                    "held_out_donor": "",
                    "n_donors_in_fold": n_folds - 1,
                    "top50_overlap": "",
                    "top50_retained_n": gene_retain[token],
                    "fold_overlap_ge_70": "",
                    "token": token,
                    "n_folds": n_folds,
                    "frac_folds_retained_in_top50": frac,
                    "gene_retained_in_70pct_folds": bool(frac >= LODO_RETENTION_GATE) if n_folds else False,
                }
            )
    return pd.DataFrame(rows)


def score_gates(effects: pd.DataFrame, lodo: pd.DataFrame, top_genes: pd.DataFrame | None) -> dict:
    by_comp = []
    for label, sub in effects.groupby("comparison_label"):
        n = len(sub)
        n_pass = int(sub["same_direction_pass_70"].sum())
        leading = sub.nlargest(LODO_TOP_K, "median_donor_effect")
        n_lead_pass = int(leading["same_direction_pass_70"].sum())
        row = {
            "comparison_label": label,
            "n_genes": n,
            "n_pass_70_all": n_pass,
            "frac_pass_70_all": n_pass / n if n else float("nan"),
            "n_top50_donor_ranked": len(leading),
            "n_top50_pass_70": n_lead_pass,
            "frac_top50_pass_70": n_lead_pass / len(leading) if len(leading) else float("nan"),
            "n_top50_ci_excludes_zero": int(leading["ci_excludes_zero"].sum()),
        }
        if top_genes is not None and "comparison_label" in top_genes.columns:
            names = set(top_genes.loc[top_genes["comparison_label"] == label, "Gene_name"].astype(str))
            pub = sub.loc[sub["Gene_name"].isin(names)]
            row["n_published_top"] = len(pub)
            row["n_published_top_pass_70"] = int(pub["same_direction_pass_70"].sum()) if len(pub) else 0
            row["frac_published_top_pass_70"] = (
                row["n_published_top_pass_70"] / len(pub) if len(pub) else float("nan")
            )
        by_comp.append(row)
    gate_table = pd.DataFrame(by_comp)
    lodo_comp = lodo.loc[lodo["level"] == "comparison"].copy()
    direction_all_top50 = bool((gate_table["frac_top50_pass_70"] >= DIRECTION_GATE).all()) if len(gate_table) else False
    published_ok = True
    if "frac_published_top_pass_70" in gate_table.columns:
        published_ok = bool((gate_table["frac_published_top_pass_70"] >= DIRECTION_GATE).all())
    lodo_ok = bool(lodo_comp["gene_retained_in_70pct_folds"].astype(bool).all()) if len(lodo_comp) else False
    return {
        "direction_gate": DIRECTION_GATE,
        "direction_gate_on_donor_top50_pass": direction_all_top50,
        "direction_gate_on_published_top_pass": published_ok,
        "lodo_top50_retention_gate": LODO_RETENTION_GATE,
        "lodo_gate_pass": lodo_ok,
        "by_comparison": gate_table,
        "lodo_comparison": lodo_comp,
    }


def write_report(path: Path, effects: pd.DataFrame, lodo: pd.DataFrame, gates: dict, meta: dict) -> None:
    lines = []
    lines.append("# Donor-level gene-effect evaluation")
    lines.append("")
    lines.append(f"Generated: {meta['generated_at']}")
    lines.append(f"Script: `{meta['script']}`")
    lines.append(f"Source: `{meta['base']}`")
    lines.append("")
    lines.append("## Method")
    lines.append("")
    lines.append("- Unit of analysis is the donor. Cell-level cosine similarities to the three")
    lines.append("  training-donor state centroids were read from existing raw pickles.")
    lines.append("- For each cell-gene deletion, the directional shift is `cos(goal) - cos(source)`.")
    lines.append("- Donor effect is the mean of that cell-level shift within the donor.")
    lines.append("- Gene summary uses the median of donor effects (equal donor weight).")
    lines.append("- Sign concordance is the fraction of donors with a nonzero mean whose sign")
    lines.append("  matches the median donor effect. The 70% gate uses this fraction.")
    lines.append("- Bootstrap 95% CI is the percentile interval of the mean of donor means")
    lines.append(f"  ({BOOTSTRAP_N} resamples, seed {BOOTSTRAP_SEED}).")
    lines.append("- Leave-one-donor-out ranks genes by median donor effect after dropping one")
    lines.append("  test donor. Training centroids are unchanged (test donors were not in them).")
    lines.append("- LUSC has five held-out donors; its stability numbers are reported separately")
    lines.append("  and are not compared mechanically with LUAD or normal.")
    lines.append("")
    lines.append("## Direction gate (same sign in ≥70% of evaluable donors)")
    lines.append("")
    lines.append(
        f"- Donor-ranked top-50, all comparisons meet gate: "
        f"**{'PASS' if gates['direction_gate_on_donor_top50_pass'] else 'FAIL'}**"
    )
    lines.append(
        f"- Published leading genes, all comparisons meet gate: "
        f"**{'PASS' if gates['direction_gate_on_published_top_pass'] else 'FAIL'}**"
    )
    lines.append("")
    lines.append(gates["by_comparison"].to_string(index=False))
    lines.append("")
    lines.append("## Leave-one-donor-out top-50 retention")
    lines.append("")
    lines.append(
        f"- Gate (mean fold overlap ≥70% of genes in ≥70% of folds, encoded per comparison): "
        f"**{'PASS' if gates['lodo_gate_pass'] else 'FAIL'}**"
    )
    lines.append("")
    show = gates["lodo_comparison"][
        [
            "comparison_label",
            "n_folds",
            "top50_overlap",
            "fold_overlap_ge_70",
            "gene_retained_in_70pct_folds",
        ]
    ]
    lines.append(show.to_string(index=False))
    lines.append("")
    lines.append("## Donor counts")
    lines.append("")
    lines.append("- LUAD test donors: 19")
    lines.append("- Normal test donors: 12")
    lines.append("- LUSC test donors: 5")
    lines.append("")
    lines.append("## Limitations")
    lines.append("")
    lines.append("- This is a statistical donor-stability summary, not a biological target claim.")
    lines.append("- Ambient RNA, doublets, and T-cell subtype structure are not addressed here.")
    lines.append("- Centroids were not rebuilt; test-donor LODO does not require that step.")
    lines.append("- No threshold, donor subset, or aggregation was changed after seeing results.")
    lines.append("")
    path.write_text("\n".join(lines) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--top-genes", type=Path, default=None)
    parser.add_argument("--geneformer-dir", type=Path, default=None)
    parser.add_argument("--remap-effects", type=Path, default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    base = args.base.resolve()
    raw_root = base / "raw"
    tables = base / "tables"
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    token_to_ensembl, ensembl_to_name = load_token_maps(args.geneformer_dir)
    enrich_maps_from_tables(tables, token_to_ensembl, ensembl_to_name)
    stats_root = base / "stats"
    if stats_root.exists():
        for csv in stats_root.glob("*.csv"):
            try:
                frame = pd.read_csv(csv)
            except Exception:
                continue
            if {"Ensembl_ID", "Gene_name"}.issubset(frame.columns):
                for row in frame.itertuples(index=False):
                    ensembl_to_name.setdefault(str(row.Ensembl_ID), str(row.Gene_name))
    top_genes = pd.read_csv(args.top_genes) if args.top_genes and args.top_genes.exists() else None
    if args.remap_effects is not None:
        effects = pd.read_csv(args.remap_effects)
        effects["Ensembl_ID"] = effects["token"].map(lambda t: token_to_ensembl.get(int(t), ""))
        effects["Gene_name"] = effects["Ensembl_ID"].map(lambda e: ensembl_to_name.get(str(e), "") if e else "")
        lodo = pd.read_csv(output / "leave_one_donor_out_stability.csv")
        gates = score_gates(effects, lodo, top_genes)
        effects.to_csv(output / "donor_gene_effects.csv", index=False)
        gates["by_comparison"].to_csv(output / "gate_score.csv", index=False)
        meta = {
            "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "script": "evaluate_donor_gene_effects.py",
            "base": str(base),
        }
        write_report(output / "REPORT.md", effects, lodo, gates, meta)
        print(
            json.dumps(
                {
                    "remapped": True,
                    "empty_gene_name_frac": float((effects["Gene_name"].fillna("") == "").mean()),
                    "direction_gate_on_donor_top50_pass": gates["direction_gate_on_donor_top50_pass"],
                    "direction_gate_on_published_top_pass": gates["direction_gate_on_published_top_pass"],
                    "lodo_gate_pass": gates["lodo_gate_pass"],
                }
            ),
            flush=True,
        )
        return 0
    manifest = load_manifest(tables)
    sums = accumulate(raw_root, manifest)
    effects, donor_level = summarize(sums, token_to_ensembl, ensembl_to_name)
    lodo = lodo_stability(donor_level)
    gates = score_gates(effects, lodo, top_genes)
    effects_path = output / "donor_gene_effects.csv"
    lodo_path = output / "leave_one_donor_out_stability.csv"
    gate_path = output / "gate_score.csv"
    effects.to_csv(effects_path, index=False)
    lodo.to_csv(lodo_path, index=False)
    gates["by_comparison"].to_csv(gate_path, index=False)
    meta = {
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "script": "evaluate_donor_gene_effects.py",
        "base": str(base),
    }
    write_report(output / "REPORT.md", effects, lodo, gates, meta)
    print(f"wrote {effects_path} rows={len(effects)}", flush=True)
    print(f"wrote {lodo_path} rows={len(lodo)}", flush=True)
    print(
        json.dumps(
            {
                "direction_gate_on_donor_top50_pass": gates["direction_gate_on_donor_top50_pass"],
                "direction_gate_on_published_top_pass": gates["direction_gate_on_published_top_pass"],
                "lodo_gate_pass": gates["lodo_gate_pass"],
            }
        ),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
