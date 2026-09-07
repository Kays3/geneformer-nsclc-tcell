#!/usr/bin/env python3
"""Donor-level Geneformer deletion-shift aggregation and leave-one-donor-out ranks."""

from __future__ import annotations

import argparse
import json
import math
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
SOURCE_N_TEST_DONORS = {"luad": 19, "lusc": 5, "normal": 12}


def min_donors_required(source: str) -> int:
    return max(3, math.ceil(SOURCE_N_TEST_DONORS[source] / 2.0))


PRE_SPECIFIED_CRITERION = """## Prespecified re-score criterion (written before the new numbers)

The original unfiltered top-50 is kept below and is **not** replaced. It is not a valid 70% same-direction test: ranking by unshrunk median donor effect selects genes evaluable in one donor, and those genes have `sign_concordance = 1` by construction.

**New scored set (fixed in advance):** a gene is donor-supported in a comparison if
`n_donors_evaluable >= max(3, ceil(n_source_test_donors / 2))`:

| Source screen | Test donors | Minimum evaluable donors |
|---|---:|---:|
| LUAD | 19 | 10 |
| Normal | 12 | 6 |
| LUSC | 5 | 3 |

**Why this floor, not a data-driven one:** two donors cannot fail a 70% same-sign test (2/2 always passes). Requiring at least half of the source test donors is the smallest rule that still measures cross-donor reproducibility rather than rarity. The cut is source-specific because LUSC has only five held-out donors; a global cut of 10 would empty the LUSC tables by arithmetic, not by biology.

The 15 published leading genes are the a-priori arm and lead the results. Unfiltered top-50 numbers are retained as an added row. No threshold is changed after seeing results. No biological-target claim.
"""


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
                    "min_donors_required": min_donors_required(source),
                    "donor_supported": n_donors >= min_donors_required(source),
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


def lodo_stability(donor_level, ranking: str = "unfiltered", min_n: int | None = None) -> pd.DataFrame:
    rows = []
    for pair in COMPARISONS:
        means = donor_level[pair]
        source, goal = pair
        label = LABELS[pair]
        floor = min_n if min_n is not None else min_donors_required(source)
        if ranking == "unfiltered":
            tokens = [token for token, donors in means.items() if len(donors) >= 2]
        else:
            tokens = [token for token, donors in means.items() if len(donors) >= floor]
        if not tokens:
            continue
        median_full = {
            token: float(np.median([v[0] for v in means[token].values()])) for token in tokens
        }
        ranked = sorted(tokens, key=lambda t: median_full[t], reverse=True)
        k = min(LODO_TOP_K, len(ranked))
        if k == 0:
            continue
        top_full = ranked[:k]
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
            top_fold = set(ranked_fold[:k])
            overlap = len(top_set & top_fold) / float(k)
            fold_overlaps.append(overlap)
            for token in top_full:
                if token in top_fold:
                    gene_retain[token] += 1
            rows.append(
                {
                    "ranking": ranking,
                    "min_donors_required": floor if ranking != "unfiltered" else 2,
                    "n_ranked": k,
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
                "ranking": ranking,
                "min_donors_required": floor if ranking != "unfiltered" else 2,
                "n_ranked": k,
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
                    "ranking": ranking,
                    "min_donors_required": floor if ranking != "unfiltered" else 2,
                    "n_ranked": k,
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


def _set_row(label: str, scored_set: str, sub: pd.DataFrame, source: str, min_req: int) -> dict:
    n = len(sub)
    n_pass = int(sub["same_direction_pass_70"].sum()) if n else 0
    n_don = sub["n_donors_evaluable"] if n else pd.Series(dtype=float)
    return {
        "comparison_label": label,
        "source": source,
        "scored_set": scored_set,
        "min_donors_required": min_req,
        "n_genes_in_set": n,
        "n_pass_70": n_pass,
        "frac_pass_70": n_pass / n if n else float("nan"),
        "n_donors_min": int(n_don.min()) if n else 0,
        "n_donors_median": float(n_don.median()) if n else float("nan"),
        "n_donors_max": int(n_don.max()) if n else 0,
        "n_ci_excludes_zero": int(sub["ci_excludes_zero"].sum()) if n else 0,
        "gate_pass": bool(n > 0 and (n_pass / n) >= DIRECTION_GATE),
    }


def published_support(effects: pd.DataFrame, top_genes: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for label, pub in top_genes.groupby("comparison_label"):
        sub = effects.loc[effects["comparison_label"] == label]
        source = str(sub["source"].iloc[0]) if len(sub) else ""
        min_req = min_donors_required(source) if source else 0
        for row in pub.itertuples(index=False):
            name = str(row.Gene_name)
            hit = sub.loc[sub["Gene_name"] == name]
            if len(hit) == 0:
                rows.append(
                    {
                        "comparison_label": label,
                        "Gene_name": name,
                        "Ensembl_ID": getattr(row, "Ensembl_ID", ""),
                        "min_donors_required": min_req,
                        "n_donors_evaluable": 0,
                        "donor_supported": False,
                        "n_cells_total": 0,
                        "median_donor_effect": float("nan"),
                        "sign_concordance": float("nan"),
                        "same_direction_pass_70": False,
                        "evaluable_for_direction_gate": False,
                    }
                )
                continue
            rec = hit.iloc[0]
            n_don = int(rec["n_donors_evaluable"])
            rows.append(
                {
                    "comparison_label": label,
                    "Gene_name": name,
                    "Ensembl_ID": rec["Ensembl_ID"],
                    "min_donors_required": min_req,
                    "n_donors_evaluable": n_don,
                    "donor_supported": bool(rec["donor_supported"]),
                    "n_cells_total": int(rec["n_cells_total"]),
                    "median_donor_effect": float(rec["median_donor_effect"]),
                    "sign_concordance": float(rec["sign_concordance"]),
                    "same_direction_pass_70": bool(rec["same_direction_pass_70"]),
                    "evaluable_for_direction_gate": bool(rec["donor_supported"]),
                }
            )
    return pd.DataFrame(rows)


def score_gates(
    effects: pd.DataFrame,
    lodo: pd.DataFrame,
    top_genes: pd.DataFrame | None,
    published: pd.DataFrame | None,
) -> dict:
    by_comp = []
    for label, sub in effects.groupby("comparison_label"):
        source = str(sub["source"].iloc[0])
        min_req = min_donors_required(source)
        unfiltered = sub.nlargest(LODO_TOP_K, "median_donor_effect")
        supported = sub.loc[sub["donor_supported"]].nlargest(LODO_TOP_K, "median_donor_effect")
        by_comp.append(_set_row(label, "unfiltered_top50", unfiltered, source, min_req))
        by_comp.append(_set_row(label, "donor_supported_top50", supported, source, min_req))
        if published is not None and len(published):
            pub = published.loc[published["comparison_label"] == label]
            eval_pub = pub.loc[pub["evaluable_for_direction_gate"]]
            n = len(eval_pub)
            n_pass = int(eval_pub["same_direction_pass_70"].sum()) if n else 0
            n_don = eval_pub["n_donors_evaluable"] if n else pd.Series(dtype=float)
            by_comp.append(
                {
                    "comparison_label": label,
                    "source": source,
                    "scored_set": "published_leading_evaluable",
                    "min_donors_required": min_req,
                    "n_genes_in_set": n,
                    "n_pass_70": n_pass,
                    "frac_pass_70": n_pass / n if n else float("nan"),
                    "n_donors_min": int(n_don.min()) if n else 0,
                    "n_donors_median": float(n_don.median()) if n else float("nan"),
                    "n_donors_max": int(n_don.max()) if n else 0,
                    "n_ci_excludes_zero": 0,
                    "gate_pass": bool(n > 0 and (n_pass / n) >= DIRECTION_GATE),
                    "n_published_total": len(pub),
                    "n_published_below_floor": int((~pub["evaluable_for_direction_gate"]).sum()) if len(pub) else 0,
                }
            )
    gate_table = pd.DataFrame(by_comp)
    lodo_comp = lodo.loc[lodo["level"] == "comparison"].copy()

    def _all_pass(set_name: str) -> bool:
        part = gate_table.loc[gate_table["scored_set"] == set_name]
        return bool(len(part) and part["gate_pass"].all())

    def _lodo_pass(ranking: str) -> bool:
        part = lodo_comp.loc[lodo_comp["ranking"] == ranking] if "ranking" in lodo_comp.columns else lodo_comp
        if not len(part):
            return False
        return bool(part["gene_retained_in_70pct_folds"].astype(bool).all())

    return {
        "direction_gate": DIRECTION_GATE,
        "published_evaluable_pass": _all_pass("published_leading_evaluable"),
        "donor_supported_top50_pass": _all_pass("donor_supported_top50"),
        "unfiltered_top50_pass": _all_pass("unfiltered_top50"),
        "lodo_unfiltered_pass": _lodo_pass("unfiltered"),
        "lodo_donor_supported_pass": _lodo_pass("donor_supported"),
        "by_comparison": gate_table,
        "lodo_comparison": lodo_comp,
        "published": published,
    }


def write_report(path: Path, effects: pd.DataFrame, lodo: pd.DataFrame, gates: dict, meta: dict) -> None:
    lines = []
    lines.append("# Donor-level gene-effect evaluation")
    lines.append("")
    lines.append(f"Generated: {meta['generated_at']}")
    lines.append(f"Script: `{meta['script']}`")
    lines.append(f"Source: `{meta['base']}`")
    lines.append("")
    lines.append(PRE_SPECIFIED_CRITERION.rstrip())
    lines.append("")
    lines.append("## Results (a-priori published genes lead)")
    lines.append("")
    lines.append(
        f"- Published leading genes, evaluable under the floor, all comparisons ≥70% same-direction: "
        f"**{'PASS' if gates['published_evaluable_pass'] else 'FAIL'}**"
    )
    lines.append(
        f"- Donor-supported top-50, all comparisons ≥70% same-direction: "
        f"**{'PASS' if gates['donor_supported_top50_pass'] else 'FAIL'}**"
    )
    lines.append(
        f"- Unfiltered top-50 (retained, not a valid test): "
        f"**{'PASS' if gates['unfiltered_top50_pass'] else 'FAIL'}**"
    )
    lines.append(
        f"- LODO on donor-supported ranking: "
        f"**{'PASS' if gates['lodo_donor_supported_pass'] else 'FAIL'}**"
    )
    lines.append(
        f"- LODO on unfiltered ranking (retained; LUSC n=5 overlap is a cohort-size artifact, not a comparison finding): "
        f"**{'PASS' if gates['lodo_unfiltered_pass'] else 'FAIL'}**"
    )
    lines.append("")
    pub = gates.get("published")
    if pub is not None and len(pub):
        lines.append("### Published leading genes (donor support)")
        lines.append("")
        cols = [
            "comparison_label",
            "Gene_name",
            "n_donors_evaluable",
            "min_donors_required",
            "donor_supported",
            "sign_concordance",
            "same_direction_pass_70",
            "n_cells_total",
        ]
        lines.append(pub[cols].to_string(index=False))
        lines.append("")
    lines.append("### Gate table (denominators beside fractions)")
    lines.append("")
    lines.append(
        "scored_set `n_genes_in_set` is the denominator. "
        "`n_donors_min/median/max` describe donor support inside that set."
    )
    lines.append("")
    lines.append(gates["by_comparison"].to_string(index=False))
    lines.append("")
    lines.append("### Leave-one-donor-out")
    lines.append("")
    show_cols = [
        c
        for c in [
            "ranking",
            "comparison_label",
            "n_folds",
            "n_ranked",
            "min_donors_required",
            "top50_overlap",
            "fold_overlap_ge_70",
            "gene_retained_in_70pct_folds",
        ]
        if c in gates["lodo_comparison"].columns
    ]
    lines.append(gates["lodo_comparison"][show_cols].to_string(index=False))
    lines.append("")
    lines.append("## Donor counts")
    lines.append("")
    lines.append("- LUAD test donors: 19 (floor 10)")
    lines.append("- Normal test donors: 12 (floor 6)")
    lines.append("- LUSC test donors: 5 (floor 3)")
    lines.append("")
    lines.append("## Limitations")
    lines.append("")
    lines.append("- This is a statistical donor-stability summary, not a biological target claim.")
    lines.append("- Unfiltered LUSC→NORMAL LODO fail tracks five-donor cohort size plus n=1 ranking, not a comparison-specific biological finding.")
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
        effects["min_donors_required"] = effects["source"].map(min_donors_required)
        effects["donor_supported"] = effects["n_donors_evaluable"] >= effects["min_donors_required"]
        lodo = pd.read_csv(output / "leave_one_donor_out_stability.csv")
        published = published_support(effects, top_genes) if top_genes is not None else None
        gates = score_gates(effects, lodo, top_genes, published)
        effects.to_csv(output / "donor_gene_effects.csv", index=False)
        gates["by_comparison"].to_csv(output / "gate_score.csv", index=False)
        if published is not None:
            published.to_csv(output / "published_gene_donor_support.csv", index=False)
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
                    "published_evaluable_pass": gates["published_evaluable_pass"],
                    "donor_supported_top50_pass": gates["donor_supported_top50_pass"],
                    "lodo_donor_supported_pass": gates["lodo_donor_supported_pass"],
                }
            ),
            flush=True,
        )
        return 0
    manifest = load_manifest(tables)
    sums = accumulate(raw_root, manifest)
    effects, donor_level = summarize(sums, token_to_ensembl, ensembl_to_name)
    lodo = pd.concat(
        [
            lodo_stability(donor_level, ranking="unfiltered"),
            lodo_stability(donor_level, ranking="donor_supported"),
        ],
        ignore_index=True,
    )
    published = published_support(effects, top_genes) if top_genes is not None else None
    gates = score_gates(effects, lodo, top_genes, published)
    if published is not None:
        published.to_csv(output / "published_gene_donor_support.csv", index=False)
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
                "published_evaluable_pass": gates["published_evaluable_pass"],
                "donor_supported_top50_pass": gates["donor_supported_top50_pass"],
                "unfiltered_top50_pass": gates["unfiltered_top50_pass"],
                "lodo_unfiltered_pass": gates["lodo_unfiltered_pass"],
                "lodo_donor_supported_pass": gates["lodo_donor_supported_pass"],
            }
        ),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
