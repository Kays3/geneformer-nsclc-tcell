#!/usr/bin/env python3
"""Ambient-RNA risk diagnostic for the published-leading candidates - PORTED, UNRUN.

STATUS (2026-09-08): drafted while the entire compute fleet
(petadimensionlab-super-server, dgx, thinkstation1/2/3) was offline in a
site-level outage. This script has NEVER BEEN EXECUTED against real data -
its H5AD default path is an unverified guess by naming-convention analogy,
not a confirmed location. Before running it:

1. Confirm the compute fleet is reachable again.
2. Confirm the actual h5ad (or equivalent AnnData object with obs columns
   `celltype` and `individual`) for the
   tcell_luad_lusc_normal_luscmax7000_heldout_allgene_perturbation cohort,
   and point NSCLC_H5AD at it. Do NOT run this against a substitute dataset -
   a diagnostic computed from the wrong cohort would answer a different
   question and look like an answer to this one.
3. Only the h5ad-dependent feature computation (`main()`) is untested. The
   candidate-loading/merge logic (`load_candidates`) IS tested against the
   real local `published_gene_donor_support.csv` in
   `test_ambient_risk_diagnostic.py`.

This is a straightforward port of
`geneformer-lung-tcell/sclc_validation/primary_test_perturbation/scripts/ambient_risk_diagnostic.py`
onto this repo's candidate table schema and this cohort's data paths. The
anchor panels, scoring model, and methodology are UNCHANGED from that script -
see its own docstring for the rationale (contamination is lineage-foreign,
unstructured across cell states, sample-driven, and library-size dependent).
Extending the anchor panel with the nine markers used in
`classify_candidate_lineage_flags.py` was deliberately NOT done here: that
would be a modeling change to an already-calibrated method, and it cannot be
validated (cross-validated AUC) without the data this script itself needs.
"""
from __future__ import annotations

import json
import os
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

# UNVERIFIED GUESS by naming-convention analogy with the sibling repo's
# `~/workspace/KD/sclc_luad_normal_htan_finetune/data/htan_sclc_luad_normal_tcells_prepared.h5ad`.
# Confirm the real path once the compute fleet returns - do not assume this is correct.
H5AD = Path(os.environ.get(
    "NSCLC_H5AD",
    Path.home() / "workspace/KD/tcell_luad_lusc_normal_luscmax7000_heldout_allgene_perturbation"
    "/data/tcell_luad_lusc_normal_prepared.h5ad",
))
GENEFORMER_DIR = Path(os.environ.get("NSCLC_GENEFORMER_DIR", Path.home() / "Geneformer/geneformer"))
CANDIDATES = Path(os.environ.get(
    "NSCLC_CANDIDATES",
    Path(__file__).resolve().parents[1] / "results" / "analysis" / "tables" / "published_gene_donor_support.csv",
))
OUT_DIR = Path(os.environ.get(
    "NSCLC_AMBIENT_OUT_DIR", Path(__file__).resolve().parents[1] / "results" / "biology" / "tables"
))
FIG_DIR = Path(__file__).resolve().parents[1] / "results" / "biology" / "figures" / "ambient_risk"

# UNCHANGED from geneformer-lung-tcell's ambient_risk_diagnostic.py - do not
# extend without being able to re-run the AUC calibration this script prints.
KNOWN_AMBIENT = [
    "HBB", "HBA1", "HBA2", "HBD", "ALAS2", "AHSP",
    "SFTPC", "SFTPB", "SFTPA1", "SFTPA2", "SCGB1A1", "SCGB3A1", "SCGB3A2", "NAPSA",
    "WFDC2", "MUC1", "AGER", "CLDN18", "EPCAM", "KRT8", "KRT18", "KRT19", "SLPI",
    "LYZ", "S100A8", "S100A9", "S100A12", "CD68", "MARCO", "FCN1", "VCAN",
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "ACTA2", "PECAM1", "VWF",
]
KNOWN_TCELL = [
    "CD3D", "CD3E", "CD3G", "CD247", "TRAC", "TRBC1", "TRBC2", "CD2", "CD5", "CD6",
    "CD7", "CD28", "LCK", "ZAP70", "LAT", "ITK", "THEMIS", "SKAP1", "CD8A", "CD8B",
    "IL7R", "CCR7", "TCF7", "LEF1", "SELL", "GZMA", "GZMK", "PRF1", "NKG7", "CTSW",
    "CD27", "ICOS", "CTLA4", "PDCD1", "FOXP3", "IKZF2", "RUNX3", "BCL11B", "CD69",
]


def load_symbol_maps(geneformer_dir: Path) -> tuple[dict[str, str], dict[str, str]]:
    """gene symbol -> Ensembl ID and back, from the local gc104M dictionaries
    used elsewhere in this repo's evaluation pipeline (evaluate_donor_gene_effects.py's
    load_token_maps), rather than the sibling repo's stats-table glob - this
    repo has no equivalent `stats/*/heldout_allgene_*.csv` tree locally."""
    name_path = geneformer_dir / "gene_name_id_dict_gc104M.pkl"
    if not name_path.exists():
        raise SystemExit(f"Missing {name_path}; point NSCLC_GENEFORMER_DIR at the gc104M dictionaries")
    obj = pickle.loads(name_path.read_bytes())
    sample = next(iter(obj.items()))
    if isinstance(sample[0], str) and str(sample[0]).startswith("ENSG"):
        ensembl_to_symbol = {str(k): str(v) for k, v in obj.items()}
    else:
        ensembl_to_symbol = {str(v): str(k) for k, v in obj.items() if str(v).startswith("ENSG")}
    symbol_to_ensembl = {v: k for k, v in ensembl_to_symbol.items()}
    return symbol_to_ensembl, ensembl_to_symbol


def load_candidates(path: Path) -> pd.DataFrame:
    """Load and collapse the published-leading candidate table to one row per
    gene, matching this repo's published_gene_donor_support.csv schema (no
    class_label/program column exists here, unlike the sibling repo's
    immune_cancer_candidates_with_donor_robustness.csv)."""
    candidates = pd.read_csv(path)
    required = {"Gene_name", "Ensembl_ID", "comparison_label", "donor_supported"}
    missing = required.difference(candidates.columns)
    if missing:
        raise SystemExit(f"Candidate table missing columns: {sorted(missing)}")
    per_gene = candidates.groupby(["Gene_name", "Ensembl_ID"], as_index=False).agg(
        n_comparisons=("comparison_label", "nunique"),
        donor_supported_in_all=("donor_supported", "all"),
    )
    return per_gene


def group_f_statistic(counts: sp.csr_matrix, labels: pd.Series) -> np.ndarray:
    """One-way F statistic per gene across the given grouping, on log1p-CPM values.
    Unchanged from the sibling script."""
    codes = pd.Categorical(labels).codes
    n_groups = int(codes.max()) + 1
    n_cells, n_genes = counts.shape
    indicator = sp.csr_matrix(
        (np.ones(n_cells), (codes, np.arange(n_cells))), shape=(n_groups, n_cells)
    )
    group_n = np.asarray(indicator.sum(axis=1)).ravel()
    group_sum = np.asarray(sp.csr_matrix(indicator @ counts).todense())
    group_sumsq = np.asarray(sp.csr_matrix(indicator @ counts.multiply(counts)).todense())

    group_mean = group_sum / group_n[:, None]
    grand_mean = np.asarray(counts.mean(axis=0)).ravel()

    between = (group_n[:, None] * (group_mean - grand_mean) ** 2).sum(axis=0) / max(n_groups - 1, 1)
    total_sumsq = group_sumsq.sum(axis=0)
    within_sumsq = total_sumsq - (group_n[:, None] * group_mean ** 2).sum(axis=0)
    within = within_sumsq / max(n_cells - n_groups, 1)
    return between / np.maximum(within, 1e-12)


def main() -> int:
    """UNRUN - see module docstring. Left as a faithful port of the
    calibration/scoring pipeline; do not execute against a substitute
    dataset."""
    raise SystemExit(
        "ambient_risk_diagnostic.py is a drafted, unrun port - the compute fleet was offline "
        "when this was written (2026-09-08) and NSCLC_H5AD is an unverified guess, not a "
        "confirmed path. Verify the real h5ad location and remove this guard before running."
    )

    import anndata as ad  # noqa: F401  (imported only once the guard above is removed)
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import cross_val_score
    from sklearn.preprocessing import StandardScaler

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Loading {H5AD}")
    adata = ad.read_h5ad(H5AD)
    raw_counts = adata.X if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    raw_counts = sp.csr_matrix(raw_counts)
    total_counts = np.asarray(raw_counts.sum(axis=1)).ravel()
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    scaling = sp.diags(1e4 / np.maximum(total_counts, 1))
    norm = sp.csr_matrix(scaling @ raw_counts)
    norm.data = np.log1p(norm.data)

    genes = pd.Index(adata.var_names)
    detect_frac = np.asarray((raw_counts > 0).sum(axis=0)).ravel() / adata.n_obs
    mean_expr = np.asarray(norm.mean(axis=0)).ravel()

    print("Computing subtype and donor structure")
    subtype_f = group_f_statistic(norm, adata.obs["celltype"])
    donor_f = group_f_statistic(norm, adata.obs["individual"])

    print("Computing library-size dependence")
    rank_lib = pd.Series(total_counts).rank().to_numpy()
    rank_lib = (rank_lib - rank_lib.mean()) / rank_lib.std()
    dense_proxy = norm.copy()
    col_mean = np.asarray(dense_proxy.mean(axis=0)).ravel()
    numerator = np.asarray(dense_proxy.T @ rank_lib).ravel() - adata.n_obs * col_mean * rank_lib.mean()
    sumsq = np.asarray(dense_proxy.multiply(dense_proxy).sum(axis=0)).ravel()
    denom = np.sqrt(np.maximum(sumsq - adata.n_obs * col_mean ** 2, 1e-12)) * np.sqrt(adata.n_obs)
    libsize_corr = numerator / np.maximum(denom, 1e-12)

    sym2ens, ens2sym = load_symbol_maps(GENEFORMER_DIR)

    features = pd.DataFrame({
        "ensembl_id": genes,
        "gene": [ens2sym.get(e, e) for e in genes],
        "detect_frac": detect_frac,
        "mean_expr": mean_expr,
        "log_subtype_f": np.log1p(subtype_f),
        "log_donor_f": np.log1p(donor_f),
        "libsize_corr": libsize_corr,
    })
    features["breadth_over_depth"] = features.detect_frac / np.maximum(features.mean_expr, 1e-6)
    features["subtype_over_donor"] = features.log_subtype_f / np.maximum(features.log_donor_f, 1e-6)

    present_symbols = set(features.gene)
    ambient_set = [g for g in KNOWN_AMBIENT if g in present_symbols]
    tcell_set = [g for g in KNOWN_TCELL if g in present_symbols]
    print(f"\nAnchors present: {len(ambient_set)} ambient, {len(tcell_set)} T-cell-intrinsic")

    scored = features[features.detect_frac >= 0.005].copy().reset_index(drop=True)
    print(f"Genes with detection >= 0.5%: {len(scored):,}")

    columns = ["log_subtype_f", "log_donor_f", "libsize_corr", "breadth_over_depth", "subtype_over_donor", "mean_expr"]
    labelled = scored[scored.gene.isin(ambient_set + tcell_set)].copy()
    labelled["is_ambient"] = labelled.gene.isin(ambient_set).astype(int)
    print(f"Anchors retained after detection filter: {int(labelled.is_ambient.sum())} ambient, "
          f"{int((1 - labelled.is_ambient).sum())} T-cell")

    scaler = StandardScaler().fit(labelled[columns])
    model = LogisticRegression(max_iter=2000, class_weight="balanced")
    auc = cross_val_score(model, scaler.transform(labelled[columns]), labelled.is_ambient,
                          cv=5, scoring="roc_auc")
    print(f"\nCross-validated AUC separating the anchor sets: {auc.mean():.3f} (+/- {auc.std():.3f})")
    if auc.mean() < 0.8:
        print("  [WARN] anchors separate poorly; treat the score as uninformative", flush=True)

    model.fit(scaler.transform(labelled[columns]), labelled.is_ambient)
    scored["ambient_risk"] = model.predict_proba(scaler.transform(scored[columns]))[:, 1]
    scored["ambient_pct"] = 100 * scored.ambient_risk.rank(pct=True)

    anchor_scores = scored[scored.gene.isin(ambient_set)].ambient_risk
    threshold = float(np.percentile(anchor_scores, 25))
    print(f"Flag threshold (25th pct of ambient anchors): {threshold:.3f}")

    scored.sort_values("ambient_risk", ascending=False).to_csv(OUT_DIR / "ambient_risk_all_genes.csv", index=False)

    per_gene = load_candidates(CANDIDATES)
    merged = per_gene.merge(
        scored[["ensembl_id", "ambient_risk", "ambient_pct", "detect_frac", "log_subtype_f", "log_donor_f"]],
        left_on="Ensembl_ID", right_on="ensembl_id", how="left",
    ).drop(columns="ensembl_id")
    merged["ambient_flag"] = np.where(
        merged.ambient_risk.isna(), "not_scored",
        np.where(merged.ambient_risk >= threshold, "AMBIENT_RISK", "ok"),
    )
    merged = merged.sort_values("ambient_risk", ascending=False)
    merged.to_csv(OUT_DIR / "cell_contamination_scores.csv", index=False)

    (OUT_DIR / "ambient_risk_manifest.json").write_text(json.dumps({
        "note": "Ambient-RISK DIAGNOSTIC, ported from geneformer-lung-tcell, run against the NSCLC cohort.",
        "h5ad": str(H5AD),
        "n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars),
        "anchors_ambient": ambient_set, "anchors_tcell": tcell_set,
        "cv_auc_mean": float(auc.mean()), "cv_auc_std": float(auc.std()),
        "flag_threshold": threshold,
    }, indent=2) + "\n")

    print(f"\nWrote outputs to {OUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
