#!/usr/bin/env python3
"""Ambient-RNA / lineage cross-expression diagnostic - Priority 1, NSCLC T-cell cohort.

WHAT THIS MEASURES, AND WHAT IT DELIBERATELY DOES NOT
-----------------------------------------------------
Priority 1 was originally scoped as "ambient AND doublet lineage flags". Two thirds of
that is **unmeasurable on the confirmed object**, established by reading the file
directly rather than by inference:

  * ``obs["doublet_status"]``  - a single distinct value (``singlet``) across all cells.
    The object is pre-filtered upstream, so a doublet flag returns zero doublets:
    true, and meaningless.
  * ``obs["pct_counts_mito"]`` / ``obs["total_counts_mito"]`` - identically zero.
    ``var["mito"]`` still flags the mitochondrial genes, so those genes were known and
    were evidently filtered out before the obs statistics were computed. **The zeros
    are an artifact of upstream filtering, not a bug in this script.**

The failure mode this module exists to prevent is precise and quiet: running the
original method against this object would **not error**. It would emit a
clean-looking *"no doublets, 0% mito"* derived from columns that are constant by
construction - the most reassuring possible wrong answer. So this module does not
compute a doublet rate or a mitochondrial fraction at all. It **states the two
retirements explicitly**, in the manifest and in the report, and it *derives* that
statement from the object at runtime rather than hard-coding it (this repo has already
been bitten once by a hard-coded literal - see the ``n_ci_excludes_zero`` fix).

``total_counts`` and ``n_genes_by_counts`` vary widely across cells, and it is tempting
to treat that variation as contamination. **It is sequencing depth.** They are carried
in the output as labelled metadata and are named as depth, never as a contamination
measure.

WHAT IS MEASURED: per-cell lineage cross-expression - ``log1p(CPM)`` module scores over
the ambient and T-cell-intrinsic anchor panels, and the ratio between them. A cell whose
ambient-lineage module is high relative to its T-cell module is carrying transcripts
foreign to its own lineage, which is what ambient contamination looks like at the level
of a single cell.

TWO TRAPS THIS MODULE IS BUILT TO FAIL LOUDLY ON
-------------------------------------------------
1. **A sibling h5ad has the identical 21,000 x 17,764 shape.**
   ``.../h5ad_for_tokenization/balanced_lusc_max_7000_per_disease_tcells_geneformer.h5ad``
   is the tokenization input: no ``doublet_status``, no mito columns, no
   ``layers["count"]``. **Any shape or cell-count assertion passes on it.** Validation
   here therefore keys on ``layers`` and ``obs`` *keys*, never on dimensions.
2. **``X`` is normalised/log; ``layers["count"]`` holds the raw integral counts.**
   Feeding ``X`` to a counts method produces numbers, not an error. The count matrix is
   selected by name and then *verified* to be integral, finite and non-negative.
   Note also that a second layer, ``counts_length_scaled``, sits next to ``count`` - it
   is length-scaled, not raw. If ``count`` is missing this module raises and names that
   neighbour rather than falling back to it. **A fallback that substitutes a different
   input and still reports success is worse than a hard failure.**

Execution is opt-in (``NSCLC_AMBIENT_ALLOW_RUN=1``). The guard has earned itself twice:
the original default path did not exist, and the correct path does not support the
original method.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

# Confirmed by direct inspection on thinkstation2 (592 MiB, 21,000 cells x 17,764 genes),
# not by naming-convention analogy. The previous default was an unverified guess.
CONFIRMED_H5AD = Path(
    "/home/kaisar/workspace/KD/tcell_luad_lusc_normal_10k_from_atlas/outputs"
    "/balanced_lusc_max_7000_per_disease_tcells.h5ad"
)
H5AD = Path(os.environ.get("NSCLC_H5AD", CONFIRMED_H5AD))
ALLOW_RUN_ENV = "NSCLC_AMBIENT_ALLOW_RUN"

COUNT_LAYER = "count"
# Present in the same object and easy to grab by accident; never a substitute for `count`.
LENGTH_SCALED_LAYER = "counts_length_scaled"

CANDIDATES = Path(os.environ.get(
    "NSCLC_CANDIDATES",
    Path(__file__).resolve().parents[1] / "results" / "analysis" / "tables" / "published_gene_donor_support.csv",
))
OUT_DIR = Path(os.environ.get(
    "NSCLC_AMBIENT_OUT_DIR", Path(__file__).resolve().parents[1] / "results" / "biology" / "tables"
))

# obs keys the measure consumes. The first group is metadata the sibling object also has;
# the second group is absent from the sibling AND is what the retirement statement is
# derived from - requiring it means the retirement claim is earned from the object rather
# than asserted about it.
REQUIRED_OBS_METADATA = ("cell_id", "celltype", "individual", "disease", "split")
REQUIRED_OBS_RETIREMENT = (
    "doublet_status", "pct_counts_mito", "total_counts_mito", "total_counts", "n_genes_by_counts",
)
# Depth, not contamination. Carried through as metadata and named as such.
DEPTH_OBS = ("total_counts", "n_genes_by_counts")

FEATURE_NAME_COL = "feature_name"

# Per-cell module scores use counts-per-million, per the Priority 1 spec.
CPM_SCALE = 1e6
# The retained gene-level port keeps its original 1e4 (counts-per-10k) scaling: changing it
# would be a modelling change to an already-calibrated method, not a port.
GENE_LEVEL_SCALE = 1e4
# Floor for the ambient:T-cell ratio denominator. A cell with a zero T-cell module score
# would otherwise divide by zero; the floor is recorded in the manifest so the resulting
# value is never mistaken for a measurement.
TCELL_SCORE_FLOOR = 1e-3

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


class ObjectContractError(RuntimeError):
    """The AnnData handed in is not the object this diagnostic is defined against."""


def require_run_opt_in(env: os._Environ | dict | None = None) -> None:
    """Execution is explicit opt-in. Raises SystemExit unless NSCLC_AMBIENT_ALLOW_RUN=1."""
    env = os.environ if env is None else env
    if env.get(ALLOW_RUN_ENV) != "1":
        raise SystemExit(
            f"Refusing to run: set {ALLOW_RUN_ENV}=1 to execute this diagnostic.\n"
            "This guard is not boilerplate - it has caught two distinct real errors:\n"
            "  (1) the original default h5ad path did not exist, and\n"
            "  (2) the correct object does not support the original doublet/mito method\n"
            "      (doublet_status is constant, mito columns are all zero), so the method\n"
            "      would have returned a clean-looking 'no doublets, 0% mito' rather than\n"
            "      failing.\n"
            f"Confirmed object: {CONFIRMED_H5AD}\n"
            f"Override the path with NSCLC_H5AD only if you have verified the replacement."
        )


def select_count_matrix(adata) -> sp.csr_matrix:
    """Return raw integral counts from ``layers['count']``, verified - never ``X``.

    ``X`` on this object is normalised/log-transformed. It has the right shape and dtype,
    so a counts method fed ``X`` returns numbers rather than raising. The matrix is
    therefore chosen by layer name and then checked to actually BE counts.
    """
    layers = getattr(adata, "layers", None)
    # anndata 0.13.x yields a spurious `None` alongside the real names in
    # Layers.keys(), even when the on-disk group holds only string keys (observed
    # 2026-09-09 on the confirmed object: h5py reports {count, counts_length_scaled},
    # anndata reports those plus None). Sorting raw keys raises TypeError there, and a
    # non-string key is never a layer worth selecting, so filter before sorting.
    layer_keys = sorted(k for k in (layers.keys() if layers is not None else [])
                        if isinstance(k, str))
    if COUNT_LAYER not in layer_keys:
        hint = ""
        if LENGTH_SCALED_LAYER in layer_keys:
            hint = (f" Layer {LENGTH_SCALED_LAYER!r} is present but is length-scaled, "
                    f"not raw counts - it is NOT a substitute.")
        raise ObjectContractError(
            f"layers[{COUNT_LAYER!r}] is missing (layers present: {layer_keys}). "
            "This is the signature of the tokenization-input sibling object, which has "
            "the IDENTICAL 21,000 x 17,764 shape - so a shape or cell-count check would "
            "have passed here. Refusing to fall back to X, which is normalised/log." + hint
        )
    matrix = layers[COUNT_LAYER]
    matrix = sp.csr_matrix(matrix) if sp.issparse(matrix) else sp.csr_matrix(np.asarray(matrix))
    if matrix.shape != adata.shape:
        raise ObjectContractError(
            f"layers[{COUNT_LAYER!r}] shape {matrix.shape} != object shape {adata.shape}"
        )
    data = matrix.data
    if data.size:
        if not np.all(np.isfinite(data)):
            raise ObjectContractError(f"layers[{COUNT_LAYER!r}] contains non-finite values")
        if np.any(data < 0):
            raise ObjectContractError(f"layers[{COUNT_LAYER!r}] contains negative values")
        if not np.array_equal(data, np.rint(data)):
            raise ObjectContractError(
                f"layers[{COUNT_LAYER!r}] is not integral - it is normalised or scaled, "
                "not raw counts. Refusing to treat it as counts."
            )
    return matrix


def build_feature_map(adata) -> dict[str, int]:
    """symbol -> column index, from the in-object ``var['feature_name']``.

    The var index is Ensembl (``ENSG00000121410``...), which is why symbol lookups
    (``EPCAM``, ``PTPRC``, ``CD3E``) all miss against it. Using the in-object map removes
    a whole class of mismatch, and needs no external gc104M dictionary.
    """
    if FEATURE_NAME_COL not in adata.var.columns:
        raise ObjectContractError(
            f"var[{FEATURE_NAME_COL!r}] is missing; symbol lookups against the Ensembl "
            "var index would silently match nothing."
        )
    names = pd.Series(adata.var[FEATURE_NAME_COL]).astype(str).to_numpy()
    if len(names) != adata.n_vars:
        raise ObjectContractError(
            f"var[{FEATURE_NAME_COL!r}] length {len(names)} != n_vars {adata.n_vars}"
        )
    blank = [i for i, n in enumerate(names) if n.strip() == "" or n.lower() == "nan"]
    if blank:
        raise ObjectContractError(
            f"var[{FEATURE_NAME_COL!r}] has {len(blank)} empty/NaN entries "
            f"(first at index {blank[0]})"
        )
    unique, counts = np.unique(names, return_counts=True)
    if len(unique) != len(names):
        duplicated = unique[counts > 1][:10].tolist()
        raise ObjectContractError(
            f"var[{FEATURE_NAME_COL!r}] is not unique: {int((counts > 1).sum())} duplicated "
            f"symbols, e.g. {duplicated}. A symbol->column map would silently drop columns."
        )
    return {name: index for index, name in enumerate(names)}


def validate_object(adata) -> tuple[sp.csr_matrix, dict[str, int]]:
    """Full contract check. Rejects the sibling object on keys, never on dimensions."""
    missing_obs = [k for k in REQUIRED_OBS_METADATA + REQUIRED_OBS_RETIREMENT
                   if k not in adata.obs.columns]
    if missing_obs:
        raise ObjectContractError(
            f"obs is missing required columns {missing_obs}. The tokenization-input sibling "
            "object has exactly this signature (obs limited to cell_id, celltype, disease, "
            "filter_pass, individual, n_counts, split) and the IDENTICAL shape, so a "
            "dimension check cannot tell the two objects apart."
        )
    counts = select_count_matrix(adata)
    feature_map = build_feature_map(adata)
    return counts, feature_map


def resolve_panel(feature_map: dict[str, int], symbols: list[str]) -> tuple[list[str], list[int]]:
    """Anchor symbols present in the object, and their column indices."""
    present = [s for s in symbols if s in feature_map]
    return present, [feature_map[s] for s in present]


def module_score(counts: sp.csr_matrix, columns: list[int],
                 library_size: np.ndarray | None = None) -> np.ndarray:
    """Mean ``log1p(CPM)`` across a panel, per cell.

    CPM is computed from the count matrix's own row sums, so the score is self-consistent
    with the matrix it is derived from rather than depending on an obs column that may
    have been produced under different upstream filtering.
    """
    if not columns:
        raise ValueError("empty anchor panel: no module score is defined")
    if library_size is None:
        library_size = np.asarray(counts.sum(axis=1)).ravel()
    scaling = sp.diags(CPM_SCALE / np.maximum(library_size.astype(float), 1.0))
    panel = sp.csr_matrix(scaling @ sp.csr_matrix(counts[:, columns], dtype=np.float64))
    panel.data = np.log1p(panel.data)
    return np.asarray(panel.sum(axis=1)).ravel() / len(columns)


def constant_component_report(adata) -> dict:
    """Derive the doublet/mito retirement statement FROM THE OBJECT.

    Deliberately not hard-coded. If a future object carries real doublet calls or
    non-zero mito fractions, the retirement does not apply and this says so loudly
    instead of repeating a stale claim.
    """
    obs = adata.obs
    doublet_values = sorted(pd.Series(obs["doublet_status"]).astype(str).unique().tolist())
    doublet_status = ("constant_by_construction" if len(doublet_values) == 1
                      else "present_and_variable__retirement_does_not_apply")

    mito_detail = {}
    for column in ("pct_counts_mito", "total_counts_mito"):
        values = np.asarray(pd.Series(obs[column]).astype(float))
        mito_detail[column] = {
            "all_zero": bool(np.all(values == 0)),
            "min": float(values.min()), "max": float(values.max()),
            "n_distinct": int(np.unique(values).size),
        }
    mito_status = ("constant_zero_by_construction"
                   if all(d["all_zero"] for d in mito_detail.values())
                   else "not_constant_zero__retirement_does_not_apply")

    var_mito_flagged = None
    if "mito" in adata.var.columns:
        var_mito_flagged = int(np.asarray(adata.var["mito"]).astype(bool).sum())

    # The mito-gene clause is derived, not asserted. On the confirmed object var["mito"]
    # exists but flags ZERO genes - the mitochondrial genes were dropped from the gene set
    # entirely, not merely flagged - so saying "var still flags the mito genes" would be
    # contradicted by the very number reported beside it.
    if var_mito_flagged is None:
        mito_gene_clause = ("var has no 'mito' flag column, so why the obs mito statistics "
                            "are zero cannot be established from this object alone.")
    elif var_mito_flagged == 0:
        mito_gene_clause = (
            "var['mito'] is present but flags ZERO genes: the mitochondrial genes are absent "
            "from this object's gene set altogether. The zero obs statistics are the residue "
            "of that upstream filtering, not a measurement taken on retained mito genes."
        )
    else:
        mito_gene_clause = (
            f"var['mito'] still flags {var_mito_flagged} genes, so the mitochondrial genes "
            "were known and were filtered before the obs statistics were taken."
        )

    return {
        "doublet_component_status": doublet_status,
        "doublet_status_values": doublet_values,
        "mito_component_status": mito_status,
        "mito_component_detail": mito_detail,
        "var_mito_flag_present": "mito" in adata.var.columns,
        "n_var_genes_flagged_mito": var_mito_flagged,
        "retirement_statement": (
            "Doublet and mitochondrial components of Priority 1 are RETIRED on this object. "
            "They are not reported as zero: doublet_status carries one distinct value and the "
            "mito obs columns are identically zero, both by upstream construction, so a "
            "doublet rate or mito fraction computed here would be vacuous rather than "
            "reassuring. " + mito_gene_clause
        ),
        "depth_columns_not_contamination": list(DEPTH_OBS),
        "depth_note": (
            "total_counts and n_genes_by_counts vary widely across cells. That variation is "
            "sequencing depth, NOT contamination, and is carried here as metadata only."
        ),
    }


def compute_cell_scores(adata, counts: sp.csr_matrix,
                        feature_map: dict[str, int]) -> tuple[pd.DataFrame, dict]:
    """Per-cell lineage cross-expression table plus its provenance block."""
    ambient_symbols, ambient_columns = resolve_panel(feature_map, KNOWN_AMBIENT)
    tcell_symbols, tcell_columns = resolve_panel(feature_map, KNOWN_TCELL)
    if not ambient_symbols or not tcell_symbols:
        raise ObjectContractError(
            f"anchor panels did not resolve against var[{FEATURE_NAME_COL!r}]: "
            f"{len(ambient_symbols)} ambient, {len(tcell_symbols)} T-cell. "
            "A symbol panel resolving to nothing is the signature of matching symbols "
            "against the Ensembl var index."
        )

    library_size = np.asarray(counts.sum(axis=1)).ravel()
    ambient_score = module_score(counts, ambient_columns, library_size)
    tcell_score = module_score(counts, tcell_columns, library_size)

    frame = pd.DataFrame({
        "cell_id": pd.Series(adata.obs["cell_id"]).astype(str).to_numpy(),
        "individual": pd.Series(adata.obs["individual"]).astype(str).to_numpy(),
        "celltype": pd.Series(adata.obs["celltype"]).astype(str).to_numpy(),
        "disease": pd.Series(adata.obs["disease"]).astype(str).to_numpy(),
        "split": pd.Series(adata.obs["split"]).astype(str).to_numpy(),
        "ambient_module_score": ambient_score,
        "tcell_module_score": tcell_score,
        # Difference of log scores: the numerically stable companion to the ratio, with
        # no divide-by-near-zero behaviour.
        "ambient_minus_tcell": ambient_score - tcell_score,
        "ambient_to_tcell_ratio": ambient_score / np.maximum(tcell_score, TCELL_SCORE_FLOOR),
        # Depth metadata. NOT a contamination measure - see depth_note in the manifest.
        "depth_total_counts": np.asarray(pd.Series(adata.obs["total_counts"]).astype(float)),
        "depth_n_genes_by_counts": np.asarray(
            pd.Series(adata.obs["n_genes_by_counts"]).astype(float)),
        "count_layer_library_size": library_size,
    })

    provenance = {
        "n_cells": int(frame.shape[0]),
        "count_layer": COUNT_LAYER,
        "count_layer_is_verified_integral": True,
        "normalisation_per_cell_scores": f"log1p(counts / library_size * {CPM_SCALE:.0f}) [CPM]",
        "library_size_source": "row sums of layers['count'] (not obs['total_counts'])",
        "feature_map_source": f"in-object var[{FEATURE_NAME_COL!r}] (no external gc104M dictionary)",
        "anchors_ambient_present": ambient_symbols,
        "anchors_ambient_missing": [s for s in KNOWN_AMBIENT if s not in feature_map],
        "anchors_tcell_present": tcell_symbols,
        "anchors_tcell_missing": [s for s in KNOWN_TCELL if s not in feature_map],
        "tcell_score_floor_for_ratio": TCELL_SCORE_FLOOR,
        "n_cells_tcell_score_below_floor": int(np.sum(tcell_score < TCELL_SCORE_FLOOR)),
    }
    return frame, provenance


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


def gene_level_candidate_merge(adata, counts: sp.csr_matrix, feature_map: dict[str, int],
                               candidates_path: Path) -> tuple[pd.DataFrame, dict]:
    """Retained gene-level ambient-risk port. Scoring model and threshold UNCHANGED.

    Two inputs changed, both mandated and both corrections rather than tuning: the matrix
    is ``layers['count']`` (it was ``X``, which is log-normalised - a counts method fed X
    returns numbers, not an error), and symbols come from the in-object feature map
    instead of the external gc104M dictionaries.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import cross_val_score
    from sklearn.preprocessing import StandardScaler

    library_size = np.asarray(counts.sum(axis=1)).ravel()
    scaling = sp.diags(GENE_LEVEL_SCALE / np.maximum(library_size, 1))
    norm = sp.csr_matrix(scaling @ counts)
    norm.data = np.log1p(norm.data)

    detect_frac = np.asarray((counts > 0).sum(axis=0)).ravel() / adata.n_obs
    mean_expr = np.asarray(norm.mean(axis=0)).ravel()
    subtype_f = group_f_statistic(norm, adata.obs["celltype"])
    donor_f = group_f_statistic(norm, adata.obs["individual"])

    rank_lib = pd.Series(library_size).rank().to_numpy()
    rank_lib = (rank_lib - rank_lib.mean()) / rank_lib.std()
    col_mean = np.asarray(norm.mean(axis=0)).ravel()
    numerator = np.asarray(norm.T @ rank_lib).ravel() - adata.n_obs * col_mean * rank_lib.mean()
    sumsq = np.asarray(norm.multiply(norm).sum(axis=0)).ravel()
    denom = np.sqrt(np.maximum(sumsq - adata.n_obs * col_mean ** 2, 1e-12)) * np.sqrt(adata.n_obs)
    libsize_corr = numerator / np.maximum(denom, 1e-12)

    index_to_symbol = {v: k for k, v in feature_map.items()}
    features = pd.DataFrame({
        "ensembl_id": pd.Index(adata.var_names),
        "gene": [index_to_symbol[i] for i in range(adata.n_vars)],
        "detect_frac": detect_frac,
        "mean_expr": mean_expr,
        "log_subtype_f": np.log1p(subtype_f),
        "log_donor_f": np.log1p(donor_f),
        "libsize_corr": libsize_corr,
    })
    features["breadth_over_depth"] = features.detect_frac / np.maximum(features.mean_expr, 1e-6)
    features["subtype_over_donor"] = features.log_subtype_f / np.maximum(features.log_donor_f, 1e-6)

    present = set(features.gene)
    ambient_set = [g for g in KNOWN_AMBIENT if g in present]
    tcell_set = [g for g in KNOWN_TCELL if g in present]

    scored = features[features.detect_frac >= 0.005].copy().reset_index(drop=True)
    columns = ["log_subtype_f", "log_donor_f", "libsize_corr",
               "breadth_over_depth", "subtype_over_donor", "mean_expr"]
    labelled = scored[scored.gene.isin(ambient_set + tcell_set)].copy()
    labelled["is_ambient"] = labelled.gene.isin(ambient_set).astype(int)

    scaler = StandardScaler().fit(labelled[columns])
    model = LogisticRegression(max_iter=2000, class_weight="balanced")
    auc = cross_val_score(model, scaler.transform(labelled[columns]), labelled.is_ambient,
                          cv=5, scoring="roc_auc")
    model.fit(scaler.transform(labelled[columns]), labelled.is_ambient)
    scored["ambient_risk"] = model.predict_proba(scaler.transform(scored[columns]))[:, 1]
    scored["ambient_pct"] = 100 * scored.ambient_risk.rank(pct=True)

    threshold = float(np.percentile(scored[scored.gene.isin(ambient_set)].ambient_risk, 25))
    merged = load_candidates(candidates_path).merge(
        scored[["ensembl_id", "ambient_risk", "ambient_pct", "detect_frac",
                "log_subtype_f", "log_donor_f"]],
        left_on="Ensembl_ID", right_on="ensembl_id", how="left",
    ).drop(columns="ensembl_id")
    merged["ambient_flag"] = np.where(
        merged.ambient_risk.isna(), "not_scored",
        np.where(merged.ambient_risk >= threshold, "AMBIENT_RISK", "ok"),
    )
    provenance = {
        "cv_auc_mean": float(auc.mean()), "cv_auc_std": float(auc.std()),
        "flag_threshold_25th_pct_of_ambient_anchors": threshold,
        "n_genes_scored": int(len(scored)),
        "normalisation": f"log1p(counts / library_size * {GENE_LEVEL_SCALE:.0f}) [CP10K, as ported]",
        "auc_warning": ("anchors separate poorly; treat the score as uninformative"
                        if auc.mean() < 0.8 else None),
    }
    return merged.sort_values("ambient_risk", ascending=False), provenance


def main() -> int:
    require_run_opt_in()

    import anndata as ad

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Loading {H5AD}")
    adata = ad.read_h5ad(H5AD)
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    counts, feature_map = validate_object(adata)
    print(f"  count layer '{COUNT_LAYER}' verified integral, finite, non-negative")

    retirements = constant_component_report(adata)
    print(f"  doublet component: {retirements['doublet_component_status']} "
          f"{retirements['doublet_status_values']}")
    print(f"  mito component:    {retirements['mito_component_status']}")

    cells, cell_provenance = compute_cell_scores(adata, counts, feature_map)
    cells.to_csv(OUT_DIR / "cell_contamination_scores.csv", index=False)
    print(f"  wrote cell_contamination_scores.csv ({len(cells):,} cells)")

    manifest = {
        "diagnostic": "Priority 1 - per-cell lineage cross-expression (ambient RNA)",
        "h5ad": str(H5AD),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "scope_note": (
            "Priority 1 on this object supports expression-based lineage cross-expression "
            "ONLY. The doublet and mito components are retired (below). Perturbation ranking "
            "is NOT rerun by this script."
        ),
        **retirements,
        "per_cell": cell_provenance,
    }

    if CANDIDATES.exists():
        merged, gene_provenance = gene_level_candidate_merge(adata, counts, feature_map, CANDIDATES)
        merged.to_csv(OUT_DIR / "candidate_gene_ambient_risk.csv", index=False)
        manifest["gene_level_candidate_merge"] = gene_provenance
        print(f"  wrote candidate_gene_ambient_risk.csv ({len(merged):,} genes), "
              f"cv AUC {gene_provenance['cv_auc_mean']:.3f}")
    else:
        raise SystemExit(
            f"Candidate table {CANDIDATES} not found. Refusing to skip the gene-level merge "
            "silently - point NSCLC_CANDIDATES at it or remove it deliberately."
        )

    (OUT_DIR / "ambient_risk_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"\nWrote outputs to {OUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
