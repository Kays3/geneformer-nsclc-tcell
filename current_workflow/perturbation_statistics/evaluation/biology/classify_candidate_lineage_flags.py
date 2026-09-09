#!/usr/bin/env python3
"""Priority 1 (ambient RNA / doublet sensitivity), identity-only first pass.

This is deliberately NOT the per-cell contamination/doublet diagnostic the
biology README asks for (`cell_contamination_scores.csv`,
`decontamination_rank_stability.csv`). That analysis needs the original
single-cell object (obs-level celltype/individual/library-size structure) for
this NSCLC/LUSC/LUAD held-out cohort, the way
`sclc_validation/primary_test_perturbation/scripts/ambient_risk_diagnostic.py`
uses `htan_sclc_luad_normal_tcells_prepared.h5ad` in the sibling SCLC repo. No
such object is present in this checkout or reachable from this sandbox.

What this script does instead: flag each published-leading candidate gene by
LINEAGE IDENTITY against a small vetted anchor panel, so the genes the biology
README already names as concerning (SFTPC, SFTPB, NAPSA, MUC1, PIGR, keratins,
FBLN1, DCN, ACKR1) get an explicit, sourced flag rather than sitting
unclassified. This is gene identity, not a per-cell measurement - it cannot
tell you whether a flagged gene's signal comes from ambient contamination in a
specific specimen, a doublet, or genuine low-level expression, and it cannot
promote anything to class 1 (T-cell-intrinsic and donor-stable). It can only
flag lineage-foreign genes for exclusion consideration (class 4) and leave
everything else explicitly pending the real per-cell diagnostic.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

PUBLISHED_SUPPORT = (
    Path(__file__).resolve().parents[1]
    / "results"
    / "analysis"
    / "tables"
    / "published_gene_donor_support.csv"
)
OUT_DIR = Path(__file__).resolve().parents[1] / "results" / "biology" / "tables"

# Ported verbatim from geneformer-lung-tcell's
# sclc_validation/primary_test_perturbation/scripts/ambient_risk_diagnostic.py,
# which calibrated these against a logistic model on real single-cell data
# (cross-validated AUC reported there). Lineage-foreign transcripts a T cell
# does not itself transcribe; detection is contamination by construction.
SIBLING_REPO_VETTED_AMBIENT = {
    "HBB", "HBA1", "HBA2", "HBD", "ALAS2", "AHSP",
    "SFTPC", "SFTPB", "SFTPA1", "SFTPA2", "SCGB1A1", "SCGB3A1", "SCGB3A2", "NAPSA",
    "WFDC2", "MUC1", "AGER", "CLDN18", "EPCAM", "KRT8", "KRT18", "KRT19", "SLPI",
    "LYZ", "S100A8", "S100A9", "S100A12", "CD68", "MARCO", "FCN1", "VCAN",
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "ACTA2", "PECAM1", "VWF",
}
SIBLING_REPO_VETTED_TCELL = {
    "CD3D", "CD3E", "CD3G", "CD247", "TRAC", "TRBC1", "TRBC2", "CD2", "CD5", "CD6",
    "CD7", "CD28", "LCK", "ZAP70", "LAT", "ITK", "THEMIS", "SKAP1", "CD8A", "CD8B",
    "IL7R", "CCR7", "TCF7", "LEF1", "SELL", "GZMA", "GZMK", "PRF1", "NKG7", "CTSW",
    "CD27", "ICOS", "CTLA4", "PDCD1", "FOXP3", "IKZF2", "RUNX3", "BCL11B", "CD69",
}
# Ported verbatim from the sibling repo's axis_consistency_sensitivity.py.
SIBLING_REPO_EXPLICIT_PATTERN = re.compile(
    r"^(?:HBA|HBB$|HBD$|HSPA1|RPL|RPS|MRPL|MRPS|S100A|SFTPA|EEF1)"
)

# Extended for genes this repo's own biology/README.md names by name as
# concerning, or that are unambiguous single-lineage markers not in the
# sibling panel. Each is a standard, single-lineage marker with no meaningful
# controversy in the literature - not a judgment call about this dataset.
EXTENDED_AMBIENT_WITH_REASON = {
    "PIGR": "polymeric immunoglobulin receptor - canonical glandular/epithelial secretory marker",
    "ACKR1": "Duffy antigen receptor - canonical venular endothelial marker",
    "FBLN1": "fibulin-1 - canonical stromal/ECM marker",
    "KRT7": "keratin 7 - canonical epithelial marker (same family as anchored KRT8/18/19)",
    "KRT17": "keratin 17 - canonical epithelial marker (same family as anchored KRT8/18/19)",
    "MGP": "matrix Gla protein - canonical vascular/stromal marker",
    "MMP12": "macrophage metalloelastase - canonical myeloid/macrophage marker",
    "TPSB2": "tryptase beta 2 - canonical mast cell marker",
    "SLC34A2": "canonical alveolar type II pneumocyte-specific transporter",
}

LINEAGE_LABEL = {
    "SFTPC": "epithelial_alveolar", "SFTPB": "epithelial_alveolar", "SFTPA1": "epithelial_alveolar",
    "SFTPA2": "epithelial_alveolar", "NAPSA": "epithelial_alveolar", "SLC34A2": "epithelial_alveolar",
    "SCGB1A1": "epithelial_alveolar", "SCGB3A1": "epithelial_alveolar", "SCGB3A2": "epithelial_alveolar",
    "WFDC2": "epithelial_alveolar", "MUC1": "epithelial_alveolar", "AGER": "epithelial_alveolar",
    "CLDN18": "epithelial_alveolar", "EPCAM": "epithelial_alveolar", "SLPI": "epithelial_alveolar",
    "KRT7": "epithelial", "KRT8": "epithelial", "KRT17": "epithelial", "KRT18": "epithelial", "KRT19": "epithelial",
    "PIGR": "epithelial",
    "LYZ": "myeloid", "S100A8": "myeloid", "S100A9": "myeloid", "S100A12": "myeloid",
    "CD68": "myeloid", "MARCO": "myeloid", "FCN1": "myeloid", "VCAN": "myeloid",
    "MMP12": "myeloid", "TPSB2": "myeloid",
    "COL1A1": "stromal_vascular", "COL1A2": "stromal_vascular", "COL3A1": "stromal_vascular",
    "DCN": "stromal_vascular", "LUM": "stromal_vascular", "ACTA2": "stromal_vascular",
    "PECAM1": "stromal_vascular", "VWF": "stromal_vascular", "FBLN1": "stromal_vascular",
    "ACKR1": "stromal_vascular", "MGP": "stromal_vascular",
    "HBB": "erythroid", "HBA1": "erythroid", "HBA2": "erythroid", "HBD": "erythroid",
    "ALAS2": "erythroid", "AHSP": "erythroid",
}


def classify_genes(gene_names: pd.Series) -> pd.DataFrame:
    ambient_reason = dict(EXTENDED_AMBIENT_WITH_REASON)
    for g in SIBLING_REPO_VETTED_AMBIENT:
        ambient_reason.setdefault(g, "sibling-repo vetted ambient anchor (ambient_risk_diagnostic.py)")

    rows = []
    for gene in gene_names:
        is_ambient_anchor = gene in ambient_reason
        is_tcell_anchor = gene in SIBLING_REPO_VETTED_TCELL
        is_pattern_hit = bool(SIBLING_REPO_EXPLICIT_PATTERN.match(gene))
        if is_ambient_anchor or is_pattern_hit:
            lineage_flag = LINEAGE_LABEL.get(gene, "ambient_pattern_match")
            interpretation_class = "4_ambient_doublet_sensitive"
            basis = ambient_reason.get(
                gene, "matches sibling-repo explicit ambient/stress/ribosomal pattern"
            )
        elif is_tcell_anchor:
            lineage_flag = "t_cell_lymphoid"
            interpretation_class = "pending_identity_consistent_with_1_t_cell_intrinsic"
            basis = "sibling-repo vetted T-cell-intrinsic anchor - identity only, not yet confirmed by per-cell burden"
        else:
            lineage_flag = "unclassified"
            interpretation_class = "pending_requires_per_cell_contamination_scoring"
            basis = "not in the vetted anchor panel or explicit pattern - identity alone is insufficient"
        rows.append(
            {
                "Gene_name": gene,
                "lineage_flag": lineage_flag,
                "interpretation_class": interpretation_class,
                "classification_basis": basis,
            }
        )
    return pd.DataFrame(rows).drop_duplicates(subset="Gene_name").reset_index(drop=True)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    published = pd.read_csv(PUBLISHED_SUPPORT)
    flags = classify_genes(pd.Series(published["Gene_name"].unique()))
    merged = published.merge(flags, on="Gene_name", how="left")
    cols = [
        "comparison_label", "Gene_name", "Ensembl_ID", "n_donors_evaluable", "donor_supported",
        "sign_concordance", "ci_excludes_zero", "same_direction_pass_70",
        "lineage_flag", "interpretation_class", "classification_basis",
    ]
    merged = merged[cols].sort_values(["comparison_label", "Gene_name"]).reset_index(drop=True)
    out_path = OUT_DIR / "candidate_interpretation_classes.csv"
    merged.to_csv(out_path, index=False)

    n_genes = len(flags)
    n_ambient = int((flags["interpretation_class"] == "4_ambient_doublet_sensitive").sum())
    n_tcell_pending = int(flags["interpretation_class"].str.startswith("pending_identity_consistent").sum())
    n_pending = int((flags["interpretation_class"] == "pending_requires_per_cell_contamination_scoring").sum())
    print(f"wrote {out_path} ({len(merged)} rows, {n_genes} unique genes)")
    print(f"  ambient/doublet-sensitive by identity: {n_ambient}/{n_genes}")
    print(f"  identity-consistent with T-cell-intrinsic (unconfirmed): {n_tcell_pending}/{n_genes}")
    print(f"  unclassified, pending per-cell scoring: {n_pending}/{n_genes}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
