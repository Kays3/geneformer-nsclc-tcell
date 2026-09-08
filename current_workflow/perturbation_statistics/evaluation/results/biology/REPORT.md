# Biological evaluation — Priority 1 (ambient RNA / doublet sensitivity), partial

**This is a partial first pass, not the Priority 1 deliverable the biology README specifies.**
It flags candidate genes by lineage **identity** against a small vetted anchor panel. It does
**not** compute per-cell contamination burden, does not compute or review doublet scores, and
does not repeat the perturbation ranking after excluding high-burden cells — those need the
original single-cell object (obs-level celltype/individual/library-size structure) for this
NSCLC/LUSC/LUAD held-out cohort, the way `ambient_risk_diagnostic.py` in the sibling SCLC repo
(`geneformer-lung-tcell`) uses `htan_sclc_luad_normal_tcells_prepared.h5ad`. No equivalent object
is present in this checkout or reachable from this sandbox (no `.h5ad`, no raw counts, no
per-cell manifest beyond donor/disease IDs).

## What this does

`classify_candidate_lineage_flags.py` reads the already-produced
`published_gene_donor_support.csv` (the 15 a-priori published leading genes per comparison, 90
rows / 74 unique genes) and flags each gene using:

- the exact ambient and T-cell-intrinsic anchor panels from
  `sclc_validation/primary_test_perturbation/scripts/ambient_risk_diagnostic.py` in the sibling
  repo (calibrated there against real single-cell data, cross-validated AUC reported in that
  script's output);
- the exact explicit ambient/stress/ribosomal regex from
  `sclc_validation/immune_axis_test/axis_consistency_sensitivity.py` in the same sibling repo;
- a small extension of nine genes this repo's own `biology/README.md` names by name as
  concerning, or that are unambiguous single-lineage markers absent from the sibling panel
  (`PIGR`, `ACKR1`, `FBLN1`, `KRT7`, `KRT17`, `MGP`, `MMP12`, `TPSB2`, `SLC34A2` — each a
  standard, uncontroversial single-lineage marker, not a judgment call about this dataset).

Output: `tables/candidate_interpretation_classes.csv`.

## Result

- **37 of 74 unique published-leading genes (50%) flag as lineage-foreign by identity alone**
  (`interpretation_class = 4_ambient_doublet_sensitive`) — epithelial/alveolar, myeloid,
  stromal/vascular, or erythroid markers, or matching the explicit stress/ribosomal pattern.
  This directly confirms the biology README's stated concern (`SFTPC`, `SFTPB`, `NAPSA`, `MUC1`,
  `PIGR`, keratins, `FBLN1`, `DCN`, `ACKR1` are all in this set) rather than leaving it as an
  unquantified worry.
- **1 gene (`GZMK`) matches the vetted T-cell-intrinsic anchor panel.** This is identity
  consistent with class 1, **not a confirmation** — it still needs per-cell/per-donor stability
  evidence before being called donor-stable and T-cell-intrinsic.
- **36 genes are not in either vetted panel or pattern** and are left explicitly
  `pending_requires_per_cell_contamination_scoring` rather than guessed at. Roughly half of the
  published-leading set cannot be classified from gene identity alone — that is the actual size
  of the gap the real Priority 1 diagnostic needs to close.

## What is still needed to finish Priority 1

1. The source single-cell object for this cohort (or a copy of its `obs` table: celltype,
   individual, library size, and a doublet score column if one already exists upstream), most
   likely on `thinkstation1`/`thinkstation2` alongside
   `~/workspace/KD/tcell_luad_lusc_normal_luscmax7000_heldout_allgene_perturbation/`. Sandbox
   policy blocks direct SSH from this session (same boundary hit repeatedly this week), so this
   needs either a compute-node run relayed back (the pattern used for the T4 GPU phases) or a
   pointed local copy.
2. With that object: port `ambient_risk_diagnostic.py`'s epithelial/alveolar/myeloid/
   stromal/erythroid burden scoring and its ambient-risk logistic model to this cohort, review or
   compute doublet scores, then produce `cell_contamination_scores.csv`.
3. Re-run the perturbation ranking excluding high-burden/high-doublet cells (reusing
   `evaluate_donor_gene_effects.py`'s `accumulate()`/`summarize()` on the reduced cell set, not a
   different method) to produce `decontamination_rank_stability.csv` and check the provisional
   gate (retain ≥60% of the top 20 candidates with the same direction after exclusion).

No re-scoring, no retuning, no threshold changes, and no target claim in this partial pass — it
only reports lineage identity, explicitly distinct from a per-cell measurement.
