# Biological evaluation

## Objective

Determine whether statistically robust perturbations represent intrinsic
T-cell programs rather than ambient lung RNA, doublets, donor composition, or
generic cellular stress.

## Priority 1: ambient RNA and doublet sensitivity

The current top rankings include alveolar and epithelial genes such as
`SFTPC`, `SFTPB`, `NAPSA`, `MUC1`, `PIGR`, and keratins, plus stromal or
vascular-associated genes such as `FBLN1`, `DCN`, and `ACKR1`. Their presence
in a T-cell-selected analysis is a major alternative explanation that must be
tested directly.

**Scope on the confirmed object, settled 2026-09-09.** Priority 1 here supports
**expression-based lineage cross-expression only**. Two of the originally-listed
components are **retired on this object**, and the ranking is **not rerun** by the
diagnostic:

- ✅ **Per-cell lineage marker burden** — done. `../ambient_risk_diagnostic.py` computes
  `log1p(CPM)` module scores over the ambient and T-cell anchor panels from
  `layers["count"]`, and writes `../results/biology/tables/cell_contamination_scores.csv`.
- ❌ **Doublet scores — RETIRED, not skipped.** `obs["doublet_status"]` is `singlet` for all
  21,000 cells: the object is pre-filtered upstream. A doublet review here would report zero
  doublets, which is true and meaningless. **The manifest states the retirement rather than
  reporting a zero**, because a clean-looking "no doublets" is the most reassuring possible
  wrong answer.
- ❌ **Mitochondrial QC — RETIRED.** `pct_counts_mito` and `total_counts_mito` are identically
  zero, and `var["mito"]` flags zero genes: the mitochondrial genes are absent from the gene
  set entirely.
- ⬜ **Ambient-RNA correction**, and **repeating the perturbation ranking** after excluding
  high-burden cells, are **not attempted** by this diagnostic. They need a separate,
  explicitly-scoped pass; the per-cell scores above are their input.
- ⬜ Whether gene effects correlate with contamination burden is likewise still open.

**Sequencing depth is not contamination.** `total_counts` and `n_genes_by_counts` vary widely
across cells; they are carried as `depth_`-prefixed metadata and must not be substituted for a
contamination measure.

See `../results/biology/REPORT.md` for results and the verification record.

Suggested provisional gate: retain at least 60% of the top 20 candidates with
the same direction after the prespecified decontamination and exclusion
analysis, while requiring individual promoted candidates to remain stable.

## Priority 2: T-cell subtype specificity

- Stratify effects across CD4, CD8, regulatory, memory-like, exhausted,
  proliferating, and other supported T-cell states.
- Rebuild centroids within subtype where sample size permits.
- Distinguish a shared T-cell disease program from changes caused by different
  subtype proportions.
- Require enough donors and cells within each subtype before interpreting a
  missing effect as biological absence.

## Priority 3: program coherence

- Compare each candidate with source-versus-goal differential expression and
  established T-cell activation, exhaustion, cytotoxicity, interferon, stress,
  and metabolic programs.
- Ask whether deletion changes a coherent multi-gene T-cell program rather
  than only moving a global embedding coordinate.
- Examine leading-edge genes rather than interpreting thousands of overlapping
  enrichment terms as independent discoveries.
- Consolidate redundant GO, Reactome, and KEGG terms into modules, then test
  module scores by donor.

## Priority 4: external and orthogonal replication

- Repeat the shortlist in a donor-balanced external NSCLC T-cell cohort.
- Test a second Geneformer checkpoint or another single-cell representation.
- Compare with published or available CRISPR, Perturb-seq, or loss-of-function
  evidence in human T cells.
- For a final experimental shortlist, use targeted knockdown or CRISPR and
  measure T-cell activation, exhaustion, cytokine secretion, proliferation,
  viability, and tumor-cell killing where appropriate.

## Candidate interpretation classes

Assign every reviewed candidate to one of these classes:

1. **T-cell-intrinsic and donor-stable** — suitable for replication.
2. **T-cell-subtype-specific** — potentially useful with a restricted claim.
3. **Tissue-environment-associated** — informative about the sample context,
   but not an intrinsic T-cell target.
4. **Ambient/doublet-sensitive** — remove from target prioritization.
5. **Generic stress or viability-related** — retain only with explicit
   functional guardrails.
6. **Unstable or underpowered** — defer pending more donors or detections.

## Minimum biological deliverables

- `cell_contamination_scores.csv` — **produced 2026-09-09**
- `decontamination_rank_stability.csv` — still outstanding; needs the ranking rerun, which
  Priority 1 as scoped does not perform
- `tcell_subtype_effects.csv`
- `candidate_interpretation_classes.csv`
- `pathway_leading_edge_modules.csv`
- `external_replication_summary.csv`
- figures showing contamination sensitivity, subtype effects, pathway modules,
  and external replication

Write future outputs under `../results/biology/`, with links to the source
cohort and exact cell-filtering definitions.
