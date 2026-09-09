# Biological evaluation — Priority 1 (ambient RNA / lineage cross-expression)

> **SUPERSEDED 2026-09-09 — the per-cell half is now DONE; see "Priority 1 completed" below.**
> The outage record that follows is kept as written, not deleted: it is the dated reason the
> earlier pass stopped where it did, and a record that quietly changes what it once said is not
> a record. Read it as history, not as the current state.

**STATUS (2026-09-08, SUPERSEDED): the per-cell half of this task is BLOCKED, not deprioritised, not skipped.**
As of 2026-09-08 ~22:13Z the entire compute fleet (`petadimensionlab-super-server`, `dgx`,
`thinkstation1`, `thinkstation2`, `thinkstation3`) went offline within about a minute of each
other - a site-level uplink/power event, confirmed independently via `tailscale status` showing
every host transmitting with zero received traffic. This is not a sandbox permission boundary:
god has direct SSH access to these hosts and could not reach any of them either. There was no
data source to point at and no machine to run the per-cell diagnostic on, for anyone, at the
time this was worked. **Parking Priority 1's per-cell half here is the correct outcome given
that outage, not a workaround chosen instead of finishing it.** Do not read the gap below as
something a more thorough pass would have closed - re-check for compute-fleet availability
before assuming the blocker still holds.

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
   `~/workspace/KD/tcell_luad_lusc_normal_luscmax7000_heldout_allgene_perturbation/`. Wait for the
   compute fleet to come back online (see STATUS above), then either run the per-cell script on
   that host and relay compact results back (the pattern used for the T4 GPU phases) or copy the
   object locally - whichever is cheaper once the hosts are reachable again.
2. With that object: run `../ambient_risk_diagnostic.py` (already ported and committed here,
   see its own docstring - it is UNRUN, its default `NSCLC_H5AD` path is an unverified guess, and
   `main()` refuses to execute until that guard is removed). It carries the same
   epithelial/alveolar/myeloid/stromal/erythroid burden scoring and ambient-risk logistic model as
   the sibling repo's script, adapted to this repo's candidate-table schema and local gc104M
   symbol dictionaries; only its candidate-loading and symbol-map helpers are tested
   (`test_ambient_risk_diagnostic.py`, 7/7 pass) since the h5ad-dependent feature computation
   cannot be exercised without the source object. Review or compute doublet scores separately if
   the h5ad has none.

   > **Annotation 2026-09-09:** items 1 and 2 are DONE. The object was located and confirmed,
   > and the script has been run against it. The description above is stale in two specifics:
   > the guard is now an opt-in (`NSCLC_AMBIENT_ALLOW_RUN=1`) rather than something to remove,
   > and symbols come from the in-object `var["feature_name"]`, not the gc104M dictionaries.
   > **Item 3 below remains outstanding and is still correct.**

3. Re-run the perturbation ranking excluding high-burden/high-doublet cells (reusing
   `evaluate_donor_gene_effects.py`'s `accumulate()`/`summarize()` on the reduced cell set, not a
   different method) to produce `decontamination_rank_stability.csv` and check the provisional
   gate (retain ≥60% of the top 20 candidates with the same direction after exclusion).

No re-scoring, no retuning, no threshold changes, and no target claim in this partial pass — it
only reports lineage identity, explicitly distinct from a per-cell measurement.

---

# Priority 1 completed — 2026-09-09

**On this object, Priority 1 supports expression-based lineage cross-expression ONLY.** The
doublet and mitochondrial components are **retired**, and the perturbation ranking is **not
rerun** here.

## The two retirements — stated, not reported as zero

| component | state in the object | disposition |
| --- | --- | --- |
| `obs["doublet_status"]` | one distinct value, `singlet`, across all 21,000 cells | **RETIRED** — `doublet_component_status: constant_by_construction` |
| `obs["pct_counts_mito"]`, `obs["total_counts_mito"]` | identically zero (min = max = 0, one distinct value) | **RETIRED** — `mito_component_status: constant_zero_by_construction` |

The object is pre-filtered upstream, so a doublet flag computed here would return zero doublets:
true, and meaningless. **This is the specific failure this section exists to prevent.** Running
the originally-scoped method would not have errored — it would have emitted a clean-looking
*"no doublets, 0 % mito"* from columns that are constant by construction, which is the most
reassuring possible wrong answer. The manifest therefore **states both retirements explicitly**
and derives them from the object at runtime rather than hard-coding them.

**One measured detail corrects the assumption behind the mito zeros.** `var["mito"]` is present
but flags **zero genes** — the mitochondrial genes are absent from this object's gene set
altogether, not merely flagged and retained. The zero obs statistics are the residue of that
upstream filtering, not a measurement taken on retained mito genes. The manifest's wording is
generated from that count, so it cannot drift from the number printed beside it.

**`total_counts` and `n_genes_by_counts` are sequencing depth, not contamination.** They are
carried in the output under `depth_`-prefixed names so they cannot be read as a contamination
measure.

## What was measured

Per-cell lineage cross-expression: `log1p(CPM)` module scores over the ambient and
T-cell-intrinsic anchor panels, their difference, and their ratio.

- Input is **`layers["count"]`** — raw integral counts — verified integral, finite and
  non-negative at load. `X` on this object is normalised/log; feeding it to a counts method
  returns numbers, not an error.
- Symbols come from the **in-object `var["feature_name"]`** (17,764 unique, none empty). The var
  index is Ensembl, which is why symbol lookups against it miss. No external gc104M dictionary
  is used.
- **All 39 ambient and all 39 T-cell anchors resolved**; none missing.

`tables/cell_contamination_scores.csv` — 21,000 rows, 176 donors, splits train 14,282 /
test 3,379 / eval 3,339. The test split's 3,379 cells match the held-out cell manifest exactly.

| | ambient module | T-cell module | ambient − T-cell |
| --- | --- | --- | --- |
| median | 0.170 | 2.013 | −1.760 |
| 95th pct | — | — | −0.518 |
| 99th pct | — | — | 0.038 |

**241 of 21,000 cells (1.15 %) carry more ambient-lineage than T-cell signal.** Contamination on
this object is a small-tail phenomenon, not a pervasive one. `ambient_minus_tcell` is the stable
statistic; `ambient_to_tcell_ratio` is reported alongside it but is floored at a T-cell score of
1e-3 (22 cells sit below that floor), so its extreme values are an artifact of the floor rather
than a measurement.

## Retained gene-level candidate merge

`tables/candidate_gene_ambient_risk.csv` — the ported gene-level ambient-risk model, scoring
model and threshold **unchanged**. Two inputs changed, both corrections rather than tuning: the
matrix is now `layers["count"]` (it was `X`), and symbols come from the in-object feature map.
Cross-validated AUC separating the anchor sets **0.975 ± 0.031**; 11,479 genes scored; flag
threshold 0.868 (25th percentile of ambient anchors, as ported). Of the 74 published-leading
genes: **36 `AMBIENT_RISK`, 36 `ok`, 2 `not_scored`** (`ACKR1`, `S100A7` fall below the 0.5 %
detection filter).

### The two flag sets are NOT the same set, despite nearly equal totals

The earlier identity-only pass flagged **37** of 74; this expression-based pass flags **36**.
**That near-match is a coincidence — the sets share only 22 genes.**

- **22 flagged by both.**
- **15 flagged by identity only**: `ACKR1`, `EEF1G`, `FBLN1`, `HSPA1B`, `KRT19`, `LYZ`, `RPL21`,
  `RPL27A`, `RPS26`, `RPS27`, `S100A2`, `S100A4`, `S100A7`, `S100A8`, **`SFTPC`** — of which two
  (`ACKR1`, `S100A7`) could not be scored at all.
- **14 flagged by expression only**: `AGR2`, `CCN1`, `FN1`, `H1-10`, `IGFBP6`, `IL1R2`, `MDK`,
  `MNDA`, `NME2`, `PABPC3`, `PLA2G2A`, `PRICKLE4`, `RAD50`, `TRAPPC5`.

`SFTPC` — the flagship concern in the biology README — is **flagged by lineage identity but not
by the expression-based ambient-risk model** on this object. That disagreement is reported, not
reconciled: it is a question for review, not something to average away. **The two passes are not
independent evidence** either — they share the same anchor panels — so agreement between them
must not be read as replication.

## Verification

- The run guard refuses without `NSCLC_AMBIENT_ALLOW_RUN=1`, confirmed on the real host.
- **The loader was pointed at the real sibling object and rejected it.** The tokenization input
  `.../h5ad_for_tokenization/balanced_lusc_max_7000_per_disease_tcells_geneformer.h5ad` loads as
  the *identical* 21,000 × 17,764 — the run log prints that shape before failing — so a shape or
  cell-count assertion would have passed on it. Rejection is on `obs`/`layers` keys.
- 39 unit tests pass under `python3 -m unittest`. The critical ones were mutation-tested:
  removing the obs-key check, hard-coding the retirement string, and adding an `X` fallback each
  turn the suite red.

## Not done here, deliberately

`decontamination_rank_stability.csv` and the ≥60 %-of-top-20 provisional gate require **rerunning
the perturbation ranking on a reduced cell set**, which is out of scope for this task and is not
attempted. The per-cell scores this pass produces are the input that work would need.

**No re-scoring, no retuning, no threshold changes, and no target claim in this pass.**
