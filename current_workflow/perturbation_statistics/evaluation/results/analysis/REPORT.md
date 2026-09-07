# Donor-level gene-effect evaluation

Generated: 2026-09-07T00:08:46Z
Script: `evaluate_donor_gene_effects.py`
Source: `/home/kaisar/workspace/KD/tcell_luad_lusc_normal_luscmax7000_heldout_allgene_perturbation`

## Method

- Unit of analysis is the donor. Cell-level cosine similarities to the three
  training-donor state centroids were read from existing raw pickles.
- For each cell-gene deletion, the directional shift is `cos(goal) - cos(source)`.
- Donor effect is the mean of that cell-level shift within the donor.
- Gene summary uses the median of donor effects (equal donor weight).
- Sign concordance is the fraction of donors with a nonzero mean whose sign
  matches the median donor effect. The 70% gate uses this fraction.
- Bootstrap 95% CI is the percentile interval of the mean of donor means
  (1000 resamples, seed 0).
- Leave-one-donor-out ranks genes by median donor effect after dropping one
  test donor. Training centroids are unchanged (test donors were not in them).
- LUSC has five held-out donors; its stability numbers are reported separately
  and are not compared mechanically with LUAD or normal.

## Direction gate (same sign in ≥70% of evaluable donors)

- Donor-ranked top-50, all comparisons meet gate: **PASS**
- Published leading genes, all comparisons meet gate: **PASS**

comparison_label  n_genes  n_pass_70_all  frac_pass_70_all  n_top50_donor_ranked  n_top50_pass_70  frac_top50_pass_70  n_top50_ci_excludes_zero  n_published_top  n_published_top_pass_70  frac_published_top_pass_70
     LUAD → LUSC    13458           4675          0.347377                    50               47                0.94                        46               15                       14                    0.933333
   LUAD → NORMAL    13458           4896          0.363798                    50               47                0.94                        47               15                       14                    0.933333
     LUSC → LUAD    11242           7422          0.660203                    50               49                0.98                        49               15                       14                    0.933333
   LUSC → NORMAL    11242           7163          0.637164                    50               48                0.96                        48               15                       12                    0.800000
   NORMAL → LUAD    14923           6096          0.408497                    50               50                1.00                        50               15                       15                    1.000000
   NORMAL → LUSC    14923           6036          0.404476                    50               44                0.88                        44               15                       14                    0.933333

## Leave-one-donor-out top-50 retention

- Gate (mean fold overlap ≥70% of genes in ≥70% of folds, encoded per comparison): **FAIL**

comparison_label  n_folds  top50_overlap fold_overlap_ge_70 gene_retained_in_70pct_folds
     LUSC → LUAD        5       0.816000                0.8                         True
   LUSC → NORMAL        5       0.784000                0.6                        False
     LUAD → LUSC       19       0.844211                1.0                         True
   LUAD → NORMAL       19       0.858947                1.0                         True
   NORMAL → LUAD       12       0.893333                1.0                         True
   NORMAL → LUSC       12       0.868333                1.0                         True

## Donor counts

- LUAD test donors: 19
- Normal test donors: 12
- LUSC test donors: 5

## Limitations

- This is a statistical donor-stability summary, not a biological target claim.
- Ambient RNA, doublets, and T-cell subtype structure are not addressed here.
- Centroids were not rebuilt; test-donor LODO does not require that step.
- No threshold, donor subset, or aggregation was changed after seeing results.

