# Donor-level gene-effect evaluation

Generated: 2026-09-08T11:53:10Z
Script: `evaluate_donor_gene_effects.py`
Source: `/home/kaisar/workspace/KD/tcell_luad_lusc_normal_luscmax7000_heldout_allgene_perturbation`

## Prespecified re-score criterion (written before the new numbers)

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

## Results (a-priori published genes lead)

- Published leading genes, evaluable under the floor, all comparisons ≥70% same-direction: **PASS**
- Donor-supported top-50, all comparisons ≥70% same-direction: **PASS**
- Unfiltered top-50 (retained, not a valid test): **PASS**
- LODO on donor-supported ranking: **FAIL**
- LODO on unfiltered ranking (retained; LUSC n=5 overlap is a cohort-size artifact, not a comparison finding): **FAIL**

### Published leading genes (donor support)

comparison_label Gene_name  n_donors_evaluable  min_donors_required  donor_supported  sign_concordance  ci_excludes_zero  same_direction_pass_70  n_cells_total
     LUAD → LUSC       HBB                  18                   10             True          0.888889              True                    True            178
     LUAD → LUSC     SFTPB                  16                   10             True          1.000000              True                    True            126
     LUAD → LUSC     SFTPC                  18                   10             True          0.944444              True                    True            314
     LUAD → LUSC     NAPSA                  15                   10             True          1.000000              True                    True             99
     LUAD → LUSC   TRAPPC5                  12                   10             True          0.750000             False                    True             28
     LUAD → LUSC    HSPA1B                  18                   10             True          0.833333              True                    True            270
     LUAD → LUSC     RAB4B                  18                   10             True          0.944444              True                    True             84
     LUAD → LUSC   SCGB3A1                  17                   10             True          1.000000              True                    True            176
     LUAD → LUSC      NME2                  16                   10             True          0.875000              True                    True             65
     LUAD → LUSC     RAD50                  17                   10             True          0.882353              True                    True             98
     LUAD → LUSC   SELENOP                  15                   10             True          0.600000             False                   False             34
     LUAD → LUSC    S100A4                  19                   10             True          0.947368              True                    True           1135
     LUAD → LUSC    GABPB1                  17                   10             True          1.000000              True                    True            248
     LUAD → LUSC       FN1                  17                   10             True          0.764706              True                    True             50
     LUAD → LUSC    IFNAR2                  18                   10             True          0.888889              True                    True             95
   LUAD → NORMAL      PIGR                  12                   10             True          0.916667              True                    True             25
   LUAD → NORMAL     EEF1G                  18                   10             True          0.611111             False                   False            196
   LUAD → NORMAL     KRT18                  16                   10             True          0.937500              True                    True             58
   LUAD → NORMAL      MUC1                  17                   10             True          0.882353              True                    True             52
   LUAD → NORMAL      AGR2                  15                   10             True          0.866667              True                    True             49
   LUAD → NORMAL     WFDC2                  15                   10             True          1.000000              True                    True             72
   LUAD → NORMAL     KRT19                  16                   10             True          0.937500              True                    True             54
   LUAD → NORMAL   SLC34A2                  14                   10             True          0.928571              True                    True             34
   LUAD → NORMAL   SCGB3A2                  17                   10             True          0.941176              True                    True            100
   LUAD → NORMAL      KRT8                  13                   10             True          1.000000              True                    True             45
   LUAD → NORMAL      KRT7                  13                   10             True          1.000000              True                    True             32
   LUAD → NORMAL      BTG1                  19                   10             True          1.000000              True                    True           1330
   LUAD → NORMAL   POLR2J3                  17                   10             True          1.000000              True                    True            335
   LUAD → NORMAL      HBA1                  16                   10             True          0.937500              True                    True             84
   LUAD → NORMAL       MDK                  12                   10             True          0.916667              True                    True             31
     LUSC → LUAD     MMP12                   2                    3            False          1.000000              True                    True             28
     LUSC → LUAD     KRT17                   3                    3             True          1.000000              True                    True             41
     LUSC → LUAD    S100A2                   3                    3             True          1.000000              True                    True             35
     LUSC → LUAD    S100A7                   3                    3             True          1.000000              True                    True             35
     LUSC → LUAD      HBA2                   3                    3             True          1.000000              True                    True            135
     LUSC → LUAD    S100A8                   3                    3             True          1.000000              True                    True            134
     LUSC → LUAD      HBA1                   3                    3             True          1.000000              True                    True            105
     LUSC → LUAD     CXCL8                   3                    3             True          1.000000              True                    True             28
     LUSC → LUAD    S100A9                   4                    3             True          1.000000              True                    True            134
     LUSC → LUAD     TPSB2                   4                    3             True          0.750000              True                    True             36
     LUSC → LUAD    GPR174                   3                    3             True          1.000000              True                    True             49
     LUSC → LUAD       HBB                   4                    3             True          1.000000              True                    True            116
     LUSC → LUAD     IL1R2                   3                    3             True          0.666667             False                   False             27
     LUSC → LUAD      GZMK                   3                    3             True          1.000000              True                    True             46
     LUSC → LUAD      MNDA                   3                    3             True          1.000000              True                    True             26
   LUSC → NORMAL     RPS27                   5                    3             True          0.800000              True                    True            560
   LUSC → NORMAL     RPL21                   5                    3             True          1.000000              True                    True            560
   LUSC → NORMAL       LYZ                   3                    3             True          1.000000              True                    True             92
   LUSC → NORMAL    UQCR11                   3                    3             True          1.000000              True                    True            156
   LUSC → NORMAL      BTG1                   5                    3             True          1.000000              True                    True            510
   LUSC → NORMAL    RPL27A                   5                    3             True          0.800000              True                    True            536
   LUSC → NORMAL    NDUFB8                   4                    3             True          0.500000             False                   False             66
   LUSC → NORMAL     PDCD7                   4                    3             True          0.500000             False                   False             26
   LUSC → NORMAL   ARL6IP1                   3                    3             True          1.000000              True                    True             74
   LUSC → NORMAL      PTMA                   5                    3             True          0.800000              True                    True            552
   LUSC → NORMAL     H1-10                   3                    3             True          1.000000              True                    True             63
   LUSC → NORMAL    FAM89B                   3                    3             True          0.666667             False                   False             30
   LUSC → NORMAL     H3-3B                   5                    3             True          1.000000              True                    True            461
   LUSC → NORMAL      CBX3                   5                    3             True          0.800000             False                    True            348
   LUSC → NORMAL     PNRC1                   4                    3             True          1.000000              True                    True            443
   NORMAL → LUAD     SFTPC                  12                    6             True          1.000000              True                    True            465
   NORMAL → LUAD     FBLN1                   6                    6             True          1.000000              True                    True             32
   NORMAL → LUAD       DCN                   9                    6             True          1.000000              True                    True             44
   NORMAL → LUAD     RPS26                  12                    6             True          1.000000              True                    True            960
   NORMAL → LUAD     TXNIP                  12                    6             True          1.000000              True                    True           1053
   NORMAL → LUAD       MGP                  10                    6             True          1.000000              True                    True             81
   NORMAL → LUAD    IGFBP6                   8                    6             True          1.000000              True                    True             28
   NORMAL → LUAD    SFTPA1                  12                    6             True          1.000000              True                    True            134
   NORMAL → LUAD     ACKR1                   8                    6             True          0.875000              True                    True             36
   NORMAL → LUAD    HSPA1B                  12                    6             True          1.000000              True                    True            641
   NORMAL → LUAD    SFTPA2                   9                    6             True          1.000000              True                    True            131
   NORMAL → LUAD      HBA2                   8                    6             True          0.750000             False                    True             34
   NORMAL → LUAD      CCN1                   8                    6             True          0.875000              True                    True             31
   NORMAL → LUAD   PLA2G2A                   7                    6             True          1.000000              True                    True             32
   NORMAL → LUAD     PSMA6                  12                    6             True          1.000000              True                    True            156
   NORMAL → LUSC     SFTPC                  12                    6             True          1.000000              True                    True            465
   NORMAL → LUSC     FBLN1                   6                    6             True          1.000000              True                    True             32
   NORMAL → LUSC    SFTPA1                  12                    6             True          0.916667              True                    True            134
   NORMAL → LUSC     MATR3                   7                    6             True          1.000000              True                    True             43
   NORMAL → LUSC    HSPA1B                  12                    6             True          1.000000              True                    True            641
   NORMAL → LUSC  PRICKLE4                   7                    6             True          1.000000              True                    True             28
   NORMAL → LUSC    SFTPA2                   9                    6             True          1.000000              True                    True            131
   NORMAL → LUSC       MGP                  10                    6             True          1.000000              True                    True             81
   NORMAL → LUSC     RAB4B                  10                    6             True          0.900000              True                    True             63
   NORMAL → LUSC   STEAP1B                  12                    6             True          1.000000              True                    True            433
   NORMAL → LUSC     RAD50                  11                    6             True          1.000000              True                    True             87
   NORMAL → LUSC     SFTPB                  12                    6             True          0.666667             False                   False            138
   NORMAL → LUSC    PABPC3                  10                    6             True          0.800000             False                    True             62
   NORMAL → LUSC    IFNAR2                  12                    6             True          0.916667              True                    True            141
   NORMAL → LUSC    DNAJB1                  12                    6             True          1.000000              True                    True            843

### Gate table (denominators beside fractions)

scored_set `n_genes_in_set` is the denominator. `n_donors_min/median/max` describe donor support inside that set. `fail_margin_genes` (only defined when `gate_pass` is True) is how many currently-passing genes would have to flip direction before this row's PASS became FAIL - a low number means the PASS is fragile, not a false statement, but one that needs the margin stated next to it.

comparison_label source                  scored_set  min_donors_required  n_genes_in_set  n_pass_70  frac_pass_70  n_donors_min  n_donors_median  n_donors_max  n_ci_excludes_zero  gate_pass  fail_margin_genes  n_published_total  n_published_below_floor
     LUAD → LUSC   luad            unfiltered_top50                   10              50         47      0.940000             1              1.0            18                  46       True               13.0                NaN                      NaN
     LUAD → LUSC   luad       donor_supported_top50                   10              50         41      0.820000            10             17.5            19                  41       True                7.0                NaN                      NaN
     LUAD → LUSC   luad published_leading_evaluable                   10              15         14      0.933333            12             17.0            19                  13       True                4.0               15.0                      0.0
   LUAD → NORMAL   luad            unfiltered_top50                   10              50         47      0.940000             1              1.0             7                  47       True               13.0                NaN                      NaN
   LUAD → NORMAL   luad       donor_supported_top50                   10              50         48      0.960000            10             15.0            19                  47       True               14.0                NaN                      NaN
   LUAD → NORMAL   luad published_leading_evaluable                   10              15         14      0.933333            12             16.0            19                  14       True                4.0               15.0                      0.0
     LUSC → LUAD   lusc            unfiltered_top50                    3              50         49      0.980000             1              1.0             4                  49       True               15.0                NaN                      NaN
     LUSC → LUAD   lusc       donor_supported_top50                    3              50         35      0.700000             3              3.0             5                  34       True                1.0                NaN                      NaN
     LUSC → LUAD   lusc published_leading_evaluable                    3              14         13      0.928571             3              3.0             4                  13       True                4.0               15.0                      1.0
   LUSC → NORMAL   lusc            unfiltered_top50                    3              50         48      0.960000             1              1.0             4                  48       True               14.0                NaN                      NaN
   LUSC → NORMAL   lusc       donor_supported_top50                    3              50         36      0.720000             3              3.0             5                  35       True                2.0                NaN                      NaN
   LUSC → NORMAL   lusc published_leading_evaluable                    3              15         12      0.800000             3              4.0             5                  11       True                2.0               15.0                      0.0
   NORMAL → LUAD normal            unfiltered_top50                    6              50         50      1.000000             1              2.0            12                  50       True               16.0                NaN                      NaN
   NORMAL → LUAD normal       donor_supported_top50                    6              50         48      0.960000             6              7.5            12                  46       True               14.0                NaN                      NaN
   NORMAL → LUAD normal published_leading_evaluable                    6              15         15      1.000000             6              9.0            12                  14       True                5.0               15.0                      0.0
   NORMAL → LUSC normal            unfiltered_top50                    6              50         44      0.880000             1              2.0            12                  44       True               10.0                NaN                      NaN
   NORMAL → LUSC normal       donor_supported_top50                    6              50         46      0.920000             6              7.0            12                  41       True               12.0                NaN                      NaN
   NORMAL → LUSC normal published_leading_evaluable                    6              15         14      0.933333             6             11.0            12                  13       True                4.0               15.0                      0.0

**PASS alone is a true statement that misleads for the rows below - each passes by a one- or two-gene margin:**

- `LUSC → LUAD` / `donor_supported_top50`: **35/50 = 0.70, AT THE THRESHOLD** - a **1-gene margin** (losing 1 more would drop it below 70%).
- `LUSC → NORMAL` / `donor_supported_top50`: **36/50 = 0.72** - a **2-gene margin** (losing 2 more would drop it below 70%).
- `LUSC → NORMAL` / `published_leading_evaluable`: **12/15 = 0.80** - a **2-gene margin** (losing 2 more would drop it below 70%).

### Leave-one-donor-out

        ranking comparison_label  n_folds  n_ranked  min_donors_required  top50_overlap fold_overlap_ge_70 gene_retained_in_70pct_folds
     unfiltered      LUSC → LUAD        5        50                    2       0.816000                0.8                         True
     unfiltered    LUSC → NORMAL        5        50                    2       0.784000                0.6                        False
     unfiltered      LUAD → LUSC       19        50                    2       0.844211                1.0                         True
     unfiltered    LUAD → NORMAL       19        50                    2       0.858947                1.0                         True
     unfiltered    NORMAL → LUAD       12        50                    2       0.893333                1.0                         True
     unfiltered    NORMAL → LUSC       12        50                    2       0.868333                1.0                         True
donor_supported      LUSC → LUAD        5        50                    3       0.660000                0.4                        False
donor_supported    LUSC → NORMAL        5        50                    3       0.712000                0.4                        False
donor_supported      LUAD → LUSC       19        50                   10       0.911579                1.0                         True
donor_supported    LUAD → NORMAL       19        50                   10       0.894737                1.0                         True
donor_supported    NORMAL → LUAD       12        50                    6       0.885000                1.0                         True
donor_supported    NORMAL → LUSC       12        50                    6       0.905000                1.0                         True

### OBSERVATION: donor-supported vs. unfiltered LODO direction split

Descriptive only - not a claim, and no threshold or ranking was changed after seeing this split:

- `LUAD → LUSC`: unfiltered 0.8442 -> donor-supported 0.9116 (improves, delta +0.0674)
- `NORMAL → LUSC`: unfiltered 0.8683 -> donor-supported 0.9050 (improves, delta +0.0367)
- `LUAD → NORMAL`: unfiltered 0.8589 -> donor-supported 0.8947 (improves, delta +0.0358)
- `NORMAL → LUAD`: unfiltered 0.8933 -> donor-supported 0.8850 (degrades, delta -0.0083)
- `LUSC → NORMAL`: unfiltered 0.7840 -> donor-supported 0.7120 (degrades, delta -0.0720)
- `LUSC → LUAD`: unfiltered 0.8160 -> donor-supported 0.6600 (degrades, delta -0.1560)

## Donor counts

- LUAD test donors: 19 (floor 10)
- Normal test donors: 12 (floor 6)
- LUSC test donors: 5 (floor 3)

## Limitations

- This is a statistical donor-stability summary, not a biological target claim.
- Unfiltered LUSC→NORMAL LODO fail tracks five-donor cohort size plus n=1 ranking, not a comparison-specific biological finding.
- Ambient RNA, doublets, and T-cell subtype structure are not addressed here.
- Centroids were not rebuilt; test-donor LODO does not require that step.
- No threshold, donor subset, or aggregation was changed after seeing results.

