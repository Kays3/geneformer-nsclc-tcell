# Relationship To KD

`KD` is the local data and generated-artifact workspace for this repository's
NSCLC Geneformer workflow. This repository contains version-controlled
methods, monitoring, validation, reporting, and migration tools; it reads
selected KD datasets, model checkpoints, statistics, and perturbation outputs.

The large files are intentionally outside Git in:

```text
/home/thinkstation2/workspace/KD/tcell_luad_lusc_normal_10k_from_atlas/
/home/thinkstation2/workspace/KD/tcell_luad_lusc_normal_luscmax7000_finetune/
/home/thinkstation2/workspace/KD/tcell_luad_lusc_normal_luscmax7000_heldout_allgene_perturbation/
/home/thinkstation2/workspace/KD/data/nsclc/nsclc_integrated.h5ad
```

(The SCLC line of work's `KD/sclc_luad_normal_htan_*` subtrees belong to the
separate [`geneformer-sclc-tcell`](https://github.com/Kays3/geneformer-sclc-tcell) repo.
The integrative epithelial/tumor-microenvironment line of work, spanning both
cancer types, belongs to
[`geneformer-epithelial-tme`](https://github.com/Kays3/geneformer-epithelial-tme)
-- it reads this atlas's malignant/epithelial cells directly, from its actual
location at
`/home/kaisar/workspace/geneformer-uv-starter/geneformer-workspace/analysis/data/nsclc/nsclc_integrated.h5ad`,
**not** the `KD/data/nsclc/nsclc_integrated.h5ad` path listed above -- that
path does not currently resolve on `thinkstation2`; this doc is stale and
should be corrected separately.)

Keep the KD directory layout stable unless all workflow scripts and path
configuration are updated together.
