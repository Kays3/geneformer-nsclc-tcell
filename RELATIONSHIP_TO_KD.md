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
separate [`geneformer-sclc-tcell`](https://github.com/Kays3/geneformer-sclc-tcell) repo.)

Keep the KD directory layout stable unless all workflow scripts and path
configuration are updated together.
