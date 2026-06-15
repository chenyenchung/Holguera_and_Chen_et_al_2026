# Source Directory Guide

This README documents top-level scripts in `src/` that are not part of the
module subdirectories.

## Module directories

The main pipeline modules are documented in their own READMEs:

- `src/preprocessing/`
- `src/stats/`
- `src/visualize/`
- `src/selector_test/`
- `src/cam_test/`

## Top-level manual scripts

These scripts are intended for manual or figure-specific analyses outside the
standard Nextflow modules.

### `src/manual_visualize.r`

- Purpose: manual scatter + depth-density visualization for one neuropil/synapse condition
- Key inputs: rotated matrix (`--synf`), annotation (`--ann`), metadata (`--meta`), preset table (`--preset`)
- Key parameters: `--np`, `--syn_type`, `--use_preset`, `--density`, `--subsample`, `--sparse_limit`
- Outputs: `<prefix>.pdf` and `<prefix>_legend.pdf` in the current working directory

Example:

```bash
Rscript src/manual_visualize.r \
  --np ME_R \
  --synf int/idv_mat/ME_R_rotated.csv.gz \
  --syn_type pre \
  --use_preset type_putative_1 \
  --density asis \
  --ann data/visual_neurons_anno.csv \
  --meta data/viz_meta.csv \
  --preset data/viz_preset.csv \
  --subsample 10000 \
  --sparse_limit 100
```

### `src/manual_per_window_enrichment.r`

- Purpose: Fisher exact-test enrichment of functional subsystem per temporal window
- Key input: `data/visual_neurons_anno.csv`
- Output: `int/Sup_tbl_3_per_window_func_enrichment.csv`

### `src/manual_lc_depth.r`

- Purpose: generate LC/LPLC partner-depth visualization figures for selected groups
- Key inputs: connections, cell type labels, rotated matrices, viz metadata, LC/LPLC origin map
- Configuration: uses hardcoded paths/parameters in the script (edit `argvs` block to customize)
- Outputs: PDFs and legend PDFs under `int/lc_lplc_dev_origin/`

### `src/manual_ds_stats.r`

- Purpose: assemble supplementary deep/superficial bias figure panels
- Key input: `int/stats/deep_superficial/combined_broad_depth_results.xlsx`
- Outputs: `int/Supp_fig_Y1.pdf`, `int/Supp_fig_Y2.pdf`, `int/Supp_fig_Y3.pdf`, `int/Supp_fig_Y4.pdf`

## Top-level helper scripts

### `src/get_selectors.R`

- Purpose: rebuild a binary selector matrix from mixture-model objects and TF list
- Note: references an external absolute path for input objects; update before use
- Output: currently writes `data/selector.csv`

### `src/convert_gene_expression_to_readable_format.R`

- Purpose: convert binary gene-expression matrices to "On/Off" human-readable tables with renamed columns
- Inputs: `data/P15_tf.csv`, `data/P15_CAM.csv`, `data/selectors.csv`, `data/visual_neurons_anno.csv`
- Outputs: `int/P15_TF_readable.csv`, `int/P15_CAM_readable.csv`, `int/Selector_readable.csv`

### `src/utils.r`

- Purpose: shared plotting and data-processing helper functions used across scripts/modules

## Execution notes

- Run these scripts from the repository root so relative paths like `data/` and `int/` resolve correctly.
- Most manual scripts assume preprocessing/statistics outputs already exist.
