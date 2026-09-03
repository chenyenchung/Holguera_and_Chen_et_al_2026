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
- `src/ds_plot/`

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

### `src/manual_c59_vs_c82.R`

- Purpose: validate the adult `toy+/ham+/D-/Vsx1-` main-OPC cluster assignment
  and visualize the distinction between clusters 59 and 82.
- Key inputs: the adult Seurat object and adult ScMarco modeled-expression
  matrix; cluster identifiers are read as numeric `FinalIdents`.
- Selection: main-OPC clusters with an adult ScMarco score greater than 0.5
  for either `toy` or `ham`.
- Normalization: selected cells are explicitly re-normalized from `RNA` counts
  with Seurat `LogNormalize` and a scale factor of 10,000.
- Outputs: stacked violin and dot-plot PDFs, an all-main-OPC ScMarco audit,
  normalized expression summaries, c59-versus-c82 Wilcoxon tests, and a text
  validation report under `int/c59_vs_c82/`.

```bash
Rscript src/manual_c59_vs_c82.R
```

The adult object is approximately 3.6 GB. Run the script in a Slurm allocation
with at least 32 GB of memory. Use `--help` to see path and analysis overrides.

## Top-level helper scripts

### `src/get_selectors.R`

- Purpose: rebuild a binary selector matrix from mixture-model objects and TF list
- Note: references an external absolute path for input objects; update before use
- Output: currently writes `data/selector.csv`

### `src/convert_gene_expression_to_readable_format.R`

- Purpose: convert binary gene-expression matrices to "On/Off" human-readable tables with renamed columns
- Inputs: `data/P15_tf.csv`, `data/P15_CAM.csv`, `data/selectors.csv`, `data/visual_neurons_anno.csv`
- Outputs: `int/P15_TF_readable.csv`, `int/P15_CAM_readable.csv`, `int/Selector_readable.csv`

### `src/export_csm_neuronal_types.R`

- Purpose: export the P15 CSM expression calls as a neuronal-type-oriented
  Excel workbook.
- Inputs: `data/P15_CAM.csv` and `data/visual_neurons_anno.csv` by default.
- Mapping rule: include only confidently annotated neuronal types and expand a
  positive Ozel cluster call to every confident neuronal type mapped to that
  cluster.
- Output: `int/P15_CSM_neuronal_types.xlsx` by default, with summary, long,
  matrix, type-metadata, mapping-audit, and provenance sheets.

```bash
Rscript src/export_csm_neuronal_types.R \
  --cam data/P15_CAM.csv \
  --ann data/visual_neurons_anno.csv \
  --output int/P15_CSM_neuronal_types.xlsx
```

### `src/utils.r`

- Purpose: shared plotting and data-processing helper functions used across scripts/modules

## Execution notes

- Run these scripts from the repository root so relative paths like `data/` and `int/` resolve correctly.
- Most manual scripts assume preprocessing/statistics outputs already exist.
