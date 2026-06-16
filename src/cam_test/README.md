# CAM Depth Bootstrap Analysis

Tests Cell Adhesion Molecule expression rows from `data/P15_CAM.csv` for
superficial vs deep depth bias and generates selector-style synapse
visualizations.

## Usage

```bash
nextflow run src/cam_test/main.nf
```

The default run processes `ME_R`, `LO_R`, and `LOP_R`, for both presynaptic and
postsynaptic synapses. All CAM rows are analyzed and visualized by default.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `camf` | `data/P15_CAM.csv` | CAM expression matrix |
| `annf` | `data/visual_neurons_anno.csv` | Cell type annotations |
| `metaf` | `data/viz_meta.csv` | Visualization metadata |
| `ref_groupsf` | `data/reference_groups.csv` | Superficial/deep reference groups |
| `sparse_limit` | `100` | Minimum synapses per CAM row/condition |
| `coefficient` | `0.5` | Threshold coefficient for depth bias |
| `n_bootstrap` | `1000` | Bootstrap iterations |
| `conf_int` | `95` | Confidence interval level |
| `genes_per_batch` | `5` | Eligible CAM rows per depth-analysis batch; `0` runs one batch |
| `subsample` | `10000` | Synapses sampled per visualization process |
| `density` | `asis` | Density mode passed to the visualization script |

## Outputs

Results are published under `int/cam_test/`.

```text
cam_test/
├── ME_R_pre/
│   └── *_cam_depth.csv
├── ME_R_post/
│   └── *_cam_depth.csv
├── LO_R_pre/
├── LO_R_post/
├── LOP_R_pre/
├── LOP_R_post/
├── combined_cam_depth.csv
├── cam_depth_results.xlsx
├── cam_depth_summary.txt
└── visualization/
    ├── ME_R_pre/
    ├── ME_R_post/
    ├── LO_R_pre/
    ├── LO_R_post/
    ├── LOP_R_pre/
    └── LOP_R_post/
```

The depth workflow reuses the selector bootstrap implementation with
`data/P15_CAM.csv` as the expression matrix. Combined FDR values are recomputed
after all batches are collected, grouped by neuropil, synapse type, and Notch
category.

Before batching, CAM rows are filtered to match the expression columns that can
be mapped through confident annotations in `visual_neurons_anno.csv`. Rows with
no `TRUE` value in those mapped columns are skipped so empty batches are not
submitted to the depth-analysis process.

CAM row names are preserved in result tables. Filenames are produced by the
existing selector-style visualization script and should be treated as generated
artifacts.
