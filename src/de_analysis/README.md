# Deep/Superficial DE Analysis

Runs P15, P30, and P50 Seurat differential expression for right-hemisphere
visual neuron types classified as superficial or deep in ME, LO, and LOP.

## Overview

The workflow has two stages:

1. Type-depth bootstrap testing for every row in
   `data/visual_neurons_anno.csv` where `Confident_annotation == "Y"` and
   `ozel2021_cluster` is not missing.
2. Seurat `FindMarkers()` contrasts using stage objects under
   `data/ozel_2021_objs/`.

Depth tests run once, separately for `pre` and `post` synapse depths in `ME_R`,
`LO_R`, and `LOP_R`. DE then fans out over `P15`, `P30`, and `P50` and uses
only types with `direction` equal to `superficial` or `deep` and
`significant_fdr == TRUE`.

## Usage

```bash
cd src/de_analysis
nextflow run main.nf
```

For a fast smoke test:

```bash
nextflow run main.nf --n_bootstrap 10
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `annf` | `data/visual_neurons_anno.csv` | Visual neuron annotations |
| `metaf` | `data/viz_meta.csv` | Neuropil depth-axis metadata |
| `ref_groupsf` | `data/reference_groups.csv` | Superficial/deep reference definitions |
| `seurat_obj_dir` | `data/ozel_2021_objs` | Directory containing stage Seurat objects |
| `stages` | `P15,P30,P50` | Stage names matched to `<stage>.rds` files |
| `sparse_limit` | `100` | Minimum synapses for a type-depth test |
| `min_neurons` | `3` | Minimum neurons for a type-depth test |
| `coefficient` | `0.5` | Depth-test threshold coefficient |
| `n_bootstrap` | `1000` | Bootstrap iterations |
| `conf_int` | `95` | Bootstrap confidence interval |
| `min_cells` | `3` | Minimum Seurat cells per DE group |

## Contrasts

For each `neuropil x syn_type`, the workflow runs:

- `all`
- `notch_on_projection`
- `notch_off_projection`
- `notch_on_intrinsic`
- `notch_off_intrinsic`

Each contrast compares superficial cells against deep cells. Ozel clusters are
included if any mapped eligible cell type matches the contrast stratum and
depth group. If a cluster is assigned to both superficial and deep within the
same contrast, it is removed from both groups and logged in the membership file.

## Outputs

Published under `int/de_analysis/`:

- `type_depth/`: per-neuropil/per-syn-type depth-test CSVs.
- `markers/<stage>/`: per-stage, per-neuropil/per-syn-type DE marker CSVs and
  membership CSVs.
- `combined_type_depth.csv`: all type-depth rows.
- `combined_de_membership.csv`: all contrast group definitions and skips.
- `combined_de_markers.csv`: all successful Seurat marker results.
- `de_analysis_summary.txt`: summary of depth tests, skipped contrasts, and DE output.

Seurat DE uses `UpdateSeuratObject()` in memory, `FinalIdents` as the Ozel
cluster field, and `FindMarkers()` on the `RNA` assay `data` slot with Wilcoxon,
`min.pct = 0`, and `logfc.threshold = 0`.
