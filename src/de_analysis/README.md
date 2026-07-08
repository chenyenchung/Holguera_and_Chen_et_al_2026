# Temporal-Cohort DE Analysis

Runs P15, P30, and P50 Seurat differential expression for right-hemisphere
visual neuron types classified into Early versus Late temporal cohorts, then
summarizes CAM candidates.

## Overview

The workflow has three stages:

1. Temporal-cohort Seurat contrasts that compare Early and Late annotated Ozel
   clusters for each stage.
2. Combined marker, membership, summary, and CAM candidate tables.
3. Early-versus-Late volcano plots.
4. A direct P15, P30, and P50 cluster-pair DE table for Ozel cluster 59 versus
   cluster 82.

Temporal-cohort DE does not depend on neuropil or synapse depth. It uses
confidently annotated types with temporal information and compares
`broad_temp == "Early"` against `broad_temp == "Late"`.

## Usage

```bash
cd src/de_analysis
nextflow run main.nf
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `annf` | `data/visual_neurons_anno.csv` | Visual neuron annotations |
| `camf` | `data/P15_CAM.csv` | CAM gene matrix used to flag candidate surface/guidance molecules |
| `seurat_obj_dir` | `data/ozel_2021_objs` | Directory containing stage Seurat objects |
| `stages` | `P15,P30,P50` | Stage names matched to `<stage>.rds` files |
| `cluster_pair_stages` | `P15,P30,P50` | Stage names used for direct cluster-pair DE |
| `cluster_pair` | `59,82` | Cluster IDs for direct cluster-pair DE; positive log2FC is enriched in the first cluster |
| `cluster_col` | `FinalIdents` | Seurat metadata column used for direct cluster-pair DE |
| `min_cells` | `3` | Minimum Seurat cells per DE group |

## Temporal-Cohort Contrasts

For each stage, the workflow runs:

- `all_projection`
- `notch_on_projection`
- `notch_off_projection`
- `all_intrinsic`

Each temporal contrast compares Early cells against Late cells, so positive
`avg_log2FC` values indicate Early-enriched genes. Ozel clusters assigned to
both cohorts within the same contrast are removed from both groups and logged.
The candidate table reports significant genes that intersect the Ozel
`P15_CAM.csv` cell-adhesion-molecule list.

## Outputs

Published under `int/de_analysis/temporal_cohort/`:

- `markers/<stage>/`: per-stage temporal-cohort marker and membership CSVs.
- `combined_temporal_cohort_de_markers.csv`: all successful temporal-cohort
  marker results.
- `combined_temporal_cohort_de_membership.csv`: temporal contrast group
  definitions and skips.
- `temporal_cohort_cam_candidates.csv`: significant temporal-cohort genes that
  are present in `P15_CAM.csv`.
- `temporal_cohort_de_summary.txt`: summary of temporal contrasts and CAM
  candidate counts.
- `volcano_plots/`: Early-vs-Late volcano plots.

Seurat DE uses `UpdateSeuratObject()` in memory, `FinalIdents` as the Ozel
cluster field, and `FindMarkers()` on the `RNA` assay `data` slot with Wilcoxon,
`min.pct = 0`, and `logfc.threshold = 0`.

## Cluster 59 versus 82

The workflow also runs a direct cluster-pair contrast for clusters 59 and 82 in
P15, P30, and P50. Positive `avg_log2FC` values indicate enrichment in cluster
59. The combined table is published under
`int/de_analysis/cluster_pair_59_vs_82/cluster_59_vs_82_de_table.csv` and sorted
globally by `q_value` ascending, then `abs_avg_log2FC` descending.
