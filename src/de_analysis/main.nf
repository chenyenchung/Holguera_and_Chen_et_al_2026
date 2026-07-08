#!/usr/bin/env nextflow

params.annf = 'data/visual_neurons_anno.csv'
params.camf = 'data/P15_CAM.csv'
params.seurat_obj_dir = 'data/ozel_2021_objs'
params.stages = ['P15', 'P30', 'P50']
params.cluster_pair_stages = ['P15', 'P30', 'P50']
params.cluster_pair = ['59', '82']
params.cluster_col = 'FinalIdents'
params.min_cells = 3

process TemporalCohortDE {
  cpus 1
  memory '32GB'
  time '6h'
  module 'r/4.5.1'

  input:
  tuple val(stage), path(seurat_obj)
  path ann
  val min_cells

  output:
  tuple val("${stage}"), path("*_temporal_cohort_de_markers.csv"), path("*_temporal_cohort_de_membership.csv")

  script:
  """
  temporal_cohort_de.r \
    --stage ${stage} \
    --ann ${ann} \
    --seurat ${seurat_obj} \
    --min_cells ${min_cells}
  """
}

process CombineTemporalCohortResults {
  cpus 1
  memory '8GB'
  time '30m'
  module 'r/4.5.1'

  input:
  path marker_csvs
  path membership_csvs
  path cam

  output:
  path "combined_temporal_cohort_de_markers.csv", emit: markers
  path "combined_temporal_cohort_de_membership.csv", emit: membership
  path "temporal_cohort_cam_candidates.csv", emit: candidates
  path "temporal_cohort_de_summary.txt", emit: summary

  script:
  """
  combine_temporal_cohort_results.r \
    --cam ${cam} \
    --markers_out combined_temporal_cohort_de_markers.csv \
    --membership_out combined_temporal_cohort_de_membership.csv \
    --candidates_out temporal_cohort_cam_candidates.csv \
    --summary temporal_cohort_de_summary.txt
  """
}

process TemporalVolcanoPlots {
  cpus 1
  memory '8GB'
  time '30m'
  module 'r/4.5.1'

  input:
  path markers_csv

  output:
  path "volcano_plots/*", emit: plots

  script:
  """
  volcano_plots.r \
    --markers ${markers_csv} \
    --out_dir volcano_plots
  """
}

process ClusterPairDE {
  cpus 1
  memory '32GB'
  time '6h'
  module 'r/4.5.1'

  input:
  tuple val(stage), path(seurat_obj)
  val cluster_1
  val cluster_2
  val cluster_col
  val min_cells

  output:
  tuple val("${stage}"), path("*_cluster_${cluster_1}_vs_${cluster_2}_de_markers.csv")

  script:
  """
  cluster_pair_de.r \
    --stage ${stage} \
    --seurat ${seurat_obj} \
    --cluster_1 ${cluster_1} \
    --cluster_2 ${cluster_2} \
    --cluster_col ${cluster_col} \
    --min_cells ${min_cells}
  """
}

process CombineClusterPairResults {
  cpus 1
  memory '8GB'
  time '30m'
  module 'r/4.5.1'

  input:
  path marker_csvs
  val cluster_1
  val cluster_2

  output:
  path "cluster_${cluster_1}_vs_${cluster_2}_de_table.csv", emit: markers

  script:
  """
  combine_cluster_pair_results.r \
    --markers_out cluster_${cluster_1}_vs_${cluster_2}_de_table.csv
  """
}

workflow {
  main:
  def STAGES = params.stages instanceof CharSequence
    ? params.stages.split(/[,;]/).collect { it.trim() }.findAll { it }
    : params.stages
  def CLUSTER_PAIR_STAGES = params.cluster_pair_stages instanceof CharSequence
    ? params.cluster_pair_stages.split(/[,;]/).collect { it.trim() }.findAll { it }
    : params.cluster_pair_stages
  def CLUSTER_PAIR = params.cluster_pair instanceof CharSequence
    ? params.cluster_pair.split(/[,;]/).collect { it.trim() }.findAll { it }
    : params.cluster_pair
  def CLUSTER_1 = CLUSTER_PAIR[0] as String
  def CLUSTER_2 = CLUSTER_PAIR[1] as String

  stage_ch = channel
    .fromList(STAGES)
    .map { stage -> [stage, file("${params.seurat_obj_dir}/${stage}.rds")] }
  cluster_pair_stage_ch = channel
    .fromList(CLUSTER_PAIR_STAGES)
    .map { stage -> [stage, file("${params.seurat_obj_dir}/${stage}.rds")] }

  temporal_de_ch = TemporalCohortDE(
    stage_ch,
    file(params.annf),
    params.min_cells
  )

  combined_temporal_ch = CombineTemporalCohortResults(
    temporal_de_ch.map { it -> it[1] }.collect(),
    temporal_de_ch.map { it -> it[2] }.collect(),
    file(params.camf)
  )

  temporal_volcano_ch = TemporalVolcanoPlots(combined_temporal_ch.markers)

  cluster_pair_de_ch = ClusterPairDE(
    cluster_pair_stage_ch,
    CLUSTER_1,
    CLUSTER_2,
    params.cluster_col,
    params.min_cells
  )

  combined_cluster_pair_ch = CombineClusterPairResults(
    cluster_pair_de_ch.map { it -> it[1] }.collect(),
    CLUSTER_1,
    CLUSTER_2
  )

  publish:
  temporal_de_results = temporal_de_ch
  temporal_combined_markers = combined_temporal_ch.markers
  temporal_combined_membership = combined_temporal_ch.membership
  temporal_cam_candidates = combined_temporal_ch.candidates
  temporal_summary = combined_temporal_ch.summary
  temporal_volcano_plots = temporal_volcano_ch.plots
  cluster_pair_de_results = cluster_pair_de_ch
  cluster_pair_combined_markers = combined_cluster_pair_ch.markers
}

output {
  temporal_de_results {
    path { input ->
      return "de_analysis/temporal_cohort/markers/${input[0]}"
    }
  }
  temporal_combined_markers { path "de_analysis/temporal_cohort/" }
  temporal_combined_membership { path "de_analysis/temporal_cohort/" }
  temporal_cam_candidates { path "de_analysis/temporal_cohort/" }
  temporal_summary { path "de_analysis/temporal_cohort/" }
  temporal_volcano_plots { path "de_analysis/temporal_cohort/volcano_plots/" }
  cluster_pair_de_results {
    path { input ->
      return "de_analysis/cluster_pair_59_vs_82/markers/${input[0]}"
    }
  }
  cluster_pair_combined_markers {
    path "de_analysis/cluster_pair_59_vs_82/"
  }
}
