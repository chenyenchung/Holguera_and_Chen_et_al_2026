#!/usr/bin/env nextflow

params.annf = 'data/visual_neurons_anno.csv'
params.metaf = 'data/viz_meta.csv'
params.ref_groupsf = 'data/reference_groups.csv'
params.seurat_obj_dir = 'data/ozel_2021_objs'
params.stages = ['P15', 'P30', 'P50']
params.broad_depth_cppf = 'src/stats/bin/broad_depth.cpp'
params.sparse_limit = 100
params.min_neurons = 3
params.coefficient = 0.5
params.n_bootstrap = 1000
params.conf_int = 95
params.min_cells = 3

process TypeDepthAnalysis {
  cpus 1
  memory '16GB'
  time '8h'
  module 'r/4.5.1'

  input:
  tuple val(np), path(syn), val(stype)
  path ann
  path meta
  path ref_groups
  path cppsrc
  val sparse_limit
  val min_neurons
  val coefficient
  val n_bootstrap
  val conf_int

  output:
  tuple val("${np}"), val("${stype}"), path("*_type_depth.csv")

  script:
  """
  type_depth.r \
    --np ${np} \
    --synf ${syn} \
    --syn_type ${stype} \
    --ann ${ann} \
    --meta ${meta} \
    --ref_groups ${ref_groups} \
    --cppsrc ${cppsrc} \
    --sparse_limit ${sparse_limit} \
    --min_neurons ${min_neurons} \
    --coefficient ${coefficient} \
    --n_bootstrap ${n_bootstrap} \
    --conf_int ${conf_int}
  """
}

process SeuratDE {
  cpus 1
  memory '32GB'
  time '6h'
  module 'r/4.5.1'

  input:
  tuple val(stage), val(np), val(stype), path(depth_csv), path(seurat_obj)
  val min_cells

  output:
  tuple val("${stage}"), val("${np}"), val("${stype}"), path("*_de_markers.csv"), path("*_de_membership.csv")

  script:
  """
  seurat_de.r \
    --stage ${stage} \
    --depth ${depth_csv} \
    --seurat ${seurat_obj} \
    --np ${np} \
    --syn_type ${stype} \
    --min_cells ${min_cells}
  """
}

process CombineDEResults {
  cpus 1
  memory '8GB'
  time '30m'
  module 'r/4.5.1'

  input:
  path depth_csvs
  path marker_csvs
  path membership_csvs

  output:
  path "combined_de_markers.csv", emit: markers
  path "combined_de_membership.csv", emit: membership
  path "combined_type_depth.csv", emit: depth
  path "de_analysis_summary.txt", emit: summary

  script:
  """
  combine_de_results.r \
    --markers_out combined_de_markers.csv \
    --membership_out combined_de_membership.csv \
    --depth_out combined_type_depth.csv \
    --summary de_analysis_summary.txt
  """
}

process VolcanoPlots {
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

workflow {
  main:
  def NP = ['ME_R', 'LO_R', 'LOP_R']
  def STYPE = ['pre', 'post']
  def MAT_PREFIX = 'int/idv_mat/'
  def STAGES = params.stages instanceof CharSequence
    ? params.stages.split(/[,;]/).collect { it.trim() }.findAll { it }
    : params.stages

  cond_ch = channel
    .fromList(NP)
    .map { np -> [np, file(MAT_PREFIX + np + '_rotated.csv.gz')] }
    .combine(channel.fromList(STYPE))
    .map { np, synfile, stype -> [np, synfile, stype] }

  depth_ch = TypeDepthAnalysis(
    cond_ch,
    file(params.annf),
    file(params.metaf),
    file(params.ref_groupsf),
    file(params.broad_depth_cppf),
    params.sparse_limit,
    params.min_neurons,
    params.coefficient,
    params.n_bootstrap,
    params.conf_int
  )

  stage_ch = channel
    .fromList(STAGES)
    .map { stage -> [stage, file("${params.seurat_obj_dir}/${stage}.rds")] }

  de_input_ch = depth_ch
    .combine(stage_ch)
    .map { np, stype, depth_csv, stage, seurat_obj -> [stage, np, stype, depth_csv, seurat_obj] }

  de_ch = SeuratDE(
    de_input_ch,
    params.min_cells
  )

  combined_ch = CombineDEResults(
    depth_ch.map { it -> it[2] }.collect(),
    de_ch.map { it -> it[3] }.collect(),
    de_ch.map { it -> it[4] }.collect()
  )

  volcano_ch = VolcanoPlots(combined_ch.markers)

  publish:
  depth_results = depth_ch
  de_results = de_ch
  combined_markers = combined_ch.markers
  combined_membership = combined_ch.membership
  combined_depth = combined_ch.depth
  summary = combined_ch.summary
  volcano_plots = volcano_ch.plots
}

output {
  depth_results {
    path { input ->
      return "de_analysis/type_depth/${input[0]}_${input[1]}"
    }
  }
  de_results {
    path { input ->
      return "de_analysis/markers/${input[0]}/${input[1]}_${input[2]}"
    }
  }
  combined_markers { path "de_analysis/" }
  combined_membership { path "de_analysis/" }
  combined_depth { path "de_analysis/" }
  summary { path "de_analysis/" }
  volcano_plots { path "de_analysis/volcano_plots/" }
}
