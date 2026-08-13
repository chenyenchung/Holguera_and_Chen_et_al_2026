#!/usr/bin/env nextflow

params.annf = "${workflow.projectDir}/../../data/visual_neurons_anno.csv"
params.utilsf = "${workflow.projectDir}/../utils.r"
params.metaf = "${workflow.projectDir}/../../data/viz_meta.csv"
params.ref_groupsf = "${workflow.projectDir}/../../data/reference_groups.csv"
params.broad_depth_cppf = "${workflow.projectDir}/../stats/bin/broad_depth.cpp"
params.broad_depth_xlsx = "${workflow.projectDir}/../../int/stats/deep_superficial/combined_broad_depth_results.xlsx"
params.mat_prefix = "${workflow.projectDir}/../../int/idv_mat"
params.neuropils = ['ME_R', 'LO_R', 'LOP_R']
params.syn_types = ['pre', 'post']
params.sparse_limit = 100
params.min_neurons = 3
params.coefficient = 0.5
params.n_bootstrap = 1000
params.conf_int = 95
params.type_chunk_size = 5
params.bootstrap_seed = 1

process PrepareTypeDepthInputs {
  cpus 1
  memory '16GB'
  time '2h'
  module 'r/4.5.1'

  input:
  tuple val(np), path(syn)
  path ann
  path meta
  path ref_groups
  val sparse_limit
  val min_neurons
  val type_chunk_size
  val syn_types

  output:
  tuple val("${np}"), path("${np}_type_depth_prepared.rds"), path("${np}_type_depth_chunks.csv"), emit: prepared

  script:
  """
  prepare_type_depth_inputs.r \
    --np ${np} \
    --synf ${syn} \
    --ann ${ann} \
    --meta ${meta} \
    --ref_groups ${ref_groups} \
    --sparse_limit ${sparse_limit} \
    --min_neurons ${min_neurons} \
    --type_chunk_size ${type_chunk_size} \
    --syn_types '${syn_types}'
  """
}

process TypeDepthAnalysis {
  cpus 1
  memory '16GB'
  time '2h'
  module 'r/4.5.1'

  input:
  tuple val(np), path(prepared), val(stype), val(shard_id), val(cell_types)
  path cppsrc
  val sparse_limit
  val min_neurons
  val coefficient
  val n_bootstrap
  val conf_int
  val bootstrap_seed

  output:
  tuple val("${np}"), val("${stype}"), val("${shard_id}"), path("*_type_depth.csv")

  script:
  """
  type_depth.r \
    --np ${np} \
    --prepared ${prepared} \
    --syn_type ${stype} \
    --shard_id ${shard_id} \
    --cell_types '${cell_types}' \
    --cppsrc ${cppsrc} \
    --sparse_limit ${sparse_limit} \
    --min_neurons ${min_neurons} \
    --coefficient ${coefficient} \
    --n_bootstrap ${n_bootstrap} \
    --conf_int ${conf_int} \
    --bootstrap_seed ${bootstrap_seed}
  """
}

process CombineTypeDepthResults {
  cpus 1
  memory '4GB'
  time '15m'
  module 'r/4.5.1'

  input:
  path depth_csvs

  output:
  path "combined_type_depth.csv", emit: depth
  path "type_depth_summary.txt", emit: summary
  path "*_*/*_type_depth.csv", emit: condition_depths

  script:
  """
  combine_type_depth_results.r \
    --depth_out combined_type_depth.csv \
    --summary type_depth_summary.txt \
    --condition_dir .
  """
}

process DsPlot {
  cpus 1
  memory '8GB'
  time '30m'
  module 'r/4.5.1'

  input:
  path ann
  path utils
  path broad_depth_xlsx
  path type_depth

  output:
  path "*.pdf"

  script:
  """
  ds_plot.r \
    --ann ${ann} \
    --utils ${utils} \
    --broad_depth_xlsx ${broad_depth_xlsx} \
    --type_depth ${type_depth} \
    --out_dir .
  """
}

process SpatialTemporalDepthFigure {
  cpus 1
  memory '8GB'
  time '30m'
  module 'r/4.5.1'

  input:
  path ann
  path utils
  path type_depth
  path figure_script

  output:
  path "Supp_Figure_18*_spatial_within_temporal*.pdf", emit: figures
  path "spatial_origin_neuronal_types.csv", emit: csv
  path "spatial_origin_neuronal_types.xlsx", emit: excel

  script:
  """
  Rscript ${figure_script} \
    --ann ${ann} \
    --utils ${utils} \
    --type_depth ${type_depth} \
    --out_dir .
  """
}

def normalizeParamList(raw) {
  if (raw instanceof List) {
    return raw.collect { it.toString().trim() }.findAll { it }
  }
  return raw.toString()
    .split(/[;,]/)
    .collect { it.trim() }
    .findAll { it }
}

workflow {
  main:
  def NP = normalizeParamList(params.neuropils)
  def STYPE = normalizeParamList(params.syn_types)
  def STYPE_ARG = STYPE.join(',')

  np_ch = channel
    .fromList(NP)
    .map { np -> [np, file("${params.mat_prefix}/${np}_rotated.csv.gz")] }

  prepared_ch = PrepareTypeDepthInputs(
    np_ch,
    file(params.annf),
    file(params.metaf),
    file(params.ref_groupsf),
    params.sparse_limit,
    params.min_neurons,
    params.type_chunk_size,
    STYPE_ARG
  )

  shard_ch = prepared_ch.prepared.flatMap { np, prepared, chunks ->
    chunks.readLines()
      .drop(1)
      .collect { line ->
        def cols = line.split(',', -1)
        [cols[0], prepared, cols[1], cols[2] as Integer, cols[3]]
      }
  }

  depth_ch = TypeDepthAnalysis(
    shard_ch,
    file(params.broad_depth_cppf),
    params.sparse_limit,
    params.min_neurons,
    params.coefficient,
    params.n_bootstrap,
    params.conf_int,
    params.bootstrap_seed
  )

  combined_depth_ch = CombineTypeDepthResults(
    depth_ch.map { it -> it[3] }.collect()
  )

  figure_ch = DsPlot(
    file(params.annf),
    file(params.utilsf),
    file(params.broad_depth_xlsx),
    combined_depth_ch.depth
  )

  spatial_figure_ch = SpatialTemporalDepthFigure(
    file(params.annf),
    file(params.utilsf),
    combined_depth_ch.depth,
    file("${workflow.projectDir}/bin/spatial_temporal_depth.r")
  )

  publish:
  type_depth_results = combined_depth_ch.condition_depths
  combined_type_depth = combined_depth_ch.depth
  type_depth_summary = combined_depth_ch.summary
  figures = figure_ch
  spatial_figures = spatial_figure_ch.figures
  spatial_origin_csv = spatial_figure_ch.csv
  spatial_origin_excel = spatial_figure_ch.excel
}

output {
  type_depth_results { path "ds_plot/type_depth/" }
  combined_type_depth { path "ds_plot/type_depth/" }
  type_depth_summary { path "ds_plot/type_depth/" }
  figures { path "./" }
  spatial_figures { path "./" }
  spatial_origin_csv { path "ds_plot/spatial_origin/" }
  spatial_origin_excel { path "ds_plot/spatial_origin/" }
}
