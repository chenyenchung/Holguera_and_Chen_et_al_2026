#!/usr/bin/env nextflow

// Default parameters
params.presetf = 'data/viz_preset.csv'
params.metaf = 'data/viz_meta.csv'
params.annf = 'data/visual_neurons_anno.csv'
params.utilsf = 'src/utils.r'
params.depth_stats_cppf = 'src/stats/bin/depth_stats.cpp'
params.broad_depth_cppf = 'src/stats/bin/broad_depth.cpp'
params.combine_scriptf = 'src/stats/bin/combine_results.r'
params.reference_depth_scriptf = 'src/stats/bin/reference_depth_distribution.r'
params.ref_groupsf = 'data/reference_groups.csv'
params.lo_l_synf = 'int/idv_mat/LO_L_rotated.csv.gz'
params.lo_r_synf = 'int/idv_mat/LO_R_rotated.csv.gz'
params.calibration_types = []
params.n_quantiles = 1000
params.n_bootstrap = 1000
params.conf_int = 95
params.ks_subsample_size = 1000
params.ks_n_iterations = 1000
params.ks_correction_method = 'fdr'
params.sparse_limit = 100
params.broad_depth_coefficient = 0.5
params.broad_depth_n_bootstrap = 1000
params.broad_depth_conf_int = 95

process DepthStatsAnalysis {
  cpus '1'
  memory '8GB'
  time '3h'
  module 'r/4.5.1'

  input:
  tuple val(np), path(syn), val(stype), val(preset)
  path ann
  path meta
  path presetf
  path utils
  path cppsrc
  val n_quantiles
  val n_bootstrap
  val conf_int
  val slimit
  val ks_subsample_size
  val ks_n_iterations
  val ks_correction_method

  output:
  tuple val("${np}"), val("${preset}"), val("${stype}"), 
        path('*.csv')

  script:
  """
  depth_stats.r \
    --np ${np} \
    --synf ${syn} \
    --syn_type ${stype} \
    --use_preset ${preset} \
    --ann ${ann} \
    --meta ${meta} \
    --preset ${presetf} \
    --utils ${utils} \
    --cppsrc ${cppsrc} \
    --sparse_limit ${slimit} \
    --n_quantiles ${n_quantiles} \
    --n_bootstrap ${n_bootstrap} \
    --conf_int ${conf_int} \
    --ks_subsample_size ${ks_subsample_size} \
    --ks_n_iterations ${ks_n_iterations} \
    --ks_correction_method ${ks_correction_method}
  """
}

process VisualizeStatSummary {
  cpus '1'
  memory '8GB'
  time '30m'
  module 'r/4.5.1'

  input:
  tuple val(np), val(preset), val(stype), path(results_csv)
  path meta
  path utils

  output:
  tuple val("${np}"), val("${preset}"), val("${stype}"), path('*.pdf'), optional: true

  script:
  def output_prefix = results_csv.baseName.replaceAll('_depth_stats', '')
  """
  visualize_stat_summary.r \
    --input_file ${results_csv} \
    --output_prefix ${output_prefix} \
    --utils ${utils}
  """
}

process CombineResults {
  cpus '1'
  memory '4GB'
  time '15m'
  module 'r/4.5.1'

  input:
  path csvs
  path combine_script

  output:
  path 'statistical_summary.txt', emit: summary
  path 'combined_depth_stats_results.xlsx', emit: excel

  script:
  """
  combine_results.r
  """
}

process FunctionalEnrichment {
  cpus '1'
  memory '4GB'
  time '30m'
  module 'r/4.5.1'

  input:
  path ann
  path utils

  output:
  path '*.pdf', emit: plots
  path '*.xlsx', emit: stats

  script:
  """
  functional_enrichment.r --anno ${ann}
  """
}

process OpcSynapseRatio {
  cpus '1'
  memory '8GB'
  time '30m'
  module 'r/4.5.1'

  input:
  path ann
  path lo_l_syn
  path lo_r_syn

  output:
  path 'LO_opc_synapse_ratio.csv', emit: csv

  script:
  """
  opc_synapse_ratio.r \
    --ann ${ann} \
    --lo_l ${lo_l_syn} \
    --lo_r ${lo_r_syn} \
    --output LO_opc_synapse_ratio.csv
  """
}

process BroadDepthAnalysis {
  cpus '1'
  memory '8GB'
  time '8h'
  module 'r/4.5.1'

  input:
  tuple val(np), path(syn), val(stype), val(preset)
  path ann
  path meta
  path presetf
  path utils
  path cppsrc
  path ref_groups
  val coefficient
  val n_bootstrap
  val conf_int
  val slimit

  output:
  tuple val("${np}"), val("${preset}"), val("${stype}"),
        path('*.csv')

  script:
  """
  broad_depth.r \
    --np ${np} \
    --synf ${syn} \
    --syn_type ${stype} \
    --use_preset ${preset} \
    --ann ${ann} \
    --meta ${meta} \
    --preset ${presetf} \
    --utils ${utils} \
    --cppsrc ${cppsrc} \
    --ref_groups ${ref_groups} \
    --sparse_limit ${slimit} \
    --coefficient ${coefficient} \
    --n_bootstrap ${n_bootstrap} \
    --conf_int ${conf_int}
  """
}

process CombineBroadDepthResults {
  cpus '1'
  memory '4GB'
  time '15m'
  module 'r/4.5.1'

  input:
  path csvs
  path combine_script

  output:
  path 'broad_depth_summary.txt', emit: summary
  path 'combined_broad_depth_results.xlsx', emit: excel

  script:
  """
  combine_results.r \
    --pattern "_broad_depth\\.csv\$" \
    --summary broad_depth_summary.txt \
    --excel combined_broad_depth_results.xlsx
  """
}

process CombineBroadDepthResultsLight {
  cpus '1'
  memory '4GB'
  time '15m'
  module 'r/4.5.1'

  input:
  path csvs
  path combine_script

  output:
  path 'broad_depth_summary_light.txt', emit: summary
  path 'combined_broad_depth_results_light.xlsx', emit: excel

  script:
  """
  combine_results_light.r \
    --pattern "_broad_depth\\.csv\$" \
    --summary broad_depth_summary_light.txt \
    --excel combined_broad_depth_results_light.xlsx
  """
}

process ReferenceDepthDistribution {
  cpus '1'
  memory '8GB'
  time '1h'
  module 'r/4.5.1'

  input:
  tuple val(np), path(syn)
  path meta
  path ref_groups
  path script_file
  val calibration_types

  output:
  tuple val("${np}"), path('*.pdf'), path('*.csv')

  script:
  def output_prefix = "${np}_reference_depth_distribution"
  """
  Rscript ${script_file} \
    --np ${np} \
    --synf ${syn} \
    --meta ${meta} \
    --ref_groups ${ref_groups} \
    --types "${calibration_types}" \
    --output_prefix ${output_prefix}
  """
}

def normalizeCalibrationTypes(raw_types) {
  if (raw_types == null) {
    return ''
  }
  if (raw_types instanceof List) {
    return raw_types.collect { it.toString().trim() }.findAll { it }.join(',')
  }
  return raw_types.toString()
    .split(/[;,]/)
    .collect { it.trim() }
    .findAll { it }
    .join(',')
}

workflow {
  main:
  // Define analysis parameters
  def NP = ['ME_L', 'ME_R', 'LOP_L', 'LOP_R', 'LO_L', 'LO_R']
  def STYPE = ['pre', 'post']
  def MAT_PREFIX = 'int/idv_mat/'
  def SPARSE_LIMIT = params.sparse_limit
  def CALIBRATION_TYPES = normalizeCalibrationTypes(params.calibration_types)
  def STATS_PRESETS = [
    'temporal_known',
    'subsystem_known',
    'temporal_new',
    'subsystem_new',
    'broad_known',
    'broad_new',
    'temporal_all',
    'type_putative_1',
    'type_putative_2',
    'type_putative_3',
    'type_putative_4',
    'type_putative_5',
    'type_putative_6',
    'type_putative_7',
    'type_putative_8',
    'type_putative_9',
    'subsystem_putative',
    'spatial_all',
    'spatial_notch'
  ]

  // Create input channel
  cond_ch = channel
    .fromList(NP)
    .map { it ->
      def synp = MAT_PREFIX + it + '_rotated.csv.gz'
      return [it, file(synp)]
    }
    .combine(channel.fromList(STYPE))
    .combine(
      channel.fromPath(file(params.presetf))
        .splitCsv(header:true)
        .map { row -> row.preset }
        .filter { preset -> preset in STATS_PRESETS }
    )

  reference_cond_ch = channel
    .fromList(NP)
    .map { it ->
      def synp = MAT_PREFIX + it + '_rotated.csv.gz'
      return [it, file(synp)]
    }
    
    
  // Perform statistical analysis
//  analysis_ch = DepthStatsAnalysis(
//    cond_ch,
//    file(params.annf),
//    file(params.metaf),
//    file(params.presetf),
//    file(params.utilsf),
//    file(params.depth_stats_cppf),
//    params.n_quantiles,
//    params.n_bootstrap,
//    params.conf_int,
//    SPARSE_LIMIT,
//    params.ks_subsample_size,
//    params.ks_n_iterations,
//    params.ks_correction_method
//  )

  // Generate visualizations
//  viz_ch = VisualizeStatSummary(
//    analysis_ch,
//    file(params.metaf),
//    file(params.utilsf)
//  )

  // Combine all results
//  combined_ch = CombineResults(
//    analysis_ch.map { it -> it[3] }.collect(),
//    file(params.combine_scriptf)
//  )

  // Run functional enrichment analysis
  functional_ch = FunctionalEnrichment(
    file(params.annf),
    file(params.utilsf)
  )

  // Calculate putative OPC synapse ratios in LO
  opc_ratio_ch = OpcSynapseRatio(
    file(params.annf),
    file(params.lo_l_synf),
    file(params.lo_r_synf)
  )

  reference_depth_ch = ReferenceDepthDistribution(
    reference_cond_ch,
    file(params.metaf),
    file(params.ref_groupsf),
    file(params.reference_depth_scriptf),
    CALIBRATION_TYPES
  )

  // Run broad depth analysis
  broad_depth_ch = BroadDepthAnalysis(
    cond_ch,
    file(params.annf),
    file(params.metaf),
    file(params.presetf),
    file(params.utilsf),
    file(params.broad_depth_cppf),
    file(params.ref_groupsf),
    params.broad_depth_coefficient,
    params.broad_depth_n_bootstrap,
    params.broad_depth_conf_int,
    SPARSE_LIMIT
  )

  // Combine broad depth results
  combined_broad_depth_ch = CombineBroadDepthResults(
    broad_depth_ch.map { it -> it[3] }.collect(),
    file(params.combine_scriptf)
  )

  combined_broad_depth_light_ch = CombineBroadDepthResultsLight(
    broad_depth_ch.map { it -> it[3] }.collect(),
    file(params.combine_scriptf)
  )

  publish:
//  stats_plots = viz_ch
//  summary = combined_ch.summary
//  excel = combined_ch.excel
  functional_plots = functional_ch.plots
  functional_stats = functional_ch.stats
  opc_synapse_ratio = opc_ratio_ch.csv
  reference_depth_distribution = reference_depth_ch
  broad_depth_results = broad_depth_ch
  broad_depth_summary = combined_broad_depth_ch.summary
  broad_depth_excel = combined_broad_depth_ch.excel
  broad_depth_excel_light = combined_broad_depth_light_ch.excel
}

output {
//  stats_plots {
//    path { input ->
//      def np = input[0]
//      def preset = input[1]
//      def stype = input[2]
//      return "stats/${preset}/${np}_${stype}"
//    }
//  }
//  summary { path "stats/" }
//  excel { path "stats/" }
  functional_plots { path "stats/functional_enrichment/" }
  functional_stats { path "stats/functional_enrichment/" }
  opc_synapse_ratio { path "stats/opc_synapse_ratio/" }
  reference_depth_distribution { path "stats/reference_depth_distribution/" }
  broad_depth_results {
    path { input ->
      def np = input[0]
      def preset = input[1]
      def stype = input[2]
      return "stats/deep_superficial/${preset}/${np}_${stype}"
    }
  }
  broad_depth_summary { path "stats/deep_superficial/" }
  broad_depth_excel { path "stats/deep_superficial/" }
  broad_depth_excel_light { path "stats/deep_superficial/" }
}
