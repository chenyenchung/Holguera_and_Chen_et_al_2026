#!/usr/bin/env nextflow

// Selector Depth Bootstrap Analysis
// Tests terminal selector genes for superficial vs deep depth bias

// Parameters
params.selectorsf = 'data/selectors.csv'
params.annf = 'data/visual_neurons_anno.csv'
params.metaf = 'data/viz_meta.csv'
params.ref_groupsf = 'data/reference_groups.csv'
params.utilsf = 'src/utils.r'
params.broad_depth_cppf = 'src/stats/bin/broad_depth.cpp'
params.selector_viz_script = 'src/visualize/bin/v_selector.r'
params.sparse_limit = 100
params.coefficient = 0.5
params.n_bootstrap = 1000
params.conf_int = 95
params.genes_per_batch = 5  // 0 = no batching, >0 = genes per batch
params.subsample = 10000
params.density = 'asis'
params.use_axis_limits = true

process SelectorDepthAnalysis {
  cpus 1
  memory '16GB'
  time '2h'
  module 'r/4.5.1'

  input:
  tuple val(np), path(syn), val(stype), val(batch_id), val(genes)
  path selectors
  path ann
  path meta
  path ref_groups
  path utils
  path cppsrc
  val coefficient
  val n_bootstrap
  val conf_int
  val slimit

  output:
  tuple val("${np}"), val("${stype}"), val("${batch_id}"), path('*.csv')

  script:
  """
  # Run analysis
  selector_depth.r \
    --np ${np} \
    --synf ${syn} \
    --syn_type ${stype} \
    --ts ${selectors} \
    --ann ${ann} \
    --meta ${meta} \
    --ref_groups ${ref_groups} \
    --cppsrc ${cppsrc} \
    --sparse_limit ${slimit} \
    --coefficient ${coefficient} \
    --n_bootstrap ${n_bootstrap} \
    --conf_int ${conf_int} \
    --batch_id ${batch_id} \
    --genes "${genes}"
  """
}

process CombineSelectorResults {
  cpus 1
  memory '4GB'
  time '15m'
  module 'r/4.5.1'

  input:
  path csvs
  path combine_script

  output:
  path 'selector_depth_summary.txt', emit: summary
  path 'combined_selector_depth.csv', emit: combined
  path 'selector_depth_results.xlsx', emit: excel

  script:
  """
  combine_selector_results.r \
    --pattern "_selector_depth\\.csv\$" \
    --summary selector_depth_summary.txt \
    --combined combined_selector_depth.csv \
    --excel selector_depth_results.xlsx
  """
}

process VisualizeSelector {
  cpus 1
  memory '14GB'
  time '2h'
  module 'r/4.5.1'

  input:
  tuple val(np), path(syn), val(stype), val(den)
  path selectors
  path ann
  path meta
  path utils
  path viz_script
  val subsample
  val slimit
  val use_axis_limits

  output:
  tuple val("${np}"), val("${stype}"), path('*.pdf'), optional: true

  script:
  """
  Rscript ${viz_script} \
    --np ${np} \
    --synf ${syn} \
    --syn_type ${stype} \
    --ts ${selectors} \
    --density ${den} \
    --ann ${ann} \
    --meta ${meta} \
    --utils ${utils} \
    --subsample ${subsample} \
    --sparse_limit ${slimit} \
    --use_axis_limits ${use_axis_limits}
  """
}

// Helper function to extract gene names from selectors.csv
def extractGeneNames(selector_file) {
  def lines = selector_file.readLines()
    if (lines.isEmpty()) return []

    def genes = lines.drop(1)
        .findAll { it.trim() } // Skip trailing empty lines
        .collect { line -> 
            line.split(',')[0].replaceAll('"', '').trim() 
        }

    println "Extracted ${genes.size()} genes from ${selector_file.name}"
    return genes
}

workflow {
  main:
  // Only right hemisphere neuropils
  def NP = ['ME_R', 'LO_R', 'LOP_R']
  def STYPE = ['pre', 'post']
  def MAT_PREFIX = 'int/idv_mat/'

  // Extract genes from selectors.csv
  selectors_file = file(params.selectorsf)
  gene_list = extractGeneNames(selectors_file)

  // Create gene batches
  if (params.genes_per_batch > 0) {
    // Split into batches
    def batches = gene_list.collate(params.genes_per_batch)
    gene_batches_ch = channel
      .fromList(batches.withIndex().collect { batch, idx ->
        [String.format("%03d", idx + 1), batch.join(',')]
      })
  } else {
    // No batching - single batch with all genes
    gene_batches_ch = channel.of(['000', gene_list.join(',')])
  }

  // Create input channel: (neuropil, synapse_file, syn_type, batch_id, genes)
  cond_ch = channel
    .fromList(NP)
    .map { np -> [np, file(MAT_PREFIX + np + '_rotated.csv.gz')] }
    .combine(channel.fromList(STYPE))
    .combine(gene_batches_ch)
    .map { np, synfile, stype, batch_id, genes ->
      [np, synfile, stype, batch_id, genes]
    }

  // Run analysis
  analysis_ch = SelectorDepthAnalysis(
    cond_ch,
    file(params.selectorsf),
    file(params.annf),
    file(params.metaf),
    file(params.ref_groupsf),
    file(params.utilsf),
    file(params.broad_depth_cppf),
    params.coefficient,
    params.n_bootstrap,
    params.conf_int,
    params.sparse_limit
  )

  // Combine results
  combined_ch = CombineSelectorResults(
    analysis_ch.map { it -> it[3] }.collect(),
    file('src/selector_test/bin/combine_selector_results.r')
  )

  viz_cond_ch = channel
    .fromList(NP)
    .map { np -> [np, file(MAT_PREFIX + np + '_rotated.csv.gz')] }
    .combine(channel.fromList(STYPE))
    .combine(channel.fromList([params.density]))
    .map { np, synfile, stype, den -> [np, synfile, stype, den] }

  viz_ch = VisualizeSelector(
    viz_cond_ch,
    file(params.selectorsf),
    file(params.annf),
    file(params.metaf),
    file(params.utilsf),
    file(params.selector_viz_script),
    params.subsample,
    params.sparse_limit,
    params.use_axis_limits
  )

  publish:
  // Publish individual results to neuropil-specific directories
  results = analysis_ch
  summary = combined_ch.summary
  combined = combined_ch.combined
  excel = combined_ch.excel
  visualizations = viz_ch
}

output {
  // Individual results go to neuropil_syntype directories
  results {
    path { input ->
      def np = input[0]
      def stype = input[1]
      return "selector_test/${np}_${stype}"
    }
  }
  // Combined results go to main selector_test directory
  summary { path "selector_test/" }
  combined { path "selector_test/" }
  excel { path "selector_test/" }
  visualizations {
    path { input ->
      def np = input[0]
      def stype = input[1]
      return "selector_test/visualization/${np}_${stype}"
    }
  }
}
