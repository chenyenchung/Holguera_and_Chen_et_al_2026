#!/usr/bin/env nextflow

// CAM Depth Bootstrap Analysis and Visualization
// Mirrors selector_test for data/P15_CAM.csv.

params.camf = 'data/P15_CAM.csv'
params.annf = 'data/visual_neurons_anno.csv'
params.metaf = 'data/viz_meta.csv'
params.ref_groupsf = 'data/reference_groups.csv'
params.utilsf = 'src/utils.r'
params.broad_depth_cppf = 'src/stats/bin/broad_depth.cpp'
params.selector_depth_script = 'src/selector_test/bin/selector_depth.r'
params.selector_viz_script = 'src/visualize/bin/v_selector.r'
params.sparse_limit = 100
params.coefficient = 0.5
params.n_bootstrap = 1000
params.conf_int = 95
params.genes_per_batch = 5
params.subsample = 10000
params.density = 'asis'

process CamDepthAnalysis {
  cpus 1
  memory '32GB'
  time '4h'
  module 'r/4.5.1'

  input:
  tuple val(np), path(syn), val(stype), val(batch_id), val(genes)
  path cam
  path ann
  path meta
  path ref_groups
  path utils
  path cppsrc
  path depth_script
  val coefficient
  val n_bootstrap
  val conf_int
  val slimit

  output:
  tuple val("${np}"), val("${stype}"), val("${batch_id}"), path('*_cam_depth.csv')

  script:
  """
  Rscript ${depth_script} \
    --np ${np} \
    --synf ${syn} \
    --syn_type ${stype} \
    --ts ${cam} \
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

  for f in *_selector_depth.csv; do
    mv "\$f" "\${f/_selector_depth.csv/_cam_depth.csv}"
  done
  """
}

process CombineCamResults {
  cpus 1
  memory '4GB'
  time '15m'
  module 'r/4.5.1'

  input:
  path csvs
  path combine_script

  output:
  path 'cam_depth_summary.txt', emit: summary
  path 'combined_cam_depth.csv', emit: combined
  path 'cam_depth_results.xlsx', emit: excel

  script:
  """
  Rscript ${combine_script} \
    --pattern "_cam_depth\\.csv\$" \
    --summary cam_depth_summary.txt \
    --combined combined_cam_depth.csv \
    --excel cam_depth_results.xlsx
  """
}

process VisualizeCam {
  cpus 1
  memory '14GB'
  time '2h'
  module 'r/4.5.1'

  input:
  tuple val(np), path(syn), val(stype), val(den)
  path cam
  path ann
  path meta
  path utils
  path viz_script
  val subsample
  val slimit

  output:
  tuple val("${np}"), val("${stype}"), path('*.pdf'), optional: true

  script:
  """
  Rscript ${viz_script} \
    --np ${np} \
    --synf ${syn} \
    --syn_type ${stype} \
    --ts ${cam} \
    --density ${den} \
    --ann ${ann} \
    --meta ${meta} \
    --utils ${utils} \
    --subsample ${subsample} \
    --sparse_limit ${slimit}
  """
}

def csvFields(line) {
  line.split(',', -1).collect { it.replaceAll('"', '').trim() }
}

// Helper function to extract CAM row names that survive filter_geneexp().
def extractEligibleCamNames(cam_file, ann_file) {
  def ann_lines = ann_file.readLines()
  if (ann_lines.isEmpty()) {
    throw new IllegalArgumentException("Annotation file is empty: ${ann_file}")
  }

  def ann_header = csvFields(ann_lines[0])
  def confident_idx = ann_header.indexOf('Confident_annotation')
  def cluster_idx = ann_header.indexOf('ozel2021_cluster')
  if (confident_idx < 0 || cluster_idx < 0) {
    throw new IllegalArgumentException(
      "Annotation file must contain Confident_annotation and ozel2021_cluster columns: ${ann_file}"
    )
  }

  def confident_clusters = ann_lines.drop(1)
    .findAll { it.trim() }
    .collect { csvFields(it) }
    .findAll { fields ->
      fields.size() > Math.max(confident_idx, cluster_idx) &&
        fields[confident_idx] == 'Y' &&
        fields[cluster_idx] &&
        fields[cluster_idx] != 'NA'
    }
    .collect { fields -> fields[cluster_idx] }
    .toSet()

  def lines = cam_file.readLines()
  if (lines.isEmpty()) return []

  def cam_header = csvFields(lines[0])
  def eligible_col_idx = (1..<cam_header.size()).findAll { idx ->
    confident_clusters.contains(cam_header[idx])
  }

  def rows = lines.drop(1).findAll { it.trim() }
  def genes = rows
    .collect { csvFields(it) }
    .findAll { fields ->
      eligible_col_idx.any { idx ->
        fields.size() > idx && fields[idx].equalsIgnoreCase('TRUE')
      }
    }
    .collect { fields -> fields[0] }
    .findAll { it.trim() }

  println "Extracted ${genes.size()} eligible CAM rows from ${cam_file.name} (${rows.size() - genes.size()} dropped)"
  return genes
}

workflow {
  main:
  def NP = ['ME_R', 'LO_R', 'LOP_R']
  def STYPE = ['pre', 'post']
  def MAT_PREFIX = 'int/idv_mat/'

  cam_file = file(params.camf)
  ann_file = file(params.annf)
  cam_list = extractEligibleCamNames(cam_file, ann_file)

  if (params.genes_per_batch > 0) {
    def batches = cam_list.collate(params.genes_per_batch)
    cam_batches_ch = channel
      .fromList(batches.withIndex().collect { batch, idx ->
        [String.format("%03d", idx + 1), batch.join(',')]
      })
  } else {
    cam_batches_ch = channel.of(['000', cam_list.join(',')])
  }

  depth_cond_ch = channel
    .fromList(NP)
    .map { np -> [np, file(MAT_PREFIX + np + '_rotated.csv.gz')] }
    .combine(channel.fromList(STYPE))
    .combine(cam_batches_ch)
    .map { np, synfile, stype, batch_id, genes ->
      [np, synfile, stype, batch_id, genes]
    }

  depth_ch = CamDepthAnalysis(
    depth_cond_ch,
    file(params.camf),
    file(params.annf),
    file(params.metaf),
    file(params.ref_groupsf),
    file(params.utilsf),
    file(params.broad_depth_cppf),
    file(params.selector_depth_script),
    params.coefficient,
    params.n_bootstrap,
    params.conf_int,
    params.sparse_limit
  )

  combined_ch = CombineCamResults(
    depth_ch.map { it -> it[3] }.collect(),
    file('src/cam_test/bin/combine_cam_results.r')
  )

  viz_cond_ch = channel
    .fromList(NP)
    .map { np -> [np, file(MAT_PREFIX + np + '_rotated.csv.gz')] }
    .combine(channel.fromList(STYPE))
    .combine(channel.fromList([params.density]))
    .map { np, synfile, stype, den -> [np, synfile, stype, den] }

  viz_ch = VisualizeCam(
    viz_cond_ch,
    file(params.camf),
    file(params.annf),
    file(params.metaf),
    file(params.utilsf),
    file(params.selector_viz_script),
    params.subsample,
    params.sparse_limit
  )

  publish:
  depth_results = depth_ch
  summary = combined_ch.summary
  combined = combined_ch.combined
  excel = combined_ch.excel
  visualizations = viz_ch
}

output {
  depth_results {
    path { input ->
      def np = input[0]
      def stype = input[1]
      return "cam_test/${np}_${stype}"
    }
  }
  summary { path "cam_test/" }
  combined { path "cam_test/" }
  excel { path "cam_test/" }
  visualizations {
    path { input ->
      def np = input[0]
      def stype = input[1]
      return "cam_test/visualization/${np}_${stype}"
    }
  }
}
