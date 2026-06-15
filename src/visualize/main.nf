#!/usr/bin/env nextflow
params.presetf = 'data/viz_preset.csv'
params.metaf = 'data/viz_meta.csv'
params.annf = 'data/visual_neurons_anno.csv'
params.tff = 'data/P15_tf.csv'
params.camf = 'data/P15_CAM.csv'
params.tsf = 'data/selectors.csv'
params.distances = 'data/TypeToTypeDistances.csv'

process Visualize {
  cpus '1'
  memory '8GB'
  time '30m'
  module 'r/gcc/4.5.0'

  input:
  tuple val(np), path(syn), val(stype), val(preset), val(den)
  path ann
  path meta
  path presetf
  path utils
  val subsample
  val slimit

  output:
  tuple val("${np}"), val("${preset}"), path('*.pdf'), optional: true

  script:
  """
  visualize.r \
    --np ${np} \
    --synf ${syn} \
    --syn_type ${stype} \
    --use_preset ${preset} \
    --density ${den} \
    --ann ${ann} \
    --meta ${meta} \
    --preset ${presetf} \
    --utils ${utils} \
    --subsample ${subsample} \
    --sparse_limit ${slimit}
  """
}

process VisualizeSelector {
  cpus '1'
  memory '14GB'
  time '2h'
  module 'r/gcc/4.5.0'

  input:
  tuple val(np), path(syn), val(stype), val(den), path(ts), val(gselection)
  path ann
  path meta
  path utils
  val subsample
  val slimit

  output:
  tuple val("${np}"), val("${gselection}"), path('*.pdf'), optional: true

  script:
  """
  v_selector.r \
    --np ${np} \
    --synf ${syn} \
    --syn_type ${stype} \
    --ts ${ts} \
    --density ${den} \
    --ann ${ann} \
    --meta ${meta} \
    --utils ${utils} \
    --subsample ${subsample} \
    --sparse_limit ${slimit}
  """
}

process SimilarityTree {
  cpus '1'
  memory '4GB'
  time '10m'
  module 'r/gcc/4.5.0'

  input:
  path distances
  path ann
  path utils

  output:
  path 'similarity_tree.pdf'

  script:
  """
  similarity_tree.r \
    --distances ${distances} \
    --annotations ${ann} \
    --utils ${utils} \
    --output similarity_tree.pdf
  """
}

process PartnerExtraction {
  cpus '2'
  memory '16GB'
  time '60m'
  module 'r/gcc/4.5.0'

  input:
  path me_l
  path me_r
  path lo_l
  path lo_r
  path lop_l
  path lop_r
  path utils

  output:
  path 'partner_summary.csv'

  script:
  """
  partner_extraction.r \
    --output partner_summary.csv \
    --me_l ${me_l} \
    --me_r ${me_r} \
    --lo_l ${lo_l} \
    --lo_r ${lo_r} \
    --lop_l ${lop_l} \
    --lop_r ${lop_r}
  """
}

process NeuropilPartnerAnalysis {
  cpus '1'
  memory '4GB'
  time '10m'
  module 'r/gcc/4.5.0'

  input:
  path partner_data
  path utils
  tuple val(neuropil), val(syn_type)
  val threshold

  output:
  tuple val("${neuropil}"), path('*.pdf')

  script:
  """
  neuropil_partner_analysis.r \
    --partner_data ${partner_data} \
    --threshold ${threshold} \
    --np ${neuropil} \
    --syn ${syn_type}
  """
}

workflow {
  main:
  def NP = ['ME_L', 'ME_R', 'LOP_L', 'LOP_R', 'LO_L', 'LO_R']
  def STYPE = ['pre', 'post']
  def DEN_ALGO = ['asis']
  def MAT_PREFIX = 'int/idv_mat/'
  def SUBSAMPLE_TO = 10000
  def SPARSE_LIMIT = 100
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
    )
    .combine(channel.fromList(DEN_ALGO))
  utils_file = file('src/utils.r')
  out_ch = Visualize(
    cond_ch,
    file(params.annf),
    file(params.metaf),
    file(params.presetf),
    utils_file,
    SUBSAMPLE_TO,
    SPARSE_LIMIT
  )
  

//  scond_ch = channel
//    .fromList(NP)
//    .map { it ->
//      def synp = MAT_PREFIX + it + '_rotated.csv.gz'
//      return [it, file(synp)]
//    }
//    .combine(channel.fromList(STYPE))
//    .combine(channel.fromList(DEN_ALGO))
//    
//  def TS_BIN = [
//    [1, file(params.tff)],
//    [2, file(params.camf)],
//    [3, file(params.tsf)]
//  ]
//  def TS_LABEL = [
//    [1, 'P15_TF'],
//    [2, 'P15_CAM'],
//    [3, "Selector"]
//  ]
//  g_ch = channel.fromList(TS_BIN)
//    .join(channel.fromList(TS_LABEL))
//    .map { v -> v[1..2] }
//  sel_ch = VisualizeSelector(
//    scond_ch.combine(g_ch),
//    file(params.annf),
//    file(params.metaf),
//    utils_file,
//    SUBSAMPLE_TO,
//    SPARSE_LIMIT
//  )
//
//  tree_ch = SimilarityTree(
//    file(params.distances),
//    file(params.annf),
//    utils_file
//  )
//
//  extraction_ch = PartnerExtraction(
//    file(MAT_PREFIX + 'ME_L_rotated.csv.gz'),
//    file(MAT_PREFIX + 'ME_R_rotated.csv.gz'),
//    file(MAT_PREFIX + 'LO_L_rotated.csv.gz'),
//    file(MAT_PREFIX + 'LO_R_rotated.csv.gz'),
//    file(MAT_PREFIX + 'LOP_L_rotated.csv.gz'),
//    file(MAT_PREFIX + 'LOP_R_rotated.csv.gz'),
//    utils_file
//  )
//
//  partner_ch = NeuropilPartnerAnalysis(
//    extraction_ch,
//    utils_file,
//    channel
//      .fromList( ["All", "ME", "LO", "LOP"] )
//      .combine(
//        channel.fromList( ["pre", "post", "both"] )
//      ),
//    0.03
//  )

  publish:
  annov = out_ch
//  selv = sel_ch
//  tree = tree_ch
//  extraction = extraction_ch
//  partners = partner_ch
}

output {
  annov {
    path { input -> 
      def subset_dict = [
        'known': 'Previously_known',
        'new': 'All_confidently_annotated_highlighting_new',
        'all': 'All_confidently_annotated',
        'spatial': 'Spatial_origin',
        'putative': 'Putative_OPC_types'
      ]
      def color_dict = [
        'temporal': 'Color_by_temporal_origin',
        'subsystem': 'Color_by_function_subsystem',
        'type': 'Color_by_cell_type',
        'broad': 'Color_by_early_late',
        'spatial': 'Color_by_spatial_origin'
      ]
      def side = input[0].split('_')[1]
      def ctype = input[1].split('_')[0]
      def subset = input[1].split('_')[1]
      def clabel = color_dict[ctype]
      def slabel = ctype == 'spatial' ? 'Spatial_origin' : subset_dict[subset]
      return "${clabel}/${slabel}_${side}"
    }
  }
//  selv {
//    path { input ->
//      def side = input[0].split('_')[1]
//      def gsel = input[1]
//      return "${gsel}/${side}"
//    }
//  }
//  tree {
//    path "similarity_trees/"
//  }
//  extraction {
//    path "partner_data/"
//  }
//  partners {
//    path { input ->
//      def neuropil = input[0]
//      return "neuropil_partners/${neuropil}/"
//    }
//  }
}
