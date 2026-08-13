# Visualization Module

## 1. Spatial Distribution Plots (`visualize.r`)

Location of synapses in 3D brain space

**Components**:
- **Scatter plots**: Each point is a synapse
- **Color coding**: By cell type, birth time, or function
- **Density profiles**: Show concentration along depth axis

**Density calculation modes**:
- `asis`: Overall density across all cells (shows general patterns)
- `pertype`: Separate density for each cell type (shows type-specific patterns)

### 2. Gene Expression Highlights (`v_selector.r`)

Synapses from cells expressing specific genes. P15 TF and CSM plots use the
P15 expression matrices; selector TF plots use stable consensus calls across
P15, P30, P40, P50, P70, and adult.

**Data source**: Özel et al. (2020) expression atlas

To run only selected spatial-distribution presets, pass a comma- or
semicolon-separated list. For example:

```bash
nextflow run src/visualize/main.nf \
  --visualization_presets spatial_notch \
  -resume
```

### 3. Cell Type Family Trees (`similarity_tree.r`)

Relationships between cell types based on connectivity

**Data source**: Connectivity dissimilarity matrix from Matsliah et al. (2024)

### 4. Partner Analysis (`partner_extraction.r` & `neuropil_partner_analysis.r`)

Donut charts showing main synaptic partners of neuronal types of interest

**Partner extraction**:
- Identifies strong connections (above threshold)
- Processes bilateral data (left/right brain)
- Outputs partner matrices for further analysis

## Dependencies
- **R packages**: ggplot2, ggrastr, dplyr, data.table, patchwork, cowplot, ggdendro
- **Data files**: Rotated matrices from preprocessing
- **Annotation files**: Cell type classifications, gene expression data
