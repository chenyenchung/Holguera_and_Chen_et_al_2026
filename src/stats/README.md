# Statistical Analysis Module

### 1. Distribution Comparison Tests

To examine if the depths of synapses differ among neurons with different
developmental origins or serving different functions, we used three
complementary methods to robustly compare meaningful differences with the
enormous number of synapses.

- Quantile Wasserstein Distance with Bootstrap Confidence Intervals
  - Measures how different two synaptic distributions are overall
- Subsampled Kolmogorov-Smirnov Test
  - Takes random samples to prevent over-sensitivity from huge sample sizes
  - Resample to assess instability from sampling
- Deep-superficial innervation test
  - Measures if each group (developmental or functional) prefers the
  superficial or deep part of the neuropil
  - Using reference cell types to mark the superficial or deep layers
  - Bootstrap resamples all three groups (superficial reference, deep reference,
    and interest group) to assess sampling variability
  - Reference cell types are configurable via `data/reference_groups.csv`
  - Default reference groups:
    - ME: Dm types (superficial) vs Pm types (deep)
    - LOP: T4/T5a,b (superficial) vs T4/T5c,d (deep)
    - LO: Tm1,2,4 (superficial) vs Tm20,5a-f (deep)

### 2. Functional Enrichment Analysis
- Tests if birth time and function are linked with Fisher's exact test
- Corrects for multiple comparisons using FDR

### 3. OPC Synapse Ratio Analysis
- Computes the fraction of lobula (LO) synapses whose presynaptic or
  postsynaptic cell type is marked as `putative_OPC` in
  `data/visual_neurons_anno.csv`
- Reports presynaptic and postsynaptic ratios separately for `LO_L` and `LO_R`
- Writes `LO_opc_synapse_ratio.csv` with columns:
  `neuropil`, `side`, `syn_type`, `opc_synapse_count`,
  `total_synapse_count`, `opc_synapse_ratio`, `n_putative_opc_types`

## Configuration Parameters

Several parameters can be tuned in the analysis, and the defaults used in
this analysis are:

```yaml
# Bootstrap settings (Wasserstein distance)
n_bootstrap: 1000          # More iterations = more accurate but slower
conf_int: 95              # Confidence level (%)

# Quantile settings
n_quantiles: 1000         # Resolution of distribution comparison (updated default)

# K-S test settings
ks_subsample_size: 1000   # Points sampled per group
ks_n_iterations: 1000     # Number of subsampling iterations
ks_correction_method: "fdr" # Multiple testing correction method (fdr/bonferroni)

# Broad depth analysis settings
broad_depth_coefficient: 0.5    # Threshold coefficient for variable threshold
broad_depth_n_bootstrap: 1000    # Bootstrap iterations for broad depth analysis
broad_depth_conf_int: 95         # Confidence level for broad depth analysis

# Data filtering
sparse_limit: 100         # Minimum neurons per group for testing
```

## Reference Group Configuration

Reference cell types for broad depth analysis can be configured via
`data/reference_groups.csv`. This file specifies which cell types mark
superficial vs deep layers in each neuropil.

**File format:**
```csv
neuropil,ref_superficial,ref_deep,pattern_type,description
ME_L,^Dm,^Pm,regex,Medulla left: Dm types mark superficial layers
...
```

**Columns:**
- `neuropil`: Neuropil name (e.g., "ME_L", "LOP_L", "LO_L")
- `ref_superficial`: Pattern to match superficial reference cell types
- `ref_deep`: Pattern to match deep reference cell types
- `pattern_type`: Either "regex" (pattern matching) or "exact" (semicolon-separated list)
- `description`: Human-readable description of the reference groups

**Configuration in Nextflow:**
```bash
nextflow run src/stats/main.nf --ref_groupsf data/reference_groups.csv
```

If the configuration file is not found, the analysis falls back to hardcoded
reference groups with a deprecation warning.

## Technical Requirements

- **R packages**: data.table, ggplot2, Rcpp, openxlsx, dplyr, patchwork
- **C++ compiler**: For optimized statistical computations
- **Memory**: ~16GB RAM for full analysis with bootstrap iterations
- **Compute time**: 1-4 hours depending on data size and bootstrap iterations

## Changelog

### Version 2.0 (Current)

**Major Improvements:**

1. **Bootstrap Methodology Fix** (broad_depth.cpp)
   - Changed bootstrap approach to resample all three groups (superficial reference,
     deep reference, and interest group) instead of only the interest group
   - This properly accounts for sampling variability in reference group medians
   - Results in more conservative (wider) confidence intervals and correct p-values
   - **Note**: Results will differ from previous versions due to this methodological correction

2. **Configurable Reference Groups** (data/reference_groups.csv)
   - Reference cell types for depth layer classification are now configurable
   - Added `data/reference_groups.csv` configuration file
   - Supports both regex patterns and exact cell type lists
   - Falls back to hardcoded values if configuration file not found
   - New parameter: `--ref_groupsf` in main.nf

3. **Robust File Path Handling**
   - Improved C++ source file path resolution in depth_stats.r and broad_depth.r
   - Scripts now try multiple potential paths before failing
   - Better error messages when files are not found
   - More reliable in different execution environments

4. **Code Cleanup**
   - Removed legacy `test_statistic_*` column naming support
   - Simplified code in combine_results.r
   - Improved code maintainability

**Breaking Changes:**
- Statistical results from broad depth analysis will differ from previous versions
  due to corrected bootstrap methodology (more conservative p-values expected)
- Old column naming (`test_statistic_*`) is no longer supported

**Migration Guide:**
- If you have custom scripts reading depth stats results, ensure they use
  `wasserstein_*` column names (not `test_statistic_*`)
- Create `data/reference_groups.csv` with your reference group definitions
  (see Reference Group Configuration section above)
