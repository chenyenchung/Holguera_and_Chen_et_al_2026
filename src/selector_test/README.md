# Selector Depth Bootstrap Analysis

Tests terminal selector genes for superficial vs deep depth bias in ME_R, LO_R, and LOP_R neuropils.

## Overview

This workflow performs neuron-level bootstrap analysis to test whether synapses from neurons expressing specific terminal selector genes show preferential localization to superficial or deep layers within visual neuropils.

**Key features:**
- Tests all 95 genes in `data/selectors.csv`
- Analyzes ME_R, LO_R, and LOP_R neuropils
- Tests both presynaptic and postsynaptic compartments separately
- Uses neuron-level bootstrap (same algorithm as `src/stats/bin/broad_depth.r`)
- Outputs comprehensive statistics with FDR correction
- Generates synapse scatter plots with overall depth-density curves

## Usage

### Basic execution

```bash
cd src/selector_test
nextflow run main.nf
```

### Custom parameters

```bash
nextflow run main.nf \
  --sparse_limit 200 \
  --n_bootstrap 10000 \
  --coefficient 0.3
```

### Gene batching

Split genes across multiple parallel jobs for faster execution:

```bash
# Process 10 genes per batch (18 batches × 3 neuropils × 2 syn_types = 108 jobs)
nextflow run main.nf --genes_per_batch 10

# Process 50 genes per batch (4 batches × 3 neuropils × 2 syn_types = 24 jobs)
nextflow run main.nf --genes_per_batch 50
```

No batching (default, backward compatible):
```bash
nextflow run main.nf --genes_per_batch 0
# OR simply
nextflow run main.nf
```

### Available parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `selectorsf` | `data/selectors.csv` | Gene expression matrix (genes × clusters) |
| `annf` | `data/visual_neurons_anno.csv` | Cell type annotations |
| `metaf` | `data/viz_meta.csv` | Visualization metadata (z-axis config) |
| `ref_groupsf` | `data/reference_groups.csv` | Reference group definitions |
| `utilsf` | `src/utils.r` | Utility functions |
| `broad_depth_cppf` | `src/stats/bin/broad_depth.cpp` | C++ bootstrap code |
| `sparse_limit` | `100` | Minimum synapses per gene to test |
| `coefficient` | `0.5` | Threshold coefficient for delta |
| `n_bootstrap` | `1000` | Number of bootstrap iterations |
| `conf_int` | `95` | Confidence interval level (%) |
| `genes_per_batch` | `0` | Genes per batch (0 = no batching, processes all genes in single job) |
| `selector_viz_script` | `src/visualize/bin/v_selector.r` | Selector synapse visualization script |
| `subsample` | `10000` | Synapses sampled per visualization process |
| `density` | `asis` | Depth-density mode used for visualization |
| `use_axis_limits` | `true` | Apply configured neuropil axis limits |

## Algorithm

### Data flow

```
1. Load synapse coordinates for neuropil (ME_R, LO_R, or LOP_R)
2. Filter sparse cell types (< sparse_limit synapses)
3. Merge gene expression data via filter_geneexp()
   → Each synapse gets TRUE/FALSE for each gene
4. Load reference groups:
   - ME_R: superficial = Dm*, deep = Pm*
   - LO_R: superficial = Tm1/2/4, deep = Tm20/5a-f
   - LOP_R: superficial = T4/T5a-b, deep = T4/T5c-d
5. For each gene and each Notch stratum (On/Off):
   a. Filter synapses to expressing neurons only within that Notch stratum
   b. Check: sufficient synapses? sufficient neurons?
   c. Prepare neuron-level depth mapping for the expressing set
   d. Run C++ bootstrap vs shared canonical reference groups
   e. Classify as superficial/deep/neither
6. Combine all workflow batches and recompute FDR within each
   neuropil × synapse-type × Notch testing family
7. Output CSV and Excel reports with the complete statistics, direction, and
   expressing neuronal-type lists
```

### Reference groups

Reference groups define canonical superficial and deep layers:

**ME_R (Medulla):**
- Superficial: Dm1-Dm30 (regex: `^Dm`)
- Deep: Pm1-Pm9 (regex: `^Pm`)

**LO_R (Lobula):**
- Superficial: Tm1, Tm2, Tm4 (exact list)
- Deep: Tm20, Tm5a, Tm5b, Tm5c, Tm5d, Tm5e, Tm5f (exact list)

**LOP_R (Lobula plate):**
- Superficial: T4a, T4b, T5a, T5b (regex: `T(4|5)(a|b)$`)
- Deep: T4c, T4d, T5c, T5d (regex: `T(4|5)(c|d)$`)

Configured in `data/reference_groups.csv`.

### Bootstrap test

Uses neuron-level resampling (from `broad_depth.cpp`):
1. Calculate observed distance to superficial vs deep reference groups
2. Bootstrap resample neurons (not synapses) 1000 times
3. Compute confidence intervals and p-values
4. Direction classification:
   - `superficial`: distance_diff > delta_threshold AND positive
   - `deep`: distance_diff > delta_threshold AND negative
   - `neither`: distance_diff within threshold

## Output

### File structure

**Without batching** (`--genes_per_batch 0` or default):
```
selector_test/
├── ME_R_pre/
│   └── ME_R_pre_selector_depth.csv
├── ME_R_post/
│   └── ME_R_post_selector_depth.csv
├── LO_R_pre/
│   └── LO_R_pre_selector_depth.csv
├── LO_R_post/
│   └── LO_R_post_selector_depth.csv
├── LOP_R_pre/
│   └── LOP_R_pre_selector_depth.csv
├── LOP_R_post/
│   └── LOP_R_post_selector_depth.csv
├── combined_selector_depth.csv
├── selector_depth_results.xlsx          # Excel with ME_R, LO_R, and LOP_R sheets
├── selector_depth_summary.txt
└── visualization/
    ├── ME_R_pre/
    ├── ME_R_post/
    ├── LO_R_pre/
    ├── LO_R_post/
    ├── LOP_R_pre/
    └── LOP_R_post/
```

**With batching** (e.g., `--genes_per_batch 10`):
```
selector_test/
├── ME_R_pre/
│   ├── ME_R_pre_batch001_selector_depth.csv
│   ├── ME_R_pre_batch002_selector_depth.csv
│   └── ... (18 batches total)
├── ME_R_post/
│   ├── ME_R_post_batch001_selector_depth.csv
│   └── ...
├── LO_R_pre/
│   └── ...
├── LO_R_post/
│   └── ...
├── LOP_R_pre/
│   └── ...
├── LOP_R_post/
│   └── ...
├── combined_selector_depth.csv          # All batches combined
├── selector_depth_results.xlsx          # Excel with ME_R, LO_R, and LOP_R sheets
├── selector_depth_summary.txt
└── visualization/
    ├── ME_R_pre/
    ├── ME_R_post/
    ├── LO_R_pre/
    ├── LO_R_post/
    ├── LOP_R_pre/
    └── LOP_R_post/
```

### Combined CSV and Excel columns

The combined report preserves the full per-test schema. Provisional FDR values
from individual workflow batches are overwritten after collection, so changing
`genes_per_batch` does not change `p_value_fdr` or `significant_fdr`.

**Gene information:**
- `types_of_interest`: Gene symbol
- `skip_reason`: Why gene was skipped (NA if tested)
- `n_neurons_expressing`: Number of neurons expressing this gene
- `n_synapses_expressing`: Number of synapses from expressing neurons
- `n_types_expressing`: Number of cell types expressing
- `types_expressing`: Semicolon-separated list of expressing types

**Reference conflict detection:**
- `conflict_with_references`: TRUE if gene expressed in both superficial and deep reference types
- `overlap_superficial`: Number of expressing types in superficial reference
- `overlap_deep`: Number of expressing types in deep reference
- `overlap_any`: TRUE if gene expressed in either superficial or deep reference types

**Reference groups:**
- `ref_superficial`: Semicolon-separated superficial types
- `ref_deep`: Semicolon-separated deep types

**Observed statistics:**
- `observed_sup_d`: Distance to superficial reference
- `observed_deep_d`: Distance to deep reference
- `observed_distance_diff`: sup_d - deep_d
- `observed_delta_thres_base`: Base threshold
- `observed_delta_thres`: Threshold for significance
- `observed_test_statistic`: Test statistic
- `observed_bias_ratio`: Bias ratio

**Bootstrap confidence intervals:**
- `bootstrap_distance_diff_{median,lower,upper}`
- `bootstrap_delta_thres_base_{median,lower,upper}`
- `bootstrap_test_stat_{median,lower,upper}`
- `bootstrap_bias_ratio_{median,lower,upper}`

**Statistical inference:**
- `p_value_exceeds_threshold`: P-value for exceeding threshold
- `p_value_bias_ratio`: P-value for bias ratio
- `p_value_fdr`: FDR-corrected p-value
- `significant_fdr`: TRUE if FDR < 0.05

**Classification:**
- `direction`: "superficial", "deep", or "neither"

**Metadata:**
- `coefficient`, `n_bootstrap`, `conf_level`: Analysis parameters
- `neuropil`, `syn_type`: Analysis condition
- `batch_id`: Batch identifier ("000" for no batching, "001", "002", etc. for batched runs)
- `sparse_limit`: Minimum synapses threshold
- `analysis_date`: Date of analysis

## Edge cases handled

1. **Genes with no expression**: Skip with `skip_reason="no_expressing_neurons"`
2. **Genes with insufficient synapses** (< sparse_limit): Skip with `skip_reason="insufficient_synapses"`
3. **Genes with insufficient neurons** (< 3): Skip with `skip_reason="insufficient_neurons"`
4. **Genes expressed in reference groups**: Flag overlap diagnostics (`overlap_*`), still test against shared references
5. **Bootstrap failures**: Skip with `skip_reason="bootstrap_failed"`
6. **NA p-values**: Exclude from FDR correction, preserve NA in output

## Integration with existing workflows

### Reused components

**From `src/utils.r`:**
- `filter_type()`: Remove sparse cell types
- `filter_geneexp()`: Merge gene expression onto synapses
- `validate_config_file()`: Configuration validation

**From `src/stats/bin/broad_depth.r`:**
- `prepare_neuron_depth_mapping()`: Create neuron-level index structures
- FDR correction pattern
- Reference group loading from CSV

**From `src/stats/bin/broad_depth.cpp`:**
- `perform_broad_depth_bootstrap_neuron_level()`: C++ bootstrap computation

**From `data/reference_groups.csv`:**
- Reference superficial/deep groups for ME_R, LO_R, and LOP_R

### Comparison to existing workflows

| Workflow | Purpose | Grouping | Output |
|----------|---------|----------|--------|
| `v_selector.r` | Visualize gene expression | Per gene | Scatter + density plots |
| `broad_depth.r` | Test cell type depth | Per cell type (preset-based) | Bootstrap statistics CSV |
| **selector_depth.r** | **Test gene expression depth** | **Per gene** | **Bootstrap statistics CSV** |

## Example results interpretation

```r
library(data.table)

# Load combined results
results <- fread("selector_test/combined_selector_depth.csv")

# Significant genes in ME_R presynaptic
me_pre <- results[neuropil == "ME_R" & syn_type == "pre" & significant_fdr == TRUE]
me_pre[order(p_value_fdr), .(gene, direction, n_neurons_expressing, p_value_fdr)]

# Genes with reference conflicts
conflicts <- results[conflict_with_references == TRUE]
conflicts[, .(gene, neuropil, syn_type, overlap_superficial, overlap_deep, direction)]

# Superficial vs deep bias summary
results[!is.na(direction), .N, by = .(neuropil, syn_type, direction)]
```

## Performance recommendations

### Choosing batch size

| Batch Size | Total Jobs | Use Case | Expected Runtime |
|------------|-----------|----------|------------------|
| 0 (default) | 6 | Small scale, backward compatibility | ~9-12 hours |
| 5 | 210 | Fast iteration during development | ~1.5-3 hours |
| **10** | **108** | **Recommended** - good balance | **~3-4 hours** |
| 50 | 24 | Limited cluster resources | ~6-8 hours |
| 1 | 1044 | Maximum parallelization (requires 300+ nodes) | ~45-70 min |

**Note**: Runtimes assume sufficient cluster resources and parallel execution. Actual time depends on cluster load.

### Edge cases

- **Batch size > total genes**: Creates single batch (equivalent to no batching)
- **Uneven division**: Last batch contains remainder (e.g., 174 ÷ 10 = 17 batches of 10 + 1 batch of 4)
- **Batch size = 0**: No batching (all genes processed in single job per neuropil/syn_type)

## Technical notes

- **Memory**: 8GB per job (handles large synapse matrices)
- **Module**: Requires `r/gcc/4.5.0` for Rcpp compilation
- **Gene extraction**: Uses logical columns added by `filter_geneexp()`, excludes annotation columns
- **Result collection**: Automatic via `combine_selector_results.r` with Excel output

## Troubleshooting

**Error: "C++ source file not found"**
- Ensure `src/stats/bin/broad_depth.cpp` exists
- Check `--broad_depth_cppf` parameter

**Error: "No reference groups found"**
- Ensure `data/reference_groups.csv` has entries for ME_R, LO_R, and LOP_R
- Check neuropil names match exactly

**All genes skipped with "insufficient_synapses"**
- Try lowering `--sparse_limit` parameter
- Check if gene expression data was merged correctly

**Bootstrap fails for specific genes**
- Check gene has sufficient neurons (>= 3)
- Verify neuron depth data is valid (no NAs)

## References

- Gene expression data: Özel et al. (2021) terminal selector atlas
- Bootstrap algorithm: Adapted from `broad_depth.r` neuron-level resampling
- Reference groups: Based on known ME/LO/LOP layer markers
