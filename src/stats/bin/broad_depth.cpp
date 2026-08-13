#include <Rcpp.h>
#include <algorithm>
#include <random>
#include <cmath>
using namespace Rcpp;

// Fast median computation using partial_sort_copy (preserves input data)
inline double fast_median(const NumericVector& data, int start, int end) {
  int n = end - start;
  if (n == 0) return NA_REAL;
  
  int mid = n / 2;
  
  if (n % 2 == 1) {
    // Odd number: need (mid+1) smallest elements, return the middle one
    std::vector<double> temp(mid + 1);
    std::partial_sort_copy(data.begin() + start, data.begin() + end,
                          temp.begin(), temp.end());
    return temp[mid];
  } else {
    // Even number: need (mid+1) smallest elements, return average of last two
    std::vector<double> temp(mid + 1);
    std::partial_sort_copy(data.begin() + start, data.begin() + end,
                          temp.begin(), temp.end());
    return (temp[mid-1] + temp[mid]) / 2.0;
  }
}

// Efficient bootstrap sampling with pre-allocated memory
inline void bootstrap_sample_inplace(const NumericVector& source,
                                   NumericVector& target,
                                   int n,
                                   std::uniform_int_distribution<int>& dist,
                                   std::mt19937& gen) {
  for (int i = 0; i < n; ++i) {
    target[i] = source[dist(gen)];
  }
}

// Neuron-level bootstrap sampling: resample neurons and collect all their depths
// Returns the actual number of synapses in the bootstrap sample
inline int bootstrap_neurons_inplace(
    const NumericVector& source_depths,
    const IntegerVector& neuron_starts,
    const IntegerVector& neuron_counts,
    NumericVector& target,
    std::uniform_int_distribution<int>& neuron_dist,
    std::mt19937& gen) {

  int n_neurons = neuron_starts.size();
  int target_idx = 0;

  // Sample n_neurons with replacement
  for (int i = 0; i < n_neurons; ++i) {
    int sampled_neuron = neuron_dist(gen);  // Which neuron to sample
    int start = neuron_starts[sampled_neuron];
    int count = neuron_counts[sampled_neuron];

    // Copy all depths from this neuron
    for (int j = 0; j < count; ++j) {
      target[target_idx++] = source_depths[start + j];
    }
  }

  return target_idx;  // Return actual number of synapses written
}

//' Optimized bootstrap analysis for broad depth testing
//' 
//' @param ref_sup_depths Depth values for superficial reference group
//' @param ref_deep_depths Depth values for deep reference group
//' @param interest_depths Depth values for types of interest
//' @param coefficient Threshold coefficient (default 0.5)
//' @param n_bootstrap Number of bootstrap iterations (default 1000)
//' @param conf_int Confidence interval percentage (default 95)
//' @param seed Random seed for reproducibility (default 1)
//' @return List with bootstrap results
// [[Rcpp::export]]
List perform_broad_depth_bootstrap(NumericVector ref_sup_depths,
                                 NumericVector ref_deep_depths,
                                 NumericVector interest_depths,
                                 double coefficient = 0.5,
                                 int n_bootstrap = 1000,
                                 double conf_int = 95.0,
                                 Rcpp::Nullable<int> seed = R_NilValue) {
  
  // Input validation
  if (ref_sup_depths.size() == 0 || ref_deep_depths.size() == 0 || interest_depths.size() == 0) {
    stop("Empty input vectors provided");
  }
  
  if (conf_int <= 0 || conf_int >= 100) {
    stop("Confidence interval must be between 0 and 100");
  }
  
  // Remove NA values using Rcpp sugar (creates new vectors)
  NumericVector sup_data = na_omit(ref_sup_depths);
  NumericVector deep_data = na_omit(ref_deep_depths);
  NumericVector interest_data = na_omit(interest_depths);
  
  int n_sup = sup_data.size();
  int n_deep = deep_data.size();
  int n_interest = interest_data.size();
  
  // Pre-allocate work vectors for bootstrap sampling
  NumericVector work_sup(n_sup);
  NumericVector work_deep(n_deep);
  NumericVector work_interest(n_interest);
  
  // Copy data to work vectors for initial median calculation
  std::copy(sup_data.begin(), sup_data.end(), work_sup.begin());
  std::copy(deep_data.begin(), deep_data.end(), work_deep.begin());
  std::copy(interest_data.begin(), interest_data.end(), work_interest.begin());
  
  // Calculate observed statistics
  double observed_median_sup = fast_median(work_sup, 0, n_sup);
  double observed_median_deep = fast_median(work_deep, 0, n_deep);
  double observed_median_interest = fast_median(work_interest, 0, n_interest);
  
  double observed_sup_d = std::abs(observed_median_interest - observed_median_sup);
  double observed_deep_d = std::abs(observed_median_interest - observed_median_deep);
  double observed_distance_diff = observed_deep_d - observed_sup_d;
  double observed_delta_thres_base = std::abs(observed_median_deep - observed_median_sup);
  double observed_delta_thres = observed_delta_thres_base * coefficient;
  
  // Pre-allocate bootstrap result vectors
  NumericVector bootstrap_sup_d(n_bootstrap);
  NumericVector bootstrap_deep_d(n_bootstrap);
  NumericVector bootstrap_delta_thres_base(n_bootstrap);
  NumericVector bootstrap_distance_diff(n_bootstrap);
  NumericVector bootstrap_test_stats(n_bootstrap);
  
  // Calculate observed test statistic
  double observed_test_stat = std::abs(observed_distance_diff) - observed_delta_thres;
  
  // Set up a deterministic random number generator by default
  std::mt19937 gen;
  int effective_seed = 1;
  if (seed.isNotNull()) {
    effective_seed = as<int>(seed);
    if (effective_seed == NA_INTEGER) {
      stop("Random seed must not be NA");
    }
  }
  gen.seed(static_cast<std::mt19937::result_type>(effective_seed));
  Rcout << "Broad-depth bootstrap seed: " << effective_seed << "\n";

  // Set up distributions for bootstrap sampling of ALL groups
  std::uniform_int_distribution<int> dist_sup(0, n_sup - 1);
  std::uniform_int_distribution<int> dist_deep(0, n_deep - 1);
  std::uniform_int_distribution<int> dist_interest(0, n_interest - 1);

  // Counter for bootstrap test statistics > 0 (one-sided p-value)
  int count_positive = 0;

  // Bootstrap approach: Resample all three groups to assess sampling variability
  // Tests: "Under repeated sampling, how often does the interest group
  // exceed the threshold distance from references?"
  // P-value represents proportion of bootstrap samples where test_stat <= 0
  // (i.e., difference does not exceed threshold). Low p-values indicate
  // stable classification as superficial or deep.

  for (int i = 0; i < n_bootstrap; ++i) {
    // Resample all three groups with replacement
    for (int j = 0; j < n_sup; ++j) {
      work_sup[j] = sup_data[dist_sup(gen)];
    }
    for (int j = 0; j < n_deep; ++j) {
      work_deep[j] = deep_data[dist_deep(gen)];
    }
    for (int j = 0; j < n_interest; ++j) {
      work_interest[j] = interest_data[dist_interest(gen)];
    }

    // Calculate bootstrap medians for all three groups
    double boot_median_sup = fast_median(work_sup, 0, n_sup);
    double boot_median_deep = fast_median(work_deep, 0, n_deep);
    double boot_median_interest = fast_median(work_interest, 0, n_interest);

    // Calculate distances using bootstrapped reference medians
    double boot_sup_d = std::abs(boot_median_interest - boot_median_sup);
    double boot_deep_d = std::abs(boot_median_interest - boot_median_deep);
    double boot_distance_diff = boot_deep_d - boot_sup_d;

    // Calculate threshold from bootstrapped reference separation
    double boot_delta_thres_base = std::abs(boot_median_deep - boot_median_sup);
    double boot_delta_thres = boot_delta_thres_base * coefficient;

    // Store bootstrap values
    bootstrap_sup_d[i] = boot_sup_d;
    bootstrap_deep_d[i] = boot_deep_d;
    bootstrap_delta_thres_base[i] = boot_delta_thres_base;
    bootstrap_distance_diff[i] = boot_distance_diff;

    // Calculate bootstrap test statistic with bootstrapped threshold
    double boot_test_stat = std::abs(boot_distance_diff) - boot_delta_thres;
    bootstrap_test_stats[i] = boot_test_stat;

    // Count positive test statistics for one-sided p-value
    if (boot_test_stat > 0) {
      count_positive++;
    }
  }
  
  // Calculate one-sided p-value for test statistic > 0
  double p_value_exceeds_threshold = 1.0 - static_cast<double>(count_positive) / n_bootstrap;
  
  // Calculate confidence intervals
  // Clone vectors for sorting (to preserve original bootstrap results)
  NumericVector sorted_distance_diff = clone(bootstrap_distance_diff);
  NumericVector sorted_delta_thres = clone(bootstrap_delta_thres_base);
  NumericVector sorted_test_stats = clone(bootstrap_test_stats);
  
  std::sort(sorted_distance_diff.begin(), sorted_distance_diff.end());
  std::sort(sorted_delta_thres.begin(), sorted_delta_thres.end());
  std::sort(sorted_test_stats.begin(), sorted_test_stats.end());
  
  double alpha = (100.0 - conf_int) / 100.0;
  int lower_idx = static_cast<int>(n_bootstrap * alpha / 2.0);
  int upper_idx = static_cast<int>(n_bootstrap * (1.0 - alpha / 2.0));
  int median_idx = static_cast<int>(n_bootstrap * 0.5);
  
  // Ensure indices are within bounds
  lower_idx = std::max(0, std::min(lower_idx, n_bootstrap - 1));
  upper_idx = std::max(0, std::min(upper_idx, n_bootstrap - 1));
  median_idx = std::max(0, std::min(median_idx, n_bootstrap - 1));
  
  List result;
  result["observed_sup_d"] = observed_sup_d;
  result["observed_deep_d"] = observed_deep_d;
  result["observed_distance_diff"] = observed_distance_diff;
  result["observed_delta_thres_base"] = observed_delta_thres_base;
  result["observed_delta_thres"] = observed_delta_thres;
  result["observed_test_statistic"] = observed_test_stat;
  result["bootstrap_distance_diff_median"] = sorted_distance_diff[median_idx];
  result["bootstrap_distance_diff_lower"] = sorted_distance_diff[lower_idx];
  result["bootstrap_distance_diff_upper"] = sorted_distance_diff[upper_idx];
  result["bootstrap_delta_thres_base_median"] = sorted_delta_thres[median_idx];
  result["bootstrap_delta_thres_base_lower"] = sorted_delta_thres[lower_idx];
  result["bootstrap_delta_thres_base_upper"] = sorted_delta_thres[upper_idx];
  result["bootstrap_test_stat_median"] = sorted_test_stats[median_idx];
  result["bootstrap_test_stat_lower"] = sorted_test_stats[lower_idx];
  result["bootstrap_test_stat_upper"] = sorted_test_stats[upper_idx];
  result["p_value_exceeds_threshold"] = p_value_exceeds_threshold;
  result["bootstrap_sup_d"] = bootstrap_sup_d;
  result["bootstrap_deep_d"] = bootstrap_deep_d;
  result["bootstrap_delta_thres_base"] = bootstrap_delta_thres_base;
  result["bootstrap_distance_diff"] = bootstrap_distance_diff;
  result["bootstrap_test_stats"] = bootstrap_test_stats;

  return result;
}

//' Neuron-level bootstrap analysis for broad depth testing
//'
//' @param ref_sup_depths Depth values for superficial reference group
//' @param ref_sup_neuron_starts Starting indices for each neuron in superficial group
//' @param ref_sup_neuron_counts Number of synapses for each neuron in superficial group
//' @param ref_deep_depths Depth values for deep reference group
//' @param ref_deep_neuron_starts Starting indices for each neuron in deep group
//' @param ref_deep_neuron_counts Number of synapses for each neuron in deep group
//' @param interest_depths Depth values for types of interest
//' @param interest_neuron_starts Starting indices for each neuron in interest group
//' @param interest_neuron_counts Number of synapses for each neuron in interest group
//' @param coefficient Threshold coefficient (default 0.5)
//' @param n_bootstrap Number of bootstrap iterations (default 1000)
//' @param conf_int Confidence interval percentage (default 95)
//' @param seed Random seed for reproducibility (default 1)
//' @return List with bootstrap results including bias_ratio statistics
// [[Rcpp::export]]
List perform_broad_depth_bootstrap_neuron_level(
    NumericVector ref_sup_depths,
    IntegerVector ref_sup_neuron_starts,
    IntegerVector ref_sup_neuron_counts,

    NumericVector ref_deep_depths,
    IntegerVector ref_deep_neuron_starts,
    IntegerVector ref_deep_neuron_counts,

    NumericVector interest_depths,
    IntegerVector interest_neuron_starts,
    IntegerVector interest_neuron_counts,

    double coefficient = 0.5,
    int n_bootstrap = 1000,
    double conf_int = 95.0,
    Rcpp::Nullable<int> seed = R_NilValue) {

  // Input validation
  if (ref_sup_depths.size() == 0 || ref_deep_depths.size() == 0 ||
      interest_depths.size() == 0) {
    stop("Empty input vectors provided");
  }

  if (ref_sup_neuron_starts.size() != ref_sup_neuron_counts.size()) {
    stop("Mismatch between ref_sup neuron starts and counts");
  }

  if (ref_deep_neuron_starts.size() != ref_deep_neuron_counts.size()) {
    stop("Mismatch between ref_deep neuron starts and counts");
  }

  if (interest_neuron_starts.size() != interest_neuron_counts.size()) {
    stop("Mismatch between interest neuron starts and counts");
  }

  if (conf_int <= 0 || conf_int >= 100) {
    stop("Confidence interval must be between 0 and 100");
  }

  // Validate that neuron indices correctly span the depth vectors
  int sum_sup = std::accumulate(ref_sup_neuron_counts.begin(), ref_sup_neuron_counts.end(), 0);
  int sum_deep = std::accumulate(ref_deep_neuron_counts.begin(), ref_deep_neuron_counts.end(), 0);
  int sum_interest = std::accumulate(interest_neuron_counts.begin(), interest_neuron_counts.end(), 0);

  if (sum_sup != ref_sup_depths.size()) {
    stop("Sum of ref_sup_neuron_counts does not match ref_sup_depths size");
  }
  if (sum_deep != ref_deep_depths.size()) {
    stop("Sum of ref_deep_neuron_counts does not match ref_deep_depths size");
  }
  if (sum_interest != interest_depths.size()) {
    stop("Sum of interest_neuron_counts does not match interest_depths size");
  }

  // Get neuron counts
  int n_neurons_sup = ref_sup_neuron_starts.size();
  int n_neurons_deep = ref_deep_neuron_starts.size();
  int n_neurons_interest = interest_neuron_starts.size();

  // Get total synapse counts (for pre-allocation)
  int n_syn_sup = ref_sup_depths.size();
  int n_syn_deep = ref_deep_depths.size();
  int n_syn_interest = interest_depths.size();

  // Calculate maximum possible bootstrap sample size
  // (worst case: sample the largest neuron n_neurons times)
  int max_count_sup = *std::max_element(ref_sup_neuron_counts.begin(), ref_sup_neuron_counts.end());
  int max_count_deep = *std::max_element(ref_deep_neuron_counts.begin(), ref_deep_neuron_counts.end());
  int max_count_interest = *std::max_element(interest_neuron_counts.begin(), interest_neuron_counts.end());

  int max_possible_sup = max_count_sup * n_neurons_sup;
  int max_possible_deep = max_count_deep * n_neurons_deep;
  int max_possible_interest = max_count_interest * n_neurons_interest;

  // Pre-allocate work vectors with maximum possible size to avoid out-of-bounds
  NumericVector work_sup(max_possible_sup);
  NumericVector work_deep(max_possible_deep);
  NumericVector work_interest(max_possible_interest);

  // Calculate observed statistics on original data
  double observed_median_sup = fast_median(ref_sup_depths, 0, n_syn_sup);
  double observed_median_deep = fast_median(ref_deep_depths, 0, n_syn_deep);
  double observed_median_interest = fast_median(interest_depths, 0, n_syn_interest);

  double observed_sup_d = std::abs(observed_median_interest - observed_median_sup);
  double observed_deep_d = std::abs(observed_median_interest - observed_median_deep);
  double observed_distance_diff = observed_deep_d - observed_sup_d;
  double observed_delta_thres_base = std::abs(observed_median_deep - observed_median_sup);
  double observed_delta_thres = observed_delta_thres_base * coefficient;
  double observed_test_stat = std::abs(observed_distance_diff) - observed_delta_thres;

  // Calculate observed bias_ratio
  double observed_bias_ratio = NA_REAL;
  if (observed_delta_thres_base > 0) {
    observed_bias_ratio = observed_distance_diff / observed_delta_thres_base;
  }

  // Pre-allocate bootstrap result vectors
  NumericVector bootstrap_sup_d(n_bootstrap);
  NumericVector bootstrap_deep_d(n_bootstrap);
  NumericVector bootstrap_delta_thres_base(n_bootstrap);
  NumericVector bootstrap_distance_diff(n_bootstrap);
  NumericVector bootstrap_test_stats(n_bootstrap);
  NumericVector bootstrap_bias_ratio(n_bootstrap);

  // Set up a deterministic random number generator by default
  std::mt19937 gen;
  int effective_seed = 1;
  if (seed.isNotNull()) {
    effective_seed = as<int>(seed);
    if (effective_seed == NA_INTEGER) {
      stop("Random seed must not be NA");
    }
  }
  gen.seed(static_cast<std::mt19937::result_type>(effective_seed));
  Rcout << "Broad-depth bootstrap seed: " << effective_seed << "\n";

  // Set up distributions for neuron-level sampling
  std::uniform_int_distribution<int> dist_neurons_sup(0, n_neurons_sup - 1);
  std::uniform_int_distribution<int> dist_neurons_deep(0, n_neurons_deep - 1);
  std::uniform_int_distribution<int> dist_neurons_interest(0, n_neurons_interest - 1);

  int count_positive = 0;
  int count_bias_exceeds = 0;

  // Bootstrap loop: resample neurons (with all their synapses)
  for (int i = 0; i < n_bootstrap; ++i) {
    // Resample neurons for all three groups and get actual sizes
    int actual_size_sup = bootstrap_neurons_inplace(
      ref_sup_depths, ref_sup_neuron_starts, ref_sup_neuron_counts,
      work_sup, dist_neurons_sup, gen
    );

    int actual_size_deep = bootstrap_neurons_inplace(
      ref_deep_depths, ref_deep_neuron_starts, ref_deep_neuron_counts,
      work_deep, dist_neurons_deep, gen
    );

    int actual_size_interest = bootstrap_neurons_inplace(
      interest_depths, interest_neuron_starts, interest_neuron_counts,
      work_interest, dist_neurons_interest, gen
    );

    // Calculate bootstrap medians using actual sizes
    double boot_median_sup = fast_median(work_sup, 0, actual_size_sup);
    double boot_median_deep = fast_median(work_deep, 0, actual_size_deep);
    double boot_median_interest = fast_median(work_interest, 0, actual_size_interest);

    // Calculate distances
    double boot_sup_d = std::abs(boot_median_interest - boot_median_sup);
    double boot_deep_d = std::abs(boot_median_interest - boot_median_deep);
    double boot_distance_diff = boot_deep_d - boot_sup_d;

    // Calculate threshold
    double boot_delta_thres_base = std::abs(boot_median_deep - boot_median_sup);
    double boot_delta_thres = boot_delta_thres_base * coefficient;

    // Store bootstrap values
    bootstrap_sup_d[i] = boot_sup_d;
    bootstrap_deep_d[i] = boot_deep_d;
    bootstrap_delta_thres_base[i] = boot_delta_thres_base;
    bootstrap_distance_diff[i] = boot_distance_diff;

    // Calculate test statistic
    double boot_test_stat = std::abs(boot_distance_diff) - boot_delta_thres;
    bootstrap_test_stats[i] = boot_test_stat;

    // Calculate bias_ratio
    double boot_bias_ratio = NA_REAL;
    if (boot_delta_thres_base > 0) {
      boot_bias_ratio = boot_distance_diff / boot_delta_thres_base;
    }
    bootstrap_bias_ratio[i] = boot_bias_ratio;

    // Count positive test statistics for one-sided p-value
    if (boot_test_stat > 0) {
      count_positive++;
    }

    // Count bias_ratio exceeding coefficient for bias_ratio p-value
    if (!NumericVector::is_na(boot_bias_ratio) && std::abs(boot_bias_ratio) > coefficient) {
      count_bias_exceeds++;
    }
  }

  // Calculate one-sided p-value for test statistic > 0
  double p_value_exceeds_threshold = 1.0 - static_cast<double>(count_positive) / n_bootstrap;

  // Calculate p-value for |bias_ratio| > coefficient
  double p_value_bias_ratio = 1.0 - static_cast<double>(count_bias_exceeds) / n_bootstrap;

  // Calculate confidence intervals
  NumericVector sorted_distance_diff = clone(bootstrap_distance_diff);
  NumericVector sorted_delta_thres = clone(bootstrap_delta_thres_base);
  NumericVector sorted_test_stats = clone(bootstrap_test_stats);
  NumericVector sorted_bias_ratio = clone(bootstrap_bias_ratio);

  std::sort(sorted_distance_diff.begin(), sorted_distance_diff.end());
  std::sort(sorted_delta_thres.begin(), sorted_delta_thres.end());
  std::sort(sorted_test_stats.begin(), sorted_test_stats.end());
  std::sort(sorted_bias_ratio.begin(), sorted_bias_ratio.end());

  double alpha = (100.0 - conf_int) / 100.0;
  int lower_idx = static_cast<int>(n_bootstrap * alpha / 2.0);
  int upper_idx = static_cast<int>(n_bootstrap * (1.0 - alpha / 2.0));
  int median_idx = static_cast<int>(n_bootstrap * 0.5);

  // Ensure indices are within bounds
  lower_idx = std::max(0, std::min(lower_idx, n_bootstrap - 1));
  upper_idx = std::max(0, std::min(upper_idx, n_bootstrap - 1));
  median_idx = std::max(0, std::min(median_idx, n_bootstrap - 1));

  // Return results
  List result;
  result["observed_sup_d"] = observed_sup_d;
  result["observed_deep_d"] = observed_deep_d;
  result["observed_distance_diff"] = observed_distance_diff;
  result["observed_delta_thres_base"] = observed_delta_thres_base;
  result["observed_delta_thres"] = observed_delta_thres;
  result["observed_test_statistic"] = observed_test_stat;
  result["observed_bias_ratio"] = observed_bias_ratio;

  result["bootstrap_distance_diff_median"] = sorted_distance_diff[median_idx];
  result["bootstrap_distance_diff_lower"] = sorted_distance_diff[lower_idx];
  result["bootstrap_distance_diff_upper"] = sorted_distance_diff[upper_idx];
  result["bootstrap_delta_thres_base_median"] = sorted_delta_thres[median_idx];
  result["bootstrap_delta_thres_base_lower"] = sorted_delta_thres[lower_idx];
  result["bootstrap_delta_thres_base_upper"] = sorted_delta_thres[upper_idx];
  result["bootstrap_test_stat_median"] = sorted_test_stats[median_idx];
  result["bootstrap_test_stat_lower"] = sorted_test_stats[lower_idx];
  result["bootstrap_test_stat_upper"] = sorted_test_stats[upper_idx];

  result["bootstrap_bias_ratio_median"] = sorted_bias_ratio[median_idx];
  result["bootstrap_bias_ratio_lower"] = sorted_bias_ratio[lower_idx];
  result["bootstrap_bias_ratio_upper"] = sorted_bias_ratio[upper_idx];

  result["p_value_exceeds_threshold"] = p_value_exceeds_threshold;
  result["p_value_bias_ratio"] = p_value_bias_ratio;

  result["bootstrap_sup_d"] = bootstrap_sup_d;
  result["bootstrap_deep_d"] = bootstrap_deep_d;
  result["bootstrap_delta_thres_base"] = bootstrap_delta_thres_base;
  result["bootstrap_distance_diff"] = bootstrap_distance_diff;
  result["bootstrap_test_stats"] = bootstrap_test_stats;
  result["bootstrap_bias_ratio"] = bootstrap_bias_ratio;

  // Add diagnostic information
  result["n_neurons_sup"] = n_neurons_sup;
  result["n_neurons_deep"] = n_neurons_deep;
  result["n_neurons_interest"] = n_neurons_interest;
  result["n_syn_sup"] = n_syn_sup;
  result["n_syn_deep"] = n_syn_deep;
  result["n_syn_interest"] = n_syn_interest;

  return result;
}
