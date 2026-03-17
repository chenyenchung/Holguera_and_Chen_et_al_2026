#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# Interactive mode for development
if (interactive()) {
  argvs$connections <- 'data/connections_princeton_no_threshold.csv.gz'
  argvs$celltypes <- 'data/consolidated_cell_types.csv.gz'
  argvs$output <- 'lc_lplc_presynaptic_partners.csv'
  argvs$threshold <- 5
} else {
  # Parse threshold as numeric
  if (!is.null(argvs$threshold)) {
    argvs$threshold <- as.numeric(argvs$threshold)
  } else {
    argvs$threshold <- 5
    cat("No threshold specified, using default: 5\n")
  }
}

# Validate input files
if (is.null(argvs$connections) || !file.exists(argvs$connections)) {
  stop("Connections file not found: ", argvs$connections)
}

if (is.null(argvs$celltypes) || !file.exists(argvs$celltypes)) {
  stop("Cell types file not found: ", argvs$celltypes)
}

if (is.null(argvs$output)) {
  argvs$output <- 'lc_lplc_presynaptic_partners.csv'
  cat("No output specified, using default: lc_lplc_presynaptic_partners.csv\n")
}

cat("========================================\n")
cat("LC/LPLC Presynaptic Partner Analysis\n")
cat("========================================\n\n")

# Step 1: Load cell type annotations
cat("Step 1: Loading cell type annotations...\n")
type_ann <- fread(argvs$celltypes)
type_ann <- type_ann[, .(root_id, primary_type)]
type_ann[, root_id := as.character(root_id)]
type_ann[, root_id := sub("^720575940", "", root_id)]

cat(sprintf("  Loaded %d neurons with type annotations\n", nrow(type_ann)))

# Identify LC/LPLC types
lc_lplc_types <- unique(type_ann[
  grepl("^LC[0-9]+", primary_type) | grepl("^LPLC[0-9]+", primary_type),
  primary_type
])

cat(sprintf("  Found %d LC/LPLC types:\n", length(lc_lplc_types)))
cat(sprintf("    %s\n", paste(sort(lc_lplc_types), collapse=", ")))

# Count total neurons per type
total_neurons <- type_ann[, .(total_neuron_count = .N), by = primary_type]
cat(sprintf("  Counted total neurons for %d cell types\n\n", nrow(total_neurons)))

# Step 2: Load connection data
cat("Step 2: Loading connection data...\n")
connections <- fread(argvs$connections)
cat(sprintf("  Loaded %d connection records\n", nrow(connections)))

connections[, pre_root_id := as.character(pre_root_id)]
connections[, post_root_id := as.character(post_root_id)]
connections[, pre_root_id := sub("^720575940", "", pre_root_id)]
connections[, post_root_id := sub("^720575940", "", post_root_id)]

# Filter for optic lobe neuropils
ol_neuropils <- c("ME_L", "ME_R", "LO_L", "LO_R", "LOP_L", "LOP_R")
connections <- connections[neuropil %in% ol_neuropils]
cat(sprintf("  After filtering for optic lobe neuropils: %d connections\n", nrow(connections)))

# Apply synapse threshold
connections <- connections[syn_count > argvs$threshold]
cat(sprintf("  After filtering for > %d synapses: %d connections\n\n",
            argvs$threshold, nrow(connections)))

# Step 3: Annotate connections with cell types
cat("Step 3: Annotating connections with cell types...\n")

# Annotate presynaptic types
connections <- merge(
  connections,
  type_ann,
  by.x = "pre_root_id",
  by.y = "root_id",
  all.x = TRUE
)
setnames(connections, "primary_type", "pre_type")

# Annotate postsynaptic types
connections <- merge(
  connections,
  type_ann,
  by.x = "post_root_id",
  by.y = "root_id",
  all.x = TRUE
)
setnames(connections, "primary_type", "post_type")

# Warn about missing annotations
missing_pre <- sum(is.na(connections$pre_type))
missing_post <- sum(is.na(connections$post_type))
if (missing_pre > 0 || missing_post > 0) {
  warning(sprintf("  Missing type annotations: %d pre, %d post",
                  missing_pre, missing_post))
}

# Filter to LC/LPLC postsynaptic only
connections <- connections[post_type %in% lc_lplc_types]

# Filter out connections where pre_type is missing/empty
connections <- connections[!is.na(pre_type) & pre_type != ""]
cat(sprintf("  Connections to LC/LPLC neurons (with valid pre_type): %d\n\n", nrow(connections)))

# Step 4: Aggregate and count
cat("Step 4: Aggregating connections and counting neurons...\n")

# Count unique presynaptic neurons connecting to each LC/LPLC type
presynaptic_summary <- connections[
  , .(connected_neuron_count = uniqueN(pre_root_id)),
  by = .(pre_type, post_type)
]

cat(sprintf("  Found %d unique (pre_type, post_type) pairs\n", nrow(presynaptic_summary)))

# Step 5: Merge with total neuron counts
cat("Step 5: Merging with total neuron counts...\n")

presynaptic_summary <- merge(
  presynaptic_summary,
  total_neurons,
  by.x = "pre_type",
  by.y = "primary_type",
  all.x = TRUE
)

# Reorder columns to match specification
presynaptic_summary <- presynaptic_summary[
  , .(pre_type, total_neuron_count, post_type, connected_neuron_count)
]

# Sort for readability: by post_type, then by connected_neuron_count (descending)
setorder(presynaptic_summary, post_type, -connected_neuron_count)

cat(sprintf("  Final dataset contains %d rows\n\n", nrow(presynaptic_summary)))

# Validation checks
cat("Validation checks:\n")
cat(sprintf("  LC/LPLC types analyzed: %d (%s)\n",
            uniqueN(presynaptic_summary$post_type),
            paste(sort(unique(presynaptic_summary$post_type)), collapse=", ")))
cat(sprintf("  Unique presynaptic types: %d\n",
            uniqueN(presynaptic_summary$pre_type)))

# Check for logical consistency
invalid_rows <- presynaptic_summary[connected_neuron_count > total_neuron_count]
if (nrow(invalid_rows) > 0) {
  warning(sprintf("  WARNING: Found %d rows where connected > total!", nrow(invalid_rows)))
  print(invalid_rows)
} else {
  cat("  All rows pass consistency check (connected <= total)\n")
}

# Check for NA values
na_rows <- presynaptic_summary[is.na(pre_type) | is.na(post_type) |
                                is.na(total_neuron_count) | is.na(connected_neuron_count)]
if (nrow(na_rows) > 0) {
  warning(sprintf("  WARNING: Found %d rows with NA values!", nrow(na_rows)))
} else {
  cat("  No NA values in output columns\n")
}

# Step 6: Write output
cat("\nStep 6: Writing output...\n")
fwrite(presynaptic_summary, argvs$output, row.names = FALSE)

cat(sprintf("\n========================================\n"))
cat(sprintf("Analysis complete!\n"))
cat(sprintf("========================================\n"))
cat(sprintf("Output written to: %s\n", argvs$output))
cat(sprintf("Total rows: %d\n", nrow(presynaptic_summary)))
cat(sprintf("Preview of first 10 rows:\n"))
print(head(presynaptic_summary, 10))
