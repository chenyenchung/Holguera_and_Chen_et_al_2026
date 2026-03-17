#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(R.utils))

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# Interactive mode for development
if (interactive()) {
  argvs$partner_summary <- 'lc_lplc_presynaptic_partners.csv'
  argvs$anno <- 'data/visual_neurons_anno.csv'
  argvs$output_dir <- 'plots/lc_lplc_heatmaps'
} else {
  if (is.null(argvs$output_dir)) {
    argvs$output_dir <- 'plots/lc_lplc_heatmaps'
    cat("No output directory specified, using default: plots/lc_lplc_heatmaps\n")
  }
}

# Validate input files
if (is.null(argvs$partner_summary) || !file.exists(argvs$partner_summary)) {
  stop("Partner summary file not found: ", argvs$partner_summary)
}

if (is.null(argvs$anno) || !file.exists(argvs$anno)) {
  stop("Annotation file not found: ", argvs$anno)
}

# Create output directory if it doesn't exist
if (!dir.exists(argvs$output_dir)) {
  dir.create(argvs$output_dir, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", argvs$output_dir))
}

cat("========================================\n")
cat("LC/LPLC Heatmap Visualization\n")
cat("========================================\n\n")

# Step 1: Load data
cat("Step 1: Loading data...\n")
partner_summary <- fread(argvs$partner_summary)
cat(sprintf("  Loaded partner summary: %d rows\n", nrow(partner_summary)))

anno <- fread(argvs$anno)
cat(sprintf("  Loaded annotations: %d cell types\n\n", nrow(anno)))

# Step 2: Filter for neurons with temporal labels
cat("Step 2: Filtering for neurons with temporal labels...\n")
# Get cell types that have temporal labels (not "unknown")
temporal_annotated <- anno[temporal_label != "unknown", cell_type]
cat(sprintf("  Found %d cell types with temporal labels\n", length(temporal_annotated)))

# Filter partner summary to only include these neurons
filtered_summary <- partner_summary[pre_type %in% temporal_annotated]
cat(sprintf("  Filtered to %d connections (from %d)\n",
            nrow(filtered_summary), nrow(partner_summary)))
cat(sprintf("  Unique OPC types: %d\n", uniqueN(filtered_summary$pre_type)))
cat(sprintf("  Unique LC/LPLC types: %d\n\n", uniqueN(filtered_summary$post_type)))

# Step 3: Calculate neuron ratio
cat("Step 3: Calculating neuron ratios...\n")
filtered_summary[, neuron_ratio := connected_neuron_count / total_neuron_count]
cat(sprintf("  Neuron ratio range: %.4f to %.4f\n\n",
            min(filtered_summary$neuron_ratio), max(filtered_summary$neuron_ratio)))

# Step 4: Create matrix for heatmap
cat("Step 4: Creating data matrix...\n")
# Create wide format matrix: rows = LC/LPLC, columns = OPC types
heatmap_matrix <- dcast(
  filtered_summary,
  post_type ~ pre_type,
  value.var = "neuron_ratio",
  fill = 0
)

# Convert to matrix format
row_names <- heatmap_matrix$post_type
heatmap_matrix <- as.matrix(heatmap_matrix[, -1])
rownames(heatmap_matrix) <- row_names

cat(sprintf("  Matrix dimensions: %d rows (LC/LPLC) × %d columns (OPC types)\n\n",
            nrow(heatmap_matrix), ncol(heatmap_matrix)))

# Step 5: Get temporal annotation for ordering
cat("Step 5: Preparing temporal annotation...\n")
# Get temporal labels for each OPC type
opc_types <- colnames(heatmap_matrix)
temporal_info <- anno[cell_type %in% opc_types, .(cell_type, temporal_label, temporal_id)]
setkey(temporal_info, cell_type)

# Define temporal label ordering (following developmental progression)
temporal_levels <- c(
  "Hth",
  "Between Hth & Hth/Opa",
  "Hth/Opa",
  "Opa/Erm",
  "Erm/Ey",
  "Ey/Hbn",
  "Hbn/Opa/Slp",
  "Slp/D",
  "D/BH-1"
)

# Convert temporal_label to ordered factor
temporal_info[, temporal_label := factor(temporal_label, levels = temporal_levels, ordered = TRUE)]

# Create annotation data frame for pheatmap
annotation_col <- data.frame(
  temporal_label = temporal_info[opc_types]$temporal_label,
  row.names = opc_types
)

cat(sprintf("  Temporal labels found:\n"))
print(table(annotation_col$temporal_label))
cat("\n")

# Step 6: Generate heatmap 1 - hierarchically clustered
cat("Step 6: Generating hierarchically clustered heatmap...\n")
output_file_1 <- file.path(argvs$output_dir, "lc_lplc_heatmap_clustered.pdf")

pdf(output_file_1, width = 12, height = 10)
pheatmap(
  heatmap_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = annotation_col,
  color = colorRampPalette(c("white", "lightblue", "blue", "darkblue"))(100),
  main = "LC/LPLC Presynaptic Partner Ratios (Hierarchically Clustered)",
  fontsize = 8,
  fontsize_row = 7,
  fontsize_col = 6
)
dev.off()
cat(sprintf("  Saved: %s\n\n", output_file_1))

# Step 7: Generate heatmap 2 - ordered by temporal labels
cat("Step 7: Generating temporal-ordered heatmap...\n")
output_file_2 <- file.path(argvs$output_dir, "lc_lplc_heatmap_temporal.pdf")

# Sort OPC types by temporal label (using ordered factor) then by cell type
temporal_info_sorted <- temporal_info[opc_types]
setorder(temporal_info_sorted, temporal_label, cell_type)
temporal_order <- temporal_info_sorted$cell_type

# Reorder matrix columns
heatmap_matrix_temporal <- heatmap_matrix[, temporal_order, drop = FALSE]

# Reorder annotation to match, preserving factor levels
annotation_col_temporal <- data.frame(
  temporal_label = temporal_info_sorted$temporal_label,
  row.names = temporal_order
)

pdf(output_file_2, width = 12, height = 10)
pheatmap(
  heatmap_matrix_temporal,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = annotation_col_temporal,
  color = colorRampPalette(c("white", "lightblue", "blue", "darkblue"))(100),
  main = "LC/LPLC Presynaptic Partner Ratios (Ordered by Temporal Label)",
  fontsize = 8,
  fontsize_row = 7,
  fontsize_col = 6
)
dev.off()
cat(sprintf("  Saved: %s\n\n", output_file_2))

# Summary
cat("========================================\n")
cat("Visualization complete!\n")
cat("========================================\n")
cat(sprintf("Output directory: %s\n", argvs$output_dir))
cat(sprintf("Files generated:\n"))
cat(sprintf("  1. %s\n", basename(output_file_1)))
cat(sprintf("  2. %s\n", basename(output_file_2)))
cat(sprintf("\nMatrix summary:\n"))
cat(sprintf("  LC/LPLC types: %d\n", nrow(heatmap_matrix)))
cat(sprintf("  OPC types: %d\n", ncol(heatmap_matrix)))
cat(sprintf("  Non-zero values: %d (%.1f%%)\n",
            sum(heatmap_matrix > 0),
            100 * sum(heatmap_matrix > 0) / length(heatmap_matrix)))
