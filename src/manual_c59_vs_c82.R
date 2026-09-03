#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) == 0L) {
  stop("Could not determine the script path from commandArgs().")
}
script_file <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)

renv::load(repo_root)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Seurat))

usage <- function() {
  cat(paste0(
    "Usage: Rscript src/manual_c59_vs_c82.R [options]\n\n",
    "Options:\n",
    "  --adult PATH          Adult Seurat object [data/ozel_2021_objs/Adult.rds]\n",
    "  --modeled PATH        Adult ScMarco modeled-expression RDS\n",
    "  --output_dir PATH     Output directory [int/c59_vs_c82]\n",
    "  --cluster_col NAME    Numeric cluster metadata column [FinalIdents]\n",
    "  --threshold NUMBER    Strict ScMarco positive-call threshold [0.5]\n",
    "  --scale_factor NUMBER LogNormalize scale factor [10000]\n",
    "  --help                Show this message\n"
  ))
}

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--help") {
      usage()
      quit(save = "no", status = 0L)
    }
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key)
    }
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      stop("Missing value for argument: ", key)
    }
    out[[sub("^--", "", key)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

resolve_repo_path <- function(path) {
  if (grepl("^/", path)) {
    return(path)
  }
  file.path(repo_root, path)
}

get_assay_layer <- function(object, assay, layer) {
  tryCatch(
    GetAssayData(object, assay = assay, layer = layer),
    error = function(layer_error) {
      tryCatch(
        suppressWarnings(GetAssayData(object, assay = assay, slot = layer)),
        error = function(slot_error) {
          stop(
            sprintf(
              "Could not read the %s layer from assay %s. Layer error: %s; slot error: %s",
              layer, assay, conditionMessage(layer_error), conditionMessage(slot_error)
            )
          )
        }
      )
    }
  )
}

argvs <- parse_args(commandArgs(trailingOnly = TRUE))
argvs$adult <- resolve_repo_path(if (is.null(argvs$adult)) {
  "data/ozel_2021_objs/Adult.rds"
} else {
  argvs$adult
})
argvs$modeled <- if (is.null(argvs$modeled)) {
  paste0(
    "/projects/rps/cgsb/desplan/File_exchange/",
    "Ozel_and_Simon_2021_related/Mixture_Modeling/Adult_modeled_expression.RDS"
  )
} else {
  resolve_repo_path(argvs$modeled)
}
argvs$output_dir <- resolve_repo_path(if (is.null(argvs$output_dir)) {
  "int/c59_vs_c82"
} else {
  argvs$output_dir
})
argvs$cluster_col <- if (is.null(argvs$cluster_col)) "FinalIdents" else argvs$cluster_col
argvs$threshold <- if (is.null(argvs$threshold)) 0.5 else as.numeric(argvs$threshold)
argvs$scale_factor <- if (is.null(argvs$scale_factor)) 10000 else as.numeric(argvs$scale_factor)

if (!is.finite(argvs$threshold)) {
  stop("--threshold must be a finite number.")
}
if (!is.finite(argvs$scale_factor) || argvs$scale_factor <= 0) {
  stop("--scale_factor must be a positive finite number.")
}
if (!file.exists(argvs$adult)) {
  stop("Adult Seurat object not found: ", argvs$adult)
}
if (!file.exists(argvs$modeled)) {
  stop("Adult ScMarco modeled-expression object not found: ", argvs$modeled)
}
dir.create(argvs$output_dir, recursive = TRUE, showWarnings = FALSE)

genes <- c("toy", "ham", "D", "Vsx1")

# The 103 adult main-OPC clusters used for the manuscript analysis. These are
# numeric FinalIdents. Cluster 36 is excluded because its main-OPC origin was
# unresolved; its exclusion does not affect the Toy-or-Ham-positive subset.
main_opc_ids <- c(
  1, 9, 11, 12, 14, 15, 17, 18, 19, 20, 21, 22, 24, 26, 29, 30, 31, 32,
  33, 34, 37, 38, 39, 40, 41, 42, 43, 45, 46, 47, 52, 53, 54, 55, 56, 59,
  61, 62, 67, 69, 70, 71, 75, 76, 80, 81, 82, 83, 84, 94, 95, 103, 104,
  105, 115, 116, 118, 119, 120, 123, 124, 125, 126, 128, 129, 136, 137,
  138, 140, 142, 144, 145, 151, 152, 155, 156, 158, 159, 160, 161, 162,
  163, 164, 166, 167, 168, 172, 173, 174, 175, 176, 177, 178, 180, 181,
  183, 225, 226, 227, 228, 229, 230, 231
)
main_opc_ids <- as.character(main_opc_ids)
if (length(main_opc_ids) != 103L || anyDuplicated(main_opc_ids)) {
  stop("The embedded main-OPC cluster whitelist must contain 103 unique IDs.")
}

expected_plot_ids <- as.character(c(
  15, 30, 31, 45, 46, 47, 52, 53, 59, 61, 69, 71, 76, 81, 82, 83, 94,
  95, 103, 104, 105, 120, 124, 128, 129, 137, 156, 161, 167, 168, 173,
  177, 181, 183, 225, 226, 231
))

message("Loading adult ScMarco modeled expression: ", argvs$modeled)
modeled <- readRDS(argvs$modeled)
if (!is.matrix(modeled) && !is.data.frame(modeled)) {
  stop("The ScMarco object must be a matrix or data.frame.")
}
missing_genes <- setdiff(genes, rownames(modeled))
if (length(missing_genes) > 0L) {
  stop("Genes missing from the ScMarco object: ", paste(missing_genes, collapse = ", "))
}

numeric_column_index <- grep("^[0-9]+$", colnames(modeled))
numeric_columns <- colnames(modeled)[numeric_column_index]
if (length(numeric_column_index) == 0L) {
  stop("No numeric cluster columns were found in the ScMarco object.")
}

# Some source files contain repeated numeric cluster names. Collapse repeats by
# taking the maximum modeled score, equivalent to an 'any positive' rule.
scmarco_long <- rbindlist(lapply(genes, function(gene) {
  score <- suppressWarnings(as.numeric(
    unlist(modeled[gene, numeric_column_index, drop = FALSE], use.names = FALSE)
  ))
  data.table(gene = gene, cluster = numeric_columns, score = score)
}))
if (anyNA(scmarco_long$score) || any(!is.finite(scmarco_long$score))) {
  stop("Non-numeric or non-finite ScMarco scores were found for the requested genes.")
}
scmarco_long <- scmarco_long[, .(score = max(score)), by = .(cluster, gene)]
scmarco_long[, positive := score > argvs$threshold]

missing_main_opc <- setdiff(main_opc_ids, unique(scmarco_long$cluster))
if (length(missing_main_opc) > 0L) {
  stop(
    "Main-OPC IDs missing from the ScMarco object: ",
    paste(missing_main_opc, collapse = ", ")
  )
}

scmarco_main <- scmarco_long[cluster %in% main_opc_ids]
score_wide <- dcast(scmarco_main, cluster ~ gene, value.var = "score")
positive_wide <- dcast(scmarco_main, cluster ~ gene, value.var = "positive")
setnames(score_wide, genes, paste0(genes, "_score"))
setnames(positive_wide, genes, paste0(genes, "_positive"))
scmarco_audit <- merge(score_wide, positive_wide, by = "cluster", sort = FALSE)
scmarco_audit[, cluster_number := as.integer(cluster)]
setorder(scmarco_audit, cluster_number)

plot_ids <- scmarco_audit[
  toy_positive | ham_positive,
  as.character(cluster)
]
plot_ids <- as.character(sort(as.integer(unique(plot_ids))))
if (!identical(plot_ids, expected_plot_ids)) {
  stop(
    paste0(
      "The adult main-OPC Toy-or-Ham-positive set differs from the manuscript set. ",
      "Observed: ", paste(plot_ids, collapse = ", "), "; expected: ",
      paste(expected_plot_ids, collapse = ", ")
    )
  )
}

signature_hits <- scmarco_audit[
  toy_positive & ham_positive & !D_positive & !Vsx1_positive,
  as.character(cluster)
]
if (!identical(signature_hits, "59")) {
  stop(
    "Expected cluster 59 to be the sole main-OPC toy+/ham+/D-/Vsx1- match; observed: ",
    paste(signature_hits, collapse = ", ")
  )
}

cluster_82_call <- scmarco_audit[cluster == "82"]
if (
  nrow(cluster_82_call) != 1L ||
    !cluster_82_call$toy_positive || cluster_82_call$ham_positive ||
    cluster_82_call$D_positive || cluster_82_call$Vsx1_positive
) {
  stop("Cluster 82 did not have the expected toy+/ham-/D-/Vsx1- ScMarco call.")
}

scmarco_audit[, `:=`(
  selected_for_plot = cluster %in% plot_ids,
  matches_toy_ham_D_Vsx1_signature = cluster %in% signature_hits
)]
setcolorder(
  scmarco_audit,
  c(
    "cluster", "cluster_number", "selected_for_plot",
    "matches_toy_ham_D_Vsx1_signature",
    paste0(genes, "_score"), paste0(genes, "_positive")
  )
)

audit_file <- file.path(argvs$output_dir, "adult_main_opc_scmarco_audit.csv")
fwrite(scmarco_audit, audit_file)

message("Loading adult Seurat object: ", argvs$adult)
adult <- readRDS(argvs$adult)
adult <- UpdateSeuratObject(adult)
assay <- "RNA"
if (!assay %in% Assays(adult)) {
  stop("RNA assay not found. Available assays: ", paste(Assays(adult), collapse = ", "))
}
if (!argvs$cluster_col %in% colnames(adult@meta.data)) {
  stop("Cluster metadata column not found: ", argvs$cluster_col)
}
missing_expression_genes <- setdiff(genes, rownames(adult[[assay]]))
if (length(missing_expression_genes) > 0L) {
  stop("Genes missing from the adult RNA assay: ", paste(missing_expression_genes, collapse = ", "))
}

adult_cluster_ids <- as.character(adult@meta.data[[argvs$cluster_col]])
missing_plot_ids <- setdiff(plot_ids, unique(adult_cluster_ids))
if (length(missing_plot_ids) > 0L) {
  stop("Selected clusters missing from the adult object: ", paste(missing_plot_ids, collapse = ", "))
}

cells_keep <- colnames(adult)[adult_cluster_ids %in% plot_ids]
adult <- subset(adult, cells = cells_keep)
DefaultAssay(adult) <- assay

counts <- get_assay_layer(adult, assay, "counts")
if (nrow(counts) == 0L || ncol(counts) == 0L) {
  stop("The adult RNA assay has no count matrix; explicit normalization cannot be verified.")
}
rm(counts)
invisible(gc())

message(
  "Normalizing selected adult cells with LogNormalize (scale.factor = ",
  format(argvs$scale_factor, scientific = FALSE), ")"
)
adult <- NormalizeData(
  adult,
  assay = assay,
  normalization.method = "LogNormalize",
  scale.factor = argvs$scale_factor,
  verbose = FALSE
)

normalized <- get_assay_layer(adult, assay, "data")[genes, , drop = FALSE]
normalized_values <- as.numeric(normalized)
if (any(!is.finite(normalized_values)) || any(normalized_values < 0)) {
  stop("Normalized expression contains non-finite or negative values.")
}

cell_clusters <- as.character(adult@meta.data[[argvs$cluster_col]])
expression_wide <- as.data.table(t(as.matrix(normalized)), keep.rownames = "cell")
expression_wide[, cluster := cell_clusters]
expression_long <- melt(
  expression_wide,
  id.vars = c("cell", "cluster"),
  measure.vars = genes,
  variable.name = "gene",
  value.name = "log_normalized_expression"
)
expression_long[, gene := factor(gene, levels = genes)]

plot_order <- c("59", "82", setdiff(plot_ids, c("59", "82")))
plot_labels <- paste0("c", plot_order)
expression_long[, plot_cluster := factor(
  paste0("c", cluster),
  levels = plot_labels
)]

expression_summary <- expression_long[, .(
  n_cells = .N,
  mean_log_normalized = mean(log_normalized_expression),
  median_log_normalized = median(log_normalized_expression),
  pct_detected = 100 * mean(log_normalized_expression > 0)
), by = .(cluster, gene)]
expression_summary[, cluster_number := as.integer(cluster)]
expression_summary[, plot_cluster := factor(paste0("c", cluster), levels = plot_labels)]
expression_summary[, plot_order_index := match(cluster, plot_order)]
setorder(expression_summary, plot_order_index, gene)

summary_file <- file.path(argvs$output_dir, "adult_main_opc_expression_summary.csv")
fwrite(
  expression_summary[, .(
    cluster, cluster_number, gene, n_cells,
    mean_log_normalized, median_log_normalized, pct_detected
  )],
  summary_file
)

cluster_colors <- setNames(rep("grey75", length(plot_labels)), plot_labels)
cluster_colors[["c59"]] <- "#0072B2"
cluster_colors[["c82"]] <- "#D55E00"

violin_plot <- ggplot(
  expression_long,
  aes(x = plot_cluster, y = log_normalized_expression, fill = plot_cluster)
) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.15) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    size = 0.7,
    stroke = 0.2,
    fill = "white"
  ) +
  facet_wrap(vars(gene), ncol = 1L, scales = "free_y") +
  scale_fill_manual(values = cluster_colors, guide = "none") +
  labs(
    x = "Adult main-OPC cluster",
    y = "Log-normalized expression",
    title = "Adult expression in ScMarco toy- or ham-positive main-OPC clusters",
    subtitle = "RNA counts re-normalized with LogNormalize; c59 and c82 highlighted"
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "italic"),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    plot.title.position = "plot"
  )

violin_file <- file.path(
  argvs$output_dir,
  "adult_main_opc_toy_or_ham_stacked_violin.pdf"
)
ggsave(violin_file, violin_plot, width = 18, height = 10, units = "in")

expression_summary[, dot_cluster := factor(
  paste0("c", cluster),
  levels = rev(plot_labels)
)]
dot_plot <- ggplot(
  expression_summary,
  aes(x = gene, y = dot_cluster, size = pct_detected, color = mean_log_normalized)
) +
  geom_point() +
  scale_size_continuous(
    name = "Cells detected (%)",
    limits = c(0, 100),
    range = c(0, 7)
  ) +
  scale_color_viridis_c(name = "Mean log-normalized\nexpression", option = "magma") +
  labs(
    x = NULL,
    y = "Adult main-OPC cluster",
    title = "Adult expression in ScMarco toy- or ham-positive main-OPC clusters",
    subtitle = "c59 and c82 are shown first in the shared cluster order"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(face = "italic"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
    plot.title.position = "plot"
  )

dotplot_file <- file.path(
  argvs$output_dir,
  "adult_main_opc_toy_or_ham_dotplot.pdf"
)
ggsave(dotplot_file, dot_plot, width = 8, height = 12, units = "in")

pair_stats <- rbindlist(lapply(genes, function(gene_name) {
  x59 <- expression_long[
    cluster == "59" & gene == gene_name,
    log_normalized_expression
  ]
  x82 <- expression_long[
    cluster == "82" & gene == gene_name,
    log_normalized_expression
  ]
  test <- wilcox.test(x59, x82, exact = FALSE)
  data.table(
    gene = gene_name,
    cluster_1 = 59L,
    cluster_2 = 82L,
    n_cells_59 = length(x59),
    n_cells_82 = length(x82),
    mean_59 = mean(x59),
    mean_82 = mean(x82),
    mean_difference_59_minus_82 = mean(x59) - mean(x82),
    median_59 = median(x59),
    median_82 = median(x82),
    pct_detected_59 = 100 * mean(x59 > 0),
    pct_detected_82 = 100 * mean(x82 > 0),
    wilcox_W = unname(test$statistic),
    p_value = test$p.value
  )
}))
pair_stats[, q_value := p.adjust(p_value, method = "BH")]

stats_file <- file.path(argvs$output_dir, "adult_c59_vs_c82_wilcox.csv")
fwrite(pair_stats, stats_file)

validation_lines <- c(
  "Adult c59 versus c82 expression validation",
  sprintf("Adult object: %s", argvs$adult),
  sprintf("Adult ScMarco object: %s", argvs$modeled),
  sprintf("Cluster metadata column: %s", argvs$cluster_col),
  "Assay: RNA",
  sprintf(
    "Normalization: LogNormalize from RNA counts, scale.factor=%s",
    format(argvs$scale_factor, scientific = FALSE)
  ),
  sprintf("ScMarco positive call: score > %s", argvs$threshold),
  sprintf("Main-OPC clusters audited: %d", length(main_opc_ids)),
  sprintf("Toy-or-Ham-positive main-OPC clusters plotted: %d", length(plot_ids)),
  paste0("Plotted cluster IDs: ", paste(plot_ids, collapse = ", ")),
  paste0(
    "toy+/ham+/D-/Vsx1- main-OPC signature matches: ",
    paste(signature_hits, collapse = ", ")
  ),
  "Cluster 59 ScMarco call: toy+, ham+, D-, Vsx1-",
  "Cluster 82 ScMarco call: toy+, ham-, D-, Vsx1-",
  "Conclusion: cluster 59 is the sole main-OPC match for the requested signature."
)
validation_file <- file.path(argvs$output_dir, "adult_c59_vs_c82_validation.txt")
writeLines(validation_lines, validation_file)

output_files <- c(
  violin_file, dotplot_file, audit_file, summary_file, stats_file, validation_file
)
missing_outputs <- output_files[!file.exists(output_files) | file.info(output_files)$size <= 0]
if (length(missing_outputs) > 0L) {
  stop("Missing or empty output files: ", paste(missing_outputs, collapse = ", "))
}

message("Analysis complete. Outputs written to: ", argvs$output_dir)
message("Validated signature match: cluster ", signature_hits)
