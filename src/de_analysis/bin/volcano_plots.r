#!/usr/bin/env Rscript
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo_root <- normalizePath(file.path(dirname(script_file), "../../.."), mustWork = TRUE)
renv::load(repo_root)
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

markers_file <- if (is.null(argvs$markers)) {
  file.path(repo_root, "int/de_analysis/combined_de_markers.csv")
} else {
  argvs$markers
}
out_dir <- if (is.null(argvs$out_dir)) {
  file.path(repo_root, "int/de_analysis/volcano_plots")
} else {
  argvs$out_dir
}
top_n <- if (is.null(argvs$top_n)) 10L else as.integer(argvs$top_n)
q_threshold <- if (is.null(argvs$q_threshold)) 0.05 else as.numeric(argvs$q_threshold)

required_cols <- c(
  "gene", "avg_log2FC", "p_val_adj", "stage", "neuropil", "syn_type",
  "contrast", "group_1", "group_2"
)

markers <- fread(markers_file)
missing_cols <- setdiff(required_cols, names(markers))
if (length(missing_cols) > 0) {
  stop("Missing marker columns: ", paste(missing_cols, collapse = ", "))
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(out_dir, "*_volcano.png"))

markers[, avg_log2FC := as.numeric(avg_log2FC)]
markers[, p_val_adj := as.numeric(p_val_adj)]
markers <- markers[!is.na(avg_log2FC) & !is.na(p_val_adj)]

positive_min_q <- min(markers[p_val_adj > 0, p_val_adj], na.rm = TRUE)
if (!is.finite(positive_min_q)) positive_min_q <- .Machine$double.xmin
markers[p_val_adj <= 0, p_val_adj := positive_min_q / 10]
markers[, neg_log10_q := -log10(p_val_adj)]

contrast_label <- function(x) {
  labels <- c(
    all = "All",
    notch_on_projection = "Notch On Projection",
    notch_off_projection = "Notch Off Projection",
    notch_on_intrinsic = "Notch On Intrinsic",
    notch_off_intrinsic = "Notch Off Intrinsic"
  )
  out <- labels[x]
  ifelse(is.na(out), gsub("_", " ", x, fixed = TRUE), out)
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

early_late_colors <- c(
  Early = "#0072B2",
  Late = "#D55E00"
)

depth_to_temporal <- function(neuropil, depth_group) {
  if (grepl("^ME", neuropil)) {
    return(ifelse(depth_group == "deep", "Early", "Late"))
  }
  ifelse(depth_group == "superficial", "Early", "Late")
}

plot_one <- function(dt, combo) {
  dt <- copy(dt)
  setorder(dt, p_val_adj)

  pos_labels <- dt[avg_log2FC > 0][order(p_val_adj, -avg_log2FC)][seq_len(min(.N, top_n))]
  neg_labels <- dt[avg_log2FC < 0][order(p_val_adj, avg_log2FC)][seq_len(min(.N, top_n))]
  labels <- unique(rbindlist(list(pos_labels, neg_labels), use.names = TRUE, fill = TRUE), by = "gene")

  group_1_temporal <- depth_to_temporal(combo$neuropil, combo$group_1)
  group_2_temporal <- depth_to_temporal(combo$neuropil, combo$group_2)

  dt[, marker_direction := fifelse(
    p_val_adj < q_threshold & avg_log2FC > 0, combo$group_1,
    fifelse(p_val_adj < q_threshold & avg_log2FC < 0, combo$group_2, "not significant")
  )]
  dt[, marker_direction := factor(
    marker_direction,
    levels = c(combo$group_2, "not significant", combo$group_1)
  )]
  direction_colors <- c(
    setNames(early_late_colors[group_2_temporal], combo$group_2),
    "not significant" = "grey75",
    setNames(early_late_colors[group_1_temporal], combo$group_1)
  )

  title <- sprintf(
    "%s %s %s: %s vs %s",
    combo$stage, combo$neuropil, combo$syn_type, combo$group_1, combo$group_2
  )
  subtitle <- sprintf(
    "%s contrast; labels: top %d positive and top %d negative markers by q-value",
    contrast_label(combo$contrast), top_n, top_n
  )

  p <- ggplot(dt, aes(x = avg_log2FC, y = neg_log10_q)) +
    geom_point(aes(color = marker_direction), alpha = 0.55, size = 0.8) +
    geom_hline(yintercept = -log10(q_threshold), linetype = "dashed", linewidth = 0.3, color = "grey45") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey55") +
    geom_text_repel(
      data = labels,
      aes(label = gene),
      size = 2.7,
      min.segment.length = 0,
      box.padding = 0.25,
      point.padding = 0.15,
      max.overlaps = Inf,
      seed = 1
    ) +
    scale_color_manual(
      values = direction_colors,
      drop = FALSE,
      name = NULL
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Average log2 fold change",
      y = "-log10(q-value)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )

  base <- paste(
    safe_name(combo$stage),
    safe_name(combo$neuropil),
    safe_name(combo$syn_type),
    safe_name(combo$contrast),
    "volcano",
    sep = "_"
  )
  pdf_file <- file.path(out_dir, paste0(base, ".pdf"))
  ggsave(pdf_file, p, width = 7, height = 5.25, device = cairo_pdf)

  data.table(
    stage = combo$stage,
    neuropil = combo$neuropil,
    syn_type = combo$syn_type,
    contrast = combo$contrast,
    group_1 = combo$group_1,
    group_2 = combo$group_2,
    group_1_temporal = group_1_temporal,
    group_2_temporal = group_2_temporal,
    group_1_color = unname(early_late_colors[group_1_temporal]),
    group_2_color = unname(early_late_colors[group_2_temporal]),
    rows = nrow(dt),
    positive_labels = sum(labels$avg_log2FC > 0),
    negative_labels = sum(labels$avg_log2FC < 0),
    pdf = pdf_file
  )
}

combo_cols <- c("stage", "neuropil", "syn_type", "contrast", "group_1", "group_2")
combos <- unique(markers[, ..combo_cols])
setorderv(combos, combo_cols)

plot_index <- rbindlist(lapply(seq_len(nrow(combos)), function(i) {
  combo <- combos[i]
  combo_dt <- markers[
    stage == combo$stage &
      neuropil == combo$neuropil &
      syn_type == combo$syn_type &
      contrast == combo$contrast &
      group_1 == combo$group_1 &
      group_2 == combo$group_2
  ]
  plot_one(combo_dt, combo)
}))

index_file <- file.path(out_dir, "volcano_plot_index.csv")
fwrite(plot_index, index_file)

cat("Wrote ", nrow(plot_index), " volcano PDF plots to ", out_dir, "\n", sep = "")
cat("Index: ", index_file, "\n", sep = "")
