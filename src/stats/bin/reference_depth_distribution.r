#!/usr/bin/env Rscript
renv::load("/scratch/ycc520/flyem")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(ggplot2))

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

required_args <- c("np", "synf", "meta", "ref_groups", "output_prefix")
missing_args <- required_args[vapply(required_args, function(x) is.null(argvs[[x]]), logical(1))]
if (length(missing_args) > 0) {
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

validate_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s file not found: %s", label, path))
  }
}

validate_columns <- function(dt, cols, label) {
  missing_cols <- setdiff(cols, colnames(dt))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "%s is missing required columns: %s",
      label, paste(missing_cols, collapse = ", ")
    ))
  }
}

parse_type_list <- function(x) {
  if (is.null(x) || !nzchar(trimws(x))) {
    return(character(0))
  }
  parts <- trimws(unlist(strsplit(x, "[,;]")))
  unique(parts[nzchar(parts)])
}

expand_reference_groups <- function(ref_config, available_groups) {
  if (ref_config$pattern_type == "regex") {
    ref_superficial <- grep(ref_config$ref_superficial, available_groups, value = TRUE)
    ref_deep <- grep(ref_config$ref_deep, available_groups, value = TRUE)
  } else if (ref_config$pattern_type == "exact") {
    ref_superficial <- intersect(
      trimws(strsplit(ref_config$ref_superficial, ";")[[1]]),
      available_groups
    )
    ref_deep <- intersect(
      trimws(strsplit(ref_config$ref_deep, ";")[[1]]),
      available_groups
    )
  } else {
    stop("Unknown pattern_type in reference_groups.csv: ", ref_config$pattern_type)
  }

  overlap <- intersect(ref_superficial, ref_deep)
  if (length(overlap) > 0) {
    stop(sprintf(
      "Reference types appear in both superficial and deep groups for %s: %s",
      argvs$np, paste(sort(overlap), collapse = ", ")
    ))
  }

  list(superficial = sort(ref_superficial), deep = sort(ref_deep))
}

collect_depths <- function(coord, type_groups, zinvert) {
  rows <- rbindlist(lapply(names(type_groups), function(group_name) {
    types <- type_groups[[group_name]]
    if (length(types) == 0) {
      return(NULL)
    }

    pre_part <- coord[pre_type %in% types, .(
      neuropil = argvs$np,
      cell_type = pre_type,
      calibration_group = group_name,
      synapse_role = "pre",
      neuron_id = pre_root_id,
      depth = pre_rz
    )]
    post_part <- coord[post_type %in% types, .(
      neuropil = argvs$np,
      cell_type = post_type,
      calibration_group = group_name,
      synapse_role = "post",
      neuron_id = post_root_id,
      depth = post_rz
    )]
    rbind(pre_part, post_part)
  }), fill = TRUE)

  if (nrow(rows) == 0) {
    return(rows)
  }

  if (zinvert) {
    rows[, depth := depth * -1]
  }

  rows[!is.na(depth)]
}

make_summary <- function(depths) {
  role_summary <- depths[, .(
    n_synapses = .N,
    n_neurons = uniqueN(neuron_id),
    median_depth = median(depth),
    q25_depth = as.numeric(quantile(depth, 0.25)),
    q75_depth = as.numeric(quantile(depth, 0.75)),
    min_depth = min(depth),
    max_depth = max(depth)
  ), by = .(neuropil, cell_type, calibration_group, synapse_role)]

  combined_summary <- depths[, .(
    n_synapses = .N,
    n_neurons = uniqueN(paste(synapse_role, neuron_id, sep = "_")),
    median_depth = median(depth),
    q25_depth = as.numeric(quantile(depth, 0.25)),
    q75_depth = as.numeric(quantile(depth, 0.75)),
    min_depth = min(depth),
    max_depth = max(depth)
  ), by = .(neuropil, cell_type, calibration_group)]
  combined_summary[, synapse_role := "combined"]

  out <- rbind(role_summary, combined_summary, fill = TRUE)
  setcolorder(out, c(
    "neuropil", "cell_type", "calibration_group", "synapse_role",
    "n_synapses", "n_neurons", "median_depth", "q25_depth", "q75_depth",
    "min_depth", "max_depth"
  ))
  setorder(out, calibration_group, cell_type, synapse_role)
  out
}

safe_density <- function(x, grid) {
  if (length(unique(x)) < 2) {
    return(rep(NA_real_, length(grid)))
  }
  den <- density(x)
  approx(den$x, den$y, xout = grid, rule = 2)$y
}

make_plot <- function(depths, summary, meta, output_file) {
  group_levels <- c("reference_superficial", "reference_deep", "extra")
  depths[, calibration_group := factor(calibration_group, levels = group_levels)]

  grid <- seq(meta$minz, meta$maxz, length.out = 1024)

  group_density <- rbindlist(lapply(levels(depths$calibration_group), function(group_name) {
    vals <- depths[calibration_group == group_name, depth]
    if (length(vals) == 0) return(NULL)
    data.table(
      depth = grid,
      density = safe_density(vals, grid),
      calibration_group = group_name
    )
  }))
  group_density <- group_density[!is.na(density)]

  combined_summary <- summary[synapse_role == "combined"]
  type_levels <- combined_summary[order(calibration_group, median_depth), cell_type]
  depths[, cell_type := factor(cell_type, levels = type_levels)]

  type_density <- rbindlist(lapply(type_levels, function(type_name) {
    vals <- depths[cell_type == type_name, depth]
    den <- safe_density(vals, grid)
    data.table(
      depth = grid,
      density = den,
      cell_type = type_name,
      calibration_group = as.character(depths[cell_type == type_name, calibration_group][1])
    )
  }))
  type_density <- type_density[!is.na(density)]
  type_density[, scaled_density := density / max(density), by = cell_type]

  group_plot <- ggplot(group_density, aes(x = depth, y = density, color = calibration_group)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(
      values = c(
        reference_superficial = "#3B7EA1",
        reference_deep = "#B85C38",
        extra = "#6B6B6B"
      ),
      drop = FALSE
    ) +
    labs(
      title = sprintf("%s reference depth calibration", argvs$np),
      x = "PC3 / z depth",
      y = "Density",
      color = "Group"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  type_plot <- ggplot(type_density, aes(x = depth, y = cell_type, fill = scaled_density)) +
    geom_raster() +
    facet_grid(calibration_group ~ ., scales = "free_y", space = "free_y", drop = TRUE) +
    scale_fill_gradient(low = "white", high = "#333333", na.value = "white") +
    labs(x = "PC3 / z depth", y = NULL, fill = "Scaled density") +
    theme_minimal() +
    theme(
      strip.text.y = element_text(angle = 0),
      legend.position = "bottom"
    )

  if (requireNamespace("patchwork", quietly = TRUE)) {
    final_plot <- group_plot / type_plot + patchwork::plot_layout(heights = c(1, 2))
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    final_plot <- gridExtra::arrangeGrob(group_plot, type_plot, ncol = 1, heights = c(1, 2))
  } else {
    warning("Neither patchwork nor gridExtra is available; writing density plot only")
    final_plot <- group_plot
  }

  height <- max(7, min(16, 5 + 0.18 * length(type_levels)))
  ggsave(output_file, plot = final_plot, width = 9, height = height, limitsize = FALSE)
}

validate_file(argvs$synf, "Synapse")
validate_file(argvs$meta, "Metadata")
validate_file(argvs$ref_groups, "Reference groups")

plot_meta_all <- fread(argvs$meta)
validate_columns(
  plot_meta_all,
  c("neuropil", "zinvert", "minz", "maxz"),
  argvs$meta
)
plot_meta <- plot_meta_all[neuropil == argvs$np]
if (nrow(plot_meta) == 0) {
  stop(sprintf("No metadata found for neuropil: %s", argvs$np))
}

ref_groups_config <- fread(argvs$ref_groups)
validate_columns(
  ref_groups_config,
  c("neuropil", "ref_superficial", "ref_deep", "pattern_type"),
  argvs$ref_groups
)
ref_config <- ref_groups_config[neuropil == argvs$np]
if (nrow(ref_config) == 0) {
  stop(sprintf("No reference groups found for neuropil: %s", argvs$np))
}
if (nrow(ref_config) > 1) {
  stop(sprintf("Multiple reference group rows found for neuropil: %s", argvs$np))
}

coord <- fread(argvs$synf, colClasses = c(pre_root_id = "character", post_root_id = "character"))
validate_columns(
  coord,
  c("pre_type", "post_type", "pre_root_id", "post_root_id", "pre_rz", "post_rz"),
  argvs$synf
)

available_groups <- unique(c(coord$pre_type, coord$post_type))
refs <- expand_reference_groups(ref_config, available_groups)

extra_types_requested <- parse_type_list(argvs$types)
reference_types <- c(refs$superficial, refs$deep)
extra_types <- setdiff(intersect(extra_types_requested, available_groups), reference_types)
missing_extra_types <- setdiff(extra_types_requested, available_groups)
if (length(missing_extra_types) > 0) {
  warning("Optional calibration types not found in ", argvs$np, ": ",
          paste(sort(missing_extra_types), collapse = ", "))
}

type_groups <- list(
  reference_superficial = refs$superficial,
  reference_deep = refs$deep,
  extra = sort(extra_types)
)

cat(sprintf("\n=== Reference depth calibration for %s ===\n", argvs$np))
cat(sprintf("Superficial reference types (%d): %s\n",
            length(type_groups$reference_superficial),
            paste(type_groups$reference_superficial, collapse = ", ")))
cat(sprintf("Deep reference types (%d): %s\n",
            length(type_groups$reference_deep),
            paste(type_groups$reference_deep, collapse = ", ")))
cat(sprintf("Extra calibration types (%d): %s\n",
            length(type_groups$extra),
            paste(type_groups$extra, collapse = ", ")))

depths <- collect_depths(coord, type_groups, isTRUE(plot_meta$zinvert))
if (nrow(depths) == 0) {
  stop(sprintf("No depth records found for calibration types in %s", argvs$np))
}

summary <- make_summary(depths)
csv_file <- paste0(argvs$output_prefix, ".csv")
pdf_file <- paste0(argvs$output_prefix, ".pdf")

fwrite(summary, csv_file)
make_plot(depths, summary, plot_meta, pdf_file)

cat("Wrote calibration summary:", csv_file, "\n")
cat("Wrote calibration plot:", pdf_file, "\n")
