#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(cowplot))

parse_args <- function(defaults = list()) {
  args <- commandArgs(trailingOnly = TRUE)
  out <- defaults
  if (length(args) == 0) return(out)

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected argument: %s", key))
    }
    key <- sub("^--", "", key)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

validate_file <- function(path, label) {
  if (is.null(path) || !file.exists(path)) {
    stop(sprintf("%s not found: %s", label, path))
  }
  invisible(path)
}

validate_columns <- function(dt, cols, label) {
  missing_cols <- setdiff(cols, colnames(dt))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "%s is missing required columns: %s",
      label,
      paste(missing_cols, collapse = ", ")
    ))
  }
  invisible(TRUE)
}

as_bool <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

reference_types <- function(ref_config, available_types) {
  if (nrow(ref_config) != 1L) {
    stop("Expected exactly one ME_R row in reference group config")
  }

  if (ref_config$pattern_type == "regex") {
    superficial <- grep(ref_config$ref_superficial, available_types, value = TRUE)
    deep <- grep(ref_config$ref_deep, available_types, value = TRUE)
  } else if (ref_config$pattern_type == "exact") {
    superficial <- intersect(
      strsplit(ref_config$ref_superficial, ";", fixed = TRUE)[[1]],
      available_types
    )
    deep <- intersect(
      strsplit(ref_config$ref_deep, ";", fixed = TRUE)[[1]],
      available_types
    )
  } else {
    stop(sprintf("Unknown pattern_type: %s", ref_config$pattern_type))
  }

  if (length(superficial) == 0 || length(deep) == 0) {
    stop("Could not identify both ME_R reference groups")
  }

  list(superficial = superficial, deep = deep)
}

reference_endpoints <- function(coord, ref_superficial, ref_deep, x_suffix, z_suffix) {
  pre_ref <- coord[
    pre_type %in% c(ref_superficial, ref_deep),
    .(
      cell_type = pre_type,
      x = get(paste0("pre_", x_suffix)),
      depth = get(paste0("pre_", z_suffix)),
      reference_group = fifelse(
        pre_type %in% ref_superficial,
        "Dm reference",
        "Pm reference"
      )
    )
  ]
  post_ref <- coord[
    post_type %in% c(ref_superficial, ref_deep),
    .(
      cell_type = post_type,
      x = get(paste0("post_", x_suffix)),
      depth = get(paste0("post_", z_suffix)),
      reference_group = fifelse(
        post_type %in% ref_superficial,
        "Dm reference",
        "Pm reference"
      )
    )
  ]

  refs <- rbind(pre_ref, post_ref)
  refs <- refs[is.finite(x) & is.finite(depth)]
  refs[, reference_group := factor(
    reference_group,
    levels = c("Pm reference", "Dm reference")
  )]
  refs
}

subsample_rows <- function(dt, n, seed) {
  n <- as.integer(n)
  if (is.na(n) || n <= 0L || nrow(dt) <= n) return(copy(dt))
  set.seed(as.integer(seed))
  dt[sample.int(nrow(dt), n)]
}

make_density <- function(refs, depth_grid) {
  split_refs <- split(refs, refs$reference_group, drop = TRUE)
  density_parts <- lapply(names(split_refs), function(group_name) {
    x <- split_refs[[group_name]]
    if (length(unique(x$depth)) < 2L) {
      return(NULL)
    }
    d <- density(x$depth, na.rm = TRUE)
    data.table(
      depth = depth_grid,
      density = approx(d$x, d$y, xout = depth_grid, rule = 2)$y,
      reference_group = factor(
        group_name,
        levels = levels(refs$reference_group)
      )
    )
  })
  rbindlist(density_parts)
}

darken_color <- function(color, factor = 0.65) {
  rgb <- grDevices::col2rgb(color) * factor
  grDevices::rgb(rgb[1, ], rgb[2, ], rgb[3, ], maxColorValue = 255)
}

defaults <- list(
  synf = "int/idv_mat/ME_R_rotated.csv.gz",
  meta = "data/viz_meta.csv",
  ref_groups = "data/reference_groups.csv",
  broad_depth_csv = paste0(
    "int/stats/deep_superficial/spatial_all/ME_R_post/",
    "ME_R_post_spatial_all_broad_depth.csv"
  ),
  utils = "src/utils.r",
  out = "int/broad_depth_schema/ME_R_broad_depth_schema.pdf",
  legend_out = NULL,
  scatter_subsample = "100000",
  seed = "1",
  width = "10",
  height = "6"
)
argvs <- parse_args(defaults)
if (is.null(argvs$legend_out)) {
  argvs$legend_out <- sub("\\.pdf$", "_legend.pdf", argvs$out)
}

validate_file(argvs$synf, "ME_R coordinate file")
validate_file(argvs$meta, "Visualization metadata")
validate_file(argvs$ref_groups, "Reference group config")
validate_file(argvs$broad_depth_csv, "Broad-depth output CSV")
validate_file(argvs$utils, "Utility file")
source(argvs$utils, chdir = FALSE)

plot_meta_all <- fread(argvs$meta)
validate_columns(
  plot_meta_all,
  c("neuropil", "x_axis", "y_axis", "axis_1_func", "axis_2_func",
    "min1", "max1", "min2", "max2", "zid", "zinvert", "minz", "maxz"),
  "Visualization metadata"
)
plot_meta <- plot_meta_all[neuropil == "ME_R"]
if (nrow(plot_meta) != 1L) {
  stop("Expected exactly one ME_R row in visualization metadata")
}
if (plot_meta$zid != "y") {
  stop("This schema script currently expects ME_R depth on the y axis")
}

coord <- fread(argvs$synf)
required_coord_cols <- c(
  "pre_type", "post_type",
  paste0("pre_", plot_meta$x_axis),
  paste0("post_", plot_meta$x_axis),
  paste0("pre_", plot_meta$y_axis),
  paste0("post_", plot_meta$y_axis)
)
validate_columns(coord, required_coord_cols, "ME_R coordinate data")

available_types <- unique(c(coord$pre_type, coord$post_type))
ref_config <- fread(argvs$ref_groups)[neuropil == "ME_R"]
refs <- reference_types(ref_config, available_types)
ref_points <- reference_endpoints(
  coord,
  ref_superficial = refs$superficial,
  ref_deep = refs$deep,
  x_suffix = plot_meta$x_axis,
  z_suffix = plot_meta$y_axis
)
if (nrow(ref_points) == 0L) {
  stop("No ME_R reference endpoints were found")
}

broad_depth <- fread(argvs$broad_depth_csv)
validate_columns(
  broad_depth,
  c("observed_delta_thres", "observed_delta_thres_base", "coefficient"),
  "Broad-depth output CSV"
)
threshold <- unique(broad_depth$observed_delta_thres)
threshold <- threshold[is.finite(threshold)]
if (length(threshold) == 0L) {
  coefficient <- unique(broad_depth$coefficient)
  coefficient <- coefficient[is.finite(coefficient)][1]
  threshold <- unique(broad_depth$observed_delta_thres_base)
  threshold <- threshold[is.finite(threshold)][1] * coefficient
} else {
  threshold <- threshold[1]
}

reference_medians <- ref_points[
  ,
  .(median_depth = median(depth, na.rm = TRUE)),
  by = reference_group
]
median_base <- diff(range(reference_medians$median_depth))
artifact_base <- unique(broad_depth$observed_delta_thres_base)
artifact_base <- artifact_base[is.finite(artifact_base)][1]
if (
  is.finite(median_base) &&
  is.finite(artifact_base) &&
  abs(median_base - artifact_base) > max(1, artifact_base * 0.02)
) {
  warning(sprintf(
    paste(
      "Reference median separation from plotted endpoints (%.3f) differs",
      "from broad-depth artifact (%.3f). The schema still uses the artifact",
      "threshold and plotted endpoint medians for placement."
    ),
    median_base,
    artifact_base
  ))
}

neutral_midpoint <- mean(reference_medians$median_depth)
neutral_lower <- neutral_midpoint - threshold / 2
neutral_upper <- neutral_midpoint + threshold / 2
axis_x_limits <- c(plot_meta$min1, plot_meta$max1)
axis_depth_limits <- c(plot_meta$min2, plot_meta$max2)
plot_x_limits <- axis_x_limits
plot_depth_limits <- axis_depth_limits
if (grepl("reverse", plot_meta$axis_1_func)) {
  plot_x_limits <- rev(plot_x_limits)
}
if (grepl("reverse", plot_meta$axis_2_func)) {
  plot_depth_limits <- rev(plot_depth_limits)
}
depth_min <- min(axis_depth_limits)
depth_max <- max(axis_depth_limits)

colors <- ih2025_colors()
reference_colors <- c(
  "Pm reference" = colors[["Early"]],
  "Dm reference" = colors[["Late"]]
)
reference_labels <- c(
  "Pm reference" = "Proximal reference (Early)",
  "Dm reference" = "Distal reference (Late)"
)
pm_median <- reference_medians[reference_group == "Pm reference", median_depth]
dm_median <- reference_medians[reference_group == "Dm reference", median_depth]
if (length(pm_median) != 1L || length(dm_median) != 1L) {
  stop("Expected one median depth for each reference group")
}

tail_regions <- data.table(
  region = c("lower", "upper"),
  xmin = -Inf,
  xmax = Inf,
  ymin = c(depth_min, neutral_upper),
  ymax = c(neutral_lower, depth_max)
)
tail_regions <- tail_regions[ymax > ymin]
tail_regions[
  ,
  zone := fifelse(
    pm_median >= ymin & pm_median <= ymax,
    "Proximal rejection zone",
    fifelse(
      dm_median >= ymin & dm_median <= ymax,
      "Distal rejection zone",
      NA_character_
    )
  )
]
if (any(is.na(tail_regions$zone))) {
  tail_regions[
    region == "lower",
    zone := if (pm_median < dm_median) {
      "Proximal rejection zone"
    } else {
      "Distal rejection zone"
    }
  ]
  tail_regions[
    region == "upper",
    zone := if (pm_median > dm_median) {
      "Proximal rejection zone"
    } else {
      "Distal rejection zone"
    }
  ]
}
tail_regions[, zone := factor(
  zone,
  levels = c("Proximal rejection zone", "Distal rejection zone")
)]

rejection_fill_colors <- c(
  "Proximal rejection zone" = colors[["Early"]],
  "Distal rejection zone" = colors[["Late"]]
)
hypothetical_colors <- c(
  "Significantly proximal" = darken_color(colors[["Early"]]),
  "Significantly distal" = darken_color(colors[["Late"]]),
  "Insignificant, proximal median" = "black",
  "Insignificant, distal median" = "black",
  "Insignificant, neutral median" = "black"
)
plot_colors <- c(reference_colors, hypothetical_colors)
plot_labels <- c(
  reference_labels,
  setNames(names(hypothetical_colors), names(hypothetical_colors))
)

region_midpoint <- function(region) {
  mean(c(region$ymin, region$ymax))
}
inward_from_boundary <- function(boundary, region, fraction = 0.25) {
  if (region$region == "lower") {
    boundary - fraction * (boundary - region$ymin)
  } else {
    boundary + fraction * (region$ymax - boundary)
  }
}
prox_region <- tail_regions[zone == "Proximal rejection zone"]
dist_region <- tail_regions[zone == "Distal rejection zone"]
if (nrow(prox_region) != 1L || nrow(dist_region) != 1L) {
  stop("Expected one proximal and one distal rejection zone")
}
prox_boundary <- if (prox_region$region == "lower") neutral_lower else neutral_upper
dist_boundary <- if (dist_region$region == "lower") neutral_lower else neutral_upper

hypo_x <- plot_x_limits[1] + diff(plot_x_limits) * seq(0.12, 0.88, length.out = 5)
hypothetical_ranges <- data.table(
  label = factor(
    c(
      "Significantly distal",
      "Significantly proximal",
      "Insignificant, proximal median",
      "Insignificant, distal median",
      "Insignificant, neutral median"
    ),
    levels = names(hypothetical_colors)
  ),
  x = hypo_x,
  median = c(
    region_midpoint(dist_region),
    region_midpoint(prox_region),
    region_midpoint(prox_region),
    region_midpoint(dist_region),
    neutral_midpoint
  )
)
hypothetical_ranges[
  label == "Significantly distal",
  `:=`(
    ymin = dist_region$ymin + 0.20 * (dist_region$ymax - dist_region$ymin),
    ymax = dist_region$ymin + 0.80 * (dist_region$ymax - dist_region$ymin)
  )
]
hypothetical_ranges[
  label == "Significantly proximal",
  `:=`(
    ymin = prox_region$ymin + 0.20 * (prox_region$ymax - prox_region$ymin),
    ymax = prox_region$ymin + 0.80 * (prox_region$ymax - prox_region$ymin)
  )
]
hypothetical_ranges[
  label == "Insignificant, proximal median",
  `:=`(
    median = inward_from_boundary(prox_boundary, prox_region, fraction = 0.25),
    ymin = min(region_midpoint(prox_region), prox_boundary) - threshold * 0.20,
    ymax = max(region_midpoint(prox_region), prox_boundary) + threshold * 0.20
  )
]
hypothetical_ranges[
  label == "Insignificant, distal median",
  `:=`(
    median = inward_from_boundary(dist_boundary, dist_region, fraction = 0.25),
    ymin = min(region_midpoint(dist_region), dist_boundary) - threshold * 0.20,
    ymax = max(region_midpoint(dist_region), dist_boundary) + threshold * 0.20
  )
]
hypothetical_ranges[
  label == "Insignificant, neutral median",
  `:=`(
    ymin = neutral_midpoint - threshold * 0.20,
    ymax = neutral_midpoint + threshold * 0.20
  )
]
hypothetical_ranges[, ymin := pmax(ymin, depth_min)]
hypothetical_ranges[, ymax := pmin(ymax, depth_max)]

scatter_points <- subsample_rows(
  ref_points,
  n = argvs$scatter_subsample,
  seed = argvs$seed
)
depth_grid <- seq(depth_min, depth_max, length.out = 1024)
density_points <- make_density(ref_points, depth_grid)

axis_scales <- list(
  scale_x_continuous(limits = plot_x_limits),
  scale_y_reverse(limits = plot_depth_limits)
)

scatter_plot <- ggplot(scatter_points, aes(x = x, y = depth)) +
  geom_rect(
    data = tail_regions,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = zone),
    inherit.aes = FALSE,
    alpha = 0.10
  ) +
  geom_point(aes(color = reference_group), size = 0.12, alpha = 0.50) +
  geom_linerange(
    data = hypothetical_ranges,
    aes(x = x, ymin = ymin, ymax = ymax, color = label),
    inherit.aes = FALSE,
    linewidth = 1.1
  ) +
  geom_point(
    data = hypothetical_ranges,
    aes(x = x, y = median, color = label),
    inherit.aes = FALSE,
    size = 3.9
  ) +
  geom_hline(
    yintercept = c(neutral_lower, neutral_upper),
    linewidth = 0.25,
    linetype = "dashed",
    color = "grey30"
  ) +
  scale_color_manual(
    values = plot_colors,
    labels = plot_labels,
    breaks = names(plot_colors),
    name = NULL
  ) +
  scale_fill_manual(
    values = rejection_fill_colors,
    name = NULL
  ) +
  axis_scales +
  theme_ih2025() +
  theme(legend.position = "bottom")

density_plot <- ggplot(density_points, aes(x = density, y = depth)) +
  geom_rect(
    data = tail_regions,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = zone),
    inherit.aes = FALSE,
    alpha = 0.10
  ) +
  geom_path(aes(color = reference_group), linewidth = 0.7) +
  geom_hline(
    yintercept = c(neutral_lower, neutral_upper),
    linewidth = 0.25,
    linetype = "dashed",
    color = "grey30"
  ) +
  scale_color_manual(values = plot_colors, guide = "none") +
  scale_fill_manual(values = rejection_fill_colors, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_y_reverse(limits = plot_depth_limits) +
  theme_ih2025() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

legend_plot <- scatter_plot +
  guides(
    fill = guide_legend(override.aes = list(alpha = 0.35), order = 1),
    color = guide_legend(
      override.aes = list(size = 3, alpha = 1, linewidth = 0.7),
      order = 2
    )
  )
legend_plot <- legend_plot + theme(legend.position = "bottom")
legendsp <- get_plot_component(legend_plot, "guide-box", return_all = TRUE)
if (!inherits(legendsp, "grob")) {
  to_keep <- sapply(legendsp, function(x) !"zeroGrob" %in% class(x))
  if (sum(to_keep) != 1) {
    stop("There should be only 1 legend but ", sum(to_keep), " was found.")
  }
  legendsp <- legendsp[to_keep][[1]]
}

scatter_plot <- scatter_plot + theme(legend.position = "none")
density_plot <- density_plot + theme(legend.position = "none")

out_plot <- scatter_plot + density_plot +
  plot_layout(widths = c(5, 1))

dir.create(dirname(argvs$out), recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = argvs$out,
  plot = out_plot,
  width = as.numeric(argvs$width),
  height = as.numeric(argvs$height)
)
ggsave(
  filename = argvs$legend_out,
  plot = legendsp
)

message(sprintf("Wrote ME_R broad-depth schema: %s", argvs$out))
message(sprintf("Wrote ME_R broad-depth schema legend: %s", argvs$legend_out))
