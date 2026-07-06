#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrastr))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))

options("ggrastr.default.dpi" = 450)

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

parse_bool <- function(x, default = FALSE) {
  if (is.null(x) || is.na(x)) return(default)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

### TODO
if (interactive()) {
  source("./src/utils.r", chdir = FALSE)
  argvs$np <- "ME_R"
  argvs$syn_type <- "post"
  argvs$use_preset <- "type_putative"
  argvs$density <- "asis"
  argvs$ann <- "data/visual_neurons_anno.csv"
  argvs$meta <- "data/viz_meta.csv"
  argvs$preset <- "data/viz_preset.csv"
  argvs$subsample <- 10000L
  argvs$sparse_limit <- 100L
  argvs$use_axis_limits <- TRUE
  syn_path <- file.path("int/idv_mat/", paste0(argvs$np, "_rotated.csv.gz"))
} else {
  argvs$subsample <- as.integer(argvs$subsample)
  argvs$sparse_limit <- as.integer(argvs$sparse_limit)
  syn_path <- argvs$synf
}
argvs$use_axis_limits <- parse_bool(argvs$use_axis_limits, default = TRUE)

# Source utility functions (soft-linked by Nextflow)
if (file.exists("./utils.r")) {
  source("./utils.r", chdir = FALSE)
} else {
  source(argvs$utils, chdir = FALSE)
}

## Load plot metadata and presets with validation
plot_meta_all <- validate_config_file(
  argvs$meta,
  c(
    "neuropil", "x_axis", "y_axis", "axis_1_func", "axis_2_func", "min1",
    "max1", "min2", "max2", "zid", "zinvert", "minz", "maxz", "outlayout",
     "outr1", "outr2", "outd1", "outd2"
  )
  )
plot_meta <- plot_meta_all[neuropil == argvs$np]
if (nrow(plot_meta) == 0) stop(sprintf("No metadata found for neuropil: %s", argvs$np))

preset_all <- validate_config_file(
  argvs$preset,
  c(
     "preset", "palette", "color_guide", "filter_func", "color_by",
      "notch_split", "do_highlight", "hl_type", "hl_col", "hl_val", "known_only"
  )
  )
preset <- preset_all[preset == argvs$use_preset]
if (nrow(preset) == 0) stop(sprintf("No preset found: %s", argvs$use_preset))

# Validate that required functions exist
validate_functions(preset, plot_meta)

## Load annotations and coordinates
opc_anno <- validate_config_file(argvs$ann, c("cell_type"))
if (!file.exists(syn_path)) stop(sprintf("Synapse coordinate file not found: %s", syn_path))
np_coord <- fread(syn_path)

## Validate required columns in coordinate data
required_coord_cols <- c(
  paste0(argvs$syn_type, "_type"),
  paste0(argvs$syn_type, "_rz"),
  paste(argvs$syn_type, plot_meta$x_axis, sep = "_"),
  paste(argvs$syn_type, plot_meta$y_axis, sep = "_")
)
missing_coord_cols <- setdiff(required_coord_cols, colnames(np_coord))
if (length(missing_coord_cols) > 0) {
  stop(sprintf("Missing required columns in coordinate data: %s", paste(missing_coord_cols, collapse=", ")))
}

## Dynamic layout generation
filter_func <- get(preset$filter_func)
color_func <- get(preset$palette)
scale_axis_1 <- get(plot_meta$axis_1_func)
scale_axis_2 <- get(plot_meta$axis_2_func)
x_axis <- paste(argvs$syn_type, plot_meta$x_axis, sep = "_")
y_axis <- paste(argvs$syn_type, plot_meta$y_axis, sep = "_")
axis_1_limits <- c(plot_meta$min1, plot_meta$max1)
axis_2_limits <- c(plot_meta$min2, plot_meta$max2)
if (grepl("reverse", plot_meta$axis_1_func)) {
  axis_1_limits <- rev(axis_1_limits)
}
if (grepl("reverse", plot_meta$axis_2_func)) {
  axis_2_limits <- rev(axis_2_limits)
}
axis_scales <- function() {
  if (argvs$use_axis_limits) {
    list(
      scale_axis_1(limits = axis_1_limits),
      scale_axis_2(limits = axis_2_limits)
    )
  } else {
    list(scale_axis_1(), scale_axis_2())
  }
}
if (preset$color_by == "cell_type") {
  preset$color_by <- paste0(argvs$syn_type, "_type")
}
spatial_notch_levels <- c("Notch On", "Notch Off")

## Filter by annotation and subsample
np_coord <- filter_type(
  np_coord, syn_type = argvs$syn_type, sparse_limit = argvs$sparse_limit
)
np_coord <- filter_func(np_coord, opc_anno, syn_type = argvs$syn_type)

## Validate that all data points fall within specified axis limits
x_values <- np_coord[[x_axis]]
y_values <- np_coord[[y_axis]]

if (preset$do_highlight && preset$hl_type == "label") {
  if (!preset$hl_col %in% colnames(np_coord)) {
    stop(sprintf("Highlight column not found in data: %s", preset$hl_col))
  }
  np_coord[, highlight := get(preset$hl_col) == preset$hl_val]
}
if (preset$do_highlight && preset$hl_type == "path") {
  if (!file.exists(preset$hl_val)) {
    stop("Cannot find highlight label file.")
  } else {
    hl_labels <- readLines(preset$hl_val)
  }
  if (!preset$hl_col %in% colnames(np_coord)) {
    stop(sprintf("Highlight column not found in data: %s", preset$hl_col))
  }
  np_coord[, highlight := get(preset$hl_col) %in% hl_labels]
}


np_raw <- copy(np_coord)
if (grepl("^ME", argvs$np)) {
  np_coord <- rsubsample(
    np_coord, n = argvs$subsample * 4, syn_type = argvs$syn_type
  )
} else {
  np_coord <- np_coord[!grepl("^(Mi|Dm|Pm)", get(paste0(argvs$syn_type, "_type")))]
  np_coord <- rsubsample(
    np_coord, n = argvs$subsample * 2, syn_type = argvs$syn_type
  )
}

if (preset$hl_col == "func") {
  np_coord <- np_coord[func != "Unannotated"]
  np_raw <- np_raw[func != "Unannotated"]
}

## Apply known_only filtering if specified (after subsampling)
if (preset$known_only) {
  np_coord <- np_coord[newly_ann == "N"]
  np_raw <- np_raw[newly_ann == "N"]
}

if (preset$notch_split) {
  required_notch_cols <- c("Notch", "ntype")
  missing_notch_cols <- setdiff(required_notch_cols, colnames(np_coord))
  if (length(missing_notch_cols) > 0) {
    stop(sprintf("Missing columns for notch_split: %s", paste(missing_notch_cols, collapse=", ")))
  }
  np_coord$notch_ntype <- paste(np_coord$Notch, np_coord$ntype, sep = "_")
  np_coord$notch_ntype <- sub("^_", "", np_coord$notch_ntype)
  np_coord <- split(np_coord, np_coord$notch_ntype, drop = TRUE)
  np_raw$notch_ntype <- paste(np_raw$Notch, np_raw$ntype, sep = "_")
  np_raw$notch_ntype <- sub("^_", "", np_raw$notch_ntype)
  np_raw <- split(np_raw, np_raw$notch_ntype, drop = TRUE)
} else {
  np_coord <- list(all = np_coord)
  np_raw <- list(all = np_raw)
}

for (i in names(np_coord)) {
  # Skip if no data after filtering
  if (nrow(np_coord[[i]]) == 0) {
    next
  }

  if (preset$color_by == "spatial_origin") {
    np_coord[[i]][, .plot_order := spatial_origin == "unknown"]
    setorder(np_coord[[i]], -.plot_order)
    np_coord[[i]][, .plot_order := NULL]
  }
  
  out_prefix <- paste(
    argvs$np, argvs$syn_type, argvs$use_preset, argvs$density, i,
    sep = "_"
  )

  ## Generate the dot plot
  if (preset$do_highlight) {
    dotp <- np_coord[[i]] |>
      ggplot(aes(x = .data[[x_axis]], y = .data[[y_axis]])) +
      rasterize(geom_point(
        aes(color = .data[[preset$color_by]], alpha = highlight)
      )) +
      guides(alpha = "none") +
      scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.05))
  } else {
    dotp <- np_coord[[i]] |>
      ggplot(aes(x = .data[[x_axis]], y = .data[[y_axis]])) +
      rasterize(geom_point(aes(color = .data[[preset$color_by]])))
  }
  dotp <- dotp +
    labs(color = preset$color_guide) +
    theme_ih2025() +
    axis_scales() +
    color_func()
  
  ## Extract legends to prevent layout fluctuation
  legendsp <- get_plot_component(dotp, "guide-box", return_all = TRUE) 
  dotp <- dotp + theme(legend.position="none")
  
 if (!inherits(legendsp, "grob")) {
   to_keep <- sapply(legendsp, function(x) !"zeroGrob" %in% class(x))
   if (sum(to_keep) != 1) {
     stop("There should be only 1 legend but ", sum(to_keep), " was found.")
   }
   legendsp <- legendsp[to_keep][[1]]
 }
  
  ## Manual calculation and interpolation of density
  if (preset$do_highlight) {
    np_den <- np_raw[[i]][highlight == TRUE, .(
      group = get(preset$color_by),
      type = get(paste(argvs$syn_type, "type", sep = "_")),
      depth = get(paste(argvs$syn_type, "rz", sep = "_"))
    )]
  } else {
    np_den <- np_raw[[i]][ , .(
      group = get(preset$color_by),
      type = get(paste(argvs$syn_type, "type", sep = "_")),
      depth = get(paste(argvs$syn_type, "rz", sep = "_"))
    )]
  }
  
  if (nrow(np_den) == 0) next
  
  
  if (plot_meta$zinvert) {
    np_den$depth <- np_den$depth * -1
  }
  
  np_den <- split(np_den, np_den$group, drop = TRUE)
  
  den_grid <- seq(plot_meta$minz, plot_meta$maxz, length.out = 1024)
  if (argvs$density == "asis") {
    np_den <- lapply(np_den, function(x) return(density(x$depth)))
    interpolated <- lapply(names(np_den), function(type) {
      d <- np_den[[type]]
      y_interp <- approx(d$x, d$y, xout = den_grid, rule = 2)$y
      return(data.frame(x = den_grid, y = y_interp, type = type))
    })
    inter <- do.call(rbind, interpolated)
  } else {
    # No need to filter here since we already filtered at the type level
    if (length(np_den) < 1) {
      next
    }
    type_avg <- lapply(names(np_den), function(tn) {
      x <- np_den[[tn]]
      id_type <- split(x, x$type, drop = TRUE)
      den_type <- lapply(id_type, function(idt) return(density(idt$depth)))
      interpolated <- lapply(names(den_type), function(type) {
        d <- den_type[[type]]
        y_interp <- approx(d$x, d$y, xout = den_grid, rule = 2)$y
        return(data.frame(x = den_grid, y_interp = y_interp, type = type))
      })
      inter <- do.call(rbind, interpolated)
      out <- data.frame(
        x = den_grid,
        y = tapply(inter$y_interp, inter$x, mean),
        type = tn
      )
      return(out)
    })
    inter <- do.call(rbind, type_avg)
  }
  
  denp <- inter |>
    ggplot(aes(x = x, y = y, color = type)) +
    geom_line() +
    theme_ih2025() +
    color_func() +
    guides(color = "none")
  
  if (plot_meta$zid == "y") denp <- denp + coord_flip()
  
  if (plot_meta$outlayout == "landscape") {
    outp <- (dotp | denp) +
      plot_layout(widths = c(plot_meta$outr1, plot_meta$outr2))
  } else {
    outp <- (dotp / denp) +
      plot_layout(heights = c(plot_meta$outr1, plot_meta$outr2))
  }
  ggsave(
    plot = outp, filename = paste0(out_prefix, ".pdf"),
    width = plot_meta$outd1, height = plot_meta$outd2
  )
  ggsave(
    plot = legendsp, paste0(out_prefix, "_legend.pdf")
  )

  if (preset$color_by == "spatial_origin" && "Notch" %in% colnames(np_coord[[i]])) {
    for (notch_status in spatial_notch_levels) {
      notch_label <- gsub(" ", "", notch_status)
      np_coord_notch <- copy(np_coord[[i]][Notch == notch_status])
      np_raw_notch <- copy(np_raw[[i]][Notch == notch_status])
      if (nrow(np_coord_notch) == 0 || nrow(np_raw_notch) == 0) next

      dotp_notch <- np_coord_notch |>
        ggplot(aes(x = .data[[x_axis]], y = .data[[y_axis]])) +
        rasterize(geom_point(aes(color = .data[[preset$color_by]]))) +
        labs(color = preset$color_guide) +
        theme_ih2025() +
        axis_scales() +
        color_func()

      legendsp_notch <- get_plot_component(dotp_notch, "guide-box", return_all = TRUE)
      dotp_notch <- dotp_notch + theme(legend.position = "none")
      if (!inherits(legendsp_notch, "grob")) {
        to_keep <- sapply(legendsp_notch, function(x) !"zeroGrob" %in% class(x))
        if (sum(to_keep) != 1) {
          stop("There should be only 1 legend but ", sum(to_keep), " was found.")
        }
        legendsp_notch <- legendsp_notch[to_keep][[1]]
      }

      np_den_notch <- np_raw_notch[, .(
        group = get(preset$color_by),
        type = get(paste(argvs$syn_type, "type", sep = "_")),
        depth = get(paste(argvs$syn_type, "rz", sep = "_"))
      )]
      if (plot_meta$zinvert) {
        np_den_notch$depth <- np_den_notch$depth * -1
      }
      np_den_notch <- split(np_den_notch, np_den_notch$group, drop = TRUE)
      np_den_notch <- Filter(
        function(x) length(unique(x$depth)) >= 2,
        np_den_notch
      )
      if (length(np_den_notch) == 0) next
      interpolated_notch <- lapply(names(np_den_notch), function(type) {
        d <- density(np_den_notch[[type]]$depth)
        y_interp <- approx(d$x, d$y, xout = den_grid, rule = 2)$y
        data.frame(x = den_grid, y = y_interp, type = type)
      })
      inter_notch <- do.call(rbind, interpolated_notch)

      denp_notch <- inter_notch |>
        ggplot(aes(x = x, y = y, color = type)) +
        geom_line() +
        theme_ih2025() +
        color_func() +
        guides(color = "none")

      if (plot_meta$zid == "y") denp_notch <- denp_notch + coord_flip()

      if (plot_meta$outlayout == "landscape") {
        outp_notch <- (dotp_notch | denp_notch) +
          plot_layout(widths = c(plot_meta$outr1, plot_meta$outr2))
      } else {
        outp_notch <- (dotp_notch / denp_notch) +
          plot_layout(heights = c(plot_meta$outr1, plot_meta$outr2))
      }

      ggsave(
        plot = outp_notch,
        filename = paste0(out_prefix, "_", notch_label, ".pdf"),
        width = plot_meta$outd1,
        height = plot_meta$outd2
      )
      ggsave(
        plot = legendsp_notch,
        paste0(out_prefix, "_", notch_label, "_legend.pdf")
      )
    }
  }
}
