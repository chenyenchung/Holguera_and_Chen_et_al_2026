#!/usr/bin/env Rscript
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo_root <- normalizePath(file.path(dirname(script_file), "../../.."), mustWork = TRUE)
renv::load(repo_root)
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Seurat))

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

required_args <- c("stage", "depth", "seurat", "np", "syn_type")
missing_args <- required_args[vapply(required_args, function(x) is.null(argvs[[x]]), logical(1))]
if (length(missing_args) > 0) {
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

argvs$cluster_col <- if (is.null(argvs$cluster_col)) "FinalIdents" else argvs$cluster_col
argvs$assay <- if (is.null(argvs$assay)) "RNA" else argvs$assay
argvs$slot <- if (is.null(argvs$slot)) "data" else argvs$slot
argvs$test_use <- if (is.null(argvs$test_use)) "wilcox" else argvs$test_use
argvs$min_cells <- if (is.null(argvs$min_cells)) 3L else as.integer(argvs$min_cells)
argvs$min_pct <- if (is.null(argvs$min_pct)) 0 else as.numeric(argvs$min_pct)
argvs$logfc_threshold <- if (is.null(argvs$logfc_threshold)) 0 else as.numeric(argvs$logfc_threshold)

contrast_specs <- data.table(
  contrast = c(
    "all",
    "notch_on_projection",
    "notch_off_projection",
    "notch_on_intrinsic",
    "notch_off_intrinsic"
  ),
  Notch = c(NA_character_, "Notch On", "Notch Off", "Notch On", "Notch Off"),
  ntype = c(NA_character_, "projection", "projection", "intrinsic", "intrinsic")
)

collapse_chr <- function(x) {
  x <- sort(unique(as.character(x[!is.na(x)])))
  if (length(x) == 0) "" else paste(x, collapse = ";")
}

make_empty_markers <- function() {
  data.table(
    gene = character(),
    p_val = numeric(),
    avg_log2FC = numeric(),
    pct.1 = numeric(),
    pct.2 = numeric(),
    p_val_adj = numeric(),
    stage = character(),
    neuropil = character(),
    syn_type = character(),
    contrast = character(),
    group_1 = character(),
    group_2 = character(),
    n_types_superficial = integer(),
    n_types_deep = integer(),
    n_clusters_superficial = integer(),
    n_clusters_deep = integer(),
    n_cells_superficial = integer(),
    n_cells_deep = integer(),
    assay = character(),
    slot = character(),
    test_use = character()
  )
}

make_membership_row <- function(contrast, status, skip_reason = NA_character_,
                                sup_types = character(), deep_types = character(),
                                sup_clusters = character(), deep_clusters = character(),
                                dropped_overlap = character(),
                                n_cells_sup = NA_integer_, n_cells_deep = NA_integer_) {
  data.table(
    stage = argvs$stage,
    neuropil = argvs$np,
    syn_type = argvs$syn_type,
    contrast = contrast,
    status = status,
    skip_reason = skip_reason,
    superficial_types = collapse_chr(sup_types),
    deep_types = collapse_chr(deep_types),
    superficial_clusters = collapse_chr(sup_clusters),
    deep_clusters = collapse_chr(deep_clusters),
    dropped_overlap_clusters = collapse_chr(dropped_overlap),
    n_types_superficial = length(unique(sup_types)),
    n_types_deep = length(unique(deep_types)),
    n_clusters_superficial = length(unique(sup_clusters)),
    n_clusters_deep = length(unique(deep_clusters)),
    n_cells_superficial = n_cells_sup,
    n_cells_deep = n_cells_deep,
    min_cells = argvs$min_cells,
    assay = argvs$assay,
    slot = argvs$slot,
    test_use = argvs$test_use,
    min_pct = argvs$min_pct,
    logfc_threshold = argvs$logfc_threshold
  )
}

depth <- fread(argvs$depth)
required_depth_cols <- c(
  "cell_type", "ozel2021_cluster", "Notch", "ntype", "direction",
  "significant_fdr", "skip_reason", "neuropil", "syn_type"
)
missing_depth_cols <- setdiff(required_depth_cols, names(depth))
if (length(missing_depth_cols) > 0) {
  stop("Missing depth columns: ", paste(missing_depth_cols, collapse = ", "))
}

depth <- depth[neuropil == argvs$np & syn_type == argvs$syn_type]
depth[, ozel2021_cluster := as.character(ozel2021_cluster)]
depth[, skip_reason := trimws(as.character(skip_reason))]
depth[skip_reason == "", skip_reason := NA_character_]

selected_depth <- depth[
  is.na(skip_reason) &
    significant_fdr == TRUE &
    direction %in% c("superficial", "deep")
]

obj <- readRDS(argvs$seurat)
obj <- UpdateSeuratObject(obj)
if (!argvs$assay %in% Assays(obj)) {
  stop(sprintf("Assay %s not found in Seurat object. Available assays: %s",
               argvs$assay, paste(Assays(obj), collapse = ", ")))
}
if (!argvs$cluster_col %in% colnames(obj@meta.data)) {
  stop("Cluster metadata column not found: ", argvs$cluster_col)
}
DefaultAssay(obj) <- argvs$assay
obj[[argvs$cluster_col]] <- as.character(obj@meta.data[[argvs$cluster_col]])

all_markers <- list()
memberships <- list()

for (i in seq_len(nrow(contrast_specs))) {
  spec <- contrast_specs[i]
  contrast_data <- copy(selected_depth)
  if (!is.na(spec$Notch)) contrast_data <- contrast_data[Notch == spec$Notch]
  if (!is.na(spec$ntype)) contrast_data <- contrast_data[ntype == spec$ntype]

  sup_types <- contrast_data[direction == "superficial", cell_type]
  deep_types <- contrast_data[direction == "deep", cell_type]
  sup_clusters <- unique(contrast_data[direction == "superficial", ozel2021_cluster])
  deep_clusters <- unique(contrast_data[direction == "deep", ozel2021_cluster])

  overlap_clusters <- intersect(sup_clusters, deep_clusters)
  if (length(overlap_clusters) > 0) {
    sup_clusters <- setdiff(sup_clusters, overlap_clusters)
    deep_clusters <- setdiff(deep_clusters, overlap_clusters)
  }

  sup_types_kept <- contrast_data[
    direction == "superficial" & ozel2021_cluster %in% sup_clusters,
    cell_type
  ]
  deep_types_kept <- contrast_data[
    direction == "deep" & ozel2021_cluster %in% deep_clusters,
    cell_type
  ]

  if (length(sup_types_kept) == 0 || length(deep_types_kept) == 0) {
    memberships[[spec$contrast]] <- make_membership_row(
      spec$contrast, "skipped", "missing_type_group",
      sup_types_kept, deep_types_kept, sup_clusters, deep_clusters, overlap_clusters
    )
    next
  }

  if (length(sup_clusters) == 0 || length(deep_clusters) == 0) {
    memberships[[spec$contrast]] <- make_membership_row(
      spec$contrast, "skipped", "missing_cluster_group",
      sup_types_kept, deep_types_kept, sup_clusters, deep_clusters, overlap_clusters
    )
    next
  }

  cells_sup <- colnames(obj)[obj@meta.data[[argvs$cluster_col]] %in% sup_clusters]
  cells_deep <- colnames(obj)[obj@meta.data[[argvs$cluster_col]] %in% deep_clusters]

  if (length(cells_sup) < argvs$min_cells || length(cells_deep) < argvs$min_cells) {
    memberships[[spec$contrast]] <- make_membership_row(
      spec$contrast, "skipped", "insufficient_cells",
      sup_types_kept, deep_types_kept, sup_clusters, deep_clusters, overlap_clusters,
      length(cells_sup), length(cells_deep)
    )
    next
  }

  cells_keep <- c(cells_sup, cells_deep)
  contrast_obj <- subset(obj, cells = cells_keep)
  contrast_obj[["de_group"]] <- ifelse(
    Cells(contrast_obj) %in% cells_sup,
    "superficial",
    "deep"
  )
  Idents(contrast_obj) <- "de_group"

  markers <- tryCatch(
    FindMarkers(
      contrast_obj,
      ident.1 = "superficial",
      ident.2 = "deep",
      assay = argvs$assay,
      slot = argvs$slot,
      test.use = argvs$test_use,
      min.pct = argvs$min_pct,
      logfc.threshold = argvs$logfc_threshold
    ),
    error = function(e) e
  )

  if (inherits(markers, "error")) {
    memberships[[spec$contrast]] <- make_membership_row(
      spec$contrast, "skipped", paste0("findmarkers_failed: ", conditionMessage(markers)),
      sup_types_kept, deep_types_kept, sup_clusters, deep_clusters, overlap_clusters,
      length(cells_sup), length(cells_deep)
    )
    next
  }

  markers_dt <- as.data.table(markers, keep.rownames = "gene")
  markers_dt[, `:=`(
    stage = argvs$stage,
    neuropil = argvs$np,
    syn_type = argvs$syn_type,
    contrast = spec$contrast,
    group_1 = "superficial",
    group_2 = "deep",
    n_types_superficial = length(unique(sup_types_kept)),
    n_types_deep = length(unique(deep_types_kept)),
    n_clusters_superficial = length(unique(sup_clusters)),
    n_clusters_deep = length(unique(deep_clusters)),
    n_cells_superficial = length(cells_sup),
    n_cells_deep = length(cells_deep),
    assay = argvs$assay,
    slot = argvs$slot,
    test_use = argvs$test_use
  )]

  all_markers[[spec$contrast]] <- markers_dt
  memberships[[spec$contrast]] <- make_membership_row(
    spec$contrast, "tested", NA_character_,
    sup_types_kept, deep_types_kept, sup_clusters, deep_clusters, overlap_clusters,
    length(cells_sup), length(cells_deep)
  )
}

membership <- rbindlist(memberships, fill = TRUE)
markers <- if (length(all_markers) > 0) {
  rbindlist(all_markers, fill = TRUE)
} else {
  make_empty_markers()
}

membership_file <- sprintf("%s_%s_%s_de_membership.csv", argvs$stage, argvs$np, argvs$syn_type)
markers_file <- sprintf("%s_%s_%s_de_markers.csv", argvs$stage, argvs$np, argvs$syn_type)
fwrite(membership, membership_file)
fwrite(markers, markers_file)

cat(sprintf(
  "Wrote %s (%d contrasts) and %s (%d marker rows)\n",
  membership_file, nrow(membership), markers_file, nrow(markers)
))
