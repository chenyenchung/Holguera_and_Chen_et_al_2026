#!/usr/bin/env Rscript
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo_root <- normalizePath(file.path(dirname(script_file), "../../.."), mustWork = TRUE)
renv::load(repo_root)
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Seurat))

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

required_args <- c("stage", "ann", "seurat")
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
    "all_projection",
    "notch_on_projection",
    "notch_off_projection",
    "all_intrinsic"
  ),
  Notch = c(NA_character_, "Notch On", "Notch Off", NA_character_),
  ntype = c("projection", "projection", "projection", "intrinsic")
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
    n_types_early = integer(),
    n_types_late = integer(),
    n_clusters_early = integer(),
    n_clusters_late = integer(),
    n_cells_early = integer(),
    n_cells_late = integer(),
    assay = character(),
    slot = character(),
    test_use = character()
  )
}

make_membership_row <- function(contrast, status, skip_reason = NA_character_,
                                early_types = character(), late_types = character(),
                                early_clusters = character(), late_clusters = character(),
                                dropped_overlap = character(),
                                n_cells_early = NA_integer_, n_cells_late = NA_integer_) {
  data.table(
    stage = argvs$stage,
    neuropil = "all",
    syn_type = "temporal_cohort",
    contrast = contrast,
    status = status,
    skip_reason = skip_reason,
    early_types = collapse_chr(early_types),
    late_types = collapse_chr(late_types),
    early_clusters = collapse_chr(early_clusters),
    late_clusters = collapse_chr(late_clusters),
    dropped_overlap_clusters = collapse_chr(dropped_overlap),
    n_types_early = length(unique(early_types)),
    n_types_late = length(unique(late_types)),
    n_clusters_early = length(unique(early_clusters)),
    n_clusters_late = length(unique(late_clusters)),
    n_cells_early = n_cells_early,
    n_cells_late = n_cells_late,
    min_cells = argvs$min_cells,
    assay = argvs$assay,
    slot = argvs$slot,
    test_use = argvs$test_use,
    min_pct = argvs$min_pct,
    logfc_threshold = argvs$logfc_threshold
  )
}

ann <- fread(argvs$ann)
required_ann_cols <- c(
  "cell_type", "Confident_annotation", "ozel2021_cluster",
  "Notch", "ntype", "broad_temp", "temporal_id"
)
missing_ann_cols <- setdiff(required_ann_cols, names(ann))
if (length(missing_ann_cols) > 0) {
  stop("Missing annotation columns: ", paste(missing_ann_cols, collapse = ", "))
}

ann <- ann[
  Confident_annotation == "Y" &
    !is.na(ozel2021_cluster) &
    temporal_id != 0 &
    broad_temp %in% c("Early", "Late")
]
ann[, ozel2021_cluster := as.character(ozel2021_cluster)]

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
  contrast_data <- copy(ann)
  if (!is.na(spec$Notch)) contrast_data <- contrast_data[Notch == spec$Notch]
  contrast_data <- contrast_data[ntype == spec$ntype]

  early_clusters <- unique(contrast_data[broad_temp == "Early", ozel2021_cluster])
  late_clusters <- unique(contrast_data[broad_temp == "Late", ozel2021_cluster])

  overlap_clusters <- intersect(early_clusters, late_clusters)
  if (length(overlap_clusters) > 0) {
    early_clusters <- setdiff(early_clusters, overlap_clusters)
    late_clusters <- setdiff(late_clusters, overlap_clusters)
  }

  early_types <- contrast_data[
    broad_temp == "Early" & ozel2021_cluster %in% early_clusters,
    cell_type
  ]
  late_types <- contrast_data[
    broad_temp == "Late" & ozel2021_cluster %in% late_clusters,
    cell_type
  ]

  if (length(early_types) == 0 || length(late_types) == 0) {
    memberships[[spec$contrast]] <- make_membership_row(
      spec$contrast, "skipped", "missing_type_group",
      early_types, late_types, early_clusters, late_clusters, overlap_clusters
    )
    next
  }

  if (length(early_clusters) == 0 || length(late_clusters) == 0) {
    memberships[[spec$contrast]] <- make_membership_row(
      spec$contrast, "skipped", "missing_cluster_group",
      early_types, late_types, early_clusters, late_clusters, overlap_clusters
    )
    next
  }

  cells_early <- colnames(obj)[obj@meta.data[[argvs$cluster_col]] %in% early_clusters]
  cells_late <- colnames(obj)[obj@meta.data[[argvs$cluster_col]] %in% late_clusters]

  if (length(cells_early) < argvs$min_cells || length(cells_late) < argvs$min_cells) {
    memberships[[spec$contrast]] <- make_membership_row(
      spec$contrast, "skipped", "insufficient_cells",
      early_types, late_types, early_clusters, late_clusters, overlap_clusters,
      length(cells_early), length(cells_late)
    )
    next
  }

  cells_keep <- c(cells_early, cells_late)
  contrast_obj <- subset(obj, cells = cells_keep)
  contrast_obj[["de_group"]] <- ifelse(
    Cells(contrast_obj) %in% cells_early,
    "Early",
    "Late"
  )
  Idents(contrast_obj) <- "de_group"

  markers <- tryCatch(
    FindMarkers(
      contrast_obj,
      ident.1 = "Early",
      ident.2 = "Late",
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
      early_types, late_types, early_clusters, late_clusters, overlap_clusters,
      length(cells_early), length(cells_late)
    )
    next
  }

  markers_dt <- as.data.table(markers, keep.rownames = "gene")
  markers_dt[, `:=`(
    stage = argvs$stage,
    neuropil = "all",
    syn_type = "temporal_cohort",
    contrast = spec$contrast,
    group_1 = "Early",
    group_2 = "Late",
    n_types_early = length(unique(early_types)),
    n_types_late = length(unique(late_types)),
    n_clusters_early = length(unique(early_clusters)),
    n_clusters_late = length(unique(late_clusters)),
    n_cells_early = length(cells_early),
    n_cells_late = length(cells_late),
    assay = argvs$assay,
    slot = argvs$slot,
    test_use = argvs$test_use
  )]

  all_markers[[spec$contrast]] <- markers_dt
  memberships[[spec$contrast]] <- make_membership_row(
    spec$contrast, "tested", NA_character_,
    early_types, late_types, early_clusters, late_clusters, overlap_clusters,
    length(cells_early), length(cells_late)
  )
}

membership <- rbindlist(memberships, fill = TRUE)
markers <- if (length(all_markers) > 0) {
  rbindlist(all_markers, fill = TRUE)
} else {
  make_empty_markers()
}

membership_file <- sprintf("%s_temporal_cohort_de_membership.csv", argvs$stage)
markers_file <- sprintf("%s_temporal_cohort_de_markers.csv", argvs$stage)
fwrite(membership, membership_file)
fwrite(markers, markers_file)
