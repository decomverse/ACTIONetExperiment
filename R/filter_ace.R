#' @export
get_mtRNA_stats <- function(ace, by = NULL, groups_use = NULL, assay = "counts", species = c("mmusculus", "hsapiens", "other"), metric = c("pct", "ratio", "counts", "frac")) {
  require(stats)
  species <- match.arg(species)
  metric <- match.arg(metric)

  mask <- grepl("^MT[:.:]|^MT-", rownames(ace), ignore.case = TRUE)

  mat <- assays(ace)[[assay]]
  cs_mat <- fastColSums(mat)
  mm <- mat[mask, , drop = F]
  cs_mm <- fastColSums(mm)

  if (!is.null(by)) {
    IDX <- .get_attr_or_split_idx(ace, by, groups_use)

    if (metric == "frac") {
      frac.list <- lapply(IDX, function(idx) {
        m <- cs_mm[idx] / cs_mat[idx]
      })
    } else if (metric == "pct") {
      frac.list <- lapply(IDX, function(idx) {
        m <- 100 * cs_mm[idx] / cs_mat[idx]
      })
    } else if (metric == "ratio") {
      frac.list <- lapply(IDX, function(idx) {
        m <- cs_mm[idx] / (cs_mat[idx] - cs_mm[idx])
      })
    } else {
      frac.list <- lapply(IDX, function(idx) {
        m <- cs_mm[idx]
      })
    }
    return(frac.list)
  } else {
    if (metric == "pct") {
      frac <- 100 * cs_mm / cs_mat
    } else if (metric == "frac") {
      frac <- cs_mm / cs_mat
    } else if (metric == "ratio") {
      frac <- cs_mm / (cs_mat - cs_mm)
    } else {
      frac <- cs_mm
    }
    return(frac)
  }
}

#' Filter columns and rows of `ACTIONetExperiment` or `SummarizedExperiment`-like object.
#' @export
filter.ace <- function(ace,
                       assay_name = "counts",
                       min_cells_per_feat = 0.003,
                       min_feats_per_cell = 1000,
                       min_umis_per_cell = NULL,
                       max_umis_per_cell = 50000,
                       max_mito_fraction = 10,
                       species = "hsapiens",
                       features_use = NULL,
                       return_fil_ace = TRUE) {
  filter_stats <- list()

  init_dim <- dim(ace)
  init_rn <- rownames(ace)
  init_cn <- colnames(ace)
  # ace = ace

  if (!is.null(max_mito_fraction)) {
    mt_frac <- get_mtRNA_stats(
      ace,
      by = NULL,
      groups_use = NULL,
      assay = assay_name,
      species = species,
      metric = "frac"
    )

    mask <- mt_frac <= max_mito_fraction
    filter_stats$mito_perc <- sum(!mask)
    ace <- ace[, mt_frac <= max_mito_fraction]
  }

  i <- 0
  repeat {
    prev_dim <- dim(ace)
    rows_mask <- rep(TRUE, NROW(ace))
    cols_mask <- rep(TRUE, NCOL(ace))
    if (!is.null(min_umis_per_cell)) {
      umi_mask <- fastColSums(SummarizedExperiment::assays(ace)[[assay_name]]) >= min_umis_per_cell
      filter_stats[[sprintf("Round%d_min_umis", i + 1)]] <- sum(!umi_mask)
      cols_mask <- cols_mask & umi_mask
    }

    if (!is.null(max_umis_per_cell)) {
      umi_mask <- fastColSums(SummarizedExperiment::assays(ace)[[assay_name]]) <= max_umis_per_cell
      filter_stats[[sprintf("Round%d_max_umis", i + 1)]] <- sum(!umi_mask)
      cols_mask <- cols_mask & umi_mask
    }

    if (!is.null(min_feats_per_cell)) {
      feature_mask <- fastColSums(SummarizedExperiment::assays(ace)[[assay_name]] > 0) >= min_feats_per_cell
      filter_stats[[sprintf("Round%d_min_feats", i + 1)]] <- sum(!feature_mask)
      cols_mask <- cols_mask & feature_mask
    }

    if (!is.null(min_cells_per_feat)) {
      if ((min_cells_per_feat < 1) & (min_cells_per_feat > 0)) {
        min_fc <- min_cells_per_feat * init_dim[2]
      } else {
        min_fc <- min_cells_per_feat
      }
      cell_count_mask <- fastRowSums(SummarizedExperiment::assays(ace)[[assay_name]] > 0) >= min_fc
      filter_stats[[sprintf("Round%d_min_cells", i + 1)]] <- sum(!cell_count_mask)
      rows_mask <- rows_mask & cell_count_mask
    }

    filter_stats[[sprintf("Round%d_combined_features", i + 1)]] <- sum(!rows_mask)
    filter_stats[[sprintf("Round%d_combined_cells", i + 1)]] <- sum(!cols_mask)
    ace <- ace[rows_mask, cols_mask]

    invisible(gc())
    i <- i + 1
    if (all(dim(ace) == prev_dim)) {
      break
    }
  }
  invisible(gc())

  if (return_fil_ace) {
    metadata(ace)$filter_stats <- filter_stats
    return(ace)
  } else {
    # fil_cols_mask = !(colnames(ace) %in% colnames(ace))
    # fil_rows_mask = !(rownames(ace) %in% rownames(ace))
    fil_cols_mask <- !(init_cn %in% colnames(ace))
    fil_rows_mask <- !(init_rn %in% rownames(ace))

    fil_cols_list <- data.frame(
      name = init_cn[fil_cols_mask],
      idx = which(fil_cols_mask)
    )

    fil_rows_list <- data.frame(
      name = init_rn[fil_rows_mask],
      idx = which(fil_rows_mask)
    )

    fil_list <- list(
      cols_filtered = fil_cols_list,
      rows_filtered = fil_rows_list
    )

    return(fil_list)
  }
}

#' Filter columns and rows of `ACTIONetExperiment` or `SummarizedExperiment` object by column attribute.
#' @export
filter.ace.by.attr <- function(ace,
                               by,
                               assay_name = "counts",
                               min_cells_per_feat = 0.003,
                               min_feats_per_cell = 1000,
                               min_umis_per_cell = NULL,
                               max_umis_per_cell = 50000,
                               max_mito_fraction = 10,
                               species = "hsapiens") {
  IDX <- .get_attr_or_split_idx(ace, by)

  if (any(duplicated(rownames(ace)))) {
    msg <- sprintf("Adding suffix to duplicate rownames.\n")
    warning(msg)
    rownames(ace) <- make.unique(rownames(ace))
  }
  if (any(duplicated(colnames(ace)))) {
    msg <- sprintf("Adding suffix to duplicate colnames.\n")
    warning(msg)
    colnames(ace) <- make.unique(colnames(ace))
  }

  fil_names <- lapply(IDX, function(idx) {
    fil_list <- filter.ace(
      ace = ace[, idx],
      assay_name = assay_name,
      min_cells_per_feat = min_cells_per_feat,
      min_umis_per_cell = min_umis_per_cell,
      max_umis_per_cell = max_umis_per_cell,
      min_feats_per_cell = min_feats_per_cell,
      max_mito_fraction = max_mito_fraction,
      species = species,
      return_fil_ace = FALSE
    )

    return(fil_list)
  })

  fil_col <- lapply(fil_names, function(i) i[["cols_filtered"]]$name)
  fil_col <- Reduce(union, fil_col)


  fil_row <- lapply(fil_names, function(i) i[["rows_filtered"]]$name)
  fil_row <- Reduce(union, fil_row)

  keep_row <- which(!(rownames(ace) %in% fil_row))
  keep_col <- which(!(colnames(ace) %in% fil_col))

  filter_stats <- list(kept_row = keep_row, kept_cols = keep_col)
  metadata(ace)$filter_stats_by_attr <- filter_stats

  # ace.fil = ace[keep_row, keep_col]
  # colData(ace.fil) <- droplevels(colData(ace.fil))
  # rowData(ace.fil) <- droplevels(rowData(ace.fil))

  ace <- ace[keep_row, keep_col]
  colData(ace) <- droplevels(colData(ace))
  rowData(ace) <- droplevels(rowData(ace))

  invisible(gc())
  return(ace)
}