#' @export
fastColSums <- function(mat) {
    if (is.sparseMatrix(mat) && class(mat) != "dMatrix") {
        mat <- as(mat, "dMatrix")
    }
    out <- fast_column_sums(mat)
    return(out)
}

#' @export
fastRowSums <- function(mat) {
    if (is.sparseMatrix(mat) && class(mat) != "dMatrix") {
        mat <- as(mat, "dMatrix")
    }
    out <- fast_row_sums(mat)
    return(out)
}

#' @export
fastColMeans <- function(mat) {
    E <- Matrix::colSums(mat) / nrow(mat)
    return(E)
}

#' @export
fastRowMeans <- function(mat) {
    E <- fastRowSums(mat) / ncol(mat)
    return(E)
}

#' @export
is.sparseMatrix <- function(A) {
    is_sparse = any(is(A) == "sparseMatrix")
    return(is_sparse)
}

#' @export
get.data.or.split <- function(
  ace,
  attr,
  groups_use = NULL,
  to_return = c("data", "levels", "split"),
  d = 2
) {

    to_return = match.arg(to_return)

    if (length(attr) == 1) {
        data_vec <- switch(d,
            SummarizedExperiment::rowData(ace)[[attr]],
            SummarizedExperiment::colData(ace)[[attr]]
        )
    } else {
        if (length(attr) != dim(ace)[d]) {
            err <- sprintf("'attr' length does not match %s of ace.\n", ifelse(d ==
                1, "NROW", "NCOL"))
            stop(err)
        }
        data_vec <- attr
    }

    if (is.null(data_vec)) {
        stop(sprintf("Invalid split conditions.\n"))
    } else {
      if (is.factor(data_vec)){
        data_vec = as.character(data_vec)
      }
        data_vec = data_vec
    }

    idx <- 1:dim(ace)[d]

    ## Ignores 'to_return'. Always returns index list.
    if (!is.null(groups_use)) {
        sub_idx <- which(data_vec %in% groups_use)
        data_vec <- data_vec[sub_idx]
        if (is.null(data_vec)) {
              stop(sprintf("Invalid split conditions.\n"))
          }
        idx_list <- split(idx[sub_idx], data_vec)
        return(idx_list)
    } else {
        idx_list <- split(idx, data_vec)
    }

    if (to_return == "data") {
        return(data_vec)
    } else if (to_return == "levels"){
        fac_vec = factor(data_vec)
        level.list = list(index = as.numeric(fac_vec), keys = levels(fac_vec))
        return(level.list)
    } else if (to_return == "split"){
        return(idx_list)
    }
}
