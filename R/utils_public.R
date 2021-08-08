#' @export
fastColSums <- function(mat) {
    if (is.sparseMatrix(mat) && class(mat) != "dgCMatrix") {
        mat <- as(mat, "dgCMatrix")
    }
    out <- fast_column_sums(mat)
    return(out)
}

#' @export
fastRowSums <- function(mat) {
    if (is.sparseMatrix(mat) && class(mat) != "dgCMatrix") {
        mat <- as(mat, "dgCMatrix")
    }
    out <- fast_row_sums(mat)
    return(out)
}

#' @export
fastColMeans <- function(mat) {
    E <- fastColSums(mat) / nrow(mat)
    return(E)
}

#' @export
fastRowMeans <- function(mat) {
    E <- fastRowSums(mat) / ncol(mat)
    return(E)
}

#' @export
fastRowVars <- function(mat) {
    mat <- as(mat, "dgTMatrix")
    E <- fastRowMeans(mat)
    V <- computeSparseRowVariances(mat@i + 1, mat@x, E, ncol(mat))
    return(V)
}

#' @export
fastBindSparse <- function(A, B, d = 1) {
    d <- d - 1
    sp_mat <- bind_sparse_mats(A, B, d)
    return(sp_mat)
}

is.sparseMatrix <- function(A) {
    return(length(which(is(A) == "sparseMatrix")) != 0)
}

#' @export
revert_ace_as_sce <- function(ace) {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = SummarizedExperiment::assays(ace),
        colData = SummarizedExperiment::colData(ace),
        rowData = SummarizedExperiment::rowData(ace)
    )

    return(sce)
}