.check_and_load_package <- function(pkg_names) {
    for (pk in pkg_names) {
        if (!require(pk, character.only = T)) {
            err <- sprintf("Package '%s' is not installed.\n", pk)
            stop(err)
        }
    }
}

.tscalet <- function(A, center = TRUE, scale = TRUE) {
    A <- Matrix::t(scale(Matrix::t(A), center = center, scale = scale))
    return(A)
}
