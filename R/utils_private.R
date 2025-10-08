## Default ACE values
.default_rowData <- function(d) {
    DF <- DataFrame(feature_name = paste0("feat_", 1:d))
    return(DF)
}

.default_colData <- function(d) {
    DF <- DataFrame(sample_name = paste0("sam_", 1:d))
    return(DF)
}

.default_rownames <- function(d) {
    n <- paste0("feat_", 1:d)
    return(n)
}

.default_colnames <- function(d) {
    n <- paste0("sam_", 1:d)
    return(n)
}

.make_chars_unique <- function(x) {
    x <- as.character(x)
    make.unique(x, sep = "_")
    return(x)
}

## Misc
.check_and_load_package <- function(pkg_names) {
    for (pk in pkg_names) {
        if (!require(pk, character.only = T)) {
            err <- sprintf("Package '%s' is not installed.\n", pk)
            stop(err)
        }
    }
}
