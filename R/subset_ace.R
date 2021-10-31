#' Subsets rows/columns of an ACTIONetExperiment (ACE) object
#'
#' @param i,j rows and columns
#'
#' @export
setMethod("[", c("ACTIONetExperiment", "ANY", "ANY"), function(x, i, j, ..., drop = TRUE) {
  
    rnets <- x@rowNets
    cnets <- x@colNets
    rmaps <- x@rowMaps
    cmaps <- x@colMaps

    if (!missing(i)) {
        if (is.character(i)) {
            fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
            i <- SummarizedExperiment:::.SummarizedExperiment.charbound(i, rownames(x),
                fmt)
        }
        i <- as.vector(i)

        if (length(rnets) > 0) {
            for (k in 1:length(rnets)) {
                tmp = rnets[[k]]
                rnets[[k]] = tmp[i, i]
            }
        }
        if (length(rmaps) > 0) {
            for (k in 1:length(rmaps)) {
                tmp = rmaps[[k]]
                rmaps[[k]] = tmp[i, , drop = FALSE]
            }
        }
    }

    if (!missing(j)) {
        if (is.character(j)) {
            fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
            j <- SummarizedExperiment:::.SummarizedExperiment.charbound(j, colnames(x),
                fmt)
        }
        j <- as.vector(j)

        if (length(cnets) > 0) {
            for (k in 1:length(cnets)) {
                tmp = cnets[[k]]
                cnets[[k]] = tmp[j, j]
            }
        }
        if (length(cmaps) > 0) {
            for (k in 1:length(cmaps)) {
                tmp = cmaps[[k]]
                cmaps[[k]] = tmp[j, , drop = FALSE]
            }
        }
    }

    out <- callNextMethod()
    BiocGenerics:::replaceSlots(out, rowNets = rnets, colNets = cnets, rowMaps = rmaps, colMaps = cmaps,
        check = FALSE)
})
