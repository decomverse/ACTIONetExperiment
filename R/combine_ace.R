.slot_is_matrix_like <- function(x) {
    inherits(x, c("matrix", "array")) || is.sparseMatrix(x) || (is.data.frame(x) && length(dim(x)) == 2)
}

.validate_bind <- function(slot_lists, mode = c("cbind", "rbind")) {
    mode <- match.arg(mode)
    for (i in seq_along(slot_lists)) {
        slot_content <- slot_lists[[i]]
        slot_class <- unique(sapply(slot_content, function(x) class(x)))
        if (length(slot_class) > 1) stop("All slots must have the same class for concatenation.")
        if (!all(sapply(slot_content, .slot_is_matrix_like))) {
            stop("All slots must be matrix-like (matrix, array, or data.frame).")
        }
        if (mode == "cbind") {
            ncol_list <- sapply(slot_content, function(x) if (is.null(dim(x))) NA else dim(x)[2])
            if (any(is.na(ncol_list))) stop("All slots must have defined number of columns.")
            if (length(unique(ncol_list)) > 1) stop("All slots must have the same number of columns for cbind.")
        } else {
            nrow_list <- sapply(slot_content, function(x) if (is.null(dim(x))) NA else dim(x)[1])
            if (any(is.na(nrow_list))) stop("All slots must have defined number of rows.")
            if (length(unique(nrow_list)) > 1) stop("All slots must have the same number of rows for rbind.")
        }
    }
}

.ace_bind_slots <- function(slot_lists, mode = c("cbind", "rbind")) {
    mode <- match.arg(mode)
    concat_fun <- if (mode == "cbind") BiocGenerics::cbind else BiocGenerics::rbind
    lapply(slot_lists, function(slot_list) {
        if (all(sapply(slot_list, function(x) length(x) == 0))) return(slot_list[[1]])
        do.call(concat_fun, slot_list)
    })
}

#' @importFrom SummarizedExperiment rbind cbind
#' @export
setMethod("cbind", "ACTIONetExperiment", function(..., deparse.level = 1) {
    args <- list(...)
    ## ace slots are dropped by default
    drop <- if ("drop" %in% names(match.call())) eval(match.call()[["drop"]], parent.frame()) else TRUE
    if (!is.logical(drop) || length(drop) != 1) stop("Argument 'drop' must be a single logical value.")
    if (!drop) {
        slot_names <- c("colMaps", "rowMaps", "colNets", "rowNets")
        slot_lists <- lapply(slot_names, function(slot) lapply(args, function(x) methods::slot(x, slot)))
        .validate_bind(slot_lists, mode = "cbind")
        concat_vals <- .ace_bind_slots(slot_lists, mode = "cbind")
        new_args <- args
        for (i in seq_along(args)) {
            for (j in seq_along(slot_names)) {
                methods::slot(new_args[[i]], slot_names[j]) <- concat_vals[[j]]
            }
        }
        out <- do.call(callNextMethod, new_args)
    } else {
        args <- .ace_clear_slots(args)
        out <- do.call(callNextMethod, args)
    }
    return(out)
})

#' @importFrom SummarizedExperiment rbind cbind
#' @export
setMethod("rbind", "ACTIONetExperiment", function(..., deparse.level = 1) {
    args <- list(...)
    ## ace slots are dropped by default
    drop <- if ("drop" %in% names(match.call())) eval(match.call()[["drop"]], parent.frame()) else TRUE
    if (!is.logical(drop) || length(drop) != 1) stop("Argument 'drop' must be a single logical value.")
    if (!drop) {
        slot_names <- c("colMaps", "rowMaps", "colNets", "rowNets")
        slot_lists <- lapply(slot_names, function(slot) lapply(args, function(x) methods::slot(x, slot)))
        .validate_bind(slot_lists, mode = "rbind")
        concat_vals <- .ace_bind_slots(slot_lists, mode = "rbind")
        new_args <- args
        for (i in seq_along(args)) {
            for (j in seq_along(slot_names)) {
                methods::slot(new_args[[i]], slot_names[j]) <- concat_vals[[j]]
            }
        }
        out <- do.call(callNextMethod, new_args)
    } else {
        args <- .ace_clear_slots(args)
        out <- do.call(callNextMethod, args)
    }
    return(out)
})


clear_ace_slots <- function(args) {
    slot_names <- c("colMaps", "rowMaps", "colNets", "rowNets")
    dropped <- character(0)
    for (obj in args) {
        for (slot in slot_names) {
            slot_val <- methods::slot(obj, slot)
            if (length(slot_val) > 0) dropped <- union(dropped, slot)
        }
    }
    if (length(dropped) > 0) {
        message(sprintf("Slots dropped: %s", paste(dropped, collapse = ", ")))
    }
    nc_rep <- S4Vectors::SimpleList()
    lapply(args, function(a) {
        BiocGenerics:::replaceSlots(a, rowNets = nc_rep, colNets = nc_rep, rowMaps = nc_rep,
            colMaps = nc_rep, check = FALSE)
    })
}
