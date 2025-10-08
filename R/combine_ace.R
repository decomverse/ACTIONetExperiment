#' @importFrom SummarizedExperiment rbind cbind
#' @export
setMethod("cbind", "ACTIONetExperiment", function(..., deparse.level = 1) {
    args <- list(...)
    ## ace slots are dropped by default
    drop <- if ("drop" %in% names(match.call())) eval(match.call()[["drop"]], parent.frame()) else TRUE
    if (!is.logical(drop) || length(drop) != 1) stop("Argument 'drop' must be a single logical value.")
    # Always clear colNets and rowNets if not empty, and colMaps/rowMaps if drop=TRUE
    args <- .ace_clear_slots(args, drop = drop)
    if (!drop) {
        slot_names <- c("colMaps", "rowMaps")
        slot_lists <- sapply(slot_names, function(slot) lapply(args, function(x) methods::slot(x, slot)), simplify = FALSE)
        .validate_bind(slot_lists, mode = "cbind")
        concat_vals <- .ace_bind_slots(slot_lists, slot_names, mode = "cbind")
        new_args <- args
        for (i in seq_along(args)) {
            for (j in seq_along(slot_names)) {
                methods::slot(new_args[[i]], slot_names[j]) <- concat_vals[[j]]
            }
        }
        out <- do.call(callNextMethod, new_args)
    } else {
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
    # Always clear colNets and rowNets if not empty, and colMaps/rowMaps if drop=TRUE
    args <- .ace_clear_slots(args, drop = drop)
    if (!drop) {
        slot_names <- c("colMaps", "rowMaps")
        slot_lists <- lapply(slot_names, function(slot) lapply(args, function(x) methods::slot(x, slot)))
        .validate_bind(slot_lists, mode = "rbind")
        concat_vals <- .ace_bind_slots(slot_lists, slot_names, mode = "rbind")
        new_args <- args
        for (i in seq_along(args)) {
            for (j in seq_along(slot_names)) {
                methods::slot(new_args[[i]], slot_names[j]) <- concat_vals[[j]]
            }
        }
        out <- do.call(callNextMethod, new_args)
    } else {
        out <- do.call(callNextMethod, args)
    }
    return(out)
})

.slot_is_matrix_like <- function(x) {
    inherits(x, c("matrix", "array")) || is.sparseMatrix(x) || (is.data.frame(x) && length(dim(x)) == 2)
}

.validate_bind <- function(slot_lists, mode = c("cbind", "rbind")) {
    mode <- match.arg(mode)
    for (i in seq_along(slot_lists)) {
        slot_content <- slot_lists[[i]]
        slot_name <- names(slot_lists)[i]
        slot_class <- unique(sapply(slot_content, function(x) class(x)))
        if (length(slot_class) > 1) stop("All objects in must have the same class for concatenation.")
        if (!all(sapply(slot_content, .slot_is_matrix_like))) {
            err <- sprintf("Objects in slot '%s' must be matrix-like (matrix, array, or data.frame).", slot_name)
            stop(err)
        }
        if (mode == "cbind") {
            if (slot_name == "colMaps") {
                # rbind: check columns match
                ncol_list <- sapply(slot_content, function(x) if (is.null(dim(x))) NA else dim(x)[2])
                if (any(is.na(ncol_list))) {
                    err <- sprintf("Objects in slot '%s' must have defined number of columns.", slot_name)
                    stop(err)
                }
                if (length(unique(ncol_list)) > 1) {
                    err <- sprintf("Objects in slot '%s' must have the same number of columns.", slot_name)
                    stop(err)
                }
            } else if (slot_name == "rowMaps") {
                # cbind: check rows match
                nrow_list <- sapply(slot_content, function(x) if (is.null(dim(x))) NA else dim(x)[1])
                if (any(is.na(nrow_list))) {
                    err <- sprintf("Objects in slot '%s' must have defined number of rows.", slot_name)
                    stop(err)
                }
                if (length(unique(nrow_list)) > 1) {
                    err <- sprintf("Objects in slot '%s' must have the same number of rows.", slot_name)
                    stop(err)
                }
            }
        } else {
            if (slot_name == "colMaps") {
                # cbind: check rows match
                nrow_list <- sapply(slot_content, function(x) if (is.null(dim(x))) NA else dim(x)[1])
                if (any(is.na(nrow_list))) {
                    err <- sprintf("Objects in slot '%s' must have defined number of rows.", slot_name)
                    stop(err)
                }
                if (length(unique(nrow_list)) > 1) {
                    err <- sprintf("Objects in slot '%s' must have the same number of rows.", slot_name)
                    stop(err)
                }
            } else if (slot_name == "rowMaps") {
                # rbind: check columns match
                ncol_list <- sapply(slot_content, function(x) if (is.null(dim(x))) NA else dim(x)[2])
                if (any(is.na(ncol_list))) {
                    err <- sprintf("Objects in slot '%s' must have defined number of columns.", slot_name)
                    stop(err)
                }
                if (length(unique(ncol_list)) > 1) {
                    err <- sprintf("Objects in slot '%s' must have the same number of columns.", slot_name)
                    stop(err)
                }
            }
        }
    }
}


.ace_bind_slots <- function(slot_lists, slot_names, mode = c("cbind", "rbind")) {
    mode <- match.arg(mode)
    out <- vector("list", length(slot_lists))
    names(out) <- slot_names
    for (i in seq_along(slot_names)) {
        slot_list <- slot_lists[[i]]
        if (all(sapply(slot_list, function(x) length(x) == 0))) {
            out[[i]] <- slot_list[[1]]
        } else {
            if (mode == "cbind") {
                # colMaps: rbind, rowMaps: cbind
                if (slot_names[i] == "colMaps") {
                    out[[i]] <- do.call(BiocGenerics::rbind, slot_list)
                } else if (slot_names[i] == "rowMaps") {
                    out[[i]] <- do.call(BiocGenerics::cbind, slot_list)
                }
            } else {
                # rbind: colMaps: cbind, rowMaps: rbind
                if (slot_names[i] == "colMaps") {
                    out[[i]] <- do.call(BiocGenerics::cbind, slot_list)
                } else if (slot_names[i] == "rowMaps") {
                    out[[i]] <- do.call(BiocGenerics::rbind, slot_list)
                }
            }
        }
    }
    out
}

.ace_clear_slots <- function(args, drop = TRUE) {
    # colNets and rowNets are always dropped if not empty
    always_drop <- c("colNets", "rowNets")
    conditional_drop <- c("colMaps", "rowMaps")
    dropped <- character(0)
    for (obj in args) {
        for (slot in always_drop) {
            slot_val <- methods::slot(obj, slot)
            if (length(slot_val) > 0) dropped <- union(dropped, slot)
        }
        if (drop) {
            for (slot in conditional_drop) {
                slot_val <- methods::slot(obj, slot)
                if (length(slot_val) > 0) dropped <- union(dropped, slot)
            }
        }
    }
    if (length(dropped) > 0) {
        message(sprintf("Slots dropped: %s", paste(dropped, collapse = ", ")))
    }
    nc_rep <- S4Vectors::SimpleList()
    lapply(args, function(a) {
        BiocGenerics:::replaceSlots(a,
            rowNets = nc_rep,
            colNets = nc_rep,
            rowMaps = if (drop) nc_rep else methods::slot(a, "rowMaps"),
            colMaps = if (drop) nc_rep else methods::slot(a, "colMaps"),
            check = FALSE)
    })
}
