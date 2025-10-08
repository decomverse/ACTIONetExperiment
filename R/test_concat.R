## Aspirations for something liek AnnData.concat()

# .ace_call_bind_maps <- function(ace_list, mode) {
#     slot_names <- c("colMaps", "rowMaps")
#     slot_lists <- sapply(slot_names, function(slot) lapply(ace_list, function(x) methods::slot(x, slot)), simplify = FALSE)
#     # Only validate and bind if at least one element is non-empty
#     for (nm in slot_names) {
#         slot_list <- slot_lists[[nm]]
#         if (!all(sapply(slot_list, function(x) length(x) == 0))) {
#             slot_lists[[nm]] <- slot_list
#         } else {
#             slot_lists[[nm]] <- slot_list # keep as is, will be ignored in binding
#         }
#     }
#     .validate_map_bind(slot_lists, mode = mode)
#     new_slots <- .ace_bind_maps(slot_lists, mode = mode)
#     new_slots
# }

# .ace_bind_maps <- function(slot_lists, mode = c("cbind", "rbind")) {
#     mode <- match.arg(mode)
#     out <- vector("list", length(slot_lists))
#     names(out) <- slot_names <- names(slot_lists)
#     for (i in seq_along(slot_names)) {
#         slot_list <- slot_lists[[i]]
#         if (all(sapply(slot_list, function(x) length(x) == 0))) {
#             out[[i]] <- slot_list[[1]]
#         } else {
#             if (mode == "cbind") {
#                 # colMaps: rbind, rowMaps: cbind
#                 if (slot_names[i] == "colMaps") {
#                     out[[i]] <- do.call(BiocGenerics::rbind, slot_list)
#                 } else if (slot_names[i] == "rowMaps") {
#                     out[[i]] <- do.call(BiocGenerics::cbind, slot_list)
#                 }
#             } else {
#                 # rbind: colMaps: cbind, rowMaps: rbind
#                 if (slot_names[i] == "colMaps") {
#                     out[[i]] <- do.call(BiocGenerics::cbind, slot_list)
#                 } else if (slot_names[i] == "rowMaps") {
#                     out[[i]] <- do.call(BiocGenerics::rbind, slot_list)
#                 }
#             }
#         }
#     }
#     out
# }

# ## Doesn't work: "slot_content" is recognized as a matrix-like, instead of a list of matrix-likes.
# .validate_map_bind <- function(slot_lists, mode = c("cbind", "rbind")) {
#     mode <- match.arg(mode)
#     for (i in seq_along(slot_lists)) {
#         slot_content <- slot_lists[[i]]
#         slot_name <- names(slot_lists)[i]
#         # If all are empty, skip validation for this slot
#         if (all(sapply(slot_content, function(x) length(x) == 0))) next
#         slot_class <- unique(sapply(slot_content, function(x) class(x)))
#         if (length(slot_class) > 1) stop("All objects in must have the same class for concatenation.")
#         if (!all(sapply(slot_content, .slot_is_matrix_like))) {
#             err <- sprintf("Objects in slot '%s' must be matrix-like (matrix, array, or data.frame).", slot_name)
#             stop(err)
#         }
#         if (mode == "cbind") {
#             if (slot_name == "colMaps") {
#                 # rbind: check columns match
#                 ncol_list <- sapply(slot_content, function(x) if (is.null(dim(x))) NA else dim(x)[2])
#                 if (any(is.na(ncol_list))) {
#                     err <- sprintf("Objects in slot '%s' must have defined number of columns.", slot_name)
#                     stop(err)
#                 }
#                 if (length(unique(ncol_list)) > 1) {
#                     err <- sprintf("Objects in slot '%s' must have the same number of columns.", slot_name)
#                     stop(err)
#                 }
#             } else if (slot_name == "rowMaps") {
#                 # cbind: check rows match
#                 nrow_list <- sapply(slot_content, function(x) if (is.null(dim(x))) NA else dim(x)[1])
#                 if (any(is.na(nrow_list))) {
#                     err <- sprintf("Objects in slot '%s' must have defined number of rows.", slot_name)
#                     stop(err)
#                 }
#                 if (length(unique(nrow_list)) > 1) {
#                     err <- sprintf("Objects in slot '%s' must have the same number of rows.", slot_name)
#                     stop(err)
#                 }
#             }
#         } else {
#             if (slot_name == "colMaps") {
#                 # cbind: check rows match
#                 nrow_list <- sapply(slot_content, function(x) if (is.null(dim(x))) NA else dim(x)[1])
#                 if (any(is.na(nrow_list))) {
#                     err <- sprintf("Objects in slot '%s' must have defined number of rows.", slot_name)
#                     stop(err)
#                 }
#                 if (length(unique(nrow_list)) > 1) {
#                     err <- sprintf("Objects in slot '%s' must have the same number of rows.", slot_name)
#                     stop(err)
#                 }
#             } else if (slot_name == "rowMaps") {
#                 # rbind: check columns match
#                 ncol_list <- sapply(slot_content, function(x) if (is.null(dim(x))) NA else dim(x)[2])
#                 if (any(is.na(ncol_list))) {
#                     err <- sprintf("Objects in slot '%s' must have defined number of columns.", slot_name)
#                     stop(err)
#                 }
#                 if (length(unique(ncol_list)) > 1) {
#                     err <- sprintf("Objects in slot '%s' must have the same number of columns.", slot_name)
#                     stop(err)
#                 }
#             }
#         }
#     }
# }


# .slot_is_matrix_like <- function(x) {
#     inherits(x, c("matrix", "array")) || is.sparseMatrix(x) || (is.data.frame(x) && length(dim(x)) == 2)
# }
