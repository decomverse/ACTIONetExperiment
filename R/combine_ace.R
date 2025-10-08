#' @importFrom SummarizedExperiment rbind cbind
#' @export
setMethod("cbind", c("ACTIONetExperiment"), function(..., deparse.level = 1) {
    args <- list(...)
    args <- .ace_clear_slots(args, slots = c("colMaps", "rowMaps", "colNets", "rowNets"))
    out <- do.call(callNextMethod, args)
    out
})

#' @importFrom SummarizedExperiment rbind cbind
#' @export
setMethod("rbind", "ACTIONetExperiment", function(..., deparse.level = 1) {
    args <- list(...)
    args <- .ace_clear_slots(args, slots = c("colMaps", "rowMaps", "colNets", "rowNets"))
    out <- do.call(callNextMethod, args)
    out
})



.ace_clear_slots <- function(ace_list, slots = c("colNets", "rowNets")) {
    # Clear specified slots if not empty
    dropped <- character(0)
    for (obj in ace_list) {
        for (slot in slots) {
            slot_val <- methods::slot(obj, slot)
            if (length(slot_val) > 0) dropped <- union(dropped, slot)
        }
    }
    if (length(dropped) > 0) {
        message(sprintf("Slots dropped: %s", paste(dropped, collapse = ", ")))
    }
    nc_rep <- S4Vectors::SimpleList()
    lapply(ace_list, function(obj) {
        slot_list <- list()
        for (slot in slots) {
            slot_list[[slot]] <- nc_rep
        }
        do.call(BiocGenerics:::replaceSlots, c(list(obj), slot_list, list(check = FALSE)))
    })
}
