#' @import hdf5r
h5addAttr.str <- function(h5group, attr.name, attr.val) {
    dtype <- H5T_STRING$new(type = "c", size = Inf)
    dtype <- dtype$set_cset(cset = "unknown")

    space <- H5S$new(type = "scalar")
    h5group$create_attr(attr_name = attr.name, dtype = dtype, space = space)
    attr <- h5group$attr_open_by_name(attr_name = attr.name, ".")
    attr$write(attr.val)
}

#' @import hdf5r
h5addAttr.str_array <- function(h5group, attr.name, attr.val) {
    dtype <- H5T_STRING$new(type = "c", size = Inf)
    dtype <- dtype$set_cset(cset = "unknown")

    space <- H5S$new(type = "simple", dims = length(attr.val), maxdims = length(attr.val))
    h5group$create_attr(attr_name = attr.name, dtype = dtype, space = space)
    attr <- h5group$attr_open_by_name(attr_name = attr.name, ".")
    attr$write(attr.val)
}

#' @import hdf5r
write.HD5DF <- function(h5file,
                        gname,
                        DF,
                        compression_level = 0) {
    string.dtype <- H5T_STRING$new(type = "c", size = Inf)
    string.dtype <- string.dtype$set_cset(cset = "unknown")

    DF <- as.data.frame(DF)

    N <- NROW(DF)

    h5group <- h5file$create_group(gname)

    h5addAttr.str(h5group, "_index", "index")
    h5addAttr.str(h5group, "encoding-version", "0.1.0")
    h5addAttr.str(h5group, "encoding-type", "dataframe")

    if (0 < NCOL(DF)) {
        # noncat.num.vars <- which(apply(DF, 2, is.numeric))
        noncat.num.vars <- which(sapply(names(DF), function(nn) is.numeric(DF[[nn]])))
        cat.vars <- which(sapply(names(DF), function(nn) length(unique(DF[[nn]])) < 128))
        # cat.vars <- which(apply(DF, 2, function(x) length(unique(x)) < 128))
        cat.vars <- setdiff(cat.vars, noncat.num.vars)
        noncat.vars <- setdiff(1:NCOL(DF), c(cat.vars, noncat.num.vars))
        cn <- colnames(DF)
        catDF <- DF[, cat.vars, drop = FALSE]
        catDF <- apply(catDF, 2, as.character)
        catDF[is.na(catDF)] <- "NA"

        numDF <- DF[, noncat.num.vars, drop = FALSE]
        numDF <- apply(numDF, 2, as.numeric)
        numDF[is.na(numDF)] <- NA

        nonNumDF <- DF[, noncat.vars, drop = FALSE]
        nonNumDF <- apply(nonNumDF, 2, as.character)
        nonNumDF[is.na(nonNumDF)] <- NA

        if (length(cn) == 0) {
            dtype <- H5T_STRING$new(type = "c", size = Inf)
            dtype <- dtype$set_cset(cset = "unknown")
            space <- H5S$new(type = "simple", dims = 0, maxdims = 10)

            h5group$create_attr(attr_name = "column-order", dtype = dtype, space = space)
        } else {
            h5addAttr.str_array(h5group, "column-order", cn)
        }

        if (length(cat.vars) > 0) {
            cat <- h5group$create_group("__categories")

            for (i in 1:length(cat.vars)) {
                x <- catDF[, i]
                l <- sort(unique(x))
                v <- match(x, l) - 1

                dtype <- H5T_STRING$new(type = "c", size = Inf)
                dtype <- dtype$set_cset(cset = "unknown")
                l.enum <- cat$create_dataset(colnames(DF)[cat.vars[i]], l,
                    gzip_level = compression_level,
                    dtype = dtype
                )


                dtype <- H5T_ENUM$new(labels = c("FALSE", "TRUE"), values = 0:1)
                space <- H5S$new(type = "scalar")
                res <- l.enum$create_attr(attr_name = "ordered", dtype = dtype, space = space)

                attr <- l.enum$attr_open_by_name(attr_name = "ordered", ".")
                attr$write(0)

                l.vec <- h5group$create_dataset(colnames(DF)[cat.vars[i]], as.integer(v),
                    gzip_level = compression_level, dtype = h5types$H5T_NATIVE_INT8
                )

                ref <- cat$create_reference(name = colnames(DF)[cat.vars[i]])

                dtype <- guess_dtype(ref)
                space <- H5S$new(type = "scalar")
                res <- l.vec$create_attr(
                    attr_name = "categories", dtype = dtype,
                    space = space
                )
                attr <- l.vec$attr_open_by_name(attr_name = "categories", ".")
                attr$write(ref)
            }
        }
        if (length(noncat.num.vars) > 0) {
            for (i in 1:NCOL(numDF)) {
                x <- numDF[, i]
                nn <- colnames(numDF)[i]
                h5group$create_dataset(nn, as.single(x),
                    gzip_level = compression_level,
                    dtype = h5types$H5T_IEEE_F32LE
                )
            }
        }

        if (length(noncat.vars) > 0) {
            for (i in 1:NCOL(nonNumDF)) {
                x <- nonNumDF[, i]
                nn <- colnames(nonNumDF)[i]
                dtype <- H5T_STRING$new(type = "c", size = Inf)
                dtype <- dtype$set_cset(cset = "unknown")
                h5group$create_dataset(nn, x,
                    gzip_level = compression_level,
                    dtype = string.dtype
                )
            }
        }
    } else {
        dtype <- H5T_STRING$new(type = "c", size = Inf)
        dtype <- dtype$set_cset(cset = "unknown")
        space <- H5S$new(type = "simple", dims = 0, maxdims = 10)

        h5group$create_attr(attr_name = "column-order", dtype = dtype, space = space)
    }

    index <- rownames(DF)
    if (length(unique(index)) < length(index)) {
        index <- make.names(index, unique = TRUE)
    }
    h5group$create_dataset("index", index, gzip_level = compression_level, dtype = string.dtype)
}

#' @import hdf5r
write.HD5SpMat <- function(h5file,
                           gname,
                           X,
                           compression_level = 0) {
    X <- Matrix::t(as(X, "dMatrix"))
    Xgroup <- h5file$create_group(gname)


    Xgroup$create_dataset("indices", X@i, gzip_level = compression_level, dtype = h5types$H5T_NATIVE_INT32)
    Xgroup$create_dataset("indptr", X@p, gzip_level = compression_level, dtype = h5types$H5T_NATIVE_INT32)
    Xgroup$create_dataset("data", as.single(X@x),
        gzip_level = compression_level,
        dtype = h5types$H5T_IEEE_F32LE
    )

    h5addAttr.str(Xgroup, "encoding-type", "csc_matrix")
    h5addAttr.str(Xgroup, "encoding-version", "0.1.0")
    h5attr(Xgroup, "shape") <- dim(X)
}

#' @import hdf5r
write.HD5List <- function(h5file,
                          gname,
                          obj_list,
                          depth = 1,
                          max_depth = 5,
                          compression_level = 0) {
    h5group <- h5file$create_group(gname)

    obj_list <- as.list(obj_list)
    obj_list <- obj_list[match(unique(names(obj_list)), names(obj_list))]

    # string.dtype <- H5T_STRING$new(type = "c", size = 100)
    # string.dtype <- string.dtype$set_cset(cset = "unknown")

    for (nn in names(obj_list)) {
        obj <- obj_list[[nn]]
        if ((sum(sapply(c("list", "SimpleList"), function(x) {
            return(length(which(is(obj) ==
                x)) != 0)
        })) != 0) & (depth < max_depth)) {
            write.HD5List(h5group, nn, obj,
                depth = depth + 1, max_depth = max_depth,
                compression_level = compression_level
            )
        } else if (sum(sapply(c("data.frame", "DataFrame", "DFrame"), function(x) {
            return(length(which(is(obj) ==
                x)) != 0)
        })) != 0) {
            write.HD5DF(h5group, nn, obj, compression_level = compression_level)
        } else if (is.sparseMatrix(obj)) {
            write.HD5SpMat(h5group, nn, obj, compression_level = compression_level)
        } else if (is.matrix(obj) | is.numeric(obj)) {
            h5group$create_dataset(nn, obj, gzip_level = compression_level, dtype = h5types$H5T_IEEE_F32LE)
        } else if (is.character(obj) & length(obj) == 1) {
            h5group$create_dataset(nn, obj, gzip_level = compression_level, dtype = H5T_STRING$new(type = "c", size = stringi::stri_length(obj)))
            # h5group$create_dataset(nn, obj, gzip_level = compression_level, dtype = string.dtype)
        } else {
            h5group[[nn]] <- obj
        }
    }
}

#' @import hdf5r
read.HD5DF <- function(h5file,
                       gname,
                       compression_level = 0) {
    h5group <- h5file[[gname]]

    if (!(h5group$attr_open_by_name("encoding-type", ".")$read() == "dataframe")) {
        err <- sprintf("%s is not a dataframe.\n", gname)
        stop(err)
    }


    key <- h5group$attr_open_by_name("_index", ".")$read()
    attr <- h5attributes(h5group[[key]])
    rn <- import.h5.data.slot(h5group, gname = key, attr)

    if (h5file$attr_exists_by_name(obj_name = gname, attr_name = "column-order")) {
        cn <- h5group$attr_open_by_name("column-order", ".")
        if (cn$get_storage_size() == 0) {
            DF <- DataFrame(row.names = rn)
            return(DF)
        }

        column.names <- cn$read()
        vars <- vector("list", length(column.names))
        names(vars) <- column.names
        for (vn in names(vars)) {
            attr <- h5attributes(h5group[[vn]])
            v <- import.h5.data.slot(h5group, vn, attr)

            # v[v == -1] <- NA
            vars[[vn]] <- v
        }

        if ("__categories" %in% names(h5group)) {
            cat <- h5group[["__categories"]]
            for (nn in names(cat)) {
                if (!(nn %in% column.names)) {
                    next
                }
                attr <- h5attributes(cat[[nn]])
                l <- import.h5.data.slot(cat, nn, attr)

                vars[[nn]] <- factor(l[vars[[nn]] + 1], l)
            }
        }
        DF <- DataFrame(vars)
        rownames(DF) <- rn
    } else {
        DF <- DataFrame(row.names = rn)
    }
    invisible(gc())
    return(DF)
}

#' @import hdf5r
read.HD5SpMat <- function(h5file,
                          gname,
                          compression_level = 0) {
    h5group <- h5file[[gname]]
    attr <- h5attributes(h5group)
    if (!(("encoding-type" %in% names(attr)) & (attr[["encoding-type"]] %in% c(
        "csc_matrix",
        "csr_matrix"
    )))) {
        err <- sprintf("%s is not a sparse matrix.\n", gname)
        stop(err)
    }

    data <- h5group[["data"]]$read()
    indices <- h5group[["indices"]]$read()
    indptr <- h5group[["indptr"]]$read()

    Dims <- attr$shape
    if (attr[["encoding-type"]] == "csc_matrix") {
        csc_sort_indices_inplace(indptr, indices, data)
        Xt <- Matrix::sparseMatrix(i = indices + 1, p = indptr, x = data, dims = Dims)
        X <- Matrix::t(Xt)
        rm(Xt)
        invisible(gc())
    } else if (attr[["encoding-type"]] == "csr_matrix") {
        csc_sort_indices_inplace(indptr, indices, data)
        Xt <- Matrix::sparseMatrix(j = indices + 1, p = indptr, x = data, dims = Dims)
        X <- Matrix::t(Xt)
        rm(Xt)
        invisible(gc())
    }

    X <- as(X, "dMatrix")
    return(X)
}


#' @import hdf5r
read.HD5Categorial <- function(h5file,
                               gname,
                               compression_level = 0) {
    h5group <- h5file[[gname]]

    codes <- h5group[["codes"]]$read()
    categories <- h5group[["categories"]]$read()

    idx <- codes + 1
    idx[idx <= 0] <- NA
    vals <- categories[idx]
    f <- factor(vals, categories)

    return(f)
}

#' @import hdf5r
read.HD5Dict <- function(h5file,
                         gname,
                         compression_level = 0) {
    ll <- read.HD5List(h5file = h5file, gname = gname, compression_level = compression_level)

    return(ll)
}

#' @import hdf5r
read.HD5StringArray <- function(h5file,
                                gname,
                                compression_level = 0) {
    arr <- h5file[[gname]]$read()

    return(arr)
}

#' @import hdf5r
read.HD5List <- function(h5file,
                         gname,
                         depth = 1,
                         max_depth = 5,
                         compression_level = 0) {
    h5group <- h5file[[gname]]

    obj_names <- names(h5group)
    L.out <- vector("list", length(obj_names))
    names(L.out) <- obj_names

    if (length(obj_names) > 0) {
        for (nn in obj_names) {
            attr <- h5attributes(h5group[[nn]])
            if (length(attr) > 0 & ("encoding-type" %in% names(attr))) {
                obj <- import.h5.data.slot(h5group, nn, attr)
            } else if (h5group[[nn]]$get_obj_type() == 2 & (depth < max_depth)) {
                obj <- read.HD5List(h5group, nn,
                    compression_level = compression_level,
                    depth = depth + 1
                )
            } else {
                obj <- h5group[[nn]]$read()
            }
            L.out[[nn]] <- obj
        }
        filter.mask <- sapply(L.out, function(x) is.null(x))
        if (sum(filter.mask) > 0) {
            L.out <- L.out[!filter.mask]
        }
    }
    invisible(gc())
    return(L.out)
}

#' @import hdf5r
#' @export
ACE2AnnData <- function(ace,
                        file,
                        main_assay = NULL,
                        full.export = TRUE,
                        compression_level = 0) {
    ACTIONetExperiment:::.check_and_load_package("hdf5r")

    # Ensure it can be case as an ACE object
    ace <- as(ace, "ACTIONetExperiment")
    if (is.null(main_assay)) {
        if ("default_assay" %in% names(metadata(ace))) {
            message(sprintf("Input main_assay is NULL. Setting main_assay to the metadata(ace)[['default_assay']]"))
            main_assay <- metadata(ace)[["default_assay"]]
        } else {
            if ("logcounts" %in% names(assays(ace))) {
                main_assay <- "logcounts"
            } else {
                main_assay <- "counts"
            }
            message(sprintf("Input main_assay is NULL. Setting main_assay to %s", main_assay))
        }
    }
    if (!(main_assay %in% names(assays(ace)))) {
        err <- sprintf("Input main_assay (%s) does not exist in assays(ace).", main_assay)
        stop(err)
    }


    if (file.exists(file)) {
        file.remove(file)
    }

    if (is.null(colnames(ace))) {
        colnames(ace) <- .default_colnames(NCOL(ace))
    }
    if (is.null(rownames(ace))) {
        rownames(ace) <- .default_rownames(NROW(ace))
    }

    colnames(ace) <- ucn <- ACTIONetExperiment:::.make_chars_unique(colnames(ace))
    rownames(ace) <- urn <- ACTIONetExperiment:::.make_chars_unique(rownames(ace))

    for (nn in names(SummarizedExperiment::assays(ace))) {
        dimnames(SummarizedExperiment::assays(ace)[[nn]]) <- list(urn, ucn)
    }

    h5file <- H5File$new(file, mode = "w")

    ## Write X (assays in ace, in either sparse or dense format)
    if (is.null(main_assay)) {
        main_mat <- Matrix::sparseMatrix(i = c(), j = c(), dims = dim(ace))
        write.HD5SpMat(h5file, gname = "X", main_mat, compression_level = compression_level)
    } else {
        if (!(main_assay %in% names(SummarizedExperiment::assays(ace)))) {
            err <- sprintf("'main_assay' is not in assays of ace'.\n")
            stop(err)
        }
        main_mat <- SummarizedExperiment::assays(ace)[[main_assay]]
        if (is.sparseMatrix(main_mat)) {
            write.HD5SpMat(h5file, gname = "X", main_mat, compression_level = compression_level)
        } else {
            h5file$create_dataset("X", main_mat, gzip_level = compression_level, dtype = h5types$H5T_IEEE_F32LE)
        }
    }

    remaining.assays <- setdiff(names(SummarizedExperiment::assays(ace)), main_assay)
    if ((full.export == T) & (0 < length(remaining.assays))) {
        layers <- h5file$create_group("layers")

        for (an in remaining.assays) {
            Xr <- SummarizedExperiment::assays(ace)[[an]]
            if (is.sparseMatrix(Xr)) {
                write.HD5SpMat(layers, gname = an, Xr, compression_level = compression_level)
            } else {
                layers$create_dataset(an, Xr, gzip_level = compression_level, dtype = h5types$H5T_IEEE_F32LE)
            }
        }
    }

    uns <- h5file$create_group("uns")
    obsm_annot <- uns$create_group("obsm_annot")
    varm_annot <- uns$create_group("varm_annot")

    metadata(ace)$main_assay <- main_assay
    obj_list <- metadata(ace)
    write.HD5List(uns, "metadata", obj_list, depth = 1, max_depth = 10, compression_level = compression_level)

    ## Write obs (colData() in ace)
    obs.DF <- as.data.frame(SummarizedExperiment::colData(ace))
    if (0 < NCOL(obs.DF)) {
        obs.DF <- as.data.frame(lapply(SummarizedExperiment::colData(ace), function(x) {
            if (is.numeric(x) & (!is.null(names(x))) & (length(unique(names(x))) < length(x) / 2)) {
                return(factor(names(x), names(x)[match(unique(x), x)]))
            } else {
                return(x)
            }
        }))
    }
    rownames(obs.DF) <- colnames(ace)

    write.HD5DF(h5file, gname = "obs", obs.DF, compression_level = compression_level)

    ## Write var (matching rowData() in ace)
    var.DF <- as.data.frame(SummarizedExperiment::rowData(ace))
    rownames(var.DF) <- rownames(ace)



    if (is(ace, "RangedSummarizedExperiment")) {
        GR <- rowRanges(ace)
        GR.df <- GenomicRanges::as.data.frame(GR)
        if (nrow(GR.df) == nrow(var.DF)) {
            BED <- data.frame(chr = as.character(GenomeInfoDb::seqnames(GR)), start = start(GR), end = end(GR))
            var.DF <- cbind(BED, elementMetadata(GR), var.DF)
            metadata(ace)$genome <- genome(ace)
        }
    }
    write.HD5DF(h5file, "var", var.DF, compression_level = compression_level)

    ## Write subset of obsm related to the cell embeddings (Dim=2 or 3)
    obsm <- h5file$create_group("obsm")
    obsm.mats <- colMaps(ace)
    obsm.public.idx <- which(colMapTypes(ace) != "internal")
    if (length(obsm.public.idx) > 0) {
        obsm.subset <- obsm.mats[obsm.public.idx]
        for (i in 1:length(obsm.subset)) {
            nn <- names(obsm.subset)[[i]]
            Y <- Matrix::t(obsm.subset[[i]])
            isDimRed <- (NROW(Y) <= 3) | (colMapTypes(ace)[[nn]] == "embedding")
            if (isDimRed) {
                # AD_nn <- paste("X", nn, sep = "_")
                AD_nn <- nn
            } else {
                AD_nn <- nn
            }

            if (is.matrix(Y)) {
                obsm$create_dataset(AD_nn, Y, gzip_level = compression_level, dtype = h5types$H5T_IEEE_F32LE)
            } else {
                write.HD5SpMat(obsm, AD_nn, Y, compression_level)
            }

            factor_info <- obsm_annot$create_group(AD_nn)
            factor_info[["type"]] <- colMapTypes(ace)[[nn]]
            factor.meta.DF <- colMapMeta(ace)[[nn]]
            if (NCOL(factor.meta.DF) > 0) {
                write.HD5DF(factor_info, "annotatation", factor.meta.DF, compression_level = 0)
            }
        }
    }

    varm <- h5file$create_group("varm")
    varm.mats <- rowMaps(ace)
    varm.public.idx <- which(rowMapTypes(ace) != "internal")
    if (length(varm.public.idx) > 0) {
        varm.subset <- varm.mats[varm.public.idx]
        for (i in 1:length(varm.subset)) {
            nn <- names(varm.subset)[[i]]
            Y <- Matrix::t(varm.subset[[i]])

            if (is.matrix(Y)) {
                varm$create_dataset(nn, Y, gzip_level = compression_level, dtype = h5types$H5T_IEEE_F32LE)
            } else {
                write.HD5SpMat(varm, nn, Y, compression_level)
            }


            factor_info <- varm_annot$create_group(nn)
            factor_info[["type"]] <- rowMapTypes(ace)[[nn]]
            factor.meta.DF <- rowMapMeta(ace)[[nn]]
            if (NCOL(factor.meta.DF) > 0) {
                write.HD5DF(factor_info, "annotatation", factor.meta.DF, compression_level = 0)
            }
        }
    }

    if (full.export) {
        print("Full export mode")

        if (length(obsm.public.idx) < length(obsm.mats)) {
            obsm.private.idx <- setdiff(1:length(obsm.mats), obsm.public.idx)
            if (length(obsm.private.idx) > 0) {
                obsm.subset <- obsm.mats[obsm.private.idx]
                for (i in 1:length(obsm.subset)) {
                    nn <- names(obsm.subset)[[i]]
                    Y <- Matrix::t(obsm.subset[[i]])
                    if (is.matrix(Y)) {
                        obsm$create_dataset(nn, Y, gzip_level = compression_level, dtype = h5types$H5T_IEEE_F32LE)
                    } else {
                        write.HD5SpMat(obsm, nn, Y, compression_level)
                    }

                    factor_info <- obsm_annot$create_group(nn)
                    factor_info[["type"]] <- colMapTypes(ace)[[nn]]
                    factor.meta.DF <- colMapMeta(ace)[[nn]]
                    if (NCOL(factor.meta.DF) > 0) {
                        write.HD5DF(factor_info, "annotation", factor.meta.DF, compression_level = 0)
                    }
                }
            }
        }

        if (length(varm.public.idx) < length(varm.mats)) {
            varm.private.idx <- setdiff(1:length(varm.mats), varm.public.idx)
            if (length(varm.private.idx) > 0) {
                varm.subset <- varm.mats[varm.private.idx]
                for (i in 1:length(varm.subset)) {
                    nn <- names(varm.subset)[[i]]
                    Y <- Matrix::t(varm.subset[[i]])

                    if (is.matrix(Y)) {
                        varm$create_dataset(nn, Y, gzip_level = compression_level, dtype = h5types$H5T_IEEE_F32LE)
                    } else {
                        write.HD5SpMat(varm, nn, Y, compression_level)
                    }

                    factor_info <- varm_annot$create_group(nn)
                    factor_info[["type"]] <- rowMapTypes(ace)[[nn]]
                    factor.meta.DF <- rowMapMeta(ace)[[nn]]
                    if (NCOL(factor.meta.DF) > 0) {
                        write.HD5DF(factor_info, "annotatation", factor.meta.DF, compression_level = 0)
                    }
                }
            }
        }

        # Export 'obsp'-associated matrices, i.e. colNets(): obs in AnnData ~ cols in SCE
        # ~ cells => cell-cell networks (such as ACTIONet)
        CN <- colNets(ace)
        if ((length(CN) > 0)) {
            obsp <- h5file$create_group("obsp")
            CN <- lapply(CN, function(x) as(x, "dMatrix"))

            for (i in 1:length(CN)) {
                write.HD5SpMat(obsp, gname = names(CN)[[i]], CN[[i]], compression_level = compression_level)
            }
        }

        # Export 'varp'-associated matrices, i.e. rowNets(): var in AnnData ~ rows in SCE
        # ~ genes => gene-gene networks (such as SCINET)
        RN <- rowNets(ace)
        if ((length(RN) > 0)) {
            varp <- h5file$create_group("varp")
            RN <- lapply(RN, function(x) as(x, "dMatrix"))

            for (i in 1:length(RN)) {
                write.HD5SpMat(varp, gname = names(RN)[[i]], RN[[i]], compression_level = compression_level)
            }
        }
    }

    h5file$close_all()
}

import.h5.data.slot <- function(h5file, gname, attr) {
    if (length(attr) == 0 || !("encoding-type" %in% names(attr))) {
        X <- h5file[[gname]]$read()
    } else {
        if (attr[["encoding-type"]] == "array") {
            X <- h5file[[gname]]$read()
        } else if (attr[["encoding-type"]] == "dataframe") {
            X <- as.matrix(read.HD5DF(h5file = h5file, gname = gname))
        } else if (attr[["encoding-type"]] %in% c("csr_matrix", "csc_matrix")) {
            X <- read.HD5SpMat(h5file = h5file, gname = gname)
        } else if (attr[["encoding-type"]] == "categorical") {
            X <- read.HD5Categorial(h5file = h5file, gname = gname)
        } else if (attr[["encoding-type"]] == "string-array") {
            X <- read.HD5StringArray(h5file = h5file, gname = gname)
        } else if (attr[["encoding-type"]] == "dict") {
            X <- read.HD5Dict(h5file = h5file, gname = gname)
        } else if (attr[["encoding-type"]] == "string") {
            X <- read.HD5StringArray(h5file = h5file, gname = gname)
        } else if (attr[["encoding-type"]] == "numeric-scalar") {
            X <- h5file[[gname]]$read()
        } else {
            stop(sprintf("Unknown format %s for X", attr[["encoding-type"]]))
        }
    }

    return(X)
}

#' @import hdf5r
#' @export
AnnData2ACE <- function(file,
                        import_X = TRUE) {
    .check_and_load_package("hdf5r")

    h5file <- H5File$new(file, mode = "r")

    objs <- names(h5file)
    input_assays <- list()

    if ("layers" %in% objs) {
        layers <- h5file[["layers"]]
        # additional_assays <- vector("list", length(names(layers)))
        # names(additional_assays) <- names(layers)
        input_assays <- vector("list", length(names(layers)))
        names(input_assays) <- names(layers)

        for (an in names(layers)) {
            attr <- h5attributes(layers[[an]])
            # additional_assays[[an]] <- import.h5.data.slot(layers, an, attr)
            input_assays[[an]] <- import.h5.data.slot(layers, an, attr)
        }
        # input_assays <- c(input_assays, additional_assays)
    }
    invisible(gc())

    if (length(input_assays) == 0 && import_X == FALSE) {
        stop("input file has no assays.")
    } else if (import_X == TRUE) {
        X.attr <- h5attributes(h5file[["X"]])
        X <- import.h5.data.slot(h5file, "X", X.attr)

        # input_assays <- list(X)
        # names(input_assays) <- main_assay

        input_assays <- c(list("X" = X), input_assays)
    }
    invisible(gc())

    default_assay_name <- names(input_assays)[1]

    if ("obs" %in% objs) {
        obs.DF <- read.HD5DF(h5file = h5file, gname = "obs")
    } else {
        obs.DF <- DataFrame(row.names = paste("Cell", 1:NCOL(X), sep = ""))
    }

    if ("var" %in% objs) {
        var.DF <- read.HD5DF(h5file = h5file, gname = "var")
    } else {
        var.DF <- DataFrame(row.names = paste("Gene", 1:NROW(X), sep = ""))
    }

    input_assays <- lapply(input_assays, function(X) {
        rownames(X) <- rownames(var.DF)
        colnames(X) <- rownames(obs.DF)
        return(X)
    })

    ace <- ACTIONetExperiment(assays = input_assays, rowData = var.DF, colData = obs.DF)
    metadata(ace)[["default_assay"]] <- default_assay_name

    rm(input_assays)
    invisible(gc())

    var.DF <- rowData(ace)
    if (ncol(var.DF) > 0) {
        if (all(colnames(var.DF) %in% c("chr", "start", "end"))) {
            GR <- GenomicRanges::makeGRangesFromDataFrame(var.DF, keep.extra.columns = T)
            SummarizedExperiment::rowRanges(ace) <- GR
        }
    }

    if ("obsm" %in% objs) {
        obsm <- h5file[["obsm"]]
        for (mn in names(obsm)) {
            attr <- h5attributes(obsm[[mn]])
            Xr <- import.h5.data.slot(obsm, mn, attr)

            if (sum(grepl(pattern = "^X_", mn))) {
                nn <- stringr::str_sub(mn, start = 3)
            } else {
                nn <- mn
            }

            if (nrow(Xr) != ncol(ace)) {
                Xr <- Matrix::t(Xr)
            }

            colMaps(ace)[[nn]] <- Xr
            rm(Xr)
            invisible(gc())
        }
    }

    if ("varm" %in% objs) {
        varm <- h5file[["varm"]]
        for (nn in names(varm)) {
            attr <- h5attributes(varm[[nn]])
            Xr <- import.h5.data.slot(varm, nn, attr)

            if (nrow(Xr) != nrow(ace)) {
                Xr <- Matrix::t(Xr)
            }

            rowMaps(ace)[[nn]] <- Xr
            rm(Xr)
            invisible(gc())
        }
    }


    if ("obsp" %in% objs) {
        obsp <- h5file[["obsp"]]
        for (pn in names(obsp)) {
            attr <- h5attributes(obsp[[pn]])
            Net <- import.h5.data.slot(obsp, pn, attr)

            colNets(ace)[[pn]] <- Net
        }
    }

    if ("varp" %in% objs) {
        varp <- h5file[["varp"]]
        for (pn in names(varp)) {
            attr <- h5attributes(varp[[pn]])
            Net <- import.h5.data.slot(varp, pn, attr)

            rowNets(ace)[[pn]] <- Net
        }
    }


    if ("uns" %in% objs) {
        uns <- h5file[["uns"]]
        if ("obsm_annot" %in% names(uns)) {
            # Import obs annotations
            obsm_annot <- uns[["obsm_annot"]]
            for (nn in names(obsm_annot)) {
                factor_annot <- obsm_annot[[nn]]
                colMapTypes(ace)[[nn]] <- factor_annot[["type"]]$read()
                if ("annotation" %in% names(factor_annot)) {
                    DF <- read.HD5DF(factor_annot, "annotation")
                    colMapMeta(ace)[[nn]] <- DF
                }
            }
        }
        if ("varm_annot" %in% names(uns)) {
            # Import obs annotations
            var_annot <- uns[["varm_annot"]]
            for (nn in names(var_annot)) {
                factor_annot <- var_annot[[nn]]
                rowMapTypes(ace)[[nn]] <- factor_annot[["type"]]$read()
                if ("annotation" %in% names(factor_annot)) {
                    DF <- read.HD5DF(factor_annot, "annotation")
                    rowMapMeta(ace)[[nn]] <- DF
                }
            }
        }
        if ("metadata" %in% names(uns)) {
            meta_data_objs <- read.HD5List(uns, "metadata",
                depth = 1, max_depth = 10,
                compression_level = compression_level
            )
        } else {
            meta_data_objs <- NULL
        }



        meta <- read.HD5List(h5file = h5file, gname = "uns")
        meta <- meta[setdiff(names(meta), c("obsm_annot", "varm_annot", "metadata"))]
        if (!is.null(meta_data_objs)) {
            if (length(meta) == 0) {
                meta <- meta_data_objs
            } else {
                meta <- c(meta_data_objs, meta)
            }
        }
        metadata(ace) <- meta
    }

    h5file$close_all()

    if (!("logcounts" %in% names(assays(ace))) & ("X" %in% names(assays(ace)))) {
        X <- assays(ace)[["X"]]
        max_samples <- min(100, min(dim(X)))
        subX <- X[1:max_samples, 1:max_samples]
        x <- as.numeric(subX)
        if (length(setdiff(unique(x), 0:max(round(x) + 1))) == 0) {
            if (!("counts" %in% names(assays(ace)))) {
                names(assays(ace))[which(names(assays(ace)) == "X")] <- "counts"
            }
        } else {
            if (!("logcounts" %in% names(assays(ace)))) {
                names(assays(ace))[which(names(assays(ace)) == "X")] <- "logcounts"
                metadata(ace)[["default_assay"]] <- "logcounts"
            }
        }
    } else {
        metadata(ace)[["default_assay"]] <- "X"
    }

    return(ace)
}