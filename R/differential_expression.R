#' Get DESeq2 Results
#'
#' This function obtains a set of comparisons from a DESeq2 analysis, given a named list of contrasts.
#' It allows additional model
#' parameters to be specified and a design matrix to be manually adjusted.
#'
#' @param dds An object of class DESeqDataSet.
#' @param res.list A list of DESeq2 result tables.
#'   Allows the function to be run multiple times if needed and append to the same list.
#' @param contrasts A named list of contrasts.
#' @param user.mat A logical indicating whether a user-specified model matrix is provided.
#'   Defaults to FALSE.
#' @param block A vector of additional terms to be considered in the model, beyond the main effect.
#'   Defaults to NULL.
#' @param design The design formula or matrix.
#'   If a matrix is provided, ensure 'user.mat' is set to TRUE.
#'   Defaults to NULL.
#' @param alpha The significance level for hypothesis testing.
#'   Defaults to 0.05.
#' @param lfc.th A numeric vector of log2 fold-change thresholds.
#'   Defaults to c(log2(1.25), log2(1.5)).
#' @param shrink.method The method used for shrinkage estimation.
#'   Defaults to "ashr".
#' @param norm.ercc A logical indicating whether to normalize to ERCC spike-ins.
#' @param BPPARAM The BiocParallelParam object specifying the parallel back-end to be used.
#'   Defaults to NULL.
#'
#' @return A named list of DESeq2 result tables for the specified contrasts.
#'
#' @import DESeq2
#' @export
#'
#' @author Jared Andrews
#'
#' @examples
#' \dontrun{
#' get_DESEQ2_res(dds, res.list, contrasts,
#'     user.mat = TRUE, block = c("term1", "term2"),
#'     design = my_design, alpha = 0.05, lfc.th = log2(2),
#'     shrink.method = "apeglm", outdir = "./my_results", BPPARAM = MulticoreParam(2)
#' )
#' }
#'
get_DESeq2_res <- function(
    dds,
    res.list,
    contrasts,
    user.mat = FALSE,
    block = NULL,
    design = NULL,
    alpha = 0.05,
    lfc.th = c(log2(1.25), log2(1.5)),
    shrink.method = "ashr",
    norm.ercc = FALSE,
    BPPARAM = NULL) {
    dir.create(file.path(outdir), showWarnings = FALSE, recursive = TRUE)

    if (shrink.method == "apeglm") {
        .package_check("apeglm")
    }

    for (i in seq_along(contrasts)) {
        rname <- names(contrasts)[i]

        # If user-supplied matrix, contrast must be in list format.
        if (user.mat) {
            con <- contrasts[[i]]
            message("Setting shrink.method to 'ashr' to work with list contrasts due to user-specified model matrix.")
            shrink.method <- "ashr"
        } else {
            con <- contrasts[[i]]
            coef <- paste(con[1], con[2], "vs", con[3], sep = "_")
            dds[[con[1]]] <- relevel(dds[[con[1]]], ref = con[3])
        }

        if (!is.null(design)) {
            desgn <- design
        } else if (!is.null(block)) {
            desgn <- as.formula(paste0("~", paste0(c(block, con[1]), collapse = "+")))
        } else {
            desgn <- as.formula(paste0("~", con[1]))
        }

        message(paste0(
            "\nDesign for ", paste(con[1], con[2], "vs", con[3], sep = "_"),
            " is ", paste0(as.character(desgn))
        ))

        dds <- DESeqDataSet(dds, design = desgn)

        # Get size factor by spike-ins if specified.
        if (norm.ercc) {
            spikes <- rownames(dds)[grep("^ERCC-", rownames(dds))]
            message(paste0("\nCalculating size factors from ", length(spikes), " ERCC spike-ins."))
            dds <- estimateSizeFactors(dds, controlGenes = rownames(dds) %in% spikes)
        }

        dds <- DESeq(dds, BPPARAM = BPPARAM)

        res1 <- results(dds, contrast = con, alpha = alpha)
        res1$ENSEMBL <- rownames(res1)
        res1$SYMBOL <- rowData(dds)$SYMBOL

        if (!is.null(shrink.method)) {
            out.name <- paste0(rname, "-shLFC")

            # ashr does not need coef, this is to ensure no error with user-supplied model matrix/list contrasts
            if (shrink.method == "ashr") {
                coef <- NULL
                shrink <- lfcShrink(dds, contrast = con, type = shrink.method)
                shrink$ENSEMBL <- rownames(shrink)
                shrink$SYMBOL <- rowData(dds)$SYMBOL
                rownames(shrink) <- shrink$SYMBOL
                shrink <- as.data.frame(shrink)
                res.list[[out.name]] <- shrink
            } else {
                shrink <- lfcShrink(dds, res = res1, coef = coef, type = shrink.method)
                shrink$ENSEMBL <- rownames(shrink)
                shrink$SYMBOL <- rowData(dds)$SYMBOL
                rownames(shrink) <- shrink$SYMBOL
                shrink <- as.data.frame(shrink)
                res.list[[out.name]] <- shrink
            }
        }

        rownames(res1) <- res1$SYMBOL
        res1 <- as.data.frame(res1)
        res.list[[rname]] <- res1

        for (l in lfc.th) {
            res <- results(dds, contrast = con, alpha = alpha, lfcThreshold = l)
            res$ENSEMBL <- rownames(res)
            res$SYMBOL <- rowData(dds)$SYMBOL

            if (!is.null(shrink.method)) {
                # ashr does not need coef, this is to ensure no error with user-supplied model matrix/list contrasts
                if (shrink.method == "ashr") {
                    coef <- NULL
                    out.name <- paste0(rname, "-shLFC", l)
                    shrink <- lfcShrink(dds, contrast = con, lfcThreshold = l, type = shrink.method)
                    shrink$ENSEMBL <- rownames(shrink)
                    shrink$SYMBOL <- rowData(dds)$SYMBOL
                    rownames(shrink) <- shrink$SYMBOL
                    shrink <- as.data.frame(shrink)
                    res.list[[out.name]] <- shrink
                } else {
                    out.name <- paste0(rname, "-shLFC", l)
                    shrink <- lfcShrink(dds, res = res, coef = coef, type = shrink.method)
                    shrink$ENSEMBL <- rownames(shrink)
                    shrink$SYMBOL <- rowData(dds)$SYMBOL
                    rownames(shrink) <- shrink$SYMBOL
                    shrink <- as.data.frame(shrink)
                    res.list[[out.name]] <- shrink
                }
            }

            rownames(res) <- res$SYMBOL
            out.name <- paste0(rname, "-LFC", l)
            res <- as.data.frame(res)
            res.list[[out.name]] <- res
        }
    }

    return(res.list)
}
