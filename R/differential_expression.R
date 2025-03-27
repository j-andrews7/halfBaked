#' Get DESeq2 Results
#'
#' This function obtains a set of comparisons from a DESeq2 analysis, given a named list of contrasts.
#' It allows additional model
#' parameters to be specified and a design matrix to be manually adjusted.
#'
#' @param dds An object of class DESeqDataSet.
#' @param contrasts A named list of contrasts.
#' @param res.list A named list to hold DESeq2 result tables.
#'   Allows the function to be run multiple times if needed and append to the same list.
#'   Defaults to an empty list.
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
#'   Defaults to `c(log2(1.25), log2(1.5))``.
#' @param shrink.method The method used for shrinkage estimation.
#'   Defaults to "ashr".
#' @param norm.ercc A logical indicating whether to normalize to ERCC spike-ins.
#' @param add.rowData A vector of column names from the rowData slot of the DESeqDataSet to be added to the results table.
#'   Defaults to NULL.
#' @param BPPARAM The BiocParallelParam object specifying the parallel back-end to be used.
#'   Defaults to NULL.
#'
#' @return A named list of \link[DESeq2:DESeqResults]{DESeq2::DESeqResults} objects for the specified contrasts.
#'   If `add.rowData` is supplied, these will be returned as \link[S4Vectors:DFrame]{S4Vectors::DFrame} objects instead.
#'
#' @import DESeq2
#' @importFrom stats as.formula relevel
#' @importFrom SummarizedExperiment rowData
#' @export
#'
#' @author Jared Andrews
#'
#' @examples
#' library(DESeq2)
#' dds_de <- makeExampleDESeqDataSet(n = 100, m = 12, betaSD = 2) # DE genes
#' rownames(dds_de) <- paste0("gene", 1001:1100)
#' dds <- makeExampleDESeqDataSet(n = 1000, m = 12) # Non-DE genes
#' dds <- rbind(dds, dds_de)
#' dds <- DESeq(dds)
#' contrasts <- list("condition" = c("condition", "A", "B"))
#' res <- get_DESeq2_res(dds, contrasts)
#'
#' names(res)
#' head(res[[1]])
#'
get_DESeq2_res <- function(
    dds,
    contrasts,
    res.list = list(),
    user.mat = FALSE,
    block = NULL,
    design = NULL,
    alpha = 0.05,
    lfc.th = c(log2(1.25), log2(1.5)),
    shrink.method = "ashr",
    norm.ercc = FALSE,
    add.rowData = NULL,
    BPPARAM = NULL) {

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
            " is ", paste0(as.character(desgn), collapse = "")
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

        if (!is.null(add.rowData) & all(add.rowData %in% colnames(rowData(dds)))) {
            message("Adding rowData columns to results table.")
            res1 <- cbind(res1, rowData(dds)[, add.rowData])
        }

        res.list[[rname]] <- res1

        if (!is.null(shrink.method)) {
            out.name <- paste0(rname, "-shLFC")

            # ashr does not need coef, this is to ensure no error with user-supplied model matrix/list contrasts
            if (shrink.method == "ashr") {
                coef <- NULL
                shrink <- lfcShrink(dds, contrast = con, type = shrink.method)
            } else {
                shrink <- lfcShrink(dds, res = res1, coef = coef, type = shrink.method)
            }

            if (!is.null(add.rowData) & all(add.rowData %in% colnames(rowData(dds)))) {
                message("Adding rowData columns to shrunken LFC results table.")
                shrink <- cbind(shrink, rowData(dds)[, add.rowData])
            }

            res.list[[out.name]] <- shrink
        }

        for (l in lfc.th) {
            res <- results(dds, contrast = con, alpha = alpha, lfcThreshold = l)
            if (!is.null(add.rowData) & all(add.rowData %in% colnames(rowData(dds)))) {
                message("Adding rowData columns to results table.")
                res <- cbind(res, rowData(dds)[, add.rowData])
            }

            if (!is.null(shrink.method)) {
                # ashr does not need coef, this is to ensure no error with user-supplied model matrix/list contrasts
                if (shrink.method == "ashr") {
                    coef <- NULL
                    out.name <- paste0(rname, "-shLFC", l)
                    shrink <- lfcShrink(dds, contrast = con, lfcThreshold = l, type = shrink.method)
                } else {
                    out.name <- paste0(rname, "-shLFC", l)
                    shrink <- lfcShrink(dds, res = res, coef = coef, type = shrink.method)
                }

                if (!is.null(add.rowData) & all(add.rowData %in% colnames(rowData(dds)))) {
                    message("Adding rowData columns to shrunken LFC results table.")
                    shrink <- cbind(shrink, rowData(dds)[, add.rowData])
                }

                res.list[[out.name]] <- shrink
            }

            out.name <- paste0(rname, "-LFC", l)
            res.list[[out.name]] <- res
        }
    }

    return(res.list)
}
