#' Get DESeq2 Results
#'
#' This function obtains a set of comparisons from a DESeq2 analysis, given a named list of contrasts.
#' It allows additional model
#' parameters to be specified and a design matrix to be manually adjusted.
#'
#' @param dds A [SummarizedExperiment::SummarizedExperiment] or [DESeq2::DESeqDataSet] object.
#' @param contrasts A named list of contrasts, e.g. `list("condition" = c("condition", "A", "B"))`.
#'   The first element is the variable of interest, the second is the test, and the third is the reference level.
#'   The name of  each element in the list will be used as a name in the results table.
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
#'   Must be one of "apeglm", "ashr", or NULL.
#'   Defaults to "ashr".
#' @param norm.ercc A logical indicating whether to normalize to ERCC spike-ins.
#' @param add.rowData A vector of column names from the rowData slot of the DESeqDataSet
#'   to be added to the results table.
#'   Defaults to NULL.
#' @param BPPARAM The BiocParallelParam object specifying the parallel back-end to be used.
#'   Defaults to NULL.
#'
#' @return A named list of [DESeq2::DESeqResults] objects for the specified contrasts.
#'   If `add.rowData` is supplied, these will be returned as [S4Vectors::DFrame] objects instead.
#'
#' @details It is important to note that LFC shrinkage is independent of the typical MLE results table.
#'   That is, if `lfc.th` is provided, the results table p-values will reflect that testing threshold.
#'   If `shrink.method` is set to `ashr` with `lfc.th`, s-values will be returned in addition to the MLE p-values.
#'
#'   It is possible to have shrunken FCs that are near 0 that still have a significant p-value.
#'   This is frustrating, as it means you end up having to do post-hoc filtering on LFC
#'   or filter with s-values, which are more difficult to interpret than the adjusted p-values.
#'
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

    if (shrink.method == "ashr") {
        .package_check("ashr")
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

        if (!is.null(shrink.method)) {
            out.name <- paste0(rname, "-shLFC")

            # ashr does not need coef, this is to ensure no error with user-supplied model matrix/list contrasts
            if (shrink.method == "ashr") {
                shrink <- lfcShrink(dds, res = res1, contrast = con, type = shrink.method)
            } else {
                shrink <- lfcShrink(dds, res = res1, coef = coef, type = shrink.method)
            }

            if (!is.null(add.rowData) & all(add.rowData %in% colnames(rowData(dds)))) {
                message("Adding rowData columns to shrunken LFC results table.")
                shrink <- cbind(shrink, rowData(dds)[, add.rowData])
            }

            res.list[[out.name]] <- shrink
        }

        # Add original results to list.
        if (!is.null(add.rowData) & all(add.rowData %in% colnames(rowData(dds)))) {
            message("Adding rowData columns to results table.")
            res1 <- cbind(res1, rowData(dds)[, add.rowData])
        }

        res.list[[rname]] <- res1

        for (l in lfc.th) {
            res <- results(dds, contrast = con, alpha = alpha, lfcThreshold = l)

            if (!is.null(shrink.method)) {
                # ashr does not need coef, this is to ensure no error with user-supplied model matrix/list contrasts
                if (shrink.method == "ashr") {
                    coef <- NULL
                    out.name <- paste0(rname, "-shLFC", round(l, 3))
                    shrink <- lfcShrink(dds, res = res, contrast = con, lfcThreshold = l, type = shrink.method)
                } else {
                    out.name <- paste0(rname, "-shLFC", round(l, 3))
                    shrink <- lfcShrink(dds, res = res, coef = coef, type = shrink.method)
                }

                if (!is.null(add.rowData) & all(add.rowData %in% colnames(rowData(dds)))) {
                    message("Adding rowData columns to shrunken LFC results table.")
                    shrink <- cbind(shrink, rowData(dds)[, add.rowData])
                }

                res.list[[out.name]] <- shrink
            }

            if (!is.null(add.rowData) & all(add.rowData %in% colnames(rowData(dds)))) {
                message("Adding rowData columns to results table.")
                res <- cbind(res, rowData(dds)[, add.rowData])
            }

            out.name <- paste0(rname, "-LFC", round(l, 3))
            res.list[[out.name]] <- res
        }
    }

    return(res.list)
}


#' Get edgeR Results
#'
#' This function obtains a set of comparisons from a edgeR analysis, given a named list of contrasts.
#' It allows additional model
#' parameters to be specified and a design matrix to be manually adjusted.
#'
#' @details
#' This function is designed to work with relatively simple designs,
#' e.g. variable of interest with a blocking factor or two.
#' If you need a bunch of interaction terms or nested designs, this function is probably not for you.
#'
#' @param se A [SummarizedExperiment::SummarizedExperiment] or [edgeR::DGEList] object.
#' @param contrasts A named list of contrasts, e.g. `list("condition_AvB" = c("condition", "A", "B"))`.
#'   The first element is the variable of interest, the second is the test, and the third is the reference level.
#'   The name of  each element in the list will be used as a name in the results table.
#' @param res.list A named list to hold edgeR result.
#'   Allows the function to be run multiple times if needed and append to the same list.
#'   Defaults to an empty list.
#' @param block A vector of additional terms to be considered in the model, beyond the main effect.
#'   Defaults to NULL.
#' @param lfc.th A numeric vector of log2 fold-change thresholds.
#'   Defaults to `c(log2(1.25), log2(1.5))`.
#' @param norm.ercc A logical indicating whether to normalize to ERCC spike-ins.
#' @param ercc.pattern A character string indicating the pattern to match ERCC spike-ins.
#'   Defaults to `^ERCC-`.
#' @param use.lrt A logical indicating whether to use the likelihood ratio test (LRT)
#'   instead of the quasi-likelihood F-test.
#'
#' @return A named list of [edgeR::TopTags-class] objects for the specified contrasts.
#'
#' @import edgeR
#' @importFrom stats as.formula relevel
#' @export
#'
#' @author Jared Andrews
#'
#' @examples
#' nlibs <- 4
#' ngenes <- 1000
#' dispersion.true <- 1/rchisq(ngenes, df=10)
#' group <- factor(c("A", "A", "B", "B"))
#' design <- model.matrix(~group)
#'
#' # Generate count data
#' y <- rnbinom(ngenes*nlibs,mu=20,size=1/dispersion.true)
#' y <- matrix(y,ngenes,nlibs)
#' d <- DGEList(y)
#'
#' res <- get_edgeR_res(d, contrasts = list("condition_AvB" = c("group", "A", "B")))
#' names(res)
#' head(res[[1]])
#'
get_edgeR_res <- function(
    se,
    contrasts,
    res.list = list(),
    block = NULL,
    lfc.th = c(log2(1.25), log2(1.5)),
    norm.ercc = FALSE,
    ercc.pattern = "^ERCC-",
    use.lrt = FALSE) {

    # Convert to DGEList if SummarizedExperiment object is provided.
    if (is(se, "SummarizedExperiment")) {
        se <- SE2DGEList(se)
    } else if (!is(se, "DGEList")) {
        stop("Input object must be a SummarizedExperiment or DGEList.")
    }

    # Normalize to ERCC spike-ins if specified.
    # As demonstrated in: https://support.bioconductor.org/p/9135179/#9135334
    if (norm.ercc) {
        spikes <- rownames(se)[grep(ercc.pattern, rownames(se))]
        message(paste0("\nCalculating size factors from ", length(spikes), " ERCC spike-ins."))
        spike_in_index <- which(rownames(se) %in% spikes)
        spike_in_factor <- as.numeric(counts(se)[spike_in_index, ]) / colSums(counts(se))
        dropped_spike_in <- se[!rownames(se) %in% spikes, ]

        # So that library size is recalculated without spike-ins.
        se <- DGEList(dropped_spike_in)

        norm.factors <- spike_in_factor / se$samples$lib.size
        norm.factors <- norm.factors / prod(norm.factors)^(1 / length(norm.factors))
        se$samples$norm.factors <- norm.factors
    } else {
        se <- calcNormFactors(se)
    }

    for (i in seq_along(contrasts)) {
        rname <- names(contrasts)[i]

        con <- contrasts[[i]]
        coef <- paste0(con[1], con[2])

        se$samples[[con[1]]] <- relevel(se$samples[[con[1]]], ref = con[3])

        if (!is.null(block)) {
            desgn <- as.formula(paste0("~", paste0(c(block, con[1]), collapse = "+")))
            mm <- model.matrix(desgn, data = se$samples)
        } else {
            desgn <- as.formula(paste0("~", con[1]))
            mm <- model.matrix(desgn, data = se$samples)
        }

        message(paste0(
            "\nDesign for ", paste(con[1], con[2], "vs", con[3], sep = "_"),
            " is ", paste0(as.character(desgn), collapse = "")
        ))

        se <- estimateDisp(se, design = mm)

        if (use.lrt) {
            message("Using likelihood ratio test (LRT) instead of quasi-likelihood F-test.")
            se_fit <- glmFit(se, design = mm)
            res <- glmLRT(se_fit, coef = coef)
        } else {
            se_fit <- glmQLFit(se, design = mm)
            res <- glmQLFTest(se_fit, coef = coef)
        }

        # Add results to list.
        res <- topTags(res, n = Inf)

        res.list[[rname]] <- res

        for (l in lfc.th) {
            message(paste0("Calculating results for LFC threshold ", round(l, 3), " using glmTreat."))
            res <- glmTreat(se_fit, coef = coef, lfc = l)
            res <- topTags(res, n = Inf)

            out.name <- paste0(rname, "-LFC", round(l, 3))
            res.list[[out.name]] <- res
        }
    }

    return(res.list)
}
