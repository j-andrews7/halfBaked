#' Run Gene Set Enrichment Analysis (GSEA) with arbitrary gene sets.
#'
#' This function performs GSEA on a named list of ranked genes.
#'
#' @param sigs A named list of gene signatures.
#' @param ranked.genes A named list of ranked genes.
#'   Using a list allows for results from many comparisons to be passed.
#' @param outdir The output directory for results.
#' @param res.name A name to be assigned to for the results added to `res.list`.
#'   Also used as a prefix for output files.
#' @param res.list A list to store results. Useful for saving lots of results to a single file,
#'   e.g. from multiple comparisons.
#'   Defaults to an empty list.
#' @param padj.th The adjusted p-value threshold for limiting which gene sets to plot.
#'   Defaults to 0.05.
#' @param min.size The minimum size of the gene sets to be considered for GSEA.
#'   Defaults to 15.
#' @param max.size The maximum size of the gene sets to be considered for GSEA.
#'   Defaults to 1000.
#' @param BPPARAM The BiocParallelParam object specifying the parallel back-end to be used.
#'   Defaults to NULL.
#' @param ... Additional arguments to pass to the [fgsea::fgsea()] function.
#'
#' @return A list of GSEA results if `res.list` is not `NULL``, otherwise, results are saved as files in the specified output directory.
#'
#' @details The function creates various output files including detailed GSEA results, enrichment plots, and tables of top enriched pathways.
#'
#' @importFrom data.table fwrite
#' @importFrom fgsea fgsea plotEnrichment plotGseaTable
#' @importFrom stringi stri_sub<-
#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 labs theme annotate element_text element_blank element_rect
#'
#' @export
#'
#' @examples
#' \dontrun{
#' run_GSEA(
#'     sigs = msigdb, ranked.genes = my_genes, outdir = "./results", outprefix = "experiment1",
#'     res.list = list()
#' )
#' }
#'
#' @author Jared Andrews
run_GSEA <- function(
    sigs,
    ranked.genes,
    outdir,
    res.name = "GSEA",
    res.list = list(),
    padj.th = 0.05,
    min.size = 15,
    max.size = 1000,
    BPPARAM = NULL,
    ...) {

    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    fgsea_res <- fgsea(
        pathways = sigs,
        stats = ranked.genes,
        eps = 1e-100,
        minSize = min.size,
        maxSize = max.size,
        BPPARAM = BPPARAM,
        ...
    )
    fgsea_res <- fgsea_res[order(padj), ]

    # Save full results.
    fwrite(fgsea_res,
        file = file.path(outdir, paste0(res.name, ".fgseaRes.txt")),
        sep = "\t", sep2 = c("", " ", "")
    )

    res.list[[res.name]] <- fgsea_res

    # Make figures
    fsig <- fgsea_res$pathway[fgsea_res$padj < padj.th]

    if (length(fsig) == 0) {
        message("No gene sets with p.adj < ", padj.th, " found, skipping plotting.")
        return(res.list)
    }

    plots <- list()
    for (f in seq_along(fsig)) {
        pathw <- fsig[f]
        if (!is.na(pathw)) {
            # reposition annotation depending on curve shape.
            if (!is.na(fgsea_res$NES[f]) && fgsea_res$NES[f] < 0) {
                maxx <- length(ranked.genes)
                maxy <- max(fgsea_res$ES[f])
                x <- round(maxx * 0.001)
                y <- maxy + (abs(maxy) * 0.001)
            } else {
                maxx <- length(ranked.genes)
                maxy <- max(fgsea_res$ES[f])
                x <- maxx - round(maxx * 0.001)
                y <- maxy - (abs(maxy) * 0.001)
            }

            # For those really long titles.
            tt <- pathw
            if (nchar(tt) > 40) {
                stri_sub(tt, 48, 47) <- "\n"
            }

            p <- plotEnrichment(
                sigs[[pathw]],
                ranked.genes
            ) + labs(title = tt) +
                theme(plot.title = element_text(size = 7),
                panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())

            # Add stats to plot.
            p <- p + annotate("text", x, y,
                label = paste0(
                    "p.val = ",
                    formatC(fgsea_res$pval[f], format = "e", digits = 2),
                    "\np.adj = ",
                    formatC(fgsea_res$padj[f], format = "e", digits = 2),
                    "\nNES = ",
                    round(fgsea_res$NES[f], digits = 2)
                ),
                vjust = "inward", hjust = "inward", size = 3
            )

            plots[[f]] <- p
        }
    }

    pdf(file.path(outdir, paste0(res.name, ".padj", padj.th, ".Swoops.pdf")), height = 10, width = 20)
    # Calculate how many pages to print assuming max 24 plots per page.
    pages <- ceiling(length(plots) / 24)
    # Print each page.
    for (i in seq(pages)) {
        end <- i * 24
        start <- end - 23
        if (end > length(plots)) {
            end <- length(plots)
        }
        grid.arrange(grobs = plots[start:end], nrow = 4, ncol = 6)
    }
    dev.off()

    pdf(file.path(outdir, paste0(res.name, ".Top10.padj", padj.th, ".pdf")), width = 12)
    top_up <- unlist(fgsea_res[ES > 0 & padj < padj.th][head(order(pval), n = 10), "pathway"])
    top_down <- unlist(fgsea_res[ES < 0 & padj < padj.th][head(order(pval), n = 10), "pathway"])
    top_pathways <- unique(c(top_up, rev(top_down)))

    p <- plotGseaTable(sigs[top_pathways], ranked.genes, fgsea_res)
    print(p)
    dev.off()

    return(res.list)
}


#' Plot top Gene Set Enrichment Analysis (GSEA) results in a barplot.
#'
#' This function condenses GSEA results by plotting the top significant gene sets.
#'
#' @param gsea.res A data.frame of GSEA results as returned by [fgsea::fgsea()].
#' @param genesets.name A name for the gene set collection or group, used to label plot.
#' @param padj.th The significance threshold (adjusted p-value) for filtering gene sets
#'   Defaults to 0.05.
#' @param top The number of top significant gene sets in each direction to consider.
#'   Defaults to 35.
#' @param color.by.direction Logical. If `TRUE`, bars are colored by up/down direction instead of significance.
#'   Defaults to FALSE.
#' @param up.color Color for upregulated pathways if `color.by.direction` is `TRUE`.
#'   Defaults to "red".
#' @param down.color Color for downregulated pathways if `color.by.direction` is `TRUE`.
#'   Defaults to "blue".
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_col coord_flip labs theme_bw ylim theme scale_x_discrete
#'   element_text scale_fill_manual element_blank geom_hline element_line
#' @importFrom viridis scale_fill_viridis
#' @importFrom stringr str_trunc
#' @importFrom tibble as_tibble
#' @importFrom dplyr arrange filter
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' \dontrun{
#' GSEA_barplot(gsea.list = my_gsea_results, outdir = "./summary", padj.th = 0.01, top = 50)
#' }
#'
#' @author Jared Andrews
GSEA_barplot <- function(gsea.res, genesets.name = NULL, padj.th = 0.05, top = 35,
                         color.by.direction = FALSE, up.color = "red", down.color = "blue") {

    df_sub <- gsea.res[gsea.res$padj < padj.th, ]

    if (nrow(df_sub) > top) {
        df_sub <- df_sub %>%
            as_tibble() %>%
            arrange(padj)

        top_up <- df_sub %>%
            filter(NES > 0) %>%
            head(top)

        top_down <- df_sub %>%
            filter(NES < 0) %>%
            head(top)

        df_sub <- rbind(top_up, top_down)
    }   

    if (nrow(df_sub) > 0) {
        df_sub <- df_sub %>%
            as_tibble() %>%
            arrange(desc(NES))

        if (color.by.direction) {
            df_sub$direction <- ifelse(df_sub$NES > 0, "Up", "Down")
            p <- ggplot(df_sub, aes(reorder(pathway, -NES), NES)) +
                geom_col(aes(fill = direction)) +
                scale_fill_manual(values = c("Up" = up.color, "Down" = down.color))
        } else {
            p <- ggplot(df_sub, aes(reorder(pathway, -NES), NES)) +
                geom_col(aes(fill = -log10(padj))) +
                scale_fill_viridis()
        }   

        p <- p +
            coord_flip() +
            labs(
                x = NULL, y = "Normalized Enrichment Score",
                title = paste0(genesets.name, " - Top ", top, "\np.adj < ", padj.th)
            ) +
            theme_bw() +
            ylim(-max(abs(df_sub$NES) + 0.1), max(abs(df_sub$NES) + 0.1)) +
            theme(axis.text.y = element_text(size = 6, color ="black"), 
            plot.title = element_text(size = 10),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(color ="black"),
            axis.ticks = element_line(color = "black")) +
            scale_x_discrete(label = function(x) str_trunc(x, 55)) + 
            geom_hline(yintercept = 0, size = 0.3)
    } else {
        stop("No pathways with p.adj < ", padj.th, " found.")
    }

    p
}