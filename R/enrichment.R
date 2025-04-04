#' Run Enrichment Analysis
#'
#' This function performs enrichment analysis using KEGG, Reactome, GO, or a custom universal method.
#' It accepts either a list of differential expression results for group enrichment or
#' a set of genes and background for simple enrichment.
#'
#' @param res.list Named list of DE results data.frames for group enrichment.
#'   If provided, `genes` and `bg` will be ignored.
#'   Default `is NULL`.
#' @param genes Character vector of gene IDs for simple enrichment.
#'   Default `is NULL`.
#' @param bg Character vector of gene IDs to be used as background for simple enrichment.
#'   Default `is NULL`.
#' @param method Enrichment method to use. Options are "KEGG", "Reactome", "GO", or "universal".
#' @param species Species to use. Options are "human" or "mouse".
#'   This determines the organism values for KEGG and Reactome.
#' @param ont Character vector of GO ontologies to test, options include "BP", "MF", "CC", and "ALL".
#'   Default is `c("BP", "MF", "CC", "ALL")`.
#' @param name Output prefix name for simple enrichment if `genes` and `bg` provided.
#' @param sig.th Significance threshold for DE genes.
#'   Default is 0.05.
#' @param sig.col Name of the column in the results data.frame containing significance values to use.
#'   Default is "padj".
#' @param lfc.th Log fold change threshold for DE genes.
#'   Default is 0.
#' @param lfc.col Name of the column in the results data.frame containing log2 fold change values.
#'   Default is "log2FoldChange".
#' @param outdir Output directory (default: "./enrichments").
#' @param OrgDb Annotation database to use (default: "org.Hs.eg.db").
#' @param id.col Name of gene ID column (default: "ENSEMBL").
#' @param id.type Type of gene ID used (default: "ENSEMBL").
#' @param ... Additional arguments passed to the enrichment functions.
#' @return Saves enrichment results and plots to the output directory; returns invisible objects.
#'
#' @importFrom clusterProfiler compareCluster enrichGO enrichKEGG enricher bitr setReadable
#'   emapplot dotplot cnetplot
#' @importFrom enrichplot pairwise_termsim treeplot
#' @importFrom ReactomePA enrichPathway
#'
#' @export
#'
#' @author Jared Andrews
run_enrichment <- function(
    res.list = NULL,
    genes = NULL,
    bg = NULL,
    method = c("KEGG", "Reactome", "GO", "universal"),
    species = c("human", "mouse"),
    ont = c("BP", "MF", "CC", "ALL"),
    name = "sample",
    sig.th = 0.05,
    sig.col = "padj",
    lfc.th = 0,
    lfc.col = "log2FoldChange",
    outdir = "./enrichments",
    OrgDb = "org.Hs.eg.db",
    id.col = "ENSEMBL",
    id.type = "ENSEMBL",
    ...) {
    species <- match.arg(species)
    method <- match.arg(method)

    # Set organism values based on species
    if (species == "human") {
        kegg_org <- "hsa"
        reactome_org <- "human"
    } else { # species == "mouse"
        kegg_org <- "mmu"
        reactome_org <- "mouse"
    }

    # Cluster enrichment mode: using a list of DE results
    if (!is.null(res.list)) {
        for (comp in names(res.list)) {
            df <- res.list[[comp]]
            comp_name <- comp
            if (lfc.th != 0) comp_name <- paste0(comp_name, "-LFC", round(lfc.th, 3), "filt")
            out_path <- file.path(outdir, comp_name)
            dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

            if (id.type == "ENSEMBL") {
                df[[id.col]] <- sapply(strsplit(as.character(df[[id.col]]), "\\."), `[`, 1)
            }
            df <- df[!is.na(df[[sig.col]]), ]

            gene_list <- bitr(df[[id.col]], fromType = id.type, toType = "ENTREZID", OrgDb = OrgDb)
            gene_list$FC <- df[[lfc.col]][match(gene_list[[id.col]], df[[id.col]])]
            gl <- sort(setNames(gene_list$FC, gene_list$ENTREZID), decreasing = TRUE)

            gene_clusters <- list(
                up = df[[id.col]][df[[sig.col]] < sig.th & df[[lfc.col]] > lfc.th],
                down = df[[id.col]][df[[sig.col]] < sig.th & df[[lfc.col]] < -lfc.th],
                all_de = df[[id.col]][df[[sig.col]] < sig.th]
            )
            skip <- FALSE
            tryCatch(
                {
                    gene_clusters$up <- bitr(gene_clusters$up, fromType = id.type, toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID
                    gene_clusters$down <- bitr(gene_clusters$down, fromType = id.type, toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID
                    gene_clusters$all_de <- bitr(gene_clusters$all_de, fromType = id.type, toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID
                },
                error = function(e) {
                    message("Error mapping IDs in ", comp_name, ": ", e)
                    skip <<- TRUE
                }
            )
            if (skip) next

            bg_ids <- bitr(df[[id.col]], fromType = id.type, toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID

            if (method == "GO") {
                for (ont_item in ont) {
                    ck <- compareCluster(
                        geneCluster = gene_clusters, fun = enrichGO,
                        universe = bg_ids, ont = ont_item, OrgDb = OrgDb, ...
                    )
                    if (!is.null(ck)) {
                        ck <- setReadable(ck, OrgDb = OrgDb, keyType = "ENTREZID")
                        ego <- pairwise_termsim(ck)
                        pdf(file.path(out_path, paste0("GO_Enrichments.Top20_", ont_item, ".pdf")),
                            width = 6, height = 4 + 0.025 * length(ego@compareClusterResult$Cluster)
                        )
                        print(dotplot(ego, showCategory = 20, font.size = 7))
                        print(dotplot(ego, size = "count", showCategory = 20, font.size = 7))
                        dev.off()
                        pdf(file.path(out_path, paste0("GO_Enrichments.termsim.Top20_", ont_item, ".pdf")),
                            width = 9, height = 9
                        )
                        print(emapplot(ego,
                            pie = "count", cex_category = 0.9,
                            cex_label_category = 0.9, layout = "kk", repel = TRUE, showCategory = 20
                        ))
                        dev.off()
                        saveRDS(ego, file = file.path(out_path, paste0("enrichGO_", ont_item, ".RDS")))
                        write.table(as.data.frame(ego),
                            file = file.path(out_path, paste0("enrichGO_", ont_item, ".txt")),
                            sep = "\t", row.names = FALSE, quote = FALSE
                        )
                    }
                }
            } else {
                enrich_fun <- switch(method,
                    KEGG = enrichKEGG,
                    Reactome = enrichPathway,
                    universal = enricher
                )
                params <- list(...)
                if (method == "KEGG") params$organism <- kegg_org
                if (method == "Reactome") params$organism <- reactome_org
                ck <- do.call(compareCluster, c(list(
                    geneCluster = gene_clusters,
                    fun = enrich_fun, universe = bg_ids,
                    keyType = if (method == "KEGG") "kegg" else "ENTREZID"
                ), params))
                if (!is.null(ck)) {
                    if (method %in% c("KEGG", "Reactome")) {
                        ck <- setReadable(ck, OrgDb = OrgDb, keyType = "ENTREZID")
                    }
                    ego <- pairwise_termsim(ck)
                    pdf(file.path(out_path, paste0(method, "_Enrichments.Top20.pdf")),
                        width = 6, height = 4 + 0.015 * length(ego@compareClusterResult$Cluster)
                    )
                    print(dotplot(ego, showCategory = 20, font.size = 7))
                    print(dotplot(ego, size = "count", showCategory = 20, font.size = 7))
                    dev.off()
                    pdf(file.path(out_path, paste0(method, "_Enrichments.termsim.Top20.pdf")),
                        width = 9, height = 9
                    )
                    print(emapplot(ego,
                        pie = "count", cex_category = 0.9,
                        cex_label_category = 0.9, layout = "kk", repel = TRUE, showCategory = 20
                    ))
                    dev.off()
                    if (nrow(as.data.frame(ego)) > 2) {
                        pdf(file.path(out_path, paste0(method, "_Enrichments.termsim.Top30_Tree.pdf")),
                            width = 17, height = 14
                        )
                        print(treeplot(ego, showCategory = 30, fontsize = 4))
                        dev.off()
                        pdf(file.path(out_path, paste0(method, "_Enrichments.termsim.Top10_FullNet.pdf")),
                            width = 15, height = 15
                        )
                        print(cnetplot(ego, showCategory = 10, layout = "kk"))
                        dev.off()
                        pdf(file.path(out_path, paste0(method, "_Enrichments.termsim.Top5_FullNet.pdf")),
                            width = 12, height = 12
                        )
                        print(cnetplot(ego, showCategory = 5, layout = "kk"))
                        dev.off()
                    }
                    saveRDS(ego, file = file.path(out_path, paste0(method, "_results.RDS")))
                    write.table(as.data.frame(ego),
                        file = file.path(out_path, paste0(method, "_results.txt")),
                        sep = "\t", row.names = FALSE, quote = FALSE
                    )
                }
            }
        }
    } else if (!is.null(genes) && !is.null(bg)) {
        # Simple enrichment mode
        out_path <- file.path(outdir, name)
        dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
        if (id.type == "ENSEMBL") {
            genes <- sapply(strsplit(as.character(genes), "\\."), `[`, 1)
            bg <- sapply(strsplit(as.character(bg), "\\."), `[`, 1)
        }
        skip <- FALSE
        tryCatch(
            {
                genes <- bitr(genes, fromType = id.type, toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID
            },
            error = function(e) {
                message("Error mapping gene IDs: ", e)
                skip <<- TRUE
            }
        )
        if (skip) {
            return(NULL)
        }
        bg <- bitr(bg, fromType = id.type, toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID
        if (method == "GO") {
            res_list <- list()
            for (ont_item in ont) {
                ego <- enrichGO(genes, OrgDb = OrgDb, universe = bg, ont = ont_item, readable = TRUE, ...)
                if (nrow(as.data.frame(ego)) > 0) {
                    ego <- pairwise_termsim(ego)
                    pdf(file.path(out_path, paste0("enrichGO_", ont_item, ".Top30.pdf")),
                        width = 6, height = 8
                    )
                    print(dotplot(ego, showCategory = 30, font.size = 7))
                    print(barplot(ego, showCategory = 30, font.size = 7))
                    dev.off()
                    if (nrow(as.data.frame(ego)) > 2) {
                        pdf(file.path(out_path, paste0("enrichGO_", ont_item, ".termsim.Top30_Tree.pdf")),
                            width = 17, height = 14
                        )
                        print(treeplot(ego, showCategory = 30, fontsize = 4))
                        dev.off()
                        pdf(file.path(out_path, paste0("enrichGO_", ont_item, ".termsim.Top10_FullNet.pdf")),
                            width = 15, height = 15
                        )
                        print(cnetplot(ego, showCategory = 10, layout = "kk"))
                        dev.off()
                        pdf(file.path(out_path, paste0("enrichGO_", ont_item, ".termsim.Top5_FullNet.pdf")),
                            width = 12, height = 12
                        )
                        print(cnetplot(ego, showCategory = 5, layout = "kk"))
                        dev.off()
                    }
                    saveRDS(ego, file = file.path(out_path, paste0("enrichGO_", ont_item, ".RDS")))
                    write.table(as.data.frame(ego),
                        file = file.path(out_path, paste0("enrichGO_", ont_item, ".txt")),
                        sep = "\t", row.names = FALSE, quote = FALSE
                    )
                    res_list[[ont_item]] <- ego
                }
            }
            return(invisible(res_list))
        } else {
            enrich_fun <- switch(method,
                KEGG = enrichKEGG,
                Reactome = enrichPathway,
                universal = enricher
            )
            params <- list(...)
            if (method == "KEGG") params$organism <- kegg_org
            if (method == "Reactome") params$organism <- reactome_org
            ego <- do.call(enrich_fun, c(list(gene = genes, universe = bg), params))
            if (nrow(as.data.frame(ego)) > 0) {
                ego <- pairwise_termsim(ego)
                pdf(file.path(out_path, paste0(method, ".Top30.pdf")), width = 6, height = 8)
                print(dotplot(ego, showCategory = 30, font.size = 7))
                print(barplot(ego, showCategory = 30, font.size = 7))
                dev.off()
                if (nrow(as.data.frame(ego)) > 2) {
                    pdf(file.path(out_path, paste0(method, ".termsim.Top30_Tree.pdf")), width = 17, height = 14)
                    print(treeplot(ego, showCategory = 30, fontsize = 4))
                    dev.off()
                    pdf(file.path(out_path, paste0(method, ".termsim.Top10_FullNet.pdf")), width = 15, height = 15)
                    print(cnetplot(ego, showCategory = 10, layout = "kk"))
                    dev.off()
                    pdf(file.path(out_path, paste0(method, ".termsim.Top5_FullNet.pdf")), width = 12, height = 12)
                    print(cnetplot(ego, showCategory = 5, layout = "kk"))
                    dev.off()
                }
                saveRDS(ego, file = file.path(out_path, paste0(method, ".RDS")))
                write.table(as.data.frame(ego),
                    file = file.path(out_path, paste0(method, ".txt")),
                    sep = "\t", row.names = FALSE, quote = FALSE
                )
                return(invisible(ego))
            }
        }
    } else {
        stop("Provide either a 'res.list' for cluster enrichment or both 'genes' and 'bg' for simple enrichment.")
    }
}
