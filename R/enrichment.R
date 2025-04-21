#' Run Enrichment Analysis
#'
#' This function performs enrichment analysis using KEGG, Reactome, GO, or a custom universal method.
#' It accepts a set of genes and background genes for enrichment.
#'
#' @details
#' This function will not return results for analyses with no significant results.
#'
#' @param genes Character vector of gene IDs for simple enrichment.
#'   Default is `NULL`.
#' @param bg Character vector of gene IDs to be used as background for simple enrichment.
#'   Default is `NULL`.
#' @param res.name Output prefix name for results.
#' @param TERM2GENE data.frame of two columns, the first for the term and the second for the gene ID.
#'   Each gene in each geneset gets its own row (long format). Gene identifiers should be ENTREZID.
#'   Required if `method` is "universal".
#' @param res.list Named list to which enrichment results should be added.
#'   Default is an empty list.
#' @param method Enrichment method to use. Options are "KEGG", "Reactome", "GO", or "universal".
#' @param species Species to use. Options are "human" or "mouse".
#'   This determines the organism values for KEGG and Reactome.
#' @param ont Character vector of GO ontologies to test, options include "BP", "MF", "CC", and "ALL".
#'   Default is `c("BP", "MF", "CC", "ALL")`.
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
#' @return Named list of enrichment results.
#'
#' @importFrom clusterProfiler compareCluster enrichGO enrichKEGG enricher bitr setReadable
#'   dotplot cnetplot
#' @importFrom enrichplot pairwise_termsim treeplot
#' @importFrom ReactomePA enrichPathway
#' @importFrom ggplot2 rel
#'
#' @export
#'
#' @author Jared Andrews
run_enrichment <- function(
    genes,
    bg,
    res.name,
    TERM2GENE = NULL,
    res.list = list(),
    method = c("KEGG", "Reactome", "GO", "universal"),
    species = c("human", "mouse"),
    ont = c("BP", "MF", "CC", "ALL"),
    sig.th = 0.05,
    sig.col = "padj",
    lfc.th = 0,
    lfc.col = "log2FoldChange",
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

    if (id.type == "ENSEMBL") {
        genes <- sapply(strsplit(as.character(genes), "\\."), `[`, 1)
        bg <- sapply(strsplit(as.character(bg), "\\."), `[`, 1)
    }

    bg <- bitr(bg, fromType = id.type, toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID
    genes <- bitr(genes, fromType = id.type, toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID

    if (method == "GO") {
        for (ont_item in ont) {
            ego <- enrichGO(genes, OrgDb = OrgDb, universe = bg, ont = ont_item, readable = TRUE, ...)

            if (nrow(as.data.frame(ego)) > 0) {
                ego <- pairwise_termsim(ego)
                res.list[[paste0(res.name, ".GO", ont_item)]] <- ego
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
        if (method == "universal") params$TERM2GENE <- TERM2GENE

        ego <- do.call(enrich_fun, c(list(gene = genes, universe = bg), params))

        if (nrow(as.data.frame(ego)) > 0) {
            ego <- pairwise_termsim(ego)
            res.list[[paste0(res.name, ".", method)]] <- ego
        }
    }

    res.list
}


#' Retrieve Genes Associated with GO Terms Containing a Specific Search Term
#'
#' This function searches for Gene Ontology (GO) terms that contain a specified search term
#' and retrieves all associated genes for the specified species and ID type.
#'
#' @param search_term A character string specifying the term to search for within GO Biological Process terms (case-insensitive).
#' @param orgdb The organism-specific database package to use for gene mapping.
#'   This should be one of the organism packages like "org.Hs.eg.db", "org.Mm.eg.db", etc.
#' @param id_type A character string specifying the type of gene identifier to return.
#'   Options include "SYMBOL", "ENTREZID", and "ENSEMBL". Default is "SYMBOL".
#'
#' @return A named list containing:
#'   - genes - A character vector of gene identifiers of the specified type associated with GO terms that contain the search term.
#'      The names of the vector are the corresponding Entrez Gene IDs (if `id_type` is not "ENTREZID").
#'   - go_terms - A character vector of GO terms that matched the search term.
#'
#' @details
#' The function performs the following steps:
#'   - Retrieves all GO terms and their descriptions.
#'   - Searches for GO terms that include the specified search term.
#'   - Retrieves all Entrez Gene IDs associated with the matching GO terms.
#'   - Maps Entrez Gene IDs to the specified type of gene identifier.
#'
#' @examples
#' # Retrieve human gene symbols associated with GO terms containing "WNT"
#' genes_wnt_human <- get_genes_by_go_term("WNT", "org.Hs.eg.db", id_type = "SYMBOL")
#' head(genes_wnt_human$genes)
#' head(genes_wnt_human$go_terms)
#'
#' @author Jared Andrews
#'
#' @export
get_genes_by_go_term <- function(search_term, orgdb, id_type = "SYMBOL") {
    # Check if GO.db and AnnotationDbi packages are installed
    for (pk in c("GO.db", "AnnotationDbi", orgdb)) {
        .package_check(pk)
    }

    # Get all GO terms
    go_terms <- as.list(GO.db::GOTERM)

    # Extract GO IDs and their associated terms
    go_ids <- names(go_terms)
    go_terms_text <- character(length(go_terms))

    for (i in seq_along(go_terms)) {
        go_terms_text[i] <- go_terms[[i]]@Term
    }

    # Search for GO terms that include the search term (case-insensitive)
    indices <- grep(search_term, go_terms_text, ignore.case = TRUE)
    matched_go_ids <- go_ids[indices]

    # Get the names of the matched GO terms
    matched_go_terms <- go_terms_text[indices]
    names(matched_go_terms) <- matched_go_ids

    # Retrieve genes associated with these GO IDs
    # Construct the name of the GO to All Genes mapping object
    suppressPackageStartupMessages(require(OrgDb, character.only = TRUE))
    org_prefix <- sub("\\.db$", "", orgdb) # Remove ".db" from package name
    go2allels_name <- paste0(org_prefix, "GO2ALLEGS")
    go2allels <- get(go2allels_name)

    genes_entrez_list <- as.list(go2allels)[matched_go_ids]

    # Flatten the list and remove NAs
    genes_entrez <- unique(unlist(genes_entrez_list))
    genes_entrez <- genes_entrez[!is.na(genes_entrez)]

    # Map Entrez Gene IDs to the specified ID type
    # Get the organism-specific database object
    org_db <- get(orgdb)

    # Check if the requested id_type is valid
    valid_id_types <- AnnotationDbi::columns(org_db)
    if (!(id_type %in% valid_id_types)) {
        stop("Invalid 'id_type'. Valid options are: ", paste(valid_id_types, collapse = ", "))
    }

    # If id_type is ENTREZID, simply return the Entrez IDs
    if (id_type == "ENTREZID") {
        genes_ids <- genes_entrez
        names(genes_ids) <- genes_entrez
    } else {
        genes_ids <- AnnotationDbi::mapIds(
            org_db,
            keys = genes_entrez,
            column = id_type,
            keytype = "ENTREZID",
            multiVals = "first"
        )
    }

    list(genes = genes_ids,
         go_terms = matched_go_terms)
}
