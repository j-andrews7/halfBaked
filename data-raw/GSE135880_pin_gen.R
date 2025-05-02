library(halfBaked)
library(DESeq2)
library(edgeR)
library(SummarizedExperiment)
library(pins)
library(msigdbr)
library(BiocParallel)
library(rrvgo)
library(org.Mm.eg.db)

# Original FASTQs were downloaded from GEO and processed with the nf-core v3.12.0 pipeline:
# nextflow run nf-core/rnaseq -r 3.12.0 -profile singularity -c "$BAKER_REF"/nf_configs/rnaseq.config -w /scratch_space/jandrews/"$LSB_JOBNAME" \
# --outdir ./nfcore_mm10 --email jared.andrews@stjude.org --input nfcore_rnaseq.samplesheet.csv --gencode --genome MM10 --aligner star_salmon \
# --pseudo_aligner salmon --max_memory 128.GB --skip_stringtie --max_multiqc_email_size 15.MB -resume

org.db <- org.Mm.eg.db

description <- "
  **Wang J. et al SciAdv 2020 - GSE135880 - Mouse OPCs with _Eed_ KO - RNA-seq Data**

  This dataset contains O4+ immunopanned oligodendrocyte precursor cells (OPCs) from cortices of P5 or P6 Eed KO or control mice.
  See the [associated publication](https://www.science.org/doi/10.1126/sciadv.aaz6477) for more details.

  Please see the Pin generation code to view how this data was processed.
  "

# Load sample metadata.
meta <- read.csv("nfcore_rnaseq.samplesheet.csv", header = TRUE,
                 stringsAsFactors = TRUE)

# Drop FASTQ file locations.
meta <- meta[, !colnames(meta) %in% c("fastq_1", "fastq_2")]

# Load counts. This object was generated using tximport via the nf-core
# RNA-seq pipeline on the salmon quants and
# is appropriate for pretty much all downstream DE packages (DESeq2, edgeR, limma).
cts <- read.table("salmon.merged.gene_counts_length_scaled.tsv", header = TRUE,
                  sep = "\t", stringsAsFactors = FALSE)

# Counts table has first two columns as gene IDs and gene symbols.
genes <- cts[, 1:2]
names(genes) <- c("ENSEMBL", "SYMBOL")
rownames(cts) <- cts[, 1]

# Remove the gene version info from the ENSEMBL IDs
genes$ENSEMBL <- gsub("\\..*", "", genes$ENSEMBL)

# 2) Using mapIds() to get a named vector of ENTREZ IDs
genes$ENTREZ <- mapIds(org.db,
                       keys=genes$ENSEMBL,
                       column="ENTREZID",
                       keytype="ENSEMBL",
                       multiVals="first")

# Set metadata rownames and ensure they match count column names.
rownames(meta) <- meta$sample

# Also carry along TPMs as an additional assay.
tpms <- read.table("salmon.merged.gene_tpm.tsv", header = TRUE,
                   sep = "\t", stringsAsFactors = FALSE)
rownames(tpms) <- tpms[, 1]
tpms <- tpms[, rownames(meta)]

cts <- cts[, rownames(meta)]

# Create a SummarizedExperiment object.
se <- SummarizedExperiment(
    assays = list(counts = round(as.matrix(cts)),
                  tpm = as.matrix(tpms),
                  log2tpm = log2(as.matrix(tpms) + 1)),
    colData = meta,
    rowData = genes
)

# Limit to reasonably expressed genes, adjust design or use `group` as needed.
design <- model.matrix(~0 + Group, data = colData(se))
keep <- filterByExpr(se, design = design)
se <- se[keep, ]

# Add experiment data to metadata.
metadata(se) <- list(description = description)

### Differential Expression Analysis ###

# Design is just being set to no design here (~ 1).
# The function will take care of setting this.
dds <- DESeqDataSet(se, design = ~1)
dds <- DESeq(dds)

# To hold the results for each comparison.
res <- list()

# Whatever can be used as name of contrast vectors.
contrasts <- list(
    "Eed_cKO.v.Control" = c("Group", "Eed_cKO", "Control")
)

res <- get_DESeq2_res(dds,
    res.list = res, contrasts = contrasts,
    add.rowData = c("ENSEMBL", "SYMBOL")
)

# Add various normalized counts
assay(se, "vst") <- vst(round(assay(se, "counts")))
assay(se, "cpm") <- cpm(se)
assay(se, "log2cpm") <- cpm(se, log = TRUE)

# Cram DE results into SE metadata.
se.meta <- metadata(se)
se.meta$DESeq2.Results <- res
metadata(se) <- se.meta

### Enrichment Analysis ###

# Set species and database stuff for enrichments
orgdb <- "org.Mm.eg.db"
msig.species <- "Mus musculus"


# Edit if you'd like
sig_th <- 0.05
lfc_th <- 0
sig_col <- "padj"
lfc_col <- "log2FoldChange"

kegg_res <- list()
rct_res <- list()
go_res <- list()

for (rez in names(res)) {
    df <- res[[rez]]

    curr_lfc_th <- lfc_th

    # Handle instances where shrunken LFCs are lower than the LFC testing threshold.
    # An annoying quirk of DESeq2 LFC shrinkage, just as stupid as it looks.
    if (grepl("-shLFC0|-LFC0", rez)) {
        curr_lfc_th <- as.numeric(unlist(strsplit(rez, "LFC"))[2])
    }

    bg <- df$ENSEMBL[!is.na(df[[sig_col]])]

    up_genes <- df$ENSEMBL[!is.na(df[[sig_col]]) &
        df[[sig_col]] < sig_th &
        df[[lfc_col]] > curr_lfc_th]

    dn_genes <- df$ENSEMBL[!is.na(df[[sig_col]]) &
        df[[sig_col]] < sig_th &
        df[[lfc_col]] < -curr_lfc_th]

    if (length(up_genes) > 0) {
        kegg_res <- run_enrichment(up_genes, bg,
            res.name = paste0(rez, "_up"),
            res.list = kegg_res, OrgDb = orgdb,
            method = "KEGG", species = "mouse"
        )

        rct_res <- run_enrichment(up_genes, bg,
            res.name = paste0(rez, "_up"),
            res.list = rct_res, OrgDb = orgdb,
            method = "Reactome", species = "mouse"
        )

        go_res <- run_enrichment(up_genes, bg,
            res.name = paste0(rez, "_up"),
            res.list = go_res, OrgDb = orgdb,
            method = "GO", species = "mouse"
        )
    }

    if (length(dn_genes) > 0) {
        kegg_res <- run_enrichment(dn_genes, bg,
            res.name = paste0(rez, "_dn"),
            res.list = kegg_res, OrgDb = orgdb,
            method = "KEGG", species = "mouse"
        )

        rct_res <- run_enrichment(dn_genes, bg,
            res.name = paste0(rez, "_dn"),
            res.list = rct_res, OrgDb = orgdb,
            method = "Reactome", species = "mouse"
        )

        go_res <- run_enrichment(dn_genes, bg,
            res.name = paste0(rez, "_dn"),
            res.list = go_res, OrgDb = orgdb,
            method = "GO", species = "mouse"
        )
    }
}

full_res <- c(kegg_res, rct_res, go_res)

# Excel sheets can only be 31 characters, so it is necessary to shorten the names
names(full_res) <- gsub("Reactome", "RCT", names(full_res))
names(full_res) <- gsub("Control", "Ctl", names(full_res))
names(full_res) <- gsub("\\.v\\.", "v", names(full_res))

# Save to SE object.
metadata(se)$DESeq2.enrichments <- full_res

### Cluster GO Term Enrichment Results ###

onts <- c("BP", "CC", "MF")

up.col <- "#56B4E9"

gobp_res <- grep("GOBP", names(full_res))
gomf_res <- grep("GOMF", names(full_res))
gocc_res <- grep("GOCC", names(full_res))

clustered_res_size <- list()
clustered_res_uniqueness <- list()
clustered_res_log10_padjscore <- list()

for (ont in onts) {
    # Get the GO terms for the current ontology
    if (ont == "BP") {
        curr_go <- gobp_res
    } else if (ont == "CC") {
        curr_go <- gocc_res
    } else {
        curr_go <- gomf_res
    }

    godat <- GOSemSim::godata(orgdb, ont = ont)

    # Ignore these words, alter as wanted
    stoppers <- tm::stopwords(kind = "en")
    stoppers <- c(
        stoppers, "regulation", "positive",
        "negative", "process", "cell", "activity"
    )

    for (i in names(full_res)[curr_go]) {
        enr.up <- full_res[[i]]

        # Collapse similar terms, score by term size, uniqueness, or significance
        if (!is.null(enr.up)) {
            if (nrow(enr.up) > 2) {
                dir.create(file.path("enrichments", "reduced"),
                    recursive = TRUE, showWarnings = FALSE
                )

                enr.up.sim <- calculateSimMatrix(enr.up$ID,
                    orgdb = orgdb,
                    ont = ont,
                    semdata = godat,
                    method = "Rel"
                )

                enr.up.scores <- -log10(enr.up$p.adjust)
                names(enr.up.scores) <- enr.up$ID

                enr.up.red.size <- reduceSimMatrix(enr.up.sim,
                    scores = "size",
                    threshold = 0.7, orgdb = orgdb
                )

                enr.up.red.unique <- reduceSimMatrix(enr.up.sim,
                    scores = "uniqueness",
                    threshold = 0.7, orgdb = orgdb
                )

                enr.up.red.score <- reduceSimMatrix(enr.up.sim,
                    scores = enr.up.scores,
                    threshold = 0.7, orgdb = orgdb
                )
                
                clustered_res_size[[i]] <- enr.up.red.size
                clustered_res_uniqueness[[i]] <- enr.up.red.unique
                clustered_res_log10_padjscore[[i]] <- enr.up.red.score
            }
        }
    }
}

# Again, tack onto SE
metadata(se)$DESeq2.clustered_GO_enrichment <- list(
    size = clustered_res_size,
    uniqueness = clustered_res_uniqueness,
    log10_padjscore = clustered_res_log10_padjscore
)

### GSEA ###

# Set categories and subcategories that we want to retrieve, see msigdbr_collections()
# These vectors must be of equal length.
cats <- c("H", "C2", "C2", "C2", "C3", "C5", "C5", "C5", "C5", "C8")
subcats <- c(
    "", "CGP", "CP:KEGG", "CP:REACTOME", "TFT:GTRD",
    "GO:MF", "GO:BP", "GO:CC", "HPO", ""
)

msig_lists <- list()

# Get gene sets for each category and subcategory and stick in the msig.lists
for (i in seq_along(cats)) {
    cat <- cats[i]
    subcat <- subcats[i]

    msig <- msigdbr(species = msig.species, category = cat, subcategory = subcat)

    # Split into a named list
    msig_ls <- msig %>% split(x = .$gene_symbol, f = .$gs_name)

    # Stick in the list
    if (subcat == "") {
        outname <- cat
    } else {
        outname <- paste0(cat, "_", subcat)
    }

    msig_lists[[outname]] <- msig_ls
}

# Add as many comparisons as wanted here.
# Should match values in `names(dds.meta$DE.Results)`.
dges <- c("Eed_cKO.v.Control")

ranked_lists <- list()

for (d in dges) {
    out.res <- metadata(se)$DESeq2.Results[[d]]
    gsea.df <- out.res[!(is.na(out.res$padj)), ]

    gsea.df.gsea <- gsea.df$stat
    names(gsea.df.gsea) <- make.names(gsea.df$SYMBOL, unique = TRUE)

    ranked_lists[[d]] <- gsea.df.gsea
}

res_list <- list()

for (i in seq_along(ranked_lists)) {
    ranked_genes <- ranked_lists[[i]]
    ranked_name <- names(ranked_lists)[i]
    message(paste0("Running pre-ranked GSEA for: ", ranked_name))

    # Loop through the gene collections
    for (j in seq_along(msig_lists)) {
        msig_ls <- msig_lists[[j]]
        msig_name <- names(msig_lists)[j]

        # Remove the colon or it'll break file paths
        msig_name <- gsub(":", "_", msig_name)

        message(paste0("Using collection: ", msig_name))

        # Run GSEA
        res_list <- run_GSEA(msig_ls, ranked_genes,
            outdir = "./GSEA/", res.name = paste0(ranked_name, ".", msig_name),
            res.list = res_list
        )
    }
}

# Excel sheets can only be 31 chars, so it is often necessary to shorten the names
names(res_list) <- gsub("REACTOME", "RCT", names(res_list))

# Add to metadata of the SummarizedExperiment object
metadata(se)$DESeq2.GSEA <- res_list

###### edgeR ######
res <- list()

res <- get_edgeR_res(se,
                     res.list = res, contrasts = contrasts,
)

# Cram DE results into SE metadata.
se.meta <- metadata(se)
se.meta$edgeR.Results <- res
metadata(se) <- se.meta

# Write to pin.
board <- board_folder("../pkgdown/assets/pins-board")
pin_write(
    board = board,
    x = se,
    type = "rds",
    name = "GSE135880_SummarizedExperiment_mm10",
    title = "GSE135880_SummarizedExperiment_mm10",
    description = "A SummarizedExperiment containing counts for O4+ immunopanned oligodendrocyte precursor cells (OPCs) from cortices of P5 or P6 Eed KO or control mice, along with DESeq2 analysis results in the metadata."
)

write_board_manifest(board)
