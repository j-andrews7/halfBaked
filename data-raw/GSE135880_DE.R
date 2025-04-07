library(halfBaked)
library(pins)

set.seed(7)

board <- board_url("https://j-andrews7.github.io/halfBaked/pins-board/")
se <- pin_read(board, "GSE135880_SummarizedExperiment_mm10")

# DEseq2 requires integer counts, so round them.
assay(se, "counts") <- round(assay(se, "counts"))

# Design is just being set to no design here (~ 1), since we'll set it later.
dds <- DESeqDataSet(se, design = ~1)

# Remove low counts genes, adjust as necessary for smallest group.
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep, ]
dds <- DESeq(dds)

# Get normalized counts and variance stabilized counts (preferred for viz).
vsd <- vst(dds)
assay(dds, "vst") <- as.matrix(assay(vsd))
lognorm <- normTransform(dds)
assay(dds, "lognorm") <- as.matrix(assay(lognorm))

res <- list()

# Whatever can be used as name of contrast vectors, though typically they should be a group vs another group.
# I use -b_OTHERVARIABLES in the name to indicate additional variables that were accounted for in the model.
# That is, fed to the `block` parameter of `get_DESEQ2_res`, if provided.
contrasts <- list(
    "Eed_cKO.v.Control" = c("Group", "Eed_cKO", "Control")
)

res <- get_DESeq2_res(dds, res.list = res, contrasts = contrasts, add.rowData = c("ENSEMBL", "SYMBOL"))

# Add to metadata of the SummarizedExperiment object.
se.meta <- metadata(se)
se.meta$DESeq2.Results <- res
metadata(se) <- se.meta

# Write to pin.
board <- board_folder(here::here("pkgdown/assets/pins-board"))
pin_write(
    board = board,
    x = se,
    type = "rds",
    name = "GSE135880_SummarizedExperiment_mm10",
    title = "GSE135880_SummarizedExperiment_mm10",
    description = "A SummarizedExperiment containing counts for O4+ immunopanned oligodendrocyte precursor cells (OPCs) from cortices of P5 or P6 Eed KO or control mice, along with DESeq2 analysis results in the metadata."
)

write_board_manifest(board)
