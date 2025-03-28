library(SummarizedExperiment)
library(readr)
library(pins)

# Original FASTQs were downloaded from GEO and processed with the nf-core v3.12.0 pipeline:
# nextflow run nf-core/rnaseq -r 3.12.0 -profile singularity -c "$BAKER_REF"/nf_configs/rnaseq.config -w /scratch_space/jandrews/"$LSB_JOBNAME" \
# --outdir ./nfcore_mm10 --email jared.andrews@stjude.org --input nfcore_rnaseq.samplesheet.csv --gencode --genome MM10 --aligner star_salmon \
# --pseudo_aligner salmon --max_memory 128.GB --skip_stringtie --max_multiqc_email_size 15.MB -resume

description <- "
  **Study - SubStudy - Experiment - RNA-seq Data**

  This dataset contains O4+ immunopanned oligodendrocyte precursor cells (OPCs) from cortices of P5 or P6 Eed KO or control mice.
  See the [associated publication](https://www.science.org/doi/10.1126/sciadv.aaz6477) for more details.

  Please see the Pin generation code to view how this data was processed.
  "

# Load sample metadata.
meta <- read.csv("nfcore_rnaseq.samplesheet.csv", header = TRUE, stringsAsFactors = TRUE)

# Drop FASTQ file locations.
meta <- meta[, !colnames(meta) %in% c("fastq_1", "fastq_2")]

# Load counts. This object was generated using tximport via the nf-core RNA-seq pipeline on the salmon quants and
# is appropriate for pretty much all downstream DE packages (DESeq2, edgeR, limma).
cts <- read.table("salmon.merged.gene_counts_length_scaled.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# nf-core RNA-seq pipeline generates a counts table with the first two columns as gene IDs and gene symbols.
genes <- cts[, 1:2]
rownames(cts) <- cts[, 1]
cts <- cts[, -c(1:2)]

# Also carry along TPMs as an additional assay.
tpms <- read.table("salmon.merged.gene_tpm.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(tpms) <- tpms[, 1]
tpms <- tpms[, -c(1:2)]

# Set metadata rownames and ensure they match count column names.
rownames(meta) <- meta$sample
all(rownames(meta) == colnames(cts))

# Create a SummarizedExperiment object.
se <- SummarizedExperiment(
    assays = list(
        counts = as.matrix(cts),
        tpm = as.matrix(tpms),
        log2tpm = log2(as.matrix(tpms) + 1)
    ),
    colData = meta,
    rowData = genes
)

# Render this notebook and add it to the metadata of the pin along with experiment metadata.
metadata(se) <- list(description = description)

# Make board and write pin to it, which will ensure it's published to the pkgdown site.
board <- board_folder(here::here("pkgdown/assets/pins-board"))
pin_write(
    board = board,
    x = se,
    type = "rds",
    name = "GSE135880_SummarizedExperiment_mm10",
    title = "GSE135880_SummarizedExperiment_mm10",
    description = "A SummarizedExperiment containing O4+ immunopanned oligodendrocyte precursor cells (OPCs) from cortices of P5 or P6 Eed KO or control mice."
)

write_board_manifest(board)
