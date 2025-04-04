% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differential_expression.R
\name{get_DESeq2_res}
\alias{get_DESeq2_res}
\title{Get DESeq2 Results}
\usage{
get_DESeq2_res(
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
  BPPARAM = NULL
)
}
\arguments{
\item{dds}{A \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment::SummarizedExperiment}  or \link[DESeq2:DESeqDataSet]{DESeq2::DESeqDataSet} object.}

\item{contrasts}{A named list of contrasts, e.g. \code{list("condition" = c("condition", "A", "B"))}.
The first element is the variable of interest, the second is the test, and the third is the reference level.
The name of  each element in the list will be used as a name in the results table.}

\item{res.list}{A named list to hold DESeq2 result tables.
Allows the function to be run multiple times if needed and append to the same list.
Defaults to an empty list.}

\item{user.mat}{A logical indicating whether a user-specified model matrix is provided.
Defaults to FALSE.}

\item{block}{A vector of additional terms to be considered in the model, beyond the main effect.
Defaults to NULL.}

\item{design}{The design formula or matrix.
If a matrix is provided, ensure 'user.mat' is set to TRUE.
Defaults to NULL.}

\item{alpha}{The significance level for hypothesis testing.
Defaults to 0.05.}

\item{lfc.th}{A numeric vector of log2 fold-change thresholds.
Defaults to `c(log2(1.25), log2(1.5))``.}

\item{shrink.method}{The method used for shrinkage estimation.
Must be one of "apeglm", "ashr", or NULL.
Defaults to "ashr".}

\item{norm.ercc}{A logical indicating whether to normalize to ERCC spike-ins.}

\item{add.rowData}{A vector of column names from the rowData slot of the DESeqDataSet
to be added to the results table.
Defaults to NULL.}

\item{BPPARAM}{The BiocParallelParam object specifying the parallel back-end to be used.
Defaults to NULL.}
}
\value{
A named list of \link[DESeq2:DESeqResults]{DESeq2::DESeqResults} objects for the specified contrasts.
If \code{add.rowData} is supplied, these will be returned as \link[S4Vectors:DFrame]{S4Vectors::DFrame} objects instead.
}
\description{
This function obtains a set of comparisons from a DESeq2 analysis, given a named list of contrasts.
It allows additional model
parameters to be specified and a design matrix to be manually adjusted.
}
\details{
It is important to note that LFC shrinkage is independent of the typical MLE results table.
That is, if \code{lfc.th} is provided, the results table p-values will reflect that testing threshold.
If \code{shrink.method} is set to \code{ashr} with \code{lfc.th}, s-values will be returned in addition to the MLE p-values.

It is possible to have shrunken FCs that are near 0 that still have a significant p-value.
This is frustrating, as it means you end up having to do post-hoc filtering on LFC
or filter with s-values, which are more difficult to interpret than the adjusted p-values.
}
\examples{
library(DESeq2)
dds_de <- makeExampleDESeqDataSet(n = 100, m = 12, betaSD = 2) # DE genes
rownames(dds_de) <- paste0("gene", 1001:1100)
dds <- makeExampleDESeqDataSet(n = 1000, m = 12) # Non-DE genes
dds <- rbind(dds, dds_de)
dds <- DESeq(dds)
contrasts <- list("condition" = c("condition", "A", "B"))
res <- get_DESeq2_res(dds, contrasts)

names(res)
head(res[[1]])

}
\author{
Jared Andrews
}
