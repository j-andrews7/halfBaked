
<!-- README.md is generated from README.Rmd. Please edit that file -->
<p align="center">
<img src="man/figures/halfBaked_hex.png" alt="halfBaked" width="330">
</p>

------------------------------------------------------------------------

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/j-andrews7/halfBaked)](https://github.com/j-andrews7/halfBaked/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/j-andrews7/halfBaked)](https://github.com/j-andrews7/halfBaked/pulls)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![check-bioc](https://github.com/j-andrews7/halfBaked/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/j-andrews7/halfBaked/actions/workflows/check-bioc.yml)
[![Codecov test
coverage](https://codecov.io/gh/j-andrews7/halfBaked/graph/badge.svg)](https://app.codecov.io/gh/j-andrews7/halfBaked)
<!-- badges: end -->

**halfBaked** is a collection of convenience functions to aid in \`omics
analyses. They span functions to automate and organize repetitive
analyses (e.g. many pairwise RNA-seq comparisons) to viz.

More specifically, they’re designed to help create template Rmd
notebooks to jumpstart exploratory data analysis and get the ball
rolling for common data modalities.

These templates are provided as vignettes and can be accessed via the
`browseVignettes("halfBaked")` function.

## Installation

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `halfBaked` from
Github:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("j-andrews7/halfBaked")
```

## Design Principles

**halfBaked** is intended to help save *experienced* bioinformaticians
some time with relatively run of the mill analyses. It is not intended
as an out of the box solution for bioinformatics newbies.

It is not recommended for use without a solid understanding of the
underlying methods and assumptions.

## Contributing

If you have a function that’d fit in here, feel free to open a PR.

If you have issues, your problems are your own.
