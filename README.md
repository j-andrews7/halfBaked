
<!-- README.md is generated from README.Rmd. Please edit that file -->

# halfBaked

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/j-andrews7/halfBaked)](https://github.com/j-andrews7/halfBaked/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/j-andrews7/halfBaked)](https://github.com/j-andrews7/halfBaked/pulls)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Bioc release
status](http://www.bioconductor.org/shields/build/release/bioc/halfBaked.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/halfBaked)
[![Bioc devel
status](http://www.bioconductor.org/shields/build/devel/bioc/halfBaked.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/halfBaked)
[![Bioc downloads
rank](https://bioconductor.org/shields/downloads/release/halfBaked.svg)](http://bioconductor.org/packages/stats/bioc/halfBaked/)
[![Bioc
support](https://bioconductor.org/shields/posts/halfBaked.svg)](https://support.bioconductor.org/tag/halfBaked)
[![Bioc
history](https://bioconductor.org/shields/years-in-bioc/halfBaked.svg)](https://bioconductor.org/packages/release/bioc/html/halfBaked.html#since)
[![Bioc last
commit](https://bioconductor.org/shields/lastcommit/devel/bioc/halfBaked.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/halfBaked/)
[![Bioc
dependencies](https://bioconductor.org/shields/dependencies/release/halfBaked.svg)](https://bioconductor.org/packages/release/bioc/html/halfBaked.html#since)
[![check-bioc](https://github.com/j-andrews7/halfBaked/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/j-andrews7/halfBaked/actions/workflows/check-bioc.yml)
[![Codecov test
coverage](https://codecov.io/gh/j-andrews7/halfBaked/graph/badge.svg)](https://app.codecov.io/gh/j-andrews7/halfBaked)
<!-- badges: end -->

The goal of `halfBaked` is to …

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `halfBaked` from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("halfBaked")
```

And the development version from
[GitHub](https://github.com/j-andrews7/halfBaked) with:

``` r
BiocManager::install("j-andrews7/halfBaked")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library("halfBaked")
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub!

## Citation

Below is the citation output from using `citation('halfBaked')` in R.
Please run this yourself to check for any updates on how to cite
**halfBaked**.

``` r
print(citation("halfBaked"), bibtex = TRUE)
#> To cite package 'halfBaked' in publications use:
#> 
#>   j-andrews7 (2025). _halfBaked_. doi:10.18129/B9.bioc.halfBaked
#>   <https://doi.org/10.18129/B9.bioc.halfBaked>,
#>   https://github.com/j-andrews7/halfBaked/halfBaked - R package version
#>   0.0.0.9000, <http://www.bioconductor.org/packages/halfBaked>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {halfBaked},
#>     author = {{j-andrews7}},
#>     year = {2025},
#>     url = {http://www.bioconductor.org/packages/halfBaked},
#>     note = {https://github.com/j-andrews7/halfBaked/halfBaked - R package version 0.0.0.9000},
#>     doi = {10.18129/B9.bioc.halfBaked},
#>   }
```

Please note that the `halfBaked` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `halfBaked` project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

- Continuous code testing is possible thanks to [GitHub
  actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
  through *[usethis](https://CRAN.R-project.org/package=usethis)*,
  *[remotes](https://CRAN.R-project.org/package=remotes)*, and
  *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized
  to use [Bioconductor’s docker
  containers](https://www.bioconductor.org/help/docker/) and
  *[BiocCheck](https://bioconductor.org/packages/3.20/BiocCheck)*.
- Code coverage assessment is possible thanks to
  [codecov](https://codecov.io/gh) and
  *[covr](https://CRAN.R-project.org/package=covr)*.
- The [documentation website](http://j-andrews7.github.io/halfBaked) is
  automatically updated thanks to
  *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
- The code is styled automatically thanks to
  *[styler](https://CRAN.R-project.org/package=styler)*.
- The documentation is formatted thanks to
  *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.20/biocthis)*.
