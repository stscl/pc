# pc

[![pc website:
https://stscl.github.io/pc/](reference/figures/pc.png)](https://stscl.github.io/pc/)

***P**attern **C**ausality Analysis*

*pc* is an R package for pattern-based causality analysis in both time
series and spatial cross-sectional data. It uses symbolic pattern
representations and cross mapping to detect directional interactions and
infer causal structure from temporal dynamics and spatial snapshots.
Built on a high-performance C++ backend with a lightweight R interface,
*pc* provides efficient and flexible tools for data-driven causality
analysis.

> *Refer to the package documentation <https://stscl.github.io/pc/> for
> more detailed information.*

## Installation

- Install from [CRAN](https://CRAN.R-project.org/package=pc) with:

``` r
install.packages("pc", dep = TRUE)
```

- Install binary version from
  [R-universe](https://stscl.r-universe.dev/pc) with:

``` r
install.packages("pc",
                 repos = c("https://stscl.r-universe.dev",
                           "https://cloud.r-project.org"),
                 dep = TRUE)
```

- Install from source code on [GitHub](https://github.com/stscl/pc)
  with:

``` r
if (!requireNamespace("devtools")) {
    install.packages("devtools")
}
devtools::install_github("stscl/pc",
                         build_vignettes = TRUE,
                         dep = TRUE)
```
