
# pc

<!-- badges: start -->

<!-- [![CRAN](https://www.r-pkg.org/badges/version/pc)](https://CRAN.R-project.org/package=pc)
[![CRAN
Release](https://www.r-pkg.org/badges/last-release/pc)](https://CRAN.R-project.org/package=pc)
[![CRAN
Checks](https://badges.cranchecks.info/worst/pc.svg)](https://cran.r-project.org/web/checks/check_results_pc.html)
[![Downloads_all](https://badgen.net/cran/dt/pc?color=orange)](https://CRAN.R-project.org/package=pc)
[![Downloads_month](https://cranlogs.r-pkg.org/badges/pc)](https://CRAN.R-project.org/package=pc)
[![License](https://img.shields.io/badge/license-GPL--3-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)
[![Lifecycle:
experimental](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental) -->
[![R-CMD-check](https://github.com/stscl/pc/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/stscl/pc/actions/workflows/R-CMD-check.yaml)
[![R-universe](https://stscl.r-universe.dev/badges/pc?color=cyan)](https://stscl.r-universe.dev/pc)

<!-- badges: end -->

<a href="https://stscl.github.io/pc/"><img src="man/figures/pc.png" align="right" hspace="5" vspace="0" width="20%" alt="pc website: https://stscl.github.io/pc/"/></a>

<p align="right"; style="font-size:11px">logo by layeyo</p>

***P**attern **C**ausality Analysis*

*pc* is an R package for pattern-based causality analysis in both time
series and spatial cross-sectional data. It uses symbolic pattern
representations and cross mapping to detect directional interactions 
and infer causal structure from temporal dynamics and spatial snapshots. 
Built on a high-performance C++ backend with a lightweight R interface, *pc*
provides efficient and flexible tools for data-driven causality analysis.

> *Refer to the package documentation <https://stscl.github.io/pc/>
> for more detailed information.*

## Installation

- Install from [CRAN](https://CRAN.R-project.org/package=pc) with:

``` r
install.packages("pc", dependencies = TRUE)
```

- Install binary version from
  [R-universe](https://stscl.r-universe.dev/pc) with:

``` r
install.packages("pc",
                 repos = c("https://stscl.r-universe.dev",
                           "https://cloud.r-project.org"),
                 dependencies = TRUE)
```

- Install from source code on [GitHub](https://github.com/stscl/pc)
  with:

``` r
if (!requireNamespace("pak")) {
    install.packages("pak")
}
pak::pak("stscl/pc", dependencies = TRUE)
```
