# pc <a href="https://stscl.github.io/pc/"><img src="man/figures/pc.png" align="right" hspace="5" vspace="0" width="15%" alt="pc website: https://stscl.github.io/pc/"/></a>

<p align="right"; style="font-size:11px">logo by layeyo</p>

<!-- badges: start -->

[![CRAN](https://www.r-pkg.org/badges/version/pc)](https://CRAN.R-project.org/package=pc)
[![CRAN
Release](https://www.r-pkg.org/badges/last-release/pc)](https://CRAN.R-project.org/package=pc)
[![CRAN
Checks](https://badges.cranchecks.info/worst/pc.svg)](https://cran.r-project.org/web/checks/check_results_pc.html)
[![Downloads_all](https://badgen.net/cran/dt/pc?color=orange)](https://CRAN.R-project.org/package=pc)
[![Downloads_month](https://cranlogs.r-pkg.org/badges/pc)](https://CRAN.R-project.org/package=pc)
[![License](https://img.shields.io/badge/license-GPL--3-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)
[![Lifecycle:
experimental](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/stscl/pc/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/stscl/pc/actions/workflows/R-CMD-check.yaml)
[![R-universe](https://stscl.r-universe.dev/badges/pc?color=cyan)](https://stscl.r-universe.dev/pc)

<!-- badges: end -->

***P**attern **C**ausality Analysis*

*pc* is an R package for pattern-based causality analysis in both time series and spatial cross-sectional data. It uses symbolic pattern representations and cross mapping to detect directional interactions and infer causal structure from temporal dynamics and spatial snapshots. Built on a high-performance C++ backend with a lightweight R interface, *pc* provides efficient and flexible tools for data-driven causality analysis.

> *Refer to the package documentation <https://stscl.github.io/pc/> for more detailed information.*

## Installation

- Install from [CRAN](https://CRAN.R-project.org/package=pc) with:

``` r
install.packages("pc", dependencies = TRUE)
```

- Install binary version from [R-universe](https://stscl.r-universe.dev/pc) with:

``` r
install.packages("pc",
                 repos = c("https://stscl.r-universe.dev",
                           "https://cloud.r-project.org"),
                 dependencies = TRUE)
```

- Install from source code on [GitHub](https://github.com/stscl/pc) with:

``` r
if (!requireNamespace("pak")) {
    install.packages("pak")
}
pak::pak("stscl/pc", dependencies = TRUE)
```

## References

Sugihara, G., May, R., Ye, H., Hsieh, C., Deyle, E., Fogarty, M., Munch, S., 2012. Detecting Causality in Complex Ecosystems. Science 338, 496–500. https://doi.org/10.1126/science.1227079.

Stavroglou, S.K., Pantelous, A.A., Stanley, H.E., Zuev, K.M., 2019. Hidden interactions in financial markets. Proceedings of the National Academy of Sciences 116, 10646–10651. https://doi.org/10.1073/pnas.1819449116.

Stavroglou, S.K., Pantelous, A.A., Stanley, H.E., Zuev, K.M., 2020. Unveiling causal interactions in complex systems. Proceedings of the National Academy of Sciences 117, 7599–7605. https://doi.org/10.1073/pnas.1918269117.

Zhang, Z., Wang, J., 2025. A model to identify causality for geographic patterns. International Journal of Geographical Information Science 1–21. https://doi.org/10.1080/13658816.2025.2581207.

Lyu, W., Lei, Y., Yi, W., Song, Y., Li, X., Dai, S., Qin, Y., Zhao, W., 2026. Causal discovery in urban data with temporal empirical dynamic modeling: The R package tEDM. Computers, Environment and Urban Systems 127, 102435. https://doi.org/10.1016/j.compenvurbsys.2026.102435.
