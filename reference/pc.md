# Pattern Causality

Pattern Causality

## Usage

``` r
# S4 method for class 'data.frame'
pc(
  data,
  target,
  source,
  libsizes = NULL,
  E = 3,
  k = E,
  tau = 1,
  style = 1,
  lib = NULL,
  pred = NULL,
  boot = 99,
  replace = FALSE,
  seed = 42L,
  dist.metric = c("euclidean", "manhattan", "maximum"),
  zero.tolerance = max(k),
  relative = TRUE,
  weighted = TRUE,
  threads = length(libsizes),
  higher.parallel = TRUE,
  verbose = TRUE,
  h = 0,
  ...
)

# S4 method for class 'sf'
pc(
  data,
  target,
  source,
  libsizes = NULL,
  E = 3,
  k = E + 1,
  tau = 1,
  style = 1,
  lib = NULL,
  pred = NULL,
  boot = 99,
  replace = FALSE,
  seed = 42L,
  dist.metric = c("euclidean", "manhattan", "maximum"),
  zero.tolerance = max(k),
  relative = TRUE,
  weighted = TRUE,
  threads = length(libsizes),
  higher.parallel = TRUE,
  verbose = TRUE,
  detrend = FALSE,
  nb = NULL,
  ...
)

# S4 method for class 'SpatRaster'
pc(
  data,
  target,
  source,
  libsizes = NULL,
  E = 3,
  k = E + 1,
  tau = 1,
  style = 1,
  lib = NULL,
  pred = NULL,
  boot = 99,
  replace = FALSE,
  seed = 42L,
  dist.metric = c("euclidean", "manhattan", "maximum"),
  zero.tolerance = max(k),
  relative = TRUE,
  weighted = TRUE,
  threads = length(libsizes),
  higher.parallel = TRUE,
  verbose = TRUE,
  detrend = FALSE,
  ...
)
```

## Arguments

- data:

  Observation data.

- target:

  Integer of column indice for the target variable.

- source:

  Integer of column indice for the source variable.

- libsizes:

  (optional) Number of observations used.

- E:

  (optional) Embedding dimensions.

- k:

  (optional) Number of nearest neighbors used for projection.

- tau:

  (optional) Step of lag.

- style:

  (optional) Embedding style (`0` includes current state, `1` excludes
  it).

- lib:

  (optional) Libraries indices.

- pred:

  (optional) Predictions indices.

- boot:

  (optional) Number of bootstraps to perform.

- replace:

  (optional) Should sampling be with replacement?

- seed:

  (optional) Random seed.

- dist.metric:

  (optional) Distance measure to be used.

- zero.tolerance:

  (optional) Maximum number of zeros tolerated in signature space.

- relative:

  (optional) Whether to calculate relative changes in embedding.

- weighted:

  (optional) Whether to weight causal strength.

- threads:

  (optional) Number of threads used.

- higher.parallel:

  (optional) Whether to use a higher level of parallelism.

- verbose:

  (optional) Whether to show the progress bar.

- h:

  (optional) Prediction horizon.

- ...:

  Additional arguments to absorb unused inputs in method dispatch.

- detrend:

  (optional) Whether to remove the linear trend.

- nb:

  (optional) Neighbours list.

## Value

A list.

- causality:

  A data.frame of causality results. When `libsizes` is `NULL`, it
  contains per-sample causality estimates; otherwise, it contains
  causality results evaluated across different library sizes.

- summary:

  A data.frame summarizing overall causality metrics. Only returned when
  `libsizes` is `NULL`.

## References

Stavroglou, S.K., Pantelous, A.A., Stanley, H.E., Zuev, K.M., 2020.
Unveiling causal interactions in complex systems. Proceedings of the
National Academy of Sciences 117, 7599–7605.

## Examples

``` r
crash = sf::read_sf(system.file("case/crash.gpkg", package = "pc"))
p1 = pc::pc(crash, 1, 2, E = 3, k = 7, threads = 1)
print(p1)
#>       type  strength
#> 1 positive 0.4951110
#> 2 negative 0.1249999
#> 3     dark 0.2307989
plot(p1)


# convergence diagnostics
p2 = pc::pc(crash, 1, 2, libsizes = seq(10,172,40), E = 3, k = 7, threads = 1)
#> Computing: [========================================] 100% (done)                         
print(p2)
#>    libsizes     type      mean        q05       q50       q95
#> 1        10 positive 0.3324836 0.12787433 0.3045097 0.5092801
#> 2        50 positive 0.3170059 0.19212321 0.3117322 0.4300975
#> 3        90 positive 0.3264837 0.22973601 0.3195043 0.4254085
#> 4       130 positive 0.3596636 0.25988689 0.3534218 0.4830695
#> 5       170 positive 0.4898578 0.44093561 0.4951110 0.5069378
#> 6        10 negative 0.1678948 0.06480651 0.1523444 0.2571607
#> 7        50 negative 0.1692400 0.03550959 0.1428533 0.3499086
#> 8        90 negative 0.1234257 0.02974710 0.1024092 0.3079307
#> 9       130 negative 0.1189617 0.03088201 0.1107992 0.2670766
#> 10      170 negative 0.1229800 0.10357129 0.1249999 0.1354733
#> 11       10     dark 0.2323689 0.10024130 0.2221773 0.3333159
#> 12       50     dark 0.2356663 0.14869228 0.2280345 0.3256868
#> 13       90     dark 0.2173155 0.14595310 0.2098594 0.3324160
#> 14      130     dark 0.2109280 0.15301645 0.2075058 0.2744680
#> 15      170     dark 0.2319855 0.22190066 0.2307989 0.2436455
plot(p2)

```
