# Pattern Causality

Pattern Causality

## Usage

``` r
# S4 method for class 'data.frame'
pc(
  data,
  source,
  target,
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
  source,
  target,
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
  source,
  target,
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

- source:

  Integer of column indice for the source variable.

- target:

  Integer of column indice for the target variable.

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
#>       type    strength
#> 1 positive 0.069745223
#> 2 negative 0.009307607
#> 3     dark 0.027701849
plot(p1)


# convergence diagnostics
p2 = pc::pc(crash, 1, 2, libsizes = seq(10,172,40), E = 3, k = 7, threads = 1)
#> Computing: [========================================] 100% (done)                         
print(p2)
#>    libsizes     type       mean         q05         q50        q95
#> 1        10 positive 0.03054525 0.008255862 0.026301738 0.07116232
#> 2        50 positive 0.04596615 0.021071881 0.041629869 0.09557344
#> 3        90 positive 0.05589145 0.028179442 0.054469716 0.08439766
#> 4       130 positive 0.06644474 0.045802551 0.065094118 0.08912378
#> 5       170 positive 0.06935769 0.066280590 0.069745223 0.07238961
#> 6        10 negative 0.02091447 0.004121190 0.017365839 0.04709580
#> 7        50 negative 0.02137358 0.003317208 0.013766347 0.04742605
#> 8        90 negative 0.01105271 0.001130310 0.008844640 0.02491195
#> 9       130 negative 0.01189137 0.003675142 0.009806801 0.02513823
#> 10      170 negative 0.00956734 0.008438496 0.009307607 0.01136995
#> 11       10     dark 0.02395373 0.007745936 0.021878864 0.04582000
#> 12       50     dark 0.02562744 0.013133235 0.023740861 0.04808648
#> 13       90     dark 0.03763731 0.019117538 0.036076026 0.06232317
#> 14      130     dark 0.03337542 0.019338071 0.028802572 0.05861903
#> 15      170     dark 0.02813939 0.025250845 0.027669575 0.02974462
plot(p2)

```
