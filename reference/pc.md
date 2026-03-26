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
  random = TRUE,
  seed = 42L,
  dist.metric = c("euclidean", "manhattan", "maximum"),
  zero.tolerance = max(k),
  relative = TRUE,
  weighted = TRUE,
  threads = length(libsizes),
  higher.parallel = TRUE,
  verbose = TRUE,
  h = 0
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
  random = TRUE,
  seed = 42L,
  dist.metric = c("euclidean", "manhattan", "maximum"),
  zero.tolerance = max(k),
  relative = TRUE,
  weighted = TRUE,
  threads = length(libsizes),
  higher.parallel = TRUE,
  verbose = TRUE,
  nb = NULL
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
  random = TRUE,
  seed = 42L,
  dist.metric = c("euclidean", "manhattan", "maximum"),
  zero.tolerance = max(k),
  relative = TRUE,
  weighted = TRUE,
  threads = length(libsizes),
  higher.parallel = TRUE,
  verbose = TRUE
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

- random:

  (optional) Whether to use random sampling.

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

- nb:

  (optional) Neighbours list.

## Value

A list (when `libsizes` is `NULL`) or data.frame. If `libsizes` is
`NULL`, a list with two components is returned:

- causality:

  A data.frame containing per-sample causality results.

- summary:

  A data.frame summarizing overall causality metrics.

## Note

`pc` only supports numeric input data.

## References

Stavroglou, S.K., Pantelous, A.A., Stanley, H.E., Zuev, K.M., 2020.
Unveiling causal interactions in complex systems. Proceedings of the
National Academy of Sciences 117, 7599–7605.

## Examples

``` r
columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
pc::pc(columbus, 1, 3, E = 5, k = 9, threads = 1)
#>       type  strength
#> 1 positive 0.0000000
#> 2 negative 0.7753134
#> 3     dark 0.6213416
```
