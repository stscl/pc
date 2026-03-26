# Optimal Parameter Search for Pattern Causality

Optimal Parameter Search for Pattern Causality

## Usage

``` r
# S4 method for class 'data.frame'
ops(
  data,
  target,
  source,
  E = 3:5,
  k = E,
  tau = 1,
  style = 1,
  lib = NULL,
  pred = NULL,
  maximize = c("positive", "negative", "dark"),
  dist.metric = c("euclidean", "manhattan", "maximum"),
  zero.tolerance = max(k),
  relative = TRUE,
  weighted = TRUE,
  threads = length(E),
  higher.parallel = TRUE,
  h = 0
)

# S4 method for class 'sf'
ops(
  data,
  target,
  source,
  E = 3:5,
  k = E + 1,
  tau = 1,
  style = 1,
  lib = NULL,
  pred = NULL,
  maximize = c("positive", "negative", "dark"),
  dist.metric = c("euclidean", "manhattan", "maximum"),
  zero.tolerance = max(k),
  relative = TRUE,
  weighted = TRUE,
  threads = length(E),
  higher.parallel = TRUE,
  nb = NULL
)

# S4 method for class 'SpatRaster'
ops(
  data,
  target,
  source,
  E = 3:5,
  k = E + 1,
  tau = 1,
  style = 1,
  lib = NULL,
  pred = NULL,
  maximize = c("positive", "negative", "dark"),
  dist.metric = c("euclidean", "manhattan", "maximum"),
  zero.tolerance = max(k),
  relative = TRUE,
  weighted = TRUE,
  threads = length(E),
  higher.parallel = TRUE
)
```

## Arguments

- data:

  Observation data.

- target:

  Integer of column indice for the target variable.

- source:

  Integer of column indice for the source variable.

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

- maximize:

  (optional) Causality metric to maximize: one of "positive",
  "negative", or "dark".

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

- h:

  (optional) Prediction horizon.

- nb:

  (optional) Neighbours list.

## Value

A list.

- `param`:

  The selected optimal parameter combination.

- `xmap`:

  A data.frame containing cross-mapping performance across parameter
  settings.

## Examples

``` r
columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
pc::ops(columbus, 1, 3, E = 2:10, maximize = "negative", threads = 1)
#> The suggested E, k, tau is 5, 9 and 1 
```
