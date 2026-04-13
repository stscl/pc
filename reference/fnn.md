# False Nearest Neighbors

False Nearest Neighbors

## Usage

``` r
# S4 method for class 'data.frame'
fnn(
  data,
  target,
  E = 10,
  k = 1,
  tau = 1,
  style = 1,
  lib = NULL,
  pred = NULL,
  dist.metric = c("euclidean", "manhattan", "maximum"),
  rt = 10,
  eps = NULL,
  threads = length(E),
  higher.parallel = TRUE,
  ...
)

# S4 method for class 'sf'
fnn(
  data,
  target,
  E = 10,
  k = 1,
  tau = 1,
  style = 1,
  lib = NULL,
  pred = NULL,
  dist.metric = c("euclidean", "manhattan", "maximum"),
  rt = 10,
  eps = NULL,
  threads = length(E),
  higher.parallel = TRUE,
  detrend = FALSE,
  nb = NULL,
  ...
)

# S4 method for class 'SpatRaster'
fnn(
  data,
  target,
  E = 10,
  k = 1,
  tau = 1,
  style = 1,
  lib = NULL,
  pred = NULL,
  dist.metric = c("euclidean", "manhattan", "maximum"),
  rt = 10,
  eps = NULL,
  threads = length(E),
  higher.parallel = TRUE,
  detrend = FALSE,
  ...
)
```

## Arguments

- data:

  Observation data.

- target:

  Integer of column indice for the target variable.

- E:

  (optional) Embedding dimensions.

- k:

  (optional) Number of nearest neighbors used for evaluation.

- tau:

  (optional) Step of lag.

- style:

  (optional) Embedding style (`0` includes current state, `1` excludes
  it).

- lib:

  (optional) Libraries indices.

- pred:

  (optional) Predictions indices.

- dist.metric:

  (optional) Distance measure to be used.

- rt:

  (optional) Relative distance threshold.

- eps:

  (optional) Absolute distance threshold.

- threads:

  (optional) Number of threads used.

- higher.parallel:

  (optional) Whether to use a higher level of parallelism.

- ...:

  Additional arguments to absorb unused inputs in method dispatch.

- detrend:

  (optional) Whether to remove the linear trend.

- nb:

  (optional) Neighbours list.

## Value

A vector.

## References

Kennel, M.B., Brown, R., Abarbanel, H.D.I., 1992. Determining embedding
dimension for phase-space reconstruction using a geometrical
construction. Physical Review A 45, 3403–3411.

## Examples

``` r
crash = sf::read_sf(system.file("case/crash.gpkg", package = "pc"))
pc::fnn(crash, 1, threads = 1)
#> [fnn] Input E values exceeding max embeddable dimension were truncated, and values < 2 were clamped to 2.
#> [fnn] Max embedding dimension E_max is auto-computed, with results returned for dimensions 1 through E_max.
#> [fnn] Output 'E:i' (where i = 1 to E_max-1) corresponds to the comparison between dimension i and i+1.
#>         E:1         E:2         E:3         E:4         E:5         E:6 
#> 0.622093023 0.209302326 0.029069767 0.005813953 0.023255814 0.005813953 
#>         E:7         E:8         E:9 
#> 0.000000000 0.005813953 0.000000000 
```
