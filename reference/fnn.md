# False Nearest Neighbors

False Nearest Neighbors

## Usage

``` r
# S4 method for class 'data.frame'
fnn(
  data,
  target,
  E = 2:10,
  tau = 1,
  style = 1,
  lib = NULL,
  pred = NULL,
  dist.metric = c("euclidean", "manhattan", "maximum"),
  rt = 10,
  eps = 2,
  threads = length(E),
  higher.parallel = TRUE,
  ...
)

# S4 method for class 'sf'
fnn(
  data,
  target,
  E = 3:10,
  tau = 1,
  style = 1,
  lib = NULL,
  pred = NULL,
  dist.metric = c("euclidean", "manhattan", "maximum"),
  rt = 10,
  eps = 2,
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
  E = 3:10,
  tau = 1,
  style = 1,
  lib = NULL,
  pred = NULL,
  dist.metric = c("euclidean", "manhattan", "maximum"),
  rt = 10,
  eps = 2,
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
columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
pc::fnn(columbus, 1, threads = 1)
#> [fnn] Output 'E:i' corresponds to the i-th valid embedding dimension.
#> [fnn] Input E values exceeding max embeddable dimension were truncated.
#> [fnn] Please map output indices to original E inputs before interpretation.
#>       E:1       E:2       E:3       E:4       E:5       E:6       E:7       E:8 
#> 0.9387755 0.7755102 0.7142857 0.5714286 0.3061224 0.1224490 0.0000000 0.0000000 
```
