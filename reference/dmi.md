# Delayed Mutual Information

Delayed Mutual Information

## Usage

``` r
# S4 method for class 'data.frame'
dmi(
  data,
  target,
  tau = 1:10,
  pred = NULL,
  k = 3,
  base = exp(1),
  normalize = FALSE,
  threads = length(tau),
  ...
)
```

## Arguments

- data:

  Observation data.

- target:

  Integer of column indice for the target variable.

- tau:

  (optional) Step of lag.

- pred:

  (optional) Predictions indices.

- k:

  (optional) Number of nearest neighbors used for evaluation.

- base:

  (optional) Logarithm base of the entropy.

- normalize:

  (optional) Whether to normalize MI values.

- threads:

  (optional) Number of threads used.

- ...:

  Additional arguments to absorb unused inputs in method dispatch.

## Value

A vector.

## References

Fraser, A.M., Swinney, H.L., 1986. Independent coordinates for strange
attractors from mutual information. Physical Review A 33, 1134–1140.

Kraskov, A., Stogbauer, H., Grassberger, P., 2004. Estimating mutual
information. Physical Review E 69, 066138.

## Examples

``` r
abun = readr::read_csv(system.file("case/abundance.csv", package = "pc"))
#> Rows: 71 Columns: 3
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> dbl (3): time, paramecium, didinium
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
pc::dmi(abun, 2)
#>      tau:1      tau:2      tau:3      tau:4      tau:5      tau:6      tau:7 
#> 0.50563138 0.09094202 0.01664535 0.11659768 0.30769606 0.26047393 0.05159755 
#>      tau:8      tau:9     tau:10 
#> 0.00000000 0.08849267 0.07227334 
```
