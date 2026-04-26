.ops_ts = \(data, target, source, E = 3, k = E+2, tau = 1, style = 1, lib = NULL, pred = NULL,
            maximize = c("positive", "negative", "dark"), dist.metric = c("euclidean", "manhattan", "maximum"), 
            zero.tolerance = max(k), relative = TRUE, weighted = TRUE, threads = length(E), lower.parallel = TRUE, h = 0) {
  maximize = match.arg(maximize)
  dist.metric = match.arg(dist.metric)
  dlist = .validate_var(data, target, source)
  tv = dlist[[1]]; sv = dlist[[2]]
  if (is.null(lib)) lib = which(!(is.na(tv) | is.na(sv)))
  if (is.null(pred)) pred = lib

  return(RcppPCops(tv, sv, lib, pred, E, tau, k, maximize, style, zero.tolerance,
                   dist.metric, relative, weighted, threads, lower.parallel, h, NULL, NULL))
}

.ops_lattice = \(data, target, source, E = 3, k = E+2, tau = 1, style = 1, lib = NULL, pred = NULL, 
                 maximize = c("positive", "negative", "dark"), dist.metric = c("euclidean", "manhattan", "maximum"), 
                 zero.tolerance = max(k), relative = TRUE, weighted = TRUE, threads = length(E), lower.parallel = TRUE, nb = NULL) {
  if (is.null(nb)) nb = sdsfun::spdep_nb(data)
  maximize = match.arg(maximize)
  dist.metric = match.arg(dist.metric)
  dlist = .validate_var(data, target, source)
  tv = dlist[[1]]; sv = dlist[[2]]
  if (is.null(lib)) lib = which(!(is.na(tv) | is.na(sv)))
  if (is.null(pred)) pred = lib

  return(RcppPCops(tv, sv, lib, pred, E, tau, k, maximize, style, zero.tolerance,
                   dist.metric, relative, weighted, threads, lower.parallel, h, nb, NULL))
}

.ops_grid = \(data, target, source, E = 3, k = E+2, tau = 1, style = 1, lib = NULL, pred = NULL,
              maximize = c("positive", "negative", "dark"), dist.metric = c("euclidean", "manhattan", "maximum"), 
              zero.tolerance = max(k), relative = TRUE, weighted = TRUE, threads = length(E), lower.parallel = TRUE) {
  maximize = match.arg(maximize)
  dist.metric = match.arg(dist.metric)
  dlist = .validate_var(data, target, source)
  tv = dlist[[1]]; sv = dlist[[2]]
  if (is.null(lib)) lib = which(!(is.na(tv) | is.na(sv)))
  if (is.null(pred)) pred = lib

  return(RcppPCops(tv, sv, lib, pred, E, tau, k, maximize, style, zero.tolerance,
                   dist.metric, relative, weighted, threads, lower.parallel, h, NULL, terra::nrow(data)))
}

#' Pattern Causality
#'
#' @note `pc` only supports numeric input data.
#'
#' @param data Observation data.
#' @param target Integer of column indice for the target variable.
#' @param source Integer of column indice for the source variable.
#' @param libsizes (optional) Number of observations used.
#' @param E (optional) Embedding dimensions.
#' @param k (optional) Number of nearest neighbors used for projection.
#' @param tau (optional) Step of lag.
#' @param style (optional) Embedding style (`0` includes current state, `1` excludes it).
#' @param lib (optional) Libraries indices.
#' @param pred (optional) Predictions indices.
#' @param boot (optional) Number of bootstraps to perform.
#' @param random (optional) Whether to use random sampling.
#' @param seed (optional) Random seed.
#' @param dist.metric (optional) Distance measure to be used.
#' @param zero.tolerance (optional) Maximum number of zeros tolerated in signature space.
#' @param relative (optional) Whether to calculate relative changes in embedding.
#' @param weighted (optional) Whether to weight causal strength.
#' @param threads (optional) Number of threads used.
#' @param lower.parallel (optional) Whether to use lower level of parallelism.
#' @param verbose (optional) Whether to show the progress bar.
#' @param nb (optional) Neighbours list.
#' @param h (optional) Prediction horizon.
#'
#' @return A list (if `libsizes` is `NULL`) or dataframe.
#'
#' @export
#' @name pc
#' @aliases pc,data.frame-method
#' @references
#' Stavroglou, S.K., Pantelous, A.A., Stanley, H.E., Zuev, K.M., 2020. Unveiling causal interactions in complex systems. Proceedings of the National Academy of Sciences 117, 7599–7605.
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' pc::pc(columbus, 1, 3, E = 6, k = 8)
#'
methods::setMethod("pc", "data.frame", .pc_ts)

#' @rdname pc
methods::setMethod("pc", "sf", .pc_lattice)

#' @rdname pc
methods::setMethod("pc", "SpatRaster", .pc_grid)
