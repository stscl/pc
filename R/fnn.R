.fnn_ts = \(data, target, E = 2:10, tau = 1, style = 1, lib = NULL, pred = NULL,
            dist.metric = c("euclidean", "manhattan", "maximum"), rt = 10, eps = 2,
            threads = length(E), higher.parallel = TRUE, ...) {
  dist.metric = match.arg(dist.metric)
  tv = .validate_var(data, target)[[1]]
  if (is.null(lib)) lib = which(!is.na(tv))
  if (is.null(pred)) pred = lib

  return(RcppFNN(tv, rt, eps, lib, pred, E, tau, style, dist.metric,
                 threads, higher.parallel, NULL, NULL))
}

.fnn_lattice = \(data, target, E = 3:10, tau = 1, style = 1, lib = NULL, pred = NULL, 
                 dist.metric = c("euclidean", "manhattan", "maximum"), rt = 10, eps = 2,
                 threads = length(E), higher.parallel = TRUE, detrend = FALSE, nb = NULL, ...) {
  if (is.null(nb)) nb = sdsfun::spdep_nb(data)
  dist.metric = match.arg(dist.metric)
  tv = .validate_var(data, target, detrend = detrend)[[1]]
  if (is.null(lib)) lib = which(!is.na(tv))
  if (is.null(pred)) pred = lib

  return(RcppFNN(tv, rt, eps, lib, pred, E, tau, style, dist.metric,
                 threads, higher.parallel, nb, NULL))
}

.fnn_grid = \(data, target, E = 3:10, tau = 1, style = 1, lib = NULL, pred = NULL, 
              dist.metric = c("euclidean", "manhattan", "maximum"), rt = 10, eps = 2, 
              threads = length(E), higher.parallel = TRUE, detrend = FALSE, ...) {
  dist.metric = match.arg(dist.metric)
  tv = .validate_var(data, target, detrend = detrend)[[1]]
  if (is.null(lib)) lib = which(!is.na(tv))
  if (is.null(pred)) pred = lib

  return(RcppFNN(tv, rt, eps, lib, pred, E, tau, style, dist.metric,
                 threads, higher.parallel, NULL, terra::nrow(data)))
}

#' False Nearest Neighbours
#'
#' @param data Observation data.
#' @param target Integer of column indice for the target variable.
#' @param E (optional) Embedding dimensions.
#' @param tau (optional) Step of lag.
#' @param style (optional) Embedding style (`0` includes current state, `1` excludes it).
#' @param lib (optional) Libraries indices.
#' @param pred (optional) Predictions indices.
#' @param dist.metric (optional) Distance measure to be used.
#' @param rt (optional) Relative distance threshold.
#' @param eps (optional) Absolute distance threshold.
#' @param threads (optional) Number of threads used.
#' @param higher.parallel (optional) Whether to use a higher level of parallelism.
#' @param detrend (optional) Whether to remove the linear trend.
#' @param nb (optional) Neighbours list.
#' @param ... Additional arguments to absorb unused inputs in method dispatch.
#'
#' @return A vector.
#' @export
#' @name fnn
#' @aliases fnn,data.frame-method
#' @references
#' Kennel, M.B., Brown, R., Abarbanel, H.D.I., 1992. Determining embedding dimension for phase-space reconstruction using a geometrical construction. Physical Review A 45, 3403–3411.
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' pc::fnn(columbus, 1, E = 3:5, threads = 1)
#'
methods::setMethod("fnn", "data.frame", .fnn_ts)

#' @rdname fnn
methods::setMethod("fnn", "sf", .fnn_lattice)

#' @rdname fnn
methods::setMethod("fnn", "SpatRaster", .fnn_grid)
