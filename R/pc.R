.pc_ts = \(data, target, source, libsizes = NULL, E = 3, k = E, tau = 1, style = 1, lib = NULL, pred = NULL, boot = 99,
           replace = TRUE, seed = 42L, dist.metric = c("euclidean", "manhattan", "maximum"), zero.tolerance = max(k),
           relative = TRUE, weighted = TRUE, threads = length(libsizes), higher.parallel = TRUE, verbose = TRUE, h = 0, ...) {
  dist.metric = match.arg(dist.metric)
  dlist = .validate_var(data, target, source)
  tv = dlist[[1]]; sv = dlist[[2]]
  if (is.null(lib)) lib = which(!(is.na(tv) | is.na(sv)))
  if (is.null(pred)) pred = lib

  if (is.null(libsizes)) {
    return(RcppPC(tv, sv, lib, pred, E, tau, style, k, zero.tolerance,
                  dist.metric, relative, weighted, threads, h, NULL, NULL))
  } else {
    return(RcppPCboot(tv, sv, libsizes, lib, pred, E, tau, style, k, zero.tolerance, dist.metric, boot,
                      replace, seed, relative, weighted, threads, higher.parallel, verbose, h, NULL, NULL))
  }
}

.pc_lattice = \(data, target, source, libsizes = NULL, E = 3, k = E+1, tau = 1, style = 1, lib = NULL, pred = NULL, boot = 99,
                replace = TRUE, seed = 42L, dist.metric = c("euclidean", "manhattan", "maximum"), zero.tolerance = max(k), 
                relative = TRUE, weighted = TRUE, threads = length(libsizes), higher.parallel = TRUE, verbose = TRUE, detrend = FALSE, nb = NULL, ...) {
  if (is.null(nb)) nb = sdsfun::spdep_nb(data)
  dist.metric = match.arg(dist.metric)
  dlist = .validate_var(data, target, source, detrend)
  tv = dlist[[1]]; sv = dlist[[2]]
  if (is.null(lib)) lib = which(!(is.na(tv) | is.na(sv)))
  if (is.null(pred)) pred = lib

  if (is.null(libsizes)) {
    return(RcppPC(tv, sv, lib, pred, E, tau, style, k, zero.tolerance,
                  dist.metric, relative, weighted, threads, 0, nb, NULL))
  } else {
    return(RcppPCboot(tv, sv, libsizes, lib, pred, E, tau, style, k, zero.tolerance, dist.metric, boot,
                      replace, seed, relative, weighted, threads, higher.parallel, verbose, 0, nb, NULL))
  }
}

.pc_grid = \(data, target, source, libsizes = NULL, E = 3, k = E+1, tau = 1, style = 1, lib = NULL, pred = NULL, boot = 99,
             replace = TRUE, seed = 42L, dist.metric = c("euclidean", "manhattan", "maximum"), zero.tolerance = max(k), 
             relative = TRUE, weighted = TRUE, threads = length(libsizes), higher.parallel = TRUE, verbose = TRUE, detrend = FALSE, ...) {
  dist.metric = match.arg(dist.metric)
  dlist = .validate_var(data, target, source, detrend)
  tv = dlist[[1]]; sv = dlist[[2]]
  if (is.null(lib)) lib = which(!(is.na(tv) | is.na(sv)))
  if (is.null(pred)) pred = lib

  if (is.null(libsizes)) {
    return(RcppPC(tv, sv, lib, pred, E, tau, style, k, zero.tolerance, dist.metric,
                  relative, weighted, threads, 0, NULL, terra::nrow(data)))
  } else {
    return(RcppPCboot(tv, sv, libsizes, lib, pred, E, tau, style, k, zero.tolerance, dist.metric, boot, replace,
                      seed, relative, weighted, threads, higher.parallel, verbose, 0, NULL, terra::nrow(data)))
  }
}

#' Pattern Causality
#'
#' @inheritParams fnn
#' @param source Integer of column indice for the source variable.
#' @param libsizes (optional) Number of observations used.
#' @param k (optional) Number of nearest neighbors used for projection.
#' @param boot (optional) Number of bootstraps to perform.
#' @param replace (optional) Should sampling be with replacement?
#' @param seed (optional) Random seed.
#' @param zero.tolerance (optional) Maximum number of zeros tolerated in signature space.
#' @param relative (optional) Whether to calculate relative changes in embedding.
#' @param weighted (optional) Whether to weight causal strength.
#' @param verbose (optional) Whether to show the progress bar.
#' @param h (optional) Prediction horizon.
#'
#' @return A list.
#' \describe{
#'   \item{causality}{
#'     A data.frame of causality results. When `libsizes` is `NULL`, it contains
#'     per-sample causality estimates; otherwise, it contains causality results
#'     evaluated across different library sizes.
#'   }
#'   \item{summary}{
#'     A data.frame summarizing overall causality metrics. Only returned when
#'     `libsizes` is `NULL`.
#'   }
#' }
#'
#' @export
#' @name pc
#' @aliases pc,data.frame-method
#' @references
#' Stavroglou, S.K., Pantelous, A.A., Stanley, H.E., Zuev, K.M., 2020. Unveiling causal interactions in complex systems. Proceedings of the National Academy of Sciences 117, 7599–7605.
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' pc::pc(columbus, 1, 3, E = 3, k = 5, threads = 1)
#'
methods::setMethod("pc", "data.frame", .pc_ts)

#' @rdname pc
methods::setMethod("pc", "sf", .pc_lattice)

#' @rdname pc
methods::setMethod("pc", "SpatRaster", .pc_grid)
