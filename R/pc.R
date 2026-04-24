
.pc_ts = \(data, target, source, libsizes = NULL, E = 3, k = E+2, tau = 1, style = 1, lib = NULL, pred = NULL, boot = 99, 
           random = TRUE, seed = 42L, dist.metric = c("euclidean", "manhattan", "maximum"), zero.tolerance = max(k),
           relative = TRUE, weighted = TRUE, threads = length(libsizes), lower.parallel = TRUE, verbose = TRUE, h = 0) {
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
                      random, seed, relative, weighted, threads, lower.parallel, verbose, h, NULL, NULL))
  }
}

.pc_lattice = \(data, target, source, libsizes = NULL, E = 3, k = E+2, tau = 1, style = 1, lib = NULL, pred = NULL, boot = 99, 
                random = TRUE, seed = 42L, dist.metric = c("euclidean", "manhattan", "maximum"), zero.tolerance = max(k),
                relative = TRUE, weighted = TRUE, threads = length(libsizes), lower.parallel = TRUE, verbose = TRUE, nb = NULL) {
  if (is.null(nb)) nb = sdsfun::spdep_nb(data)
  dist.metric = match.arg(dist.metric)
  dlist = .validate_var(data, target, source)
  tv = dlist[[1]]; sv = dlist[[2]]
  if (is.null(lib)) lib = which(!(is.na(tv) | is.na(sv)))
  if (is.null(pred)) pred = lib
  
  if (is.null(libsizes)) {
    return(RcppPC(tv, sv, lib, pred, E, tau, style, k, zero.tolerance, 
                  dist.metric, relative, weighted, threads, 0, nb, NULL))
  } else {
    return(RcppPCboot(tv, sv, libsizes, lib, pred, E, tau, style, k, zero.tolerance, dist.metric, boot, 
                      random, seed, relative, weighted, threads, lower.parallel, verbose, 0, nb, NULL))
  }
}

.pc_grid = \(data, target, source, libsizes = NULL, E = 3, k = E+2, tau = 1, style = 1, lib = NULL, pred = NULL, boot = 99, 
             random = TRUE, seed = 42L, dist.metric = c("euclidean", "manhattan", "maximum"), zero.tolerance = max(k),
             relative = TRUE, weighted = TRUE, threads = length(libsizes), lower.parallel = TRUE, verbose = TRUE) {
  dist.metric = match.arg(dist.metric)
  dlist = .validate_var(data, target, source)
  tv = dlist[[1]]; sv = dlist[[2]]
  if (is.null(lib)) lib = which(!(is.na(tv) | is.na(sv)))
  if (is.null(pred)) pred = lib
  
  if (is.null(libsizes)) {
    return(RcppPC(tv, sv, lib, pred, E, tau, style, k, zero.tolerance, dist.metric,
                  relative, weighted, threads, 0, NULL, terra::nrow(data)))
  } else {
    return(RcppPCboot(tv, sv, libsizes, lib, pred, E, tau, style, k, zero.tolerance, dist.metric, boot, random, 
                      seed, relative, weighted, threads, lower.parallel, verbose, 0, NULL, terra::nrow(data)))
  }
}

#' Pattern Causality
#' 
#' @note `pc` only supports numeric input data.
#'
#' @inheritParams te
#' @param lag (optional) Lag of the agent variables.
#' @param bin (optional) Number of discretization bins.
#' @param method (optional) Discretization method. One of
#'   `"sd"`, `"equal"`, `"geometric"`, `"quantile"`,
#'   `"natural("jenks")"`, or `"headtail"("headtails")`.
#' @param max.order (optional) Maximum combination order.
#' @param threads (optional) Number of threads used.
#' @param nb (optional) Neighbours list.
#'
#' @return A list.
#' \describe{
#'   \item{vars}{Character vector indicating the variable combination associated with each information component.}
#'   \item{types}{Character vector indicating the information type of each component.}
#'   \item{values}{Numeric vector giving the magnitude of each information component.}
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
#' pc::pc(columbus, 1, 3, E = 6, k = 8)
#'
methods::setMethod("pc", "data.frame", .pc_ts)

#' @rdname pc
methods::setMethod("pc", "sf", .pc_lattice)

#' @rdname pc
methods::setMethod("pc", "SpatRaster", .pc_grid)
