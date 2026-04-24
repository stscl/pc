
.pc_ts = \(data, target, source, libsizes = NULL, E = 3, k = E+2, tau = 1, style = 1, lib = NULL, pred = NULL, boot = 99, 
           random = TRUE, seed = 42L, dist.metric = c("euclidean", "manhattan", "maximum"), zero.tolerance = max(k),
           relative = TRUE, weighted = TRUE, threads = detectThreads(), lower.parallel = TRUE, verbose = TRUE, h = 0) {
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
                relative = TRUE, weighted = TRUE, threads = detectThreads(), lower.parallel = TRUE, verbose = TRUE, nb = NULL) {
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
             relative = TRUE, weighted = TRUE, threads = detectThreads(), lower.parallel = TRUE, verbose = TRUE) {
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

#' SURD
#' 
#' Synergistic-Unique-Redundant Decomposition
#' 
#' @note `pc` only supports numeric input data. Both `bin` and `method`
#'   support variable-specific settings using R-style recycling:
#'   \itemize{
#'     \item length 1: applied to the target and all agent variables
#'     \item length 2: first for the target, second for all agents
#'     \item length > 2: first for the target, remaining values are recycled across agents
#'   }
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
#' Martinez-Sanchez, A., Arranz, G., Lozano-Duran, A., 2024. Decomposing causality into its synergistic, unique, and redundant components. Nature Communications 15.
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' pc::pc(columbus, 1, 2:3)
#'
methods::setMethod("pc", "data.frame", .pc_ts)

#' @rdname pc
methods::setMethod("pc", "sf", .pc_lattice)

#' @rdname pc
methods::setMethod("pc", "SpatRaster", .pc_grid)
