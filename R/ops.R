.ops_ts = \(data, target, source, E = 3:5, k = E, tau = 1, style = 1, lib = NULL, pred = NULL,
            maximize = c("positive", "negative", "dark"), dist.metric = c("euclidean", "manhattan", "maximum"), 
            zero.tolerance = max(k), relative = TRUE, weighted = TRUE, threads = length(E), higher.parallel = TRUE, h = 0) {
  maximize = match.arg(maximize)
  dist.metric = match.arg(dist.metric)
  dlist = .validate_var(data, target, source)
  tv = dlist[[1]]; sv = dlist[[2]]
  if (is.null(lib)) lib = which(!(is.na(tv) | is.na(sv)))
  if (is.null(pred)) pred = lib

  return(RcppPCops(tv, sv, lib, pred, E, tau, k, maximize, style, zero.tolerance,
                   dist.metric, relative, weighted, threads, higher.parallel, h, NULL, NULL))
}

.ops_lattice = \(data, target, source, E = 3:5, k = E+1, tau = 1, style = 1, lib = NULL, pred = NULL, 
                 maximize = c("positive", "negative", "dark"), dist.metric = c("euclidean", "manhattan", "maximum"), 
                 zero.tolerance = max(k), relative = TRUE, weighted = TRUE, threads = length(E), higher.parallel = TRUE, nb = NULL) {
  if (is.null(nb)) nb = sdsfun::spdep_nb(data)
  maximize = match.arg(maximize)
  dist.metric = match.arg(dist.metric)
  dlist = .validate_var(data, target, source)
  tv = dlist[[1]]; sv = dlist[[2]]
  if (is.null(lib)) lib = which(!(is.na(tv) | is.na(sv)))
  if (is.null(pred)) pred = lib

  return(RcppPCops(tv, sv, lib, pred, E, tau, k, maximize, style, zero.tolerance, 
                   dist.metric, relative, weighted, threads, higher.parallel, 0, nb, NULL))
}

.ops_grid = \(data, target, source, E = 3:5, k = E+1, tau = 1, style = 1, lib = NULL, pred = NULL,
              maximize = c("positive", "negative", "dark"), dist.metric = c("euclidean", "manhattan", "maximum"), 
              zero.tolerance = max(k), relative = TRUE, weighted = TRUE, threads = length(E), higher.parallel = TRUE) {
  maximize = match.arg(maximize)
  dist.metric = match.arg(dist.metric)
  dlist = .validate_var(data, target, source)
  tv = dlist[[1]]; sv = dlist[[2]]
  if (is.null(lib)) lib = which(!(is.na(tv) | is.na(sv)))
  if (is.null(pred)) pred = lib

  return(RcppPCops(tv, sv, lib, pred, E, tau, k, maximize, style, zero.tolerance, dist.metric, 
                   relative, weighted, threads, higher.parallel, 0, NULL, terra::nrow(data)))
}

#' Optimal Parameter Search for Pattern Causality
#'
#' @inheritParams pc
#' @param maximize (optional) Causality metric to maximize: one of "positive", "negative", or "dark".
#'
#' @return A list.
#' \describe{
#'   \item{\code{param}}{The selected optimal parameter combination.}
#'   \item{\code{xmap}}{A data.frame containing cross-mapping performance across parameter settings.}
#' }
#'
#' @export
#' @name ops
#' @aliases ops,data.frame-method
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' pc::ops(columbus, 1, 3, E = 5:10, maximize = "negative", threads = 1)
#'
methods::setMethod("ops", "data.frame", .ops_ts)

#' @rdname ops
methods::setMethod("ops", "sf", .ops_lattice)

#' @rdname ops
methods::setMethod("ops", "SpatRaster", .ops_grid)
