
.pc_ts = \(data, target, agent, lag = 1, bin = 5, method = "equal",
             max.order = 10, threads = 1, base = 2.0, normalize = TRUE) {
  mat = .convert2mat(data, contain_type = FALSE)
  return(RcppSURD(mat, abs(target), abs(agent), lag, max.order, 
                  threads, base, normalize, abs(bin), method))
}

.pc_lattice = \(data, target, agent, lag = 1, bin = 5, method = "equal", 
                  max.order = 10, threads = 1, base = 2.0, normalize = TRUE, nb = NULL) {
  if (is.null(nb)) nb = sdsfun::spdep_nb(data)
  mat = .convert2mat(data, contain_type = FALSE)
  return(RcppSURD(mat, abs(target), abs(agent), lag, max.order, 
                  threads, base, normalize, abs(bin), method, nb))
}

.pc_grid = \(data, target, agent, lag = 1, bin = 5, method = "equal",
               max.order = 10, threads = 1, base = 2.0, normalize = TRUE) {
  mat = .convert2mat(data, contain_type = FALSE)
  return(RcppSURD(mat, abs(target), abs(agent), lag, max.order, threads, base, 
                  normalize, abs(bin), method, NULL, terra::nrow(data[[1]])))
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
