.dmi_ts = \(data, target, tau = 1:10, pred = NULL, k = 3,
            base = exp(1), normalize = FALSE, threads = length(tau), ...) {
  tv = .validate_var(data, target)[[1]]
  if (is.null(pred)) pred = which(!is.na(tv))
  return(RcppDMI(tv, tau, pred, k, 0, base, normalize, threads))
}

#' Delayed Mutual Information
#'
#' @param data Observation data.
#' @param target Integer of column indice for the target variable.
#' @param E (optional) Embedding dimensions.
#' @param k (optional) Number of nearest neighbors used for evaluation.
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
#' @name dmi
#' @aliases dmi,data.frame-method
#' @references
#' Kennel, M.B., Brown, R., Abarbanel, H.D.I., 1992. Determining embedding dimension for phase-space reconstruction using a geometrical construction. Physical Review A 45, 3403–3411.
#'
#' @examples
#' abun = readr::read_csv(system.file("case/abundance.csv", package = "pc"))
#' pc::dmi(abun, 2)
#'
methods::setMethod("dmi", "data.frame", .dmi_ts)
