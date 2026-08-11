.dmi_ts = \(data, target, tau = 1:10, pred = NULL, k = 3, base = exp(1), 
            normalize = FALSE, threads = length(tau), ...) {
  tv = .validate_var(data, target)[[1]]
  if (is.null(pred)) pred = which(!is.na(tv))
  return(RcppDMI(tv, tau, pred, k, 0, base, normalize, threads))
}

#' Delayed Mutual Information
#'
#' @inheritParams fnn
#' @param base (optional) Logarithm base of the entropy.
#' @param normalize (optional) Whether to normalize MI values.
#'
#' @return A vector.
#' @export
#' @name dmi
#' @aliases dmi,data.frame-method
#' @references
#' Fraser, A.M., Swinney, H.L., 1986. Independent coordinates for strange attractors from mutual information. Physical Review A 33, 1134–1140.
#' Kraskov, A., Stogbauer, H., Grassberger, P., 2004. Estimating mutual information. Physical Review E 69, 066138.
#'
#' @examples
#' abun = readr::read_csv(system.file("case/abundance.csv", package = "pc"))
#' pc::dmi(abun, 2)
#'
methods::setMethod("dmi", "data.frame", .dmi_ts)
