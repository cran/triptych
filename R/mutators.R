#' Adding confidence regions
#' 
#' Confidence regions are supposed to contain the true "parameter" with a given degree of confidence.
#' Here, "parameter" refers to a murphy curve, a reliability curve, or a ROC curve, respectively.
#'
#' @param x An object to which a confidence region should be added.
#' @param level A single value for the level of confidence.
#' @param method A string that gives the name of method to generate the confidence regions. Currently, one of: "resampling_cases", "resampling_Bernoulli".
#' @param ... Additional arguments passed to methods.
#' 
#' @return 
#' The object given to `x`, but with information about the confidence regions.
#' This information can be accessed conveniently by using [regions()] on the
#' curve component of interest.
#' 
#' @examples
#' data(ex_binary, package = "triptych")
#' 
#' tr <- triptych(ex_binary) |>
#'   dplyr::slice(1, 9)
#' 
#' # Bootstrap resampling is expensive
#' # (the number of bootstrap samples is small to keep execution times short)
#' 
#' tr <- add_confidence(tr, level = 0.9, method = "resampling_cases", n_boot = 20)
#' regions(tr$murphy)
#' regions(tr$reliability)
#' regions(tr$roc)
#' 
#' @seealso [resampling_cases()], [resampling_Bernoulli()]
#' 
#' @export
add_confidence <- function(x, level = 0.9, method = "resampling_cases", ...) {
  UseMethod("add_confidence")
}

#' Adding consistency regions for reliability curves
#' 
#' Consistency regions are created under the assumption that the forecasts are calibrated.
#' A reliability curve that significantly violates the consistency region indicates
#' a miscalibrated forecast.
#'
#' @param x An object to which a consistency region should be added.
#' @param level A single value for the level of confidence.
#' @param method A string that gives the name of method to generate the consistency regions. Currently, only: "resampling_Bernoulli".
#' @param ... Additional arguments passed to methods.
#'
#' @return 
#' The object given to `x`, but with information about the consistency regions.
#' This information can be accessed conveniently by using [regions()] on the
#' reliability curve component.
#' 
#' @examples
#' data(ex_binary, package = "triptych")
#' 
#' tr <- triptych(ex_binary) |>
#'   dplyr::slice(1, 9)
#' 
#' # Bootstrap resampling is expensive
#' # (the number of bootstrap samples is small to keep execution times short)
#' 
#' tr <- add_consistency(tr, level = 0.9, method = "resampling_Bernoulli", n_boot = 20)
#' regions(tr$reliability)
#' 
#' @seealso [resampling_Bernoulli()]
#' 
#' @export
add_consistency <- function(x, level = 0.9, method = "resampling_Bernoulli", ...) {
  UseMethod("add_consistency")
}
