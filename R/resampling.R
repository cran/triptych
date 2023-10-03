#' Bootstrap case resampling for triptych objects
#' 
#' This function is intended to be called from [add_confidence()],
#' by specifying `"resampling_cases"` in the `method` argument.
#' 
#' @details
#' Case resampling assumes independent and identically distributed forecast-observation pairs.
#' A given number of bootstrap samples are the basis for pointwise computed confidence intervals.
#' For every bootstrap sample, we draw forecast-observations pairs with replacement until the size of the original data set is reached.
#'
#' @param x One of the triptych objects.
#' @param level A single value that determines which quantiles of
#'   the bootstrap sample to return. These quantiles envelop `level * n_boot`
#'   bootstrap draws.
#' @param n_boot The number of bootstrap samples.
#' @param ... Additional arguments passed to other methods.
#'
#' @return
#' A list of tibbles that contain the information to draw confidence regions.
#' The length is equal to the number of forecasting methods in `x`.
#'
#' @examples
#' data(ex_binary, package = "triptych")
#' 
#' # Bootstrap resampling is expensive
#' # (the number of bootstrap samples is small to keep execution times short)
#'
#' tr <- triptych(ex_binary) |>
#'   dplyr::slice(1, 9) |>
#'   add_confidence(level = 0.9, method = "resampling_cases", n_boot = 20)
#'
#' @export
resampling_cases <- function(x, level = 0.9, n_boot = 1000, ...) {
  UseMethod("resampling_cases")
}

#' Bootstrap (binary) observation resampling for triptych objects
#' 
#' This function is intended to be called from [add_consistency()] or [add_confidence()],
#' by specifying `"resampling_Bernoulli"` in the respective `method` argument.
#'
#' @details
#' Bootstrap (binary) observation resampling assumes conditionally independent observations given the forecast value.
#' A given number of bootstrap samples are the basis for pointwise computed confidence/consistency intervals.
#' For every bootstrap sample, we sample observations from a Bernoulli distribution conditional on (recalibrated) forecast values.
#'
#' @param x One of the triptych objects.
#' @param level A single value that determines which quantiles of
#'   the bootstrap sample to return. These quantiles envelop `level * n_boot`
#'   bootstrap draws.
#' @param n_boot The number of bootstrap samples.
#' @param ... Additional arguments passed to other methods.
#'
#' @return
#' A list of tibbles that contain the information to draw confidence regions.
#' The length is equal to the number of forecasting methods in `x`.
#'
#' @examples
#' data(ex_binary, package = "triptych")
#' 
#' # Bootstrap resampling is expensive
#' # (the number of bootstrap samples is small to keep execution times short)
#' 
#' tr_consistency <- triptych(ex_binary) |>
#'   dplyr::slice(1, 9) |>
#'   add_consistency(level = 0.9, method = "resampling_Bernoulli", n_boot = 20)
#'
#' tr_confidence <- triptych(ex_binary) |>
#'   dplyr::slice(1, 9) |>
#'   add_confidence(level = 0.9, method = "resampling_Bernoulli", n_boot = 20)
#' 
#' @export
resampling_Bernoulli <- function(x, level = 0.9, n_boot = 1000, ...) {
  UseMethod("resampling_Bernoulli")
}

# # todo
# resampling_residuals <- function(x, ...) {
#   UseMethod("resampling_residuals")
# }

bootstrap_sample_cases <- function(x, y, n_boot, diagram, at, ...) {
  stopifnot(length(x) == length(y))
  n <- length(x)
  evl <- \(x) eval_diag(x, at = at, ...)
  replicate(n_boot, simplify = FALSE, {
    s <- sample(n, replace = TRUE)
    diagram(x[s], y = y[s]) |> evl()
  })
}

#' @importFrom stats rbinom
bootstrap_sample_Bernoulli <- function(x, prob, n_boot, diagram, at, ..., resample_x = TRUE) {
  stopifnot(length(x) == length(prob))
  n <- length(x)
  evl <- \(x) eval_diag(x, at = at, ...)
  if (isTRUE(resample_x)) {
    return(replicate(n_boot, simplify = FALSE, {
      s <- sample(n, replace = TRUE)
      diagram(x[s], y = stats::rbinom(n, 1, prob[s])) |> evl()
    }))
  }
  replicate(n_boot, simplify = FALSE, {
    diagram(x, y = stats::rbinom(n, 1, prob)) |> evl()
  })
}

bootstrap_quantile <- function(l, probs) {
  # l is a list from 'replicate' with 'simplify = FALSE'
  # each list element comes from 'eval_diag' which is usually also a list
  unlist(l, recursive = TRUE) |>
    matrix(ncol = length(l)) |>
    apply(MARGIN = 1L, FUN = stats::quantile, probs = probs, na.rm = TRUE, type = 1)
}
