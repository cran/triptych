#' Evaluation of forecasts using a Triptych
#'
#' A triptych visualizes three important aspects of predictive performance:
#' Economic utility via Murphy curves, miscalibration via reliability curves,
#' and discrimination ability via ROC curves.
#' The `triptych` S3 class has plotting methods for `ggplot2`.
#'
#' @param x A data frame, list, matrix, or other object that can be coerced to a tibble. Contains numeric forecasts, and observations (optional).
#' @param y_var A variable in `x` that contains observations. Specified as the argument `var`in [dplyr::pull()].
#' @param ... Additional arguments passed to [murphy()], [reliability()], [roc()], and [mcbdsc()].
#' @param y A numeric vector of observations. If supplied, overrides `y_var`. Otherwise, defaults to `dplyr::pull(x, y_var)`.
#'
#' @return A `triptych` object, that is a tibble subclass, and contains five columns:
#'   * `forecast`: Contains the names.
#'   * `murphy`: Contains a `vctrs_vctr` subclass of Murphy curves.
#'   * `reliability`: Contains a `vctrs_vctr` subclass of reliability curves.
#'   * `roc`: Contains a `vctrs_vctr` subclass of ROC curves.
#'   * `mcbdsc`: Contains a `vctrs_vctr` subclass of score decompositions.
#'
#' @seealso Vector class constructors: [murphy()], [reliability()], [roc()], [mcbdsc()]
#'
#'   Adding uncertainty quantification: [add_consistency()], [add_confidence()]
#'
#'   Visualization: [plot.triptych()], [autoplot.triptych()]
#'
#' @examples
#' data(ex_binary, package = "triptych")
#' 
#' tr <- triptych(ex_binary)
#' identical(tr, triptych(ex_binary, y))
#' identical(tr, triptych(ex_binary, 1))
#' tr
#' 
#' # 1. Choose 4 predictions
#' # 2. Add consistency bands (for reliability curves)
#' #    (Bootstrap resampling is expensive, the number of bootstrap samples
#' #     is small to keep execution times short)
#' # 3. Create patchwork object
#' # 4. Adjust the title of the legend
#' dplyr::slice(tr, 1, 3, 6, 9) |>
#'   add_consistency(level = 0.9, method = "resampling_Bernoulli", n_boot = 20) |>
#'   autoplot() &
#'   ggplot2::guides(colour = ggplot2::guide_legend("Forecast"))
#'
#' @export
triptych <- function(x, y_var = "y", ..., y = NULL) {
  x <- tibble::as_tibble(x)
  if (is.null(y)) {
    y_var <- tidyselect::vars_pull(names(x), !!rlang::enquo(y_var))
    y <- x[[y_var]]
    x <- dplyr::select(x, !y_var)
  }
  y <- vec_cast(y, to = double())
  new_triptych(tibble::tibble(
    forecast = names(x),
    murphy = murphy(x, y = y, ...),
    reliability = reliability(x, y = y, ...),
    roc = roc(x, y = y, ...),
    mcbdsc = mcbdsc(x, y = y, ...)
  ), y)
}

new_triptych <- function(x, y) {
  stopifnot(is.data.frame(x))
  tibble::new_tibble(
    x,
    y = y,
    class = "triptych",
    nrow = nrow(x)
  )
}

# vec_cast.triptych <- function(x, to, ...) {
#   UseMethod("vec_cast.triptych")
# }
vec_cast.triptych.triptych <- function(x, to, ...) {
  triptych_cast(x, to, ...)
}
triptych_cast <- function(x, to, ..., x_arg = "", to_arg = "") {
  out <- tib_cast(x, to, ..., x_arg = x_arg, to_arg = to_arg)
  new_triptych(out, y = attr(to, "y"))
}
vec_cast.triptych.double <- function(x, to, ...) {
  triptych(x, y = attr(to, "y"))
}
# vec_cast.triptych.data.frame <- function(x, to, ...) {
#   triptych(x, y = attr(to, "y"))
# }

vec_ptype2.triptych.triptych <- function(x, y, ...) {
  triptych_ptype2(x, y, ...)
}
triptych_ptype2 <- function(x, y, ...) {
  out <- tib_ptype2(x, y)
  new_triptych(out, y = attr(x, "y"))
}

#' @export
observations.triptych <- function(x, ...) {
  attr(x, "y")
}

#' @export
forecasts.triptych <- function(x, ...) {
  f <- function(o) tibble::tibble(x = o$x)
  g <- function(...) vec_rbind(..., .names_to = "forecast")
  purrr::map(x$murphy, f) |>
    do.call(g, args = _)
}

#' @export
add_confidence.triptych <- function(x, level = 0.9, ...) {
  x$murphy <- add_confidence(x$murphy, level = level, ...)
  x$reliability <- add_confidence(x$reliability, level = level, ...)
  x$roc <- add_confidence(x$roc, level = level, ...)
  x
}

#' @export
add_consistency.triptych <- function(x, level = 0.9, ...) {
  x$reliability <- add_consistency(x$reliability, level = level, ...)
  x
}
