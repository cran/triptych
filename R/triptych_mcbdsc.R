#' Evaluation of forecasts using score decompositions
#'
#' A score decomposition splits the mean score into the three components of
#' miscalibration (MCB), discrimination (DSC), and uncertainty (UNC). Plotting
#' the DSC component against the MCB component allows for a quick visual
#' inspection of predictive performance for many forecasting methods.
#'
#' @param score A string specifying the score function.
#'   One of: `"Brier_score"` (default), `"log_score"`, `"MR_score"`.
#' @param ... Unused.
#' @inheritParams triptych
#'
#' @return A `triptych_mcbdsc` object, that is a `vctrs_vctr` subclass, and has
#'   a length equal to number of forecasting methods supplied in `x`. Each entry
#'   is named according to the corresponding forecasting method,
#'   and contains a list of named objects:
#'   \itemize{
#'     \item `estimate`: A data frame of the score decomposition.
#'     \item `region`: An empty list. Adding confidence regions is not yet supported.
#'     \item `x`: The numeric vector of original forecasts.
#'   }
#'   Access is most convenient through [estimates()], [regions()], and [forecasts()].
#'
#' @seealso Accessors: [estimates()], [regions()], [forecasts()], [observations()]
#'
#'   Visualization: [plot.triptych_mcbdsc()], [autoplot.triptych_mcbdsc()]
#'
#' @examples
#' data(ex_binary, package = "triptych")
#'
#' md <- mcbdsc(ex_binary)
#' md
#' 
#' autoplot(md)
#' estimates(md)
#'
#' @name mcbdsc
NULL

#' @rdname mcbdsc
#' @export
mcbdsc <- function(x, y_var = "y", ..., y = NULL, score = "Brier_score") {
  x <- tibble::as_tibble(x)
  if (is.null(y)) {
    y_var <- tidyselect::vars_pull(names(x), !!rlang::enquo(y_var))
    y <- x[[y_var]]
    x <- dplyr::select(x, !y_var)
  }
  y <- vec_cast(y, to = double())
  score <- vec_cast(score, to = character())
  stopifnot(identical(length(score), 1L))
  x <- dplyr::mutate_all(x, vec_cast, to = double())
  vec_cast(x, to = new_triptych_mcbdsc(y = y, score = score))
}

new_triptych_mcbdsc <- function(x = list(), y = double(), score = character()) {
  new_vctr(x, y = y, score = score, class = "triptych_mcbdsc")
}

# formatting
#' @export
format.triptych_mcbdsc <- function(x, ...) {
  sprintf("<named list[%i]>", sapply(x, length))
}
#' @export
vec_ptype_abbr.triptych_mcbdsc <- function(x, ..., prefix_named = FALSE, suffix_shape = TRUE) {
  "trpt_mcbdsc"
}

# coercion

vec_ptype2.triptych_mcbdsc <- function(x, y, ..., x_arg = "", y_arg = "") {
  UseMethod("vec_ptype2.triptych_mcbdsc")
}
#' @export
vec_ptype2.triptych_mcbdsc.triptych_mcbdsc <- function(x, y, ..., x_arg = "", y_arg = "") {
  if (!has_compatible_observations(x, y)) {
    stop_incompatible_type(
      x,
      y,
      x_arg = x_arg,
      y_arg = y_arg,
      details = "Observations are not compatible."
    )
  }
  new_triptych_mcbdsc(list(), observations(x))
}

# casting

#' @param r A reference triptych_mcbdsc object whose attributes are used for casting.
#'
#' @rdname mcbdsc
#' @export
as_mcbdsc <- function(x, r) {
  stopifnot(inherits(r, "triptych_mcbdsc"))
  x <- tibble::as_tibble(x)
  x <- dplyr::mutate_all(x, vec_cast, to = double())
  vec_cast(x, to = r)
}


vec_cast.triptych_mcbdsc <- function(x, to, ...) {
  UseMethod("vec_cast.triptych_mcbdsc")
}
#' @export
vec_cast.triptych_mcbdsc.triptych_mcbdsc <- function(x, to, ..., x_arg = "", to_arg = "") {
  if (!has_compatible_observations(x, to)) {
    stop_incompatible_cast(
      x,
      to,
      x_arg = x_arg,
      to_arg = to_arg,
      details = "Observations are not compatible."
    )
  }
  x
}
#' @export
vec_cast.triptych_mcbdsc.list <- function(x, to, ...) {
  x <- lapply(x, vec_cast, to = to)
  f <- \(...) vec_c(..., .name_spec = "{outer}_{inner}")
  do.call(f, x)
}
#' @export
vec_cast.triptych_mcbdsc.data.frame <- function(x, to, ...) {
  x <- lapply(x, vec_cast, to = to)
  f <- \(...) vec_c(..., .name_spec = "{outer}_{inner}")
  do.call(f, x)
}
#' @export
vec_cast.triptych_mcbdsc.tbl_df <- function(x, to, ...) {
  x <- lapply(x, vec_cast, to = to)
  f <- \(...) vec_c(..., .name_spec = "{outer}_{inner}")
  do.call(f, x)
}
#' @export
vec_cast.triptych_mcbdsc.double <- function(x, to, ...) {
  y <- observations(to)
  score <- attr(to, "score")
  f <- get(score)
  mscore <- mean(f(x, y))
  mscore_r <- mean(f(recalibrate_mean(x, y), y))
  mscore_u <- mean(f(mean(y), y))

  list(
    estimate = tibble::tibble(
      mean_score = mscore,
      MCB = mscore - mscore_r,
      DSC = mscore_u - mscore_r,
      UNC = mscore_u
    ),
    region = list(),
    x = x
  ) |>
    list() |>
    new_triptych_mcbdsc(y = y, score = score)
}

Brier_score <- function(x, y) (x - y)^2
log_score <- function(x, y) ifelse(y > 0.5, -log(x), -log(1 - x))
MR_score <- function(x, y) {
  dplyr::case_when(
    x < 0.5 & y > 0.5 ~ 1.0,
    x > 0.5 & y < 0.5 ~ 1.0,
    x == 0.5          ~ 0.5,
    TRUE              ~ 0.0
  )
}


eval_diag.triptych_mcbdsc <- function(x, ...) {
  purrr::map(x, .f = \(o) with(o$estimate, c(mean_score, MCB, DSC, UNC)))
}

#' @export
observations.triptych_mcbdsc <- function(x, ...) {
  attr(x, "y")
}

#' @export
forecasts.triptych_mcbdsc <- function(x, ...) {
  f <- function(o) tibble::tibble(x = o$x)
  g <- function(...) vec_rbind(..., .names_to = "forecast")
  purrr::map(x, f) |>
    do.call(g, args = _)
}

#' @rdname estimates
#' @export
estimates.triptych_mcbdsc <- function(x, ...) {
  f <- function(o) o$estimate
  g <- function(...) vec_rbind(..., .names_to = "forecast")
  purrr::map(x, f) |>
    do.call(g, args = _)
}

#' @rdname regions
#' @export
regions.triptych_mcbdsc <- function(x, ...) {
  f <- function(o) o$region
  g <- function(...) vec_rbind(..., .names_to = "forecast")
  if (!any(sapply(x, \(o) tibble::is_tibble(o$region)))) return(NULL)
  purrr::map(x, f) |>
    do.call(g, args = _)
}


# add_confidence.triptych_mcbdsc <- function(x, level = 0.9, method = "resampling_cases", ...) {
#   m <- get(method)(x, level, ...)
#   for (i in seq_along(x)) {
#     x[[i]]$region <- m[[i]]
#   }
#   x
# }
# 
# resampling_cases.triptych_mcbdsc <- function(x, level = 0.9, n_boot = 1000, ...) {
#   #saved_seed <- .Random.seed
#   y <- observations(x)
#   warning("Experimental: Bootstrapping may be unreliable for the MCB component.")
#   purrr::map(
#     .x = x,
#     level = level,
#     n_boot = n_boot,
#     .f = function(o, level, n_boot) {
#       bootstrap_sample_cases(o$x, y, n_boot, mcbdsc) |>
#         unlist() |>
#         matrix(nrow = n_boot, byrow = TRUE) |>
#         list(s = _)
#     }
#   )
# }
# 
# resampling_Bernoulli.triptych_mcbdsc <- function(x, level = 0.9, n_boot = 1000, resample_x = TRUE, ...) {
#   #saved_seed <- .Random.seed
#   y <- observations(x)
#   n_obs <- length(y)
#   warning("Experimental: Bootstrapping may be unreliable for the MCB component.")
# 
#   purrr::map(
#     .x = x,
#     level = level,
#     n_boot = n_boot,
#     .f = function(o, level, n_boot) {
#       xr <- recalibrate_mean(o$x, y)
#       bootstrap_sample_Bernoulli(o$x, xr, n_boot, mcbdsc, resample_x = resample_x) |>
#         unlist() |>
#         matrix(nrow = n_boot, byrow = TRUE) |>
#         list(s = _)
#     }
#   )
# }
