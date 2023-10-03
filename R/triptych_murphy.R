#' Evaluation of forecasts using Murphy curves
#'
#' A Murphy curve visualizes economic utility by displaying the mean elementary
#' scores across all threshold values.
#'
#' @param ref_var A variable in `x` that contains reference forecasts. Specified as the argument `var` in [dplyr::pull()]). Ignored, if the default name is not present in `x`. 
#' @param ... Unused.
#' @param ref A numeric vector of reference forecasts. If supplied, overrides `ref_var`. Otherwise, ignored (see `ref_var`) or `dplyr::pull(x, ref_var)`.
#' @inheritParams triptych
#'
#' @return A `triptych_murphy` object, that is a `vctrs_vctr` subclass, and has
#'   a length equal to number of forecasting methods supplied in `x`. Each entry
#'   is named according to the corresponding forecasting method,
#'   and contains a list of named objects:
#'   \itemize{
#'     \item `estimate`: A data frame with the threshold and corresponding mean score values.
#'     \item `region`: Either an empty list, or a data frame of point confidence intervals
#'       added by [add_confidence()].
#'     \item `x`: The numeric vector of original forecasts.
#'   }
#'   Access is most convenient through [estimates()], [regions()], and [forecasts()].
#'
#' @seealso Accessors: [estimates()], [regions()], [forecasts()], [observations()]
#'
#'   Adding uncertainty quantification: [add_confidence()]
#'
#'   Visualization: [plot.triptych_murphy()], [autoplot.triptych_murphy()]
#'
#' @examples
#' data(ex_binary, package = "triptych")
#'
#' mr <- murphy(ex_binary)
#' mr
#' 
#' # 1. Choose 4 predictions
#' # 2. Visualize
#' # 3. Adjust the title of the legend
#' mr[c(1, 3, 6, 9)] |>
#'   autoplot() +
#'   ggplot2::guides(colour = ggplot2::guide_legend("Forecast"))
#'   
#' # Build yourself using accessors
#' library(ggplot2)
#' df_est <- estimates(mr[c(1, 3, 6, 9)])
#' ggplot(df_est) +
#'   geom_path(aes(x = knot, y = mean_score, col = forecast))
#'
#' @name murphy
NULL

#' @rdname murphy
#' @export
murphy <- function(x, y_var = "y", ref_var = "ref", ..., y = NULL, ref = NULL) {
  x <- tibble::as_tibble(x)
  if (is.null(y)) {
    y_var <- tidyselect::vars_pull(names(x), !!rlang::enquo(y_var))
    y <- x[[y_var]]
    x <- dplyr::select(x, !y_var)
  }
  if (is.null(ref)) {
    # ref_var must be present in x, if ref_var is user-supplied...
    # ...otherwise, it's okay if ref_var is not present in x.
    if (!missing(ref_var)) {
      ref_var <- tidyselect::vars_pull(names(x), !!rlang::enquo(ref_var))
      ref <- x[[ref_var]]
      x <- dplyr::select(x, !ref_var)
    } else if (ref_var %in% names(x)) {
      ref <- x[[ref_var]]
      x <- dplyr::select(x, !ref_var)
    }
  }
  y <- vec_cast(y, to = double())
  ref <- vec_cast(ref, to = double())
  x <- dplyr::mutate_all(x, vec_cast, to = double())
  vec_cast(x, to = new_triptych_murphy(y = y, ref = ref))
}

new_triptych_murphy <- function(x = list(), y = double(), ref = double()) {
  new_vctr(x, y = y, ref = ref, class = "triptych_murphy")
}

# formatting
#' @export
format.triptych_murphy <- function(x, ...) {
  sprintf("<named list[%i]>", sapply(x, length))
}
#' @export
vec_ptype_abbr.triptych_murphy <- function(x, ..., prefix_named = FALSE, suffix_shape = TRUE) {
  "trpt_mur"
}

# coercion

vec_ptype2.triptych_murphy <- function(x, y, ..., x_arg = "", y_arg = "") {
  UseMethod("vec_ptype2.triptych_murphy")
}
#' @export
vec_ptype2.triptych_murphy.triptych_murphy <- function(x, y, ..., x_arg = "", y_arg = "") {
  if (!has_compatible_observations(x, y)) {
    stop_incompatible_type(
      x,
      y,
      x_arg = x_arg,
      y_arg = y_arg,
      details = "Observations are not compatible."
    )
  }
  new_triptych_murphy(list(), observations(x), attr(x, "ref"))
}

# casting

#' @param r A reference triptych_murphy object whose attributes are used for casting.
#'
#' @rdname murphy
#' @export
as_murphy <- function(x, r) {
  stopifnot(inherits(r, "triptych_murphy"))
  x <- tibble::as_tibble(x)
  x <- dplyr::mutate_all(x, vec_cast, to = double())
  vec_cast(x, to = r)
}

vec_cast.triptych_murphy <- function(x, to, ...) {
  UseMethod("vec_cast.triptych_murphy")
}
#' @export
vec_cast.triptych_murphy.triptych_murphy <- function(x, to, ..., x_arg = "", to_arg = "") {
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
vec_cast.triptych_murphy.list <- function(x, to, ...) {
  x <- lapply(x, vec_cast, to = to)
  f <- \(...) vec_c(..., .name_spec = "{outer}_{inner}")
  do.call(f, x)
}
#' @export
vec_cast.triptych_murphy.data.frame <- function(x, to, ...) {
  x <- lapply(x, vec_cast, to = to)
  f <- \(...) vec_c(..., .name_spec = "{outer}_{inner}")
  do.call(f, x)
}
#' @export
vec_cast.triptych_murphy.tbl_df <- function(x, to, ...) {
  x <- lapply(x, vec_cast, to = to)
  f <- \(...) vec_c(..., .name_spec = "{outer}_{inner}")
  do.call(f, x)
}
#' @export
vec_cast.triptych_murphy.double <- function(x, to, ...) {
  y <- observations(to)
  ref <- attr(to, "ref")
  list(
    estimate = if (!length(ref)) C_murphydiag_prob(x, y) else C_murphydiag_prob_ref(x, y, ref),
    region = list(),
    x = x
  ) |>
    list() |>
    new_triptych_murphy(y = y, ref = ref)
}


eval_diag.triptych_murphy <- function(x, at, ...) {
  purrr::map(x, at = at, .f = \(o, at) {
    ind <- findInterval(at, o$estimate$knot, all.inside = TRUE)
    w <- with(o$estimate, (at - knot[ind]) / (knot[ind + 1] - knot[ind]))
    with(o$estimate, {
      (1 - w) * val_right[ind] + w * val_left[ind + 1]
    })
  })
}

#' @export
observations.triptych_murphy <- function(x, ...) {
  attr(x, "y")
}

#' @export
forecasts.triptych_murphy <- function(x, ...) {
  f <- function(o) tibble::tibble(x = o$x)
  g <- function(...) vec_rbind(..., .names_to = "forecast")
  purrr::map(x, f) |>
    do.call(g, args = _)
}

#' @rdname estimates
#' @export
estimates.triptych_murphy <- function(x, at = NULL, ...) {
  f <- if (is.null(at)) {
    function(o) {
      tidyr::pivot_longer(
        data = o$estimate,
        cols = dplyr::starts_with("val"),
        names_to = "limit",
        names_prefix = "^val_",
        values_to = "mean_score")
    }
  } else {
    function(o) {
      mean_score <- with(o$estimate, {
        ind <- findInterval(at, knot, all.inside = TRUE)
        w <- (at - knot[ind]) / (knot[ind + 1] - knot[ind])
        (1 - w) * val_right[ind] + w * val_left[ind + 1]
      })
      tibble::tibble(
        knot = at,
        mean_score = mean_score
      )
    }
  }
  g <- function(...) vec_rbind(..., .names_to = "forecast")
  purrr::map(x, f) |>
    do.call(g, args = _)
}

#' @rdname regions
#' @export
regions.triptych_murphy <- function(x, ...) {
  if (!has_regions(x)) return(NULL)
  f <- function(o) o$region
  g <- function(...) vec_rbind(..., .names_to = "forecast")
  purrr::map(x, f) |>
    do.call(g, args = _)
}

#' @export
has_regions.triptych_murphy <- function(x, ...) {
  any(sapply(x, \(o) tibble::is_tibble(o$region)))
}

#' @export
add_confidence.triptych_murphy <- function(x, level = 0.9, method = "resampling_cases", ...) {
  stopifnot(method %in% c("resampling_cases", "resampling_Bernoulli"))
  m <- get(method)(x, level = level, ...)
  for (i in seq_along(x)) {
    x[[i]]$region <- m[[i]]
  }
  x
}

#' @rdname resampling_cases
#' @export
resampling_cases.triptych_murphy <- function(x, level = 0.9, n_boot = 1000, ...) {
  #saved_seed <- .Random.seed
  y <- observations(x)
  thresholds <- seq(0, 1, length.out = 1000)
  purrr::map(
    .x = x,
    level = level,
    n_boot = n_boot,
    .f = function(o, level, n_boot) {
      bounds <- bootstrap_sample_cases(o$x, y, n_boot, murphy, thresholds) |>
        bootstrap_quantile(probs = 0.5 + c(-0.5, 0.5) * level)
      tibble::tibble(
        threshold = thresholds,
        lower = bounds[1L, ],
        upper = bounds[2L, ],
        method = paste0("resampling_cases_", n_boot),
        level = level
      )
    }
  )
}

#' @rdname resampling_Bernoulli
#' @export
resampling_Bernoulli.triptych_murphy <- function(x, level = 0.9, n_boot = 1000, ...) {
  #saved_seed <- .Random.seed
  y <- observations(x)
  thresholds <- seq(0, 1, length.out = 1000)
  purrr::map(
    .x = x,
    level = level,
    n_boot = n_boot,
    .f = function(o, level, n_boot) {
      xr <- recalibrate_mean2(o$x, y)
      bounds <- bootstrap_sample_Bernoulli(o$x, xr, n_boot, murphy, thresholds) |>
        bootstrap_quantile(probs = 0.5 + c(-0.5, 0.5) * level)
      tibble::tibble(
        threshold = thresholds,
        lower = bounds[1L, ],
        upper = bounds[2L, ],
        method = paste0("resampling_Bernoulli_", n_boot),
        level = level
      )
    }
  )
}
