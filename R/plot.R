#' Plot methods for the triptych classes
#'
#' @param x An object that inherits from one of the triptych classes.
#' @param object An object that inherits from one of the triptych classes.
#' @param ... Arguments passed from `autoplot.triptych()` to the other methods for triptych classes.
#'
#'
#' @return
#' For an object of class `'triptych'`: A patchwork object (invisibly).
#'
#' For all other triptych objects: A ggplot object (invisibly).
#'
#' Every `plot()` method wraps the corresponding `autoplot()` method,
#' followed by an explicit `print()` call.
#' That is, it always draws a plot, even during assignment or within a loop.
#'
#' @examples
#' data(ex_binary, package = "triptych")
#' tr <- triptych(ex_binary)
#'
#' dplyr::slice(tr, 1, 3, 6, 9) |> autoplot()
#' autoplot(tr$murphy)
#' autoplot(tr$reliability)
#' autoplot(tr$roc)
#' autoplot(tr$mcbdsc)
#'
#' @name plot.triptych
NULL

#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot

#' @rdname plot.triptych
#' @export
plot.triptych <- function(x, ...) {
  p <- ggplot2::autoplot(x, ...)
  print(p)
}

#' @rdname plot.triptych
#' @export
autoplot.triptych <- function(object, ...) {
  diags <- c("murphy", "reliability", "roc")
  plots <- dplyr::select(object, dplyr::any_of(diags))
  if (!length(plots)) {
    diags_string <- paste(diags, collapse = ", ")
    stop(paste(
      sprintf("No diagrams matching any of: %s", diags_string),
      "Did you try to plot the column of 'forecast' names?",
      sep = "\n"
    ))
  }

  pp <- purrr::map(plots, autoplot_w_facet, ...) |>
    patchwork::wrap_plots(guides = "collect") &
    ggplot2::theme(
      legend.position = "bottom",
      aspect.ratio = 1,
      plot.title = ggplot2::element_text(hjust = 0.5)
    )

  pp
}

autoplot_w_facet <- function(x, ...) {
  if (!has_regions(x) && !inherits(x, "triptych_reliability")) {
    return(autoplot(x, ...))
  }
  autoplot(x, ...) +
    ggplot2::facet_wrap(
      facets = ggplot2::vars(.data$forecast),
      nrow = triptych_facet_wrap_nrow(length(x))
    ) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_blank()
    ) +
    ggplot2::scale_x_continuous(n.breaks = 3) +
    ggplot2::scale_y_continuous(n.breaks = 3)
}

triptych_facet_wrap_nrow <- function(n) ceiling(sqrt(n))

#' @rdname plot.triptych
#' @export
plot.triptych_murphy <- function(x, ...) {
  p <- ggplot2::autoplot(x, ...)
  print(p)
}

#' @rdname plot.triptych
#' @export
autoplot.triptych_murphy <- function(object, ...) {
  nm <- names(object)
  df_e <- estimates(object)
  df_r <- regions(object)
  df_e$forecast <- factor(df_e$forecast, levels = nm)
  if (!is.null(df_r)) df_r$forecast <- factor(df_r$forecast, levels = nm)
  ggplot2::ggplot(df_e) +
    ggplot2::ggtitle("Murphy") +
    ggplot2::ylab("Mean elementary score") +
    ggplot2::xlab("Threshold value") +
    ggplot2::theme_bw() +
    {
      if (!is.null(df_r)) {
        list(
          ggplot2::geom_ribbon(
            mapping = ggplot2::aes(
              x = .data$threshold,
              ymin = .data$lower,
              ymax = .data$upper
            ),
            data = df_r,
            fill = "grey50",
            alpha = 0.3
          ),
          ggplot2::facet_wrap(
            facets = ggplot2::vars(.data$forecast)
          )
        )
      }
    } +
    ggplot2::geom_path(
      mapping = ggplot2::aes(
        x = .data$knot,
        y = .data$mean_score,
        col = .data$forecast,
      )
    )
}


autolayer.triptych_murphy <- function(object, ...) {
  list(
    if (has_regions(object)) {
      ggplot2::geom_ribbon(
        mapping = ggplot2::aes(
          x = .data$threshold,
          ymin = .data$lower,
          ymax = .data$upper,
          fill = .data$forecast
        ),
        data = regions(object),
        alpha = 0.3
      )
    },
    ggplot2::geom_path(
      mapping = ggplot2::aes(
        x = .data$knot,
        y = .data$mean_score,
        col = .data$forecast,
      ),
      data = estimates(object)
    )
  )
}

#' @rdname plot.triptych
#' @export
plot.triptych_reliability <- function(x, ...) {
  p <- ggplot2::autoplot(x, ...)
  print(p)
}

#' @param breaks A vector of bin boundaries for the `geom_histogram()` layer. Set to `NA` to disable.
#'
#' @rdname plot.triptych
#' @export
autoplot.triptych_reliability <- function(object, ..., breaks = seq(0, 1, length.out = 11)) {
  nm <- names(object)
  df_e <- estimates(object)
  df_f <- forecasts(object)
  df_r <- regions(object)
  df_e$forecast <- factor(df_e$forecast, levels = nm)
  df_f$forecast <- factor(df_f$forecast, levels = nm)
  if (!is.null(df_r)) df_r$forecast <- factor(df_r$forecast, levels = nm)
  if (!isTRUE(is.na(breaks))) mdens <- max_density(df_f, breaks)
  ggplot2::ggplot(df_e) +
    ggplot2::ggtitle("Reliability") +
    ggplot2::ylab("Conditional event probability") +
    ggplot2::xlab("Forecast value") +
    ggplot2::theme_bw() +
    ggplot2::theme(aspect.ratio = 1) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, y = 0, xend = 1, yend = 1),
      colour = "black"
    ) +
    {
      if (!isTRUE(is.na(breaks))) {
        ggplot2::geom_histogram(
          mapping = ggplot2::aes(
            x = .data$x,
            y = 0.2 * ggplot2::after_stat(.data$density) / mdens
          ),
          data = df_f,
          fill = NA,
          colour = "black",
          breaks = breaks
        )
      }
    } +
    {
      if (length(object) > 1L) {
        ggplot2::facet_wrap(facets = ggplot2::vars(.data$forecast))
      }
    } +
    {
      if (!is.null(df_r)) {
        ggplot2::geom_ribbon(
          mapping = ggplot2::aes(
            x = .data$x,
            ymin = .data$lower,
            ymax = .data$upper
          ),
          data = df_r,
          fill = "grey50",
          alpha = 0.3
        )
      }
    } +
    ggplot2::geom_path(
      mapping = ggplot2::aes(
        x = .data$x,
        y = .data$CEP,
        col = .data$forecast
      )
    )
}

#' @importFrom graphics hist
max_density <- function(fcs, breaks) {
  dplyr::group_by(fcs, .data$forecast) |>
    dplyr::reframe(density = graphics::hist(.data$x, breaks, plot = FALSE)$density) |>
    dplyr::pull(.data$density) |>
    max()
}

autolayer.triptych_reliability <- function(object, ...) {
  list(
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, y = 0, xend = 1, yend = 1),
      colour = "black"
    ),
    if (has_regions(object)) {
      ggplot2::geom_ribbon(
        mapping = ggplot2::aes(
          x = .data$x,
          ymin = .data$lower,
          ymax = .data$upper,
          fill = .data$forecast
        ),
        data = regions(object),
        alpha = 0.3
      )
    },
    ggplot2::geom_path(
      mapping = ggplot2::aes(
        x = .data$x,
        y = .data$CEP,
        col = .data$forecast
      ),
      data = estimates(object)
    )
  )
}

#' @rdname plot.triptych
#' @export
plot.triptych_roc <- function(x, ...) {
  p <- ggplot2::autoplot(x, ...)
  print(p)
}

#' @rdname plot.triptych
#' @export
autoplot.triptych_roc <- function(object, ...) {
  nm <- names(object)
  df_e <- estimates(object)
  df_r <- regions(object)
  df_e$forecast <- factor(df_e$forecast, levels = nm)
  if (!is.null(df_r)) df_r$forecast <- factor(df_r$forecast, levels = nm)
  ggplot2::ggplot(df_e) +
    ggplot2::ggtitle("ROC") +
    ggplot2::ylab("Hit rate") +
    ggplot2::xlab("False alarm rate") +
    ggplot2::theme_bw() +
    ggplot2::theme(aspect.ratio = 1) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, y = 0, xend = 1, yend = 1),
      colour = "black"
    ) +
    {
      if (!is.null(df_r)) {
        list(
          ggplot2::geom_polygon(
            mapping = ggplot2::aes(
              x = .data$FAR,
              y = .data$HR
            ),
            data = df_r,
            fill = "grey50",
            alpha = 0.3
          ),
          ggplot2::facet_wrap(facets = ggplot2::vars(.data$forecast))
        )
      }
    } +
    ggplot2::geom_path(
      mapping = ggplot2::aes(
        x = .data$FAR,
        y = .data$HR,
        col = .data$forecast
      )
    )
}

autolayer.triptych_roc <- function(object, ...) {
  list(
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, y = 0, xend = 1, yend = 1),
      colour = "black"
    ),
    if (has_regions(object)) {
      ggplot2::geom_polygon(
        mapping = ggplot2::aes(
          x = .data$FAR,
          y = .data$HR,
          fill = .data$forecast
        ),
        data = regions(object),
        alpha = 0.3
      )
    },
    ggplot2::geom_path(
      mapping = ggplot2::aes(
        x = .data$FAR,
        y = .data$HR,
        col = .data$forecast
      ),
      data = estimates(object)
    )
  )
}

#' @rdname plot.triptych
#' @export
plot.triptych_mcbdsc <- function(x, ...) {
  p <- ggplot2::autoplot(x, ...)
  print(p)
}

#' @importFrom scales oob_squish_infinite
#' @importFrom geomtextpath geom_labelabline
#' @importFrom ggrepel geom_text_repel
#'
#' @param n_isolines The number of isolines showing mean scores.
#' @param colour_values A colour specification passed to the `values` argument
#'   of `scale_colour_manual()`. Recycled if length 1.
#' @param colour_unc A colour specification highlighting the UNC component layers.
#' @param MCBDSC_repel A boolean value indicating whether labels should be placed
#'   by the `ggrepel` package.
#' @param MCB_lim The plot limits for the x-axis (the MCB component).
#' @param DSC_lim The plot limits for the y-axis (the DSC component).
#'
#' @rdname plot.triptych
#' @export
autoplot.triptych_mcbdsc <- function(
    object,
    ...,
    n_isolines = 10,
    colour_values = "black",
    colour_unc = "#00BF7D",
    MCBDSC_repel = FALSE,
    MCB_lim = NA,
    DSC_lim = NA
) {
  df_e <- estimates(object)

  # Limits and out-of-bound identification
  default_lims <- function(x) c(0, 1.1 * max(x[is.finite(x)]))
  is_in_range <- function(x, xrange) x >= xrange[1] & x <= xrange[2]
  split_by_out_of_bounds <- function(df, MCB_lim, DSC_lim) {
    oob_state <- factor(
      x = with(df, dplyr::case_when(
        is_in_range(DSC, DSC_lim) & is_in_range(MCB, MCB_lim) ~ "within",
        is_in_range(DSC, DSC_lim) & !is.finite(MCB)           ~ "infty",
        TRUE                                                  ~ "oob"
      )),
      levels = c("within", "infty", "oob")
    )
    res <- split(df, oob_state)
    if (nrow(res$infty)) {
      res$infty$x_geom_text <- MCB_lim[2]
    }
    res
  }
  if (anyNA(MCB_lim)) {
    MCB_lim <- default_lims(df_e$MCB)
  }
  if (anyNA(DSC_lim)) {
    DSC_lim <- default_lims(df_e$DSC)
  }
  df_e_by_state <- split_by_out_of_bounds(df_e, MCB_lim, DSC_lim)
  # Check that the plot is not empty of points!
  if (!nrow(df_e_by_state$within)) {
    warning(paste(
      "The given limits for the MCB-DSC plot exclude all forecasts.",
      "The default choices are used instead."
    ))
    MCB_lim <- default_lims(df_e$MCB)
    DSC_lim <- default_lims(df_e$DSC)
    df_e_by_state <- split_by_out_of_bounds(df_e, MCB_lim, DSC_lim)
  }
  if (nrow(df_e_by_state$oob)) {
    message(paste(
      "The following forecasts are not included in the MCB-DSC plot as their",
      "miscalibration measure is outside the plot limits:",
      paste(df_e_by_state$oob$forecast, collapse = ", ")
    ))
  }

  # Reasonable score values for isolines
  choose_isolines <- function(unc, MCB_lim, DSC_lim) {
    scores <- pretty(x = unc - c(-1.1 * MCB_lim[2], DSC_lim[2]), n = n_isolines)
    tibble::tibble(
      slope = 1,
      intercept = unc - scores,
      label = scores
    ) |>
      # Remove a line if its intercept is too close to zero (less than 1/5 times the line distance).
      dplyr::filter(abs(.data$intercept) > abs(diff(.data$intercept)[1]) / 5)
  }
  df_iso_abline <- choose_isolines(df_e$UNC[1], MCB_lim, DSC_lim)

  # replicate colour_values if it has length 1
  if (length(colour_values) == 1L & nrow(df_e) > 1L) {
    colour_values <- rep(colour_values, nrow(df_e))
  }

  ggplot2::ggplot(df_e) +
    ggplot2::geom_segment(
      data = tibble::tibble(max_val = 2 * max(MCB_lim, DSC_lim)),
      mapping = ggplot2::aes(x = 0, y = 0, xend = .data$max_val, yend = .data$max_val),
      colour = colour_unc,
      linewidth = 1
    ) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(x = 0, y = 0),
      colour = colour_unc,
      fill = colour_unc,
      size = 2,
      shape = 15
    ) +
    geomtextpath::geom_labelabline(
      mapping = ggplot2::aes(
        intercept = 0,
        slope = 1,
        label = paste("UNC:", prettyNum(.data$UNC[1], digits = 3))
      ),
      colour = colour_unc,
      hjust = 0.85,
      size = 7 * 0.36,
      text_only = TRUE,
      boxcolour = NA,
      straight = TRUE
    ) +
    ggplot2::geom_abline(
      data = df_iso_abline,
      mapping = ggplot2::aes(intercept = .data$intercept, slope = .data$slope),
      colour = "gray50"
    ) +
    geomtextpath::geom_labelabline(
      data = df_iso_abline,
      mapping = ggplot2::aes(intercept = .data$intercept, slope = .data$slope, label = .data$label),
      colour = "gray50",
      hjust = 0.85,
      size = 7 * 0.36,
      text_only = TRUE,
      boxcolour = NA,
      straight = TRUE
    ) +
    ggplot2::geom_point(
      data = df_e_by_state$within,
      mapping = ggplot2::aes(x = .data$MCB, y = .data$DSC, colour = .data$forecast)
    ) +
    {
      if (!isTRUE(MCBDSC_repel)) {
        ggplot2::geom_text(
          data = df_e_by_state$within,
          mapping = ggplot2::aes(x = .data$MCB, y = .data$DSC, label = .data$forecast, colour = .data$forecast),
          size = 3,
          vjust = 0,
          hjust = 0,
          check_overlap = TRUE,
          position = ggplot2::position_nudge(
            x = diff(MCB_lim) / 80,
            y = -diff(DSC_lim) / 40
          )
        )
      } else if (isTRUE(MCBDSC_repel)) {
        ggrepel::geom_text_repel(
          data = df_e_by_state$within,
          mapping = ggplot2::aes(x = .data$MCB, y = .data$DSC, label = .data$forecast, colour = .data$forecast),
          size = 3
        )
      }
    } +
    {
      if (nrow(df_e_by_state$infty)) {
        list(
          ggplot2::geom_rug(
            data = df_e_by_state$infty,
            mapping = ggplot2::aes(x = .data$MCB, y = .data$DSC, colour = .data$forecast),
            sides = "r",
            linewidth = 2
          ),
          ggplot2::geom_text(
            data = df_e_by_state$infty,
            mapping = ggplot2::aes(x = .data$x_geom_text, y = .data$DSC, label = .data$forecast, colour = .data$forecast),
            size = 3,
            hjust = 1,
            check_overlap = TRUE
          )
        )
      }
    } +
    ggplot2::scale_colour_manual(values = colour_values) +
    ggplot2::scale_x_continuous(oob = scales::oob_squish_infinite) +
    ggplot2::coord_cartesian(xlim = MCB_lim, ylim = DSC_lim) +
    ggplot2::labs(x = "MCB", y = "DSC") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      aspect.ratio = 1,
      panel.border = ggplot2::element_rect(
        colour = "black",
        fill = NA,
        linewidth = 1
      )
    )
}
