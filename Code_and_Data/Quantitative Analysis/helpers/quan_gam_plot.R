
# ------------------------------------------------------------
# Plot a GAM smooth with cluster-robust confidence intervals
#
# This function:
#   - evaluates a specified smooth term on a grid of values
#   - computes pointwise cluster-robust confidence intervals
#     using a user-supplied covariance matrix
#   - produces and saves a ggplot figure
#
# The smooth is shown on the linear predictor (link) scale.
#
# Args:
#   model        : Fitted mgcv::gam object.
#   data         : Original data used for model fitting.
#   var          : Name of the smooth variable (character).
#   smooth_label : Label of the smooth term (e.g. "s(x)").
#   V            : Cluster-robust covariance matrix.
#   x_label      : Label for the x-axis.
#   y_label      : Label for the y-axis.
#   file_name    : Output file name (PDF).
#   n            : Number of grid points for smooth evaluation.
#   n_ticks      : Number of x-axis ticks.
#   is_time      : Logical; whether the x-axis represents calendar time.
#   time_map     : Optional data frame mapping var to calendar dates.
#   date_var     : Name of the date variable in time_map (character).
#   date_fmt     : Date format for x-axis labels.
#
# Returns:
#   A ggplot object representing the smooth with confidence intervals.
#   The figure is also saved to disk as a PDF file.
# ------------------------------------------------------------
plot_smooth_clusterCI <- function(
    model,
    data,
    var,
    smooth_label,
    V,
    x_label,
    y_label,
    file_name,
    n = 300,
    n_ticks = 6,
    is_time = FALSE,
    time_map = NULL,
    date_var = NULL,
    date_fmt = "%Y-%m"
) {
  
  # Construct grid for smooth evaluation
  x_seq <- seq(
    min(data[[var]], na.rm = TRUE),
    max(data[[var]], na.rm = TRUE),
    length.out = n
  )
  
  nd <- data[rep(1, n), ]
  nd[[var]] <- x_seq
  
  # Compute smooth estimates with cluster-robust confidence intervals
  sm <- smooth_ci(
    model        = model,
    newdata      = nd,
    smooth_label = smooth_label,
    V            = V
  ) %>%
    dplyr::rename(
      x     = !!rlang::sym(var),
      est   = fit,
      lower = lwr,
      upper = upr
    )
  
  # Construct x-axis scale
  if (!is_time) {
    
    unique_x <- sm %>% dplyr::select(x) %>% dplyr::distinct() %>% dplyr::arrange(x)
    n_ticks  <- min(n_ticks, nrow(unique_x))
    tick_ids <- round(seq(1, nrow(unique_x), length.out = n_ticks))
    tick_df  <- unique_x[tick_ids, , drop = FALSE]
    
    scale_x <- ggplot2::scale_x_continuous(breaks = tick_df$x)
    
  } else {
    
    if (is.null(time_map) || is.null(date_var)) {
      stop("For is_time = TRUE, both time_map and date_var must be provided.")
    }
    
    time_map2 <- time_map %>%
      dplyr::select(dplyr::all_of(c(var, date_var))) %>%
      dplyr::distinct() %>%
      dplyr::arrange(.data[[var]])
    
    n_ticks  <- min(n_ticks, nrow(time_map2))
    tick_ids <- round(seq(1, nrow(time_map2), length.out = n_ticks))
    tick_df  <- time_map2[tick_ids, ]
    
    scale_x <- ggplot2::scale_x_continuous(
      breaks = tick_df[[var]],
      labels = format(tick_df[[date_var]], date_fmt)
    )
  }
  
  # Create smooth plot
  y_range <- diff(range(sm$est))
  
  p <- ggplot2::ggplot(sm, ggplot2::aes(x = x, y = est)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.25
    ) +
    ggplot2::geom_line(linewidth = 0.4) +
    ggplot2::labs(x = x_label, y = y_label) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::coord_cartesian(
      ylim = c(
        min(sm$est) - 0.5 * y_range,
        max(sm$est) + 0.5 * y_range
      )
    ) +
    scale_x
  
  # Save figure as vector PDF
  ggplot2::ggsave(
    filename = file_name,
    plot     = p,
    device   = "pdf",
    width    = 8,
    height   = 4.5
  )
  
  return(p)
}


# ------------------------------------------------------------
# Forest plot of parametric GAM effects with cluster-robust CIs
#
# This function:
#   - extracts parametric coefficients from a fitted GAM
#   - computes odds ratios and cluster-robust confidence intervals
#   - produces a forest plot on the log-odds-ratio scale
#
# The function is intended for presenting discrete covariate
# effects (e.g. region, sector, month indicators, state).
#
# Args:
#   model          : Fitted mgcv::gam object.
#   V              : Cluster-robust covariance matrix.
#   file_name      : Output file name (PDF).
#   keep_regex     : Optional regular expression used to select which
#                    parametric coefficients are displayed in the forest plot.
#                    For example:
#                      - keep_regex = "vg_westeast" keeps only region effects
#                      - keep_regex = "sector_total|month" keeps sector and
#                        month indicators.
#                    If NULL (default), all parametric terms are included.
#   drop_intercept : Logical; whether to drop the intercept term from the plot.
#                    Defaults to TRUE.
#   level          : Confidence level (default: 0.95).
#   x_label        : Label for the x-axis.
#   title          : Optional plot title.
#   width          : Width of the saved figure.
#   height         : Height of the saved figure.
#
# Returns:
#   A ggplot object representing the forest plot.
#   The figure is also saved to disk as a PDF file.
# ------------------------------------------------------------
forestplot_parametric_clusterCI <- function(
    model,
    V,
    file_name,
    keep_regex = NULL,
    drop_intercept = TRUE,
    level = 0.95,
    x_label = "Odds Ratio (log scale)",
    title = NULL,
    width = 8,
    height = 6
) {
  
  # Cluster-robust coefficient table
  ct_raw <- lmtest::coeftest(model, vcov. = V)
  
  ct <- tibble::tibble(
    term      = rownames(ct_raw),
    Estimate  = ct_raw[, "Estimate"],
    Std_Error = ct_raw[, "Std. Error"],
    p_value   = ct_raw[, "Pr(>|z|)", drop = TRUE]
  )
  
  # Drop intercept if requested
  if (drop_intercept) {
    ct <- dplyr::filter(ct, term != "(Intercept)")
  }
  
  # Keep selected parametric terms only
  if (!is.null(keep_regex)) {
    ct <- dplyr::filter(ct, grepl(keep_regex, term))
  }
  
  if (nrow(ct) == 0) {
    stop("No parametric terms left after filtering.")
  }
  
  # Odds ratios and cluster-robust confidence intervals
  alpha <- 1 - level
  crit  <- qnorm(1 - alpha / 2)
  
  ct <- ct %>%
    dplyr::mutate(
      OR  = exp(Estimate),
      lwr = exp(Estimate - crit * Std_Error),
      upr = exp(Estimate + crit * Std_Error)
    )
  
  # Significance symbols (lmtest convention)
  ct <- ct %>%
    dplyr::mutate(
      stars = dplyr::case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        p_value < 0.10  ~ ".",
        TRUE            ~ ""
      )
    )
  
  # Ordering and labeling
  month_levels <- c(
    "Jan","Feb","Mar","Apr","May","Jun",
    "Jul","Aug","Sep","Oct","Nov","Dec"
  )
  
  ct <- ct %>%
    dplyr::mutate(
      group = dplyr::case_when(
        grepl("^vg_westeast", term) ~ "Region",
        grepl("^fedstaifo", term)   ~ "Region",
        grepl("^sector_total", term)~ "Sector",
        grepl("^month", term)       ~ "Month",
        TRUE                        ~ "Other"
      ),
      label = dplyr::case_when(
        grepl("^vg_westeast", term) ~ sub("^vg_westeast", "", term),
        grepl("^fedstaifo", term)   ~ sub("^fedstaifo", "", term),
        grepl("^sector_total", term)~ sub("^sector_total", "", term),
        grepl("^month", term)       ~ sub("^month", "", term),
        TRUE                        ~ term
      ),
      order_in_group = dplyr::case_when(
        group == "Month" ~ match(label, month_levels),
        TRUE             ~ seq_along(label)
      )
    ) %>%
    dplyr::arrange(
      factor(group, levels = c("Region", "Sector", "Month")),
      order_in_group
    ) %>%
    dplyr::mutate(
      label = factor(label, levels = rev(unique(label)))
    )
  
  # Forest plot
  x_star <- max(ct$upr, na.rm = TRUE) * 1.15
  
  p <- ggplot2::ggplot(ct, ggplot2::aes(x = OR, y = label)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = lwr, xmax = upr),
      height = 0.2
    ) +
    ggplot2::geom_point(size = 1) +
    ggplot2::geom_text(
      ggplot2::aes(x = x_star, label = stars),
      hjust = 0,
      size = 4
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::coord_cartesian(
      xlim = c(min(ct$lwr, na.rm = TRUE), x_star * 1.05)
    ) +
    ggplot2::labs(
      x = x_label,
      y = NULL,
      title = title
    ) +
    ggplot2::theme_minimal(base_size = 12)
  
  # Save figure as vector PDF
  ggplot2::ggsave(
    filename = file_name,
    plot     = p,
    device   = "pdf",
    width    = width,
    height   = height
  )
  
  return(p)
}
