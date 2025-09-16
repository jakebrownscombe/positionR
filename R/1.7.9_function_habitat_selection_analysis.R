# Habitat Selection Analysis Functions for Depth

library(dplyr)
library(ggplot2)

#' Analyze depth selection patterns
#'
#' @param habitat_data Data frame with presence/absence points and depth_m
#' @param fish_ids Vector of fish IDs to analyze (default: all)
#' @param time_periods Vector of time periods to analyze (default: all)
#' @param depth_bins Number of depth bins for analysis (default: 20)
#' @return List with depth selection results
analyze_depth_selection <- function(habitat_data,
                                    fish_ids = NULL,
                                    time_periods = NULL,
                                    depth_bins = 20) {

  cat("=== DEPTH SELECTION ANALYSIS ===\n")

  # Filter data if specified
  if (!is.null(fish_ids)) {
    habitat_data <- habitat_data %>% filter(fish_id %in% fish_ids)
  }
  if (!is.null(time_periods)) {
    habitat_data <- habitat_data %>% filter(time_period %in% time_periods)
  }

  # Get data summary
  presence_data <- habitat_data %>% filter(type == "presence")
  absence_data <- habitat_data %>% filter(type == "absence")

  cat("Total points:", nrow(habitat_data), "\n")
  cat("Presence points:", nrow(presence_data), "\n")
  cat("Absence points:", nrow(absence_data), "\n")

  # Depth range analysis
  depth_range <- range(habitat_data$depth_m, na.rm = TRUE)
  cat("Depth range:", round(depth_range, 2), "m\n")

  # Create depth bins
  depth_breaks <- seq(depth_range[1], depth_range[2], length.out = depth_bins + 1)
  habitat_data$depth_bin <- cut(habitat_data$depth_m, breaks = depth_breaks, include.lowest = TRUE)

  # Calculate selection ratios by depth bin
  selection_analysis <- habitat_data %>%
    group_by(depth_bin) %>%
    summarise(
      n_presence = sum(type == "presence"),
      n_absence = sum(type == "absence"),
      total_points = n(),
      presence_prop = n_presence / total_points,
      depth_mid = mean(depth_m, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      # Calculate selection ratio (use vs availability)
      availability_prop = n_absence / sum(n_absence),
      use_prop = n_presence / sum(n_presence),
      selection_ratio = ifelse(availability_prop > 0, use_prop / availability_prop, NA),
      log_selection_ratio = log(selection_ratio)
    )

  # Summary statistics by group
  group_summary <- habitat_data %>%
    group_by(fish_id, time_period, type) %>%
    summarise(
      n_points = n(),
      mean_depth = mean(depth_m, na.rm = TRUE),
      median_depth = median(depth_m, na.rm = TRUE),
      sd_depth = sd(depth_m, na.rm = TRUE),
      min_depth = min(depth_m, na.rm = TRUE),
      max_depth = max(depth_m, na.rm = TRUE),
      q25_depth = quantile(depth_m, 0.25, na.rm = TRUE),
      q75_depth = quantile(depth_m, 0.75, na.rm = TRUE),
      .groups = 'drop'
    )

  # Overall comparison
  overall_comparison <- habitat_data %>%
    group_by(type) %>%
    summarise(
      n_points = n(),
      mean_depth = mean(depth_m, na.rm = TRUE),
      median_depth = median(depth_m, na.rm = TRUE),
      sd_depth = sd(depth_m, na.rm = TRUE),
      .groups = 'drop'
    )

  cat("\nOverall depth comparison:\n")
  print(overall_comparison)

  # Statistical test (if both presence and absence exist)
  stat_test <- NULL
  if (nrow(presence_data) > 0 && nrow(absence_data) > 0) {
    stat_test <- tryCatch({
      test_result <- wilcox.test(presence_data$depth_m, absence_data$depth_m)
      list(
        test = "Wilcoxon rank-sum test",
        p_value = test_result$p.value,
        statistic = test_result$statistic,
        significant = test_result$p.value < 0.05
      )
    }, error = function(e) {
      list(test = "Failed", error = e$message)
    })
  }

  cat("=== ANALYSIS COMPLETE ===\n")

  return(list(
    selection_by_depth = selection_analysis,
    group_summary = group_summary,
    overall_comparison = overall_comparison,
    statistical_test = stat_test,
    data_summary = list(
      total_points = nrow(habitat_data),
      n_presence = nrow(presence_data),
      n_absence = nrow(absence_data),
      depth_range = depth_range,
      fish_ids = sort(unique(habitat_data$fish_id)),
      time_periods = sort(unique(habitat_data$time_period))
    )
  ))
}

#' Plot depth selection patterns
#'
#' @param habitat_data Data frame with presence/absence points and depth_m
#' @param fish_select Vector of fish IDs to plot (default: all)
#' @param time_select Vector of time periods to plot (default: all)
#' @param plot_type Type of plot: "histogram", "boxplot", "density", "selection_ratio"
#' @param depth_bins Number of bins for histogram/selection ratio plots (default: 20)
#' @return ggplot object
#' @export
plot_depth_selection <- function(habitat_data,
                                 fish_select = NULL,
                                 time_select = NULL,
                                 plot_type = "histogram",
                                 depth_bins = 20) {

  # Filter data if specified
  if (!is.null(fish_select)) {
    habitat_data <- habitat_data %>% filter(fish_id %in% fish_select)
  }
  if (!is.null(time_select)) {
    habitat_data <- habitat_data %>% filter(time_period %in% time_select)
  }

  cat("Plotting depth selection for", nrow(habitat_data), "points\n")

  if (plot_type == "histogram") {
    p <- ggplot(habitat_data, aes(x = depth_m, fill = type)) +
      geom_histogram(bins = depth_bins, alpha = 0.7, position = "identity") +
      scale_fill_manual(values = c("presence" = "blue", "absence" = "red")) +
      labs(
        title = "Depth Distribution: Presence vs Absence",
        x = "Depth (m)",
        y = "Number of Points",
        fill = "Point Type"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")

  } else if (plot_type == "boxplot") {
    p <- ggplot(habitat_data, aes(x = type, y = depth_m, fill = type)) +
      geom_boxplot(alpha = 0.7) +
      scale_fill_manual(values = c("presence" = "blue", "absence" = "red")) +
      labs(
        title = "Depth Distribution: Presence vs Absence",
        x = "Point Type",
        y = "Depth (m)",
        fill = "Point Type"
      ) +
      theme_minimal() +
      theme(legend.position = "none")

  } else if (plot_type == "density") {
    p <- ggplot(habitat_data, aes(x = depth_m, fill = type)) +
      geom_density(alpha = 0.7) +
      scale_fill_manual(values = c("presence" = "blue", "absence" = "red")) +
      labs(
        title = "Depth Density: Presence vs Absence",
        x = "Depth (m)",
        y = "Density",
        fill = "Point Type"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")

  } else if (plot_type == "selection_ratio") {
    # Calculate selection ratios
    depth_range <- range(habitat_data$depth_m, na.rm = TRUE)
    depth_breaks <- seq(depth_range[1], depth_range[2], length.out = depth_bins + 1)
    habitat_data$depth_bin <- cut(habitat_data$depth_m, breaks = depth_breaks, include.lowest = TRUE)

    # Check if we need to group by fish_id and time_period
    grouping_vars <- c("depth_bin")
    if (length(unique(habitat_data$fish_id)) > 1) {
      grouping_vars <- c(grouping_vars, "fish_id")
    }
    if (length(unique(habitat_data$time_period)) > 1) {
      grouping_vars <- c(grouping_vars, "time_period")
    }

    selection_data <- habitat_data %>%
      group_by(across(all_of(grouping_vars))) %>%
      summarise(
        n_presence = sum(type == "presence"),
        n_absence = sum(type == "absence"),
        depth_mid = mean(depth_m, na.rm = TRUE),
        .groups = 'drop'
      )

    # Calculate selection ratios within each group
    if ("fish_id" %in% grouping_vars || "time_period" %in% grouping_vars) {
      # Group-specific selection ratios
      group_vars_for_totals <- grouping_vars[grouping_vars != "depth_bin"]
      selection_data <- selection_data %>%
        group_by(across(all_of(group_vars_for_totals))) %>%
        mutate(
          availability_prop = n_absence / sum(n_absence),
          use_prop = n_presence / sum(n_presence),
          selection_ratio = ifelse(availability_prop > 0, use_prop / availability_prop, NA)
        ) %>%
        ungroup()
    } else {
      # Overall selection ratios
      selection_data <- selection_data %>%
        mutate(
          availability_prop = n_absence / sum(n_absence),
          use_prop = n_presence / sum(n_presence),
          selection_ratio = ifelse(availability_prop > 0, use_prop / availability_prop, NA)
        )
    }

    selection_data <- selection_data %>%
      filter(!is.na(selection_ratio) & is.finite(selection_ratio))

    p <- ggplot(selection_data, aes(x = depth_mid, y = selection_ratio)) +
      geom_col(fill = "steelblue", alpha = 0.7) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +
      labs(
        title = "Depth Selection Ratios",
        subtitle = "Values > 1 indicate selection, < 1 indicate avoidance",
        x = "Depth (m)",
        y = "Selection Ratio (Use/Availability)",
        caption = "Red line shows neutral selection (ratio = 1)"
      ) +
      theme_minimal()

    # Add faceting only if we have the grouping variables in the data
    if ("fish_id" %in% names(selection_data) && "time_period" %in% names(selection_data)) {
      if (length(unique(selection_data$fish_id)) > 1 && length(unique(selection_data$time_period)) > 1) {
        p <- p + facet_grid(fish_id ~ time_period, labeller = label_both)
      }
    } else if ("fish_id" %in% names(selection_data)) {
      if (length(unique(selection_data$fish_id)) > 1) {
        p <- p + facet_wrap(~fish_id, labeller = label_both)
      }
    } else if ("time_period" %in% names(selection_data)) {
      if (length(unique(selection_data$time_period)) > 1) {
        p <- p + facet_wrap(~time_period, labeller = label_both)
      }
    }
  } else {
    # Add faceting for non-selection-ratio plots
    if (length(unique(habitat_data$fish_id)) > 1 && length(unique(habitat_data$time_period)) > 1) {
      p <- p + facet_grid(fish_id ~ time_period, labeller = label_both)
    } else if (length(unique(habitat_data$fish_id)) > 1) {
      p <- p + facet_wrap(~fish_id, labeller = label_both)
    } else if (length(unique(habitat_data$time_period)) > 1) {
      p <- p + facet_wrap(~time_period, labeller = label_both)
    }
  }

  return(p)
}

#' Create comparison plots for multiple fish or time periods
#'
#' @param habitat_data Data frame with presence/absence points and depth_m
#' @param comparison_var Variable to compare across: "fish_id", "time_period", or "time_period_posix"
#' @param plot_type Type of plot: "boxplot", "violin", "density_ridges", "smooth"
#' @return ggplot object
#' @export
plot_depth_comparison <- function(habitat_data,
                                  comparison_var = "fish_id",
                                  plot_type = "boxplot") {

  if (!comparison_var %in% c("fish_id", "time_period", "time_period_posix")) {
    stop("comparison_var must be 'fish_id', 'time_period', or 'time_period_posix'")
  }

  # Create comparison variable as factor
  habitat_data[[paste0(comparison_var, "_factor")]] <- as.factor(habitat_data[[comparison_var]])

  if (plot_type == "boxplot") {
    p <- ggplot(habitat_data, aes_string(x = paste0(comparison_var, "_factor"),
                                         y = "depth_m",
                                         fill = "type")) +
      geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
      scale_fill_manual(values = c("presence" = "blue", "absence" = "red")) +
      labs(
        title = paste("Depth Selection by", stringr::str_to_title(gsub("_", " ", comparison_var))),
        x = stringr::str_to_title(gsub("_", " ", comparison_var)),
        y = "Depth (m)",
        fill = "Point Type"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")

  } else if (plot_type == "violin") {
    p <- ggplot(habitat_data, aes_string(x = paste0(comparison_var, "_factor"),
                                         y = "depth_m",
                                         fill = "type")) +
      geom_violin(alpha = 0.7, position = position_dodge(width = 0.8)) +
      scale_fill_manual(values = c("presence" = "blue", "absence" = "red")) +
      labs(
        title = paste("Depth Distribution by", stringr::str_to_title(gsub("_", " ", comparison_var))),
        x = stringr::str_to_title(gsub("_", " ", comparison_var)),
        y = "Depth (m)",
        fill = "Point Type"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")

  } else if (plot_type == "density_ridges") {
    if (requireNamespace("ggridges", quietly = TRUE)) {
      p <- ggplot(habitat_data, aes_string(x = "depth_m",
                                           y = paste0(comparison_var, "_factor"),
                                           fill = "type")) +
        ggridges::geom_density_ridges(alpha = 0.7) +
        scale_fill_manual(values = c("presence" = "blue", "absence" = "red")) +
        labs(
          title = paste("Depth Distribution by", stringr::str_to_title(gsub("_", " ", comparison_var))),
          x = "Depth (m)",
          y = stringr::str_to_title(gsub("_", " ", comparison_var)),
          fill = "Point Type"
        ) +
        theme_minimal() +
        theme(legend.position = "bottom")
    } else {
      warning("ggridges package not available, using boxplot instead")
      p <- plot_depth_comparison(habitat_data, comparison_var, "boxplot")
    }
    
  } else if (plot_type == "smooth") {
    # Treat comparison variable as continuous (especially useful for time_period_posix)
    p <- ggplot(habitat_data, aes_string(x = comparison_var, y = "depth_m", color = "type")) +
      geom_point(alpha = 0.4, size = 0.8) +
      geom_smooth(method = "loess", se = TRUE, alpha = 0.2, size = 1.2) +
      scale_color_manual(values = c("presence" = "blue", "absence" = "red")) +
      labs(
        title = paste("Depth Selection Trends by", stringr::str_to_title(gsub("_", " ", comparison_var))),
        x = stringr::str_to_title(gsub("_", " ", comparison_var)),
        y = "Depth (m)",
        color = "Point Type"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    # Add better time axis formatting if using POSIXct
    if (comparison_var == "time_period_posix") {
      p <- p + scale_x_datetime(date_labels = "%H:%M", date_breaks = "2 hours")
    }
  }

  return(p)
}

#' Create spatial plot of presence/absence points colored by depth
#'
#' @param habitat_data Data frame with x, y, depth_m, type columns
#' @param fish_select Fish ID to plot (default: 1)
#' @param time_select Time period to plot (default: 0)
#' @param background_raster Optional background raster
#' @param point_size Size of points (default: 1)
#' @return ggplot object
plot_depth_spatial <- function(habitat_data,
                               fish_select = 1,
                               time_select = 0,
                               background_raster = NULL,
                               point_size = 1) {

  # Filter data
  plot_data <- habitat_data %>%
    filter(fish_id == fish_select & time_period == time_select)

  cat("Plotting", nrow(plot_data), "points for fish", fish_select, "time", time_select, "\n")

  if (nrow(plot_data) == 0) {
    stop("No data found for fish ", fish_select, " and time period ", time_select)
  }

  # Create base plot
  p <- ggplot()

  # Add background raster if provided
  if (!is.null(background_raster)) {
    tryCatch({
      if (requireNamespace("raster", quietly = TRUE)) {
        x_range <- range(plot_data$x)
        y_range <- range(plot_data$y)
        x_buffer <- diff(x_range) * 0.1
        y_buffer <- diff(y_range) * 0.1

        plot_extent <- raster::extent(
          x_range[1] - x_buffer, x_range[2] + x_buffer,
          y_range[1] - y_buffer, y_range[2] + y_buffer
        )

        cropped_raster <- raster::crop(background_raster, plot_extent)
        raster_df <- raster::as.data.frame(cropped_raster, xy = TRUE)

        if (ncol(raster_df) >= 3) {
          names(raster_df)[3] <- "values"
          raster_df <- raster_df[!is.na(raster_df$values), ]

          p <- p +
            geom_raster(data = raster_df, aes(x = x, y = y, alpha = values), fill = "lightgray") +
            scale_alpha_continuous(range = c(0.1, 0.4), guide = "none")
        }
      }
    }, error = function(e) {
      warning("Could not add background raster: ", e$message)
    })
  }

  # Add points colored by depth and shaped by type
  p <- p +
    geom_point(data = plot_data,
               aes(x = x, y = y, color = depth_m, shape = type),
               size = point_size,
               alpha = 0.8) +
    scale_color_viridis_c(name = "Depth (m)", option = "viridis") +
    scale_shape_manual(name = "Type", values = c("presence" = 16, "absence" = 1)) +
    coord_fixed() +
    labs(
      title = paste("Spatial Distribution of Points by Depth - Fish", fish_select, "Time", time_select),
      subtitle = paste("Presence points (filled) vs Absence points (hollow)"),
      x = "X Coordinate (m)",
      y = "Y Coordinate (m)"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  return(p)
}
