# Comparative Habitat Selection Analysis Functions
# Compare habitat selection using positioning-based vs space use-based presence-absence data

library(randomForest)
library(dplyr)
library(ggplot2)
library(patchwork)

#' Comparative habitat selection analysis using positioning and/or space use data
#'
#' @param positioning_data Presence-absence data from positioning analysis (optional)
#' @param space_use_data Presence-absence data from space use analysis (optional)
#' @param analysis_type Type of analysis ("positioning", "space_use", or "both")
#' @param formula Model formula (default: presence ~ depth_m + fish_id + time_period)
#' @param sample_size Number of points to sample per dataset (default: 10000)
#' @param ntree Number of trees in random forest (default: 500)
#' @param depth_variable Name of depth variable column (default: "depth_m")
#' @param fish_id_variable Name of fish ID variable column (default: "fish_id")
#' @param time_variable Name of time period variable column (default: "time_period")
#' @param type_variable Name of presence/absence type variable column (default: "type")
#' @param n_depth_points Number of depth points for prediction curves (default: 100)
#' @param create_plots Whether to create plots (default: TRUE)
#' @param create_comparison Whether to create comparison plots when both datasets provided (default: TRUE)
#' @return List containing models, partial dependence data, and plots for each analysis type
#' @export
analyze_comparative_habitat_selection <- function(positioning_data = NULL,
                                                  space_use_data = NULL,
                                                  analysis_type = "both",
                                                  formula = presence ~ depth_m + fish_id + time_period,
                                                  sample_size = 10000,
                                                  ntree = 500,
                                                  depth_variable = "depth_m",
                                                  fish_id_variable = "fish_id", 
                                                  time_variable = "time_period",
                                                  type_variable = "type",
                                                  n_depth_points = 100,
                                                  create_plots = TRUE,
                                                  create_comparison = TRUE) {
  
  cat("=== COMPARATIVE HABITAT SELECTION ANALYSIS ===\n")
  
  # Validate inputs
  if (analysis_type == "both" && (is.null(positioning_data) || is.null(space_use_data))) {
    stop("Both positioning_data and space_use_data required when analysis_type = 'both'")
  }
  if (analysis_type == "positioning" && is.null(positioning_data)) {
    stop("positioning_data required when analysis_type = 'positioning'")
  }
  if (analysis_type == "space_use" && is.null(space_use_data)) {
    stop("space_use_data required when analysis_type = 'space_use'")
  }
  
  results <- list()
  
  # Analyze positioning data
  if (analysis_type %in% c("positioning", "both")) {
    cat("\n--- ANALYZING POSITIONING-BASED DATA ---\n")
    results$positioning <- analyze_single_dataset(
      data = positioning_data,
      dataset_name = "Positioning",
      formula = formula,
      sample_size = sample_size,
      ntree = ntree,
      depth_variable = depth_variable,
      fish_id_variable = fish_id_variable,
      time_variable = time_variable,
      type_variable = type_variable,
      n_depth_points = n_depth_points,
      create_plots = create_plots
    )
  }
  
  # Analyze space use data
  if (analysis_type %in% c("space_use", "both")) {
    cat("\n--- ANALYZING SPACE USE-BASED DATA ---\n")
    results$space_use <- analyze_single_dataset(
      data = space_use_data,
      dataset_name = "Actual Track",
      formula = formula,
      sample_size = sample_size,
      ntree = ntree,
      depth_variable = depth_variable,
      fish_id_variable = fish_id_variable,
      time_variable = time_variable,
      type_variable = type_variable,
      n_depth_points = n_depth_points,
      create_plots = create_plots
    )
  }
  
  # Create comparison plots if both datasets analyzed
  if (analysis_type == "both" && create_comparison) {
    cat("\n--- CREATING COMPARISON PLOTS ---\n")
    results$comparison <- create_comparison_analysis(
      positioning_results = results$positioning,
      space_use_results = results$space_use,
      depth_variable = depth_variable
    )
  }
  
  # Add analysis metadata
  results$metadata <- list(
    analysis_type = analysis_type,
    formula = formula,
    sample_size = sample_size,
    ntree = ntree,
    datasets_analyzed = names(results)[names(results) != "metadata"]
  )
  
  cat("\n=== COMPARATIVE ANALYSIS COMPLETE ===\n")
  
  return(results)
}

#' Analyze a single presence-absence dataset
analyze_single_dataset <- function(data, dataset_name, formula, sample_size, ntree,
                                   depth_variable, fish_id_variable, time_variable, 
                                   type_variable, n_depth_points, create_plots) {
  
  # Prepare data
  # Drop geometry if this is an sf object
  if (inherits(data, "sf")) {
    rf_data <- sf::st_drop_geometry(data)
  } else {
    rf_data <- data
  }
  
  rf_data$presence <- as.factor(ifelse(rf_data[[type_variable]] == "presence", 1, 0))
  
  # Remove rows with missing values
  required_cols <- c(depth_variable, fish_id_variable, "presence")
  
  # Only include time_variable if it exists in the data
  if (time_variable %in% names(rf_data)) {
    required_cols <- c(required_cols, time_variable)
  }
  
  complete_rows <- complete.cases(rf_data[, required_cols])
  rf_data_complete <- rf_data[complete_rows, ]
  
  cat(dataset_name, "data - Original:", nrow(rf_data), "Complete:", nrow(rf_data_complete), "\n")
  
  # Sample data if needed
  if (nrow(rf_data_complete) > sample_size) {
    rf_data_sample <- rf_data_complete %>% sample_n(sample_size)
    cat("Sampled", sample_size, "points for", dataset_name, "model\n")
  } else {
    rf_data_sample <- rf_data_complete
    cat("Using all", nrow(rf_data_sample), "points for", dataset_name, "model\n")
  }
  
  # Fit random forest model
  cat("Fitting", dataset_name, "random forest model...\n")
  rf_model <- randomForest(
    formula,
    data = rf_data_sample,
    ntree = ntree,
    importance = TRUE
  )
  
  cat(dataset_name, "Model OOB Error:", round(rf_model$err.rate[ntree, "OOB"] * 100, 2), "%\n")
  
  # Calculate partial dependence
  depth_range <- range(rf_data_complete[[depth_variable]], na.rm = TRUE)
  depth_seq <- seq(depth_range[1], depth_range[2], length.out = n_depth_points)
  
  fish_ids <- unique(rf_data_complete[[fish_id_variable]])
  
  # Only get time periods if time variable exists
  if (time_variable %in% names(rf_data_complete)) {
    time_periods <- unique(rf_data_complete[[time_variable]])
  } else {
    time_periods <- c("no_time")  # Dummy value when no time aggregation
  }
  
  # Calculate overall partial dependence
  pred_data_overall <- calculate_overall_partial_dependence_single(
    rf_model, depth_seq, fish_ids, time_periods,
    depth_variable, fish_id_variable, time_variable, rf_data_complete
  )
  
  # Calculate by fish and time if plots requested
  pd_data_fish <- NULL
  pd_data_time <- NULL
  
  if (create_plots) {
    pd_data_fish <- calculate_partial_dependence_fish_single(
      rf_model, depth_seq, fish_ids, time_periods,
      depth_variable, fish_id_variable, time_variable, rf_data_complete
    )
    
    pd_data_time <- calculate_partial_dependence_time_single(
      rf_model, depth_seq, fish_ids, time_periods,
      depth_variable, fish_id_variable, time_variable, rf_data_complete
    )
  }
  
  # Create plots
  plots <- list()
  if (create_plots) {
    plots$overall <- create_single_overall_plot(pred_data_overall, depth_variable, dataset_name)
    plots$fish <- create_single_fish_plot(pd_data_fish, pred_data_overall, depth_variable, dataset_name)
    
    # Only create time plot if we have time data
    if (time_variable %in% names(rf_data_complete) && !is.null(pd_data_time) && nrow(pd_data_time) > 0) {
      plots$time <- create_single_time_plot(pd_data_time, pred_data_overall, depth_variable, time_variable, dataset_name)
    }
  }
  
  return(list(
    model = rf_model,
    importance = importance(rf_model),
    partial_dependence = list(
      overall = pred_data_overall,
      by_fish = pd_data_fish,
      by_time = pd_data_time
    ),
    plots = plots,
    dataset_name = dataset_name,
    sample_size_used = nrow(rf_data_sample),
    depth_range = depth_range
  ))
}

#' Create comparison analysis between positioning and space use results
create_comparison_analysis <- function(positioning_results, space_use_results, depth_variable) {
  
  # Combine overall predictions for comparison
  pos_overall <- positioning_results$partial_dependence$overall
  pos_overall$method <- "Positioning"
  
  space_overall <- space_use_results$partial_dependence$overall
  space_overall$method <- "Actual Track"
  
  combined_overall <- rbind(pos_overall, space_overall)
  
  # Create comparison plots
  plots <- list()
  
  # Overall comparison
  plots$overall_comparison <- ggplot(combined_overall, aes(x = !!sym(depth_variable), y = prob_presence, color = method)) +
    geom_line(size = 2, alpha = 0.8) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray", alpha = 0.7) +
    labs(
      title = "Habitat Selection Comparison: Positioning vs Actual Track",
      subtitle = "Overall patterns across all fish and time periods",
      x = paste(tools::toTitleCase(gsub("_", " ", depth_variable))),
      y = "Probability of Presence (Selection)",
      color = "Analysis Method"
    ) +
    scale_color_manual(values = c("Positioning" = "blue", "Actual Track" = "red")) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    )
  
  # Side-by-side comparison removed - using faceted plots instead
  
  # Difference plot
  # Interpolate to common depth grid for difference calculation
  common_depths <- seq(
    max(min(pos_overall[[depth_variable]]), min(space_overall[[depth_variable]])),
    min(max(pos_overall[[depth_variable]]), max(space_overall[[depth_variable]])),
    length.out = 100
  )
  
  pos_interp <- approx(pos_overall[[depth_variable]], pos_overall$prob_presence, xout = common_depths)$y
  space_interp <- approx(space_overall[[depth_variable]], space_overall$prob_presence, xout = common_depths)$y
  
  difference_data <- data.frame(
    depth = common_depths,
    difference = pos_interp - space_interp,
    abs_difference = abs(pos_interp - space_interp)
  )
  names(difference_data)[1] <- depth_variable
  
  plots$difference <- ggplot(difference_data, aes(x = !!sym(depth_variable), y = difference)) +
    geom_line(size = 1.5, color = "purple") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_ribbon(aes(ymin = pmin(difference, 0), ymax = 0), fill = "red", alpha = 0.3) +
    geom_ribbon(aes(ymin = 0, ymax = pmax(difference, 0)), fill = "blue", alpha = 0.3) +
    labs(
      title = "Selection Probability Difference: Positioning - Actual Track",
      subtitle = "Blue = Positioning higher | Red = Actual Track higher",
      x = paste(tools::toTitleCase(gsub("_", " ", depth_variable))),
      y = "Difference in Probability"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"))
  
  # Calculate comparison statistics
  stats <- list(
    mean_difference = mean(difference_data$difference),
    mean_abs_difference = mean(difference_data$abs_difference),
    max_difference = max(difference_data$abs_difference),
    correlation = cor(pos_interp, space_interp),
    positioning_range = range(pos_overall$prob_presence),
    space_use_range = range(space_overall$prob_presence)
  )
  
  return(list(
    plots = plots,
    combined_data = combined_overall,
    difference_data = difference_data,
    statistics = stats
  ))
}

# Helper functions for single dataset analysis
calculate_partial_dependence_fish_single <- function(model, depth_seq, fish_ids, time_periods,
                                                     depth_var, fish_var, time_var, data) {
  pd_data_fish <- data.frame()
  
  for(fish in fish_ids) {
    fish_predictions <- data.frame()
    
    for(time in time_periods) {
      # Create prediction data frame
      if (time_var %in% names(data)) {
        pred_data <- data.frame(depth_seq, fish, time)
        names(pred_data) <- c(depth_var, fish_var, time_var)
      } else {
        # No time variable - just depth and fish
        pred_data <- data.frame(depth_seq, fish)
        names(pred_data) <- c(depth_var, fish_var)
      }
      
      pred_data$prob_presence <- predict(model, pred_data, type = "prob")[,2]
      fish_predictions <- rbind(fish_predictions, pred_data)
    }
    
    fish_avg <- fish_predictions %>%
      group_by(!!sym(depth_var)) %>%
      summarise(prob_presence = mean(prob_presence), .groups = 'drop') %>%
      mutate(fish_id_factor = as.factor(fish))
    
    pd_data_fish <- rbind(pd_data_fish, fish_avg)
  }
  
  return(pd_data_fish)
}

calculate_partial_dependence_time_single <- function(model, depth_seq, fish_ids, time_periods,
                                                     depth_var, fish_var, time_var, data) {
  pd_data_time <- data.frame()
  
  # Skip time-based partial dependence if no time variable
  if (!time_var %in% names(data)) {
    return(pd_data_time)
  }
  
  for(time in time_periods) {
    time_predictions <- data.frame()
    
    for(fish in fish_ids) {
      pred_data <- data.frame(depth_seq, fish, time)
      names(pred_data) <- c(depth_var, fish_var, time_var)
      pred_data$prob_presence <- predict(model, pred_data, type = "prob")[,2]
      time_predictions <- rbind(time_predictions, pred_data)
    }
    
    time_avg <- time_predictions %>%
      group_by(!!sym(depth_var)) %>%
      summarise(prob_presence = mean(prob_presence), .groups = 'drop') %>%
      mutate(time_period_factor = as.factor(time))
    
    pd_data_time <- rbind(pd_data_time, time_avg)
  }
  
  return(pd_data_time)
}

calculate_overall_partial_dependence_single <- function(model, depth_seq, fish_ids, time_periods,
                                                       depth_var, fish_var, time_var, data) {
  all_predictions <- data.frame()
  
  for(fish in fish_ids) {
    for(time in time_periods) {
      # Create prediction data frame
      if (time_var %in% names(data)) {
        pred_data <- data.frame(depth_seq, fish, time)
        names(pred_data) <- c(depth_var, fish_var, time_var)
      } else {
        # No time variable - just depth and fish
        pred_data <- data.frame(depth_seq, fish)
        names(pred_data) <- c(depth_var, fish_var)
      }
      
      pred_data$prob_presence <- predict(model, pred_data, type = "prob")[,2]
      all_predictions <- rbind(all_predictions, pred_data)
    }
  }
  
  pred_data_overall <- all_predictions %>%
    group_by(!!sym(depth_var)) %>%
    summarise(prob_presence = mean(prob_presence), .groups = 'drop')
  
  return(pred_data_overall)
}

# Plotting functions for single datasets
create_single_overall_plot <- function(overall_data, depth_var, dataset_name) {
  ggplot(overall_data, aes(x = !!sym(depth_var), y = prob_presence)) +
    geom_line(size = 2, color = "darkblue") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", alpha = 0.7) +
    labs(
      title = paste(dataset_name, "Habitat Selection - Overall Pattern"),
      x = paste(tools::toTitleCase(gsub("_", " ", depth_var))),
      y = "Probability of Presence"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"))
}

create_single_fish_plot <- function(pd_data_fish, overall_data, depth_var, dataset_name) {
  ggplot() +
    geom_line(data = pd_data_fish, 
              aes(x = !!sym(depth_var), y = prob_presence, group = fish_id_factor), 
              alpha = 0.3, color = "lightblue", size = 0.8) +
    geom_line(data = overall_data, 
              aes(x = !!sym(depth_var), y = prob_presence), 
              color = "darkblue", size = 2) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", alpha = 0.7) +
    labs(
      title = paste(dataset_name, "- Individual Fish Patterns"),
      x = paste(tools::toTitleCase(gsub("_", " ", depth_var))),
      y = "Probability of Presence"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"))
}

create_single_time_plot <- function(pd_data_time, overall_data, depth_var, time_var, dataset_name) {
  ggplot() +
    geom_line(data = pd_data_time, 
              aes(x = !!sym(depth_var), y = prob_presence, group = time_period_factor, color = time_period_factor), 
              alpha = 0.7, size = 1) +
    geom_line(data = overall_data, 
              aes(x = !!sym(depth_var), y = prob_presence), 
              color = "black", size = 2) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", alpha = 0.7) +
    labs(
      title = paste(dataset_name, "- Time Period Patterns"),
      x = paste(tools::toTitleCase(gsub("_", " ", depth_var))),
      y = "Probability of Presence",
      color = tools::toTitleCase(gsub("_", " ", time_var))
    ) +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          legend.position = "bottom") +
    scale_color_viridis_d()
}

#' Print comparative habitat selection summary
#'
#' @param comp_results Output from analyze_comparative_habitat_selection
#' @export
print_comparative_summary <- function(comp_results) {
  
  cat("=== COMPARATIVE HABITAT SELECTION SUMMARY ===\n")
  
  # Analysis type
  cat("Analysis Type:", comp_results$metadata$analysis_type, "\n")
  cat("Datasets Analyzed:", paste(comp_results$metadata$datasets_analyzed, collapse = ", "), "\n\n")
  
  # Individual dataset summaries
  for (dataset in comp_results$metadata$datasets_analyzed) {
    if (dataset == "comparison") next
    
    results <- comp_results[[dataset]]
    cat("--- ", toupper(dataset), " RESULTS ---\n")
    cat("OOB Error Rate:", round(results$model$err.rate[nrow(results$model$err.rate), "OOB"] * 100, 2), "%\n")
    cat("Sample Size Used:", results$sample_size_used, "\n")
    cat("Depth Range:", round(results$depth_range, 2), "\n")
    
    # Top variable importance
    imp_sorted <- results$importance[order(results$importance[, "MeanDecreaseGini"], decreasing = TRUE), ]
    cat("Top Variable (MeanDecreaseGini):", rownames(imp_sorted)[1], "(", round(imp_sorted[1, "MeanDecreaseGini"], 3), ")\n")
    
    # Prediction range
    pred_range <- range(results$partial_dependence$overall$prob_presence)
    cat("Probability Range:", round(pred_range, 3), "\n")
    cat("Mean Probability:", round(mean(results$partial_dependence$overall$prob_presence), 3), "\n\n")
  }
  
  # Comparison statistics if available
  if ("comparison" %in% names(comp_results)) {
    cat("--- COMPARISON STATISTICS ---\n")
    stats <- comp_results$comparison$statistics
    cat("Mean Difference (Pos - Space):", round(stats$mean_difference, 4), "\n")
    cat("Mean Absolute Difference:", round(stats$mean_abs_difference, 4), "\n")
    cat("Maximum Difference:", round(stats$max_difference, 4), "\n")
    cat("Correlation:", round(stats$correlation, 4), "\n")
  }
  
  cat("\n=== SUMMARY COMPLETE ===\n")
}

#' Plot comparative habitat selection results
#'
#' @param comp_results Output from analyze_comparative_habitat_selection
#' @param plot_type Type of plots to show ("individual", "comparison", "all")
#' @export
plot_comparative_results <- function(comp_results, plot_type = "all") {
  
  if (plot_type %in% c("individual", "all")) {
    # Show individual dataset plots
    for (dataset in comp_results$metadata$datasets_analyzed) {
      if (dataset == "comparison") next
      
      cat("\nPlots for", dataset, "analysis:\n")
      
      if (!is.null(comp_results[[dataset]]$plots$overall)) {
        print(comp_results[[dataset]]$plots$overall)
      }
      if (!is.null(comp_results[[dataset]]$plots$fish)) {
        print(comp_results[[dataset]]$plots$fish)
      }
      if (!is.null(comp_results[[dataset]]$plots$time)) {
        print(comp_results[[dataset]]$plots$time)
      }
    }
  }
  
  if (plot_type %in% c("comparison", "all") && "comparison" %in% names(comp_results)) {
    cat("\nComparison plots:\n")
    
    # Create faceted comparison plot
    comparison_plot <- comp_results$comparison$plots$overall_comparison
    difference_plot <- comp_results$comparison$plots$difference
    
    # Combine plots using patchwork for side-by-side layout
    combined_comparison <- comparison_plot | difference_plot + 
      plot_layout(widths = c(5, 1)) +
      plot_annotation(
        title = "Comparative Habitat Selection Analysis",
        subtitle = "Positioning vs Actual Track",
        theme = theme(
          plot.title = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 14)
        )
      )
    
    print(combined_comparison)
  }
}
