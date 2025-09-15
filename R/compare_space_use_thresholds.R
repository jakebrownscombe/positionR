#' Compare Space Use Results Across Multiple Probability Thresholds
#'
#' Analyzes how space use areas change across different probability thresholds
#' and compares them to reference track-based space use estimates.
#'
#' @param multi_threshold_results Output from calculate_space_use with multiple prob_thresholds
#' @param reference_tracks Track data (e.g., fish_simulation$tracks) for overlap analysis
#' @param reference_space_use Track-based space use results for area comparison
#' @param methods Character vector of methods to compare (default: all available)
#' @param overlap_method Method for calculating track overlap (default: "intersection")
#' @param create_plots Whether to generate visualization plots (default: TRUE)
#' @param plot_by_fish Whether to include individual fish patterns in plots (default: TRUE)
#'
#' @return List containing comparison results and plots:
#' \itemize{
#'   \item area_comparison: Data frame with area ratios by threshold
#'   \item track_overlap: Data frame with track overlap analysis
#'   \item summary_stats: Summary statistics and optimal thresholds
#'   \item plots: List of ggplot objects for visualization
#' }
#'
#' @details
#' This function compares positioning-based space use areas (calculated with
#' different probability thresholds) against reference track-based space use.
#' 
#' Key analyses:
#' \itemize{
#'   \item Area ratio: positioning_area / track_area for each threshold
#'   \item Track overlap: percentage of positioning area containing track points
#'   \item Threshold sensitivity: how areas change across thresholds
#'   \item Optimal threshold: which threshold produces areas closest to track-based
#' }
#'
#' @examples
#' \dontrun{
#' # Calculate positioning results with multiple thresholds
#' multi_results <- calculate_space_use(
#'   track_data = position_points_coords,
#'   prob_column = "probability",
#'   prob_thresholds = c(0.1, 0.2, 0.3, 0.4, 0.5),
#'   methods = c("constrained_convex_hull", "grid_cell_count")
#' )
#' 
#' # Calculate reference track-based space use
#' track_results <- calculate_space_use(
#'   track_data = fish_simulation$tracks,
#'   methods = c("constrained_convex_hull", "grid_cell_count")
#' )
#' 
#' # Compare results
#' comparison <- compare_space_use_thresholds(
#'   multi_threshold_results = multi_results,
#'   reference_tracks = fish_simulation$tracks,
#'   reference_space_use = track_results
#' )
#' 
#' # View area ratio plots
#' print(comparison$plots$area_ratio_boxplot)
#' print(comparison$plots$area_ratio_by_fish)
#' }
#'
#' @export
compare_space_use_thresholds <- function(multi_threshold_results,
                                       reference_tracks,
                                       reference_space_use,
                                       methods = NULL,
                                       overlap_method = "intersection",
                                       create_plots = TRUE,
                                       plot_by_fish = TRUE) {
  
  # Load required packages
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' needed for data manipulation")
  }
  if (create_plots && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for plotting")
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' needed for spatial operations")
  }
  
  cat("=== SPACE USE THRESHOLD COMPARISON ===\n")
  
  # Validate inputs
  if (!is.list(multi_threshold_results) || 
      !("results_by_threshold" %in% names(multi_threshold_results))) {
    stop("multi_threshold_results must be output from calculate_space_use with multiple prob_thresholds")
  }
  
  if (!is.data.frame(reference_tracks)) {
    stop("reference_tracks must be a data frame with track data")
  }
  
  if (!is.list(reference_space_use) || 
      !("space_use_estimates" %in% names(reference_space_use))) {
    stop("reference_space_use must be output from calculate_space_use")
  }
  
  # Extract thresholds and methods
  thresholds <- multi_threshold_results$thresholds
  threshold_results <- multi_threshold_results$results_by_threshold
  
  cat("Analyzing", length(thresholds), "thresholds:", paste(thresholds, collapse = ", "), "\n")
  
  # Determine methods to analyze
  first_result <- threshold_results[[1]]
  
  # Get spatial data methods, but exclude outline/auxiliary data
  threshold_methods <- names(first_result$spatial_data)
  reference_methods <- names(reference_space_use$spatial_data)
  
  # Filter out outline and auxiliary data that are not actual space use methods
  threshold_methods <- threshold_methods[!grepl("_outline$", threshold_methods)]
  reference_methods <- reference_methods[!grepl("_outline$", reference_methods)]
  
  available_methods <- intersect(threshold_methods, reference_methods)
  
  if (is.null(methods)) {
    methods <- available_methods
  } else {
    methods <- intersect(methods, available_methods)
  }
  
  if (length(methods) == 0) {
    stop("No matching methods found between multi_threshold_results and reference_space_use")
  }
  
  cat("Analyzing methods:", paste(methods, collapse = ", "), "\n")
  
  # Initialize results
  area_comparison_list <- list()
  track_overlap_list <- list()
  
  # Process each threshold
  for (i in seq_along(thresholds)) {
    threshold <- thresholds[i]
    threshold_name <- names(threshold_results)[i]
    
    cat("\nProcessing threshold:", threshold, "\n")
    
    current_results <- threshold_results[[threshold_name]]
    
    # Check if we have valid results for this threshold
    if (is.null(current_results) || !("space_use_estimates" %in% names(current_results))) {
      cat("  No valid results for threshold", threshold, "- skipping\n")
      next
    }
    
    # Extract area data for this threshold
    threshold_areas <- current_results$space_use_estimates
    
    # Check if threshold areas exist
    if (is.null(threshold_areas) || nrow(threshold_areas) == 0) {
      cat("  No space use estimates for threshold", threshold, "- skipping\n")
      next
    }
    
    # Process each method
    for (method in methods) {
      cat("  Method:", method, "\n")
      
      # Get reference areas for this method
      method_col <- paste0(method, "_area_hectares")
      if (!method_col %in% names(reference_space_use$space_use_estimates)) {
        warning("Method column not found in reference data: ", method_col)
        next
      }
      
      # Merge positioning and reference areas
      positioning_areas <- threshold_areas %>%
        dplyr::select(fish_id, time_period_label, !!dplyr::sym(method_col)) %>%
        dplyr::rename(positioning_area_hectares = !!dplyr::sym(method_col))
      
      reference_areas <- reference_space_use$space_use_estimates %>%
        dplyr::select(fish_id, time_period_label, !!dplyr::sym(method_col)) %>%
        dplyr::rename(track_area_hectares = !!dplyr::sym(method_col))
      
      # Check if positioning areas exist
      if (nrow(positioning_areas) == 0) {
        cat("    No data found for threshold", threshold, "method", method, "- skipping\n")
        next
      }
      
      # Combine areas
      combined_areas <- positioning_areas %>%
        dplyr::left_join(reference_areas, by = c("fish_id", "time_period_label"))
      
      # Check if we have any matches after the join
      if (nrow(combined_areas) == 0 || all(is.na(combined_areas$track_area_hectares))) {
        cat("    No matching reference data for threshold", threshold, "method", method, "- skipping\n")
        next
      }
      
      # Add calculated columns
      combined_areas <- combined_areas %>%
        dplyr::mutate(
          threshold = threshold,
          method = method,
          positioning_area_m2 = positioning_area_hectares * 10000,
          track_area_m2 = track_area_hectares * 10000,
          area_ratio = track_area_hectares / positioning_area_hectares,
          area_difference_hectares = positioning_area_hectares - track_area_hectares,
          area_difference_percent = (area_difference_hectares / track_area_hectares) * 100
        )
      
      # Skip if no valid combinations found
      if (nrow(combined_areas) == 0) {
        cat("    No matching fish-time combinations found for threshold", threshold, "method", method, "- skipping\n")
        next
      }
      
      area_comparison_list[[paste(threshold, method, sep = "_")]] <- combined_areas
      
      # Calculate track overlap for this threshold and method
      if (method %in% names(current_results$spatial_data) && overlap_method != "none") {
        spatial_data <- current_results$spatial_data[[method]]
        
        if (!is.null(spatial_data) && nrow(spatial_data) > 0) {
          overlap_results <- calculate_track_overlap(
            spatial_data = spatial_data,
            reference_tracks = reference_tracks,
            method = overlap_method
          )
          
          # Only add threshold/method if we have results
          if (nrow(overlap_results) > 0) {
            overlap_results$threshold <- threshold
            overlap_results$method <- method
            track_overlap_list[[paste(threshold, method, sep = "_")]] <- overlap_results
          }
        }
      }
    }
  }
  
  # Combine results
  if (length(area_comparison_list) == 0) {
    stop("No valid area comparison results found. Check that your thresholds produce data and that fish_id/time_period_label columns match between datasets.")
  }
  
  area_comparison <- dplyr::bind_rows(area_comparison_list)
  track_overlap <- dplyr::bind_rows(track_overlap_list)
  
  # Calculate summary statistics
  summary_stats <- calculate_threshold_summary_stats(area_comparison, track_overlap, methods)
  
  # Create plots if requested
  plots <- list()
  if (create_plots) {
    cat("\nCreating visualization plots...\n")
    plots <- create_threshold_comparison_plots(
      area_comparison, 
      track_overlap, 
      plot_by_fish = plot_by_fish
    )
  }
  
  cat("\n=== THRESHOLD COMPARISON COMPLETE ===\n")
  cat("Thresholds analyzed:", length(thresholds), "\n")
  cat("Methods analyzed:", length(methods), "\n")
  cat("Fish-time combinations:", nrow(area_comparison) / (length(thresholds) * length(methods)), "\n")
  
  return(list(
    area_comparison = area_comparison,
    track_overlap = track_overlap,
    summary_stats = summary_stats,
    plots = plots,
    parameters = list(
      thresholds = thresholds,
      methods = methods,
      overlap_method = overlap_method,
      prob_column = multi_threshold_results$prob_column
    )
  ))
}

#' Calculate track overlap for spatial data
calculate_track_overlap <- function(spatial_data, reference_tracks, method = "intersection") {
  
  if (method == "none") {
    # Skip overlap calculation entirely
    return(data.frame(
      fish_id = numeric(0),
      time_period_label = character(0),
      positioning_area_m2 = numeric(0),
      track_points_total = numeric(0),
      track_points_in_area = numeric(0),
      overlap_percentage = numeric(0)
    ))
  }
  
  # Standardize track data column names
  tracks <- reference_tracks
  if ("path_id" %in% names(tracks) && !("fish_id" %in% names(tracks))) {
    tracks$fish_id <- tracks$path_id
  }
  
  # Group spatial data by fish and time
  overlap_results <- spatial_data %>%
    dplyr::group_by(fish_id, time_period_label) %>%
    dplyr::group_modify(~ {
      current_spatial <- .x
      current_fish <- .y$fish_id
      current_time <- .y$time_period_label
      
      # Filter tracks to match fish and time
      fish_tracks <- tracks %>%
        dplyr::filter(fish_id == current_fish)
      
      # Add time period to tracks if needed
      if ("datetime" %in% names(fish_tracks) && !("time_period_label" %in% names(fish_tracks))) {
        fish_tracks$time_period_label <- as.character(as.Date(fish_tracks$datetime))
      }
      
      # Filter by time period
      if ("time_period_label" %in% names(fish_tracks)) {
        fish_tracks <- fish_tracks %>%
          dplyr::filter(time_period_label == current_time)
      }
      
      if (nrow(fish_tracks) == 0) {
        return(data.frame(
          positioning_area_m2 = 0,
          track_points_total = 0,
          track_points_in_area = 0,
          overlap_percentage = 0
        ))
      }
      
      # Fast bounding box method - similar to calculate_contour_coverage efficiency
      tryCatch({
        # Get bounds of positioning area
        x_min <- min(current_spatial$x, na.rm = TRUE)
        x_max <- max(current_spatial$x, na.rm = TRUE)
        y_min <- min(current_spatial$y, na.rm = TRUE)
        y_max <- max(current_spatial$y, na.rm = TRUE)
        
        # Simple point-in-bounds check (much faster than sf operations)
        points_in_area <- (fish_tracks$x >= x_min & fish_tracks$x <= x_max & 
                          fish_tracks$y >= y_min & fish_tracks$y <= y_max)
        n_points_in_area <- sum(points_in_area, na.rm = TRUE)
        
        # Estimate area from spatial data points (simple approximation)
        area_width <- x_max - x_min
        area_height <- y_max - y_min
        positioning_area_m2 <- area_width * area_height
        
        data.frame(
          positioning_area_m2 = positioning_area_m2,
          track_points_total = nrow(fish_tracks),
          track_points_in_area = n_points_in_area,
          overlap_percentage = (n_points_in_area / nrow(fish_tracks)) * 100
        )
        
      }, error = function(e) {
        warning("Error calculating overlap for fish ", current_fish, ", time ", current_time, ": ", e$message)
        data.frame(
          positioning_area_m2 = 0,
          track_points_total = nrow(fish_tracks),
          track_points_in_area = 0,
          overlap_percentage = 0
        )
      })
    }) %>%
    dplyr::ungroup()
  
  return(overlap_results)
}

#' Calculate summary statistics for threshold comparison
calculate_threshold_summary_stats <- function(area_comparison, track_overlap, methods) {
  
  # Handle empty area_comparison data
  if (nrow(area_comparison) == 0) {
    cat("Warning: No area comparison data available - returning empty summary statistics\n")
    return(list(
      optimal_thresholds = data.frame(
        method = character(0),
        optimal_threshold = numeric(0),
        min_area_ratio_diff = numeric(0)
      ),
      mean_area_ratios = data.frame(
        threshold = numeric(0),
        method = character(0),
        mean_area_ratio = numeric(0),
        median_area_ratio = numeric(0),
        sd_area_ratio = numeric(0),
        n_observations = numeric(0)
      ),
      mean_overlap_percentages = data.frame(
        threshold = numeric(0),
        method = character(0),
        mean_overlap_percentage = numeric(0),
        median_overlap_percentage = numeric(0),
        sd_overlap_percentage = numeric(0),
        n_observations = numeric(0)
      )
    ))
  }
  
  # Calculate optimal thresholds (closest area ratio to 1.0)
  optimal_thresholds <- area_comparison %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(
      optimal_threshold = threshold[which.min(abs(area_ratio - 1.0))],
      min_area_ratio_diff = min(abs(area_ratio - 1.0)),
      .groups = 'drop'
    )
  
  # Calculate mean area ratios by threshold and method
  mean_area_ratios <- area_comparison %>%
    dplyr::group_by(threshold, method) %>%
    dplyr::summarise(
      mean_area_ratio = mean(area_ratio, na.rm = TRUE),
      median_area_ratio = median(area_ratio, na.rm = TRUE),
      sd_area_ratio = sd(area_ratio, na.rm = TRUE),
      n_observations = dplyr::n(),
      .groups = 'drop'
    )
  
  # Calculate mean overlap percentages
  if (nrow(track_overlap) == 0) {
    mean_overlap <- data.frame(
      threshold = numeric(0),
      method = character(0),
      mean_overlap_percentage = numeric(0),
      median_overlap_percentage = numeric(0),
      sd_overlap_percentage = numeric(0),
      n_observations = numeric(0)
    )
  } else {
    mean_overlap <- track_overlap %>%
      dplyr::group_by(threshold, method) %>%
      dplyr::summarise(
        mean_overlap_percentage = mean(overlap_percentage, na.rm = TRUE),
        median_overlap_percentage = median(overlap_percentage, na.rm = TRUE),
        sd_overlap_percentage = sd(overlap_percentage, na.rm = TRUE),
        n_observations = dplyr::n(),
        .groups = 'drop'
      )
  }
  
  return(list(
    optimal_thresholds = optimal_thresholds,
    mean_area_ratios = mean_area_ratios,
    mean_overlap_percentages = mean_overlap
  ))
}

#' Create visualization plots for threshold comparison
create_threshold_comparison_plots <- function(area_comparison, track_overlap, plot_by_fish = TRUE) {
  
  plots <- list()
  
  # Handle empty area_comparison data
  if (nrow(area_comparison) == 0) {
    cat("Warning: No area comparison data available - cannot create plots\n")
    return(plots)
  }
  
  # Area ratio boxplot (Threshold vs Size Ratio)
  plots$area_ratio_boxplot <- ggplot2::ggplot(area_comparison, 
                                             ggplot2::aes(x = factor(threshold), y = area_ratio)) +
    ggplot2::geom_boxplot(alpha = 0.7) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
    ggplot2::facet_wrap(~ method, scales = "free_y") +
    ggplot2::labs(
      title = "Area Ratio vs Probability Threshold",
      subtitle = "Track-based Area / Positioning Area (Red line = equal areas)",
      x = "Probability Threshold",
      y = "Area Ratio (Track-based / Positioning)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  # Area ratio by individual fish (if requested)
  if (plot_by_fish) {
    plots$area_ratio_by_fish <- ggplot2::ggplot(area_comparison, 
                                               ggplot2::aes(x = threshold, y = area_ratio, 
                                                           color = factor(fish_id), 
                                                           group = fish_id)) +
      ggplot2::geom_line(alpha = 0.8, size = 1) +
      ggplot2::geom_point(alpha = 0.8, size = 2) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
      ggplot2::facet_wrap(~ method, scales = "free_y") +
      ggplot2::labs(
        title = "Area Ratio by Individual Fish",
        subtitle = "How each fish responds to probability thresholds",
        x = "Probability Threshold",
        y = "Area Ratio (Track-based / Positioning)",
        color = "Fish ID"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::scale_color_viridis_d()
  }
  
  # Track overlap boxplot
  if (nrow(track_overlap) > 0) {
    plots$track_overlap_boxplot <- ggplot2::ggplot(track_overlap, 
                                                  ggplot2::aes(x = factor(threshold), y = overlap_percentage)) +
      ggplot2::geom_boxplot(alpha = 0.7) +
      ggplot2::facet_wrap(~ method, scales = "free_y") +
      ggplot2::labs(
        title = "Track Overlap vs Probability Threshold",
        subtitle = "Percentage of positioning area containing actual track points",
        x = "Probability Threshold",
        y = "Track Overlap Percentage (%)"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    
    # Track overlap by individual fish (if requested)
    if (plot_by_fish) {
      plots$track_overlap_by_fish <- ggplot2::ggplot(track_overlap, 
                                                    ggplot2::aes(x = threshold, y = overlap_percentage, 
                                                                color = factor(fish_id), 
                                                                group = fish_id)) +
        ggplot2::geom_line(alpha = 0.8, size = 1) +
        ggplot2::geom_point(alpha = 0.8, size = 2) +
        ggplot2::facet_wrap(~ method, scales = "free_y") +
        ggplot2::labs(
          title = "Track Overlap by Individual Fish",
          subtitle = "How track overlap changes with probability thresholds",
          x = "Probability Threshold",
          y = "Track Overlap Percentage (%)",
          color = "Fish ID"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::scale_color_viridis_d()
    }
  }
  
  return(plots)
}