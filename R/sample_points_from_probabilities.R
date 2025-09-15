#' Sample points from position probability surfaces
#'
#' Generates random points based on probability surfaces from fish positioning
#' results. Points are sampled with probability proportional to the specified
#' probability column, allowing for Monte Carlo sampling from position estimates.
#' Can process multiple fish and time periods simultaneously.
#'
#' @param positioning_results A list returned by \code{\link{calculate_fish_positions}}
#'   containing position probabilities and summary statistics.
#' @param prob_column Character. Name of the probability column to use for sampling.
#'   Default is "weighted_mean_DE_normalized_scaled". Available options include:
#'   \itemize{
#'     \item weighted_mean_DE_normalized - Raw detection probability
#'     \item weighted_mean_DE_normalized_scaled - Scaled detection probability (0-1)
#'     \item non_det_DE_normalized - Non-detection probability
#'     \item non_det_DE_normalized_scaled - Scaled non-detection probability
#'     \item integrated_prob - Integrated position probability
#'   }
#' @param n_points Integer. Number of points to sample per fish-time combination. 
#'   Default is 100.
#' @param fish_select Integer, character vector, or NULL. Fish ID(s) to sample from. 
#'   If NULL, samples from all fish. Default is NULL.
#' @param time_select Numeric vector, character vector, POSIXct vector, or NULL. 
#'   Time period(s) to sample from. If NULL, samples from all time periods. Default is NULL.
#' @param min_prob_threshold Numeric. Minimum probability threshold (0-1).
#'   Only cells with probability above this threshold are eligible for sampling.
#'   Default is 0.001 to exclude zero-probability cells.
#' @param max_prob_threshold Numeric. Maximum probability threshold (0-1).
#'   Only cells with probability below this threshold are eligible for sampling.
#'   Default is 1.0 (no upper limit). Set to lower values (e.g., 0.05) to exclude
#'   high-probability cells and focus on uncertainty regions.
#' @param seed Integer. Random seed for reproducible sampling. Default is NULL.
#' @param by_group Logical. If TRUE, samples n_points for each fish-time combination.
#'   If FALSE, samples n_points total distributed across all combinations.
#'   Default is TRUE.
#' @param crs Coordinate reference system for the output sf object. Can be:
#'   \itemize{
#'     \item NULL (default) - uses CRS from positioning_results
#'     \item Numeric EPSG code (e.g., 4326 for WGS84, 32618 for UTM Zone 18N)
#'     \item Character proj4 string
#'     \item An sf/sfc object from which to extract CRS
#'   }
#'
#' @return An sf object containing the sampled points with columns:
#'   \item{fish_id}{Fish identifier}
#'   \item{time_period}{Time period identifier}
#'   \item{time_period_posix}{POSIXct datetime for the time period}
#'   \item{time_period_label}{Human-readable time period label}
#'   \item{x}{X coordinates}
#'   \item{y}{Y coordinates}
#'   \item{probability}{The probability value used for sampling}
#'   \item{sample_id}{Sequential sample identifier}
#'   \item{group_id}{Unique identifier for each fish-time combination}
#'   \item{geometry}{sf point geometry}
#'
#' @details
#' This function performs weighted random sampling where each spatial cell has
#' a probability of being selected proportional to its probability value in
#' the specified column. This is useful for:
#' \itemize{
#'   \item Monte Carlo analysis of position uncertainty
#'   \item Generating representative sample locations for further analysis
#'   \item Creating random tracks based on position probability surfaces
#'   \item Uncertainty propagation in downstream analyses
#' }
#'
#' The sampling process:
#' \enumerate{
#'   \item Filters data by fish_select and time_select if specified
#'   \item Removes cells below min_prob_threshold
#'   \item Normalizes probabilities to sum to 1 within each time period
#'   \item Performs weighted sampling without replacement (or with if n_points > available cells)
#' }
#'
#' @examples
#' \dontrun{
#' # Calculate fish positions first
#' positioning_results <- calculate_fish_positions(station_detections, distances, stations)
#'
#' # Sample 100 points per fish-time combination for all data
#' sampled_points <- sample_points_from_probabilities(
#'   positioning_results,
#'   prob_column = "weighted_mean_DE_normalized_scaled",
#'   n_points = 100,
#'   seed = 123
#' )
#'
#' # Sample from multiple specific fish and time periods
#' multi_samples <- sample_points_from_probabilities(
#'   positioning_results,
#'   prob_column = "integrated_prob",
#'   n_points = 50,
#'   fish_select = c(1, 2, 3),
#'   time_select = c("2025-07-15", "2025-07-16", "2025-07-17"),
#'   by_group = TRUE  # 50 points per fish-time combination
#' )
#'
#' # Sample 1000 points total distributed across all groups
#' distributed_samples <- sample_points_from_probabilities(
#'   positioning_results,
#'   n_points = 1000,
#'   by_group = FALSE  # Distribute 1000 points across all combinations
#' )
#'
#' # Sample only from low-probability uncertainty regions
#' uncertainty_samples <- sample_points_from_probabilities(
#'   positioning_results,
#'   prob_column = "integrated_prob",
#'   n_points = 200,
#'   min_prob_threshold = 0.001,
#'   max_prob_threshold = 0.05,  # Exclude high-probability areas
#'   seed = 456
#' )
#'
#' # Sample with specific CRS (UTM Zone 18N)
#' utm_samples <- sample_points_from_probabilities(
#'   positioning_results,
#'   n_points = 100,
#'   crs = 32618  # EPSG code for UTM Zone 18N
#' )
#'
#' # Plot sampled points colored by time
#' library(ggplot2)
#' ggplot() +
#'   geom_sf(data = multi_samples, aes(color = time_period_posix)) +
#'   scale_color_viridis_c() +
#'   facet_wrap(~fish_id) +
#'   theme_minimal() +
#'   labs(title = "Sampled Position Points by Fish and Time")
#'
#' # Temporal centroid analysis
#' library(dplyr)
#' temporal_centroids <- multi_samples %>%
#'   group_by(fish_id, time_period_posix, time_period_label) %>%
#'   summarise(
#'     n_samples = n(),
#'     mean_prob = mean(probability),
#'     centroid = sf::st_union(geometry) %>% sf::st_centroid(),
#'     .groups = 'drop'
#'   ) %>%
#'   sf::st_as_sf()
#'
#' # Create animated track from sampled points
#' track_samples <- multi_samples %>%
#'   arrange(fish_id, time_period_posix) %>%
#'   group_by(fish_id) %>%
#'   summarise(
#'     track = sf::st_cast(sf::st_union(geometry), "LINESTRING")
#'   )
#' }
#'
#' @seealso \code{\link{calculate_fish_positions}}, \code{\link{plot_fish_positions}}
#'
#' @export
sample_points_from_probabilities <- function(positioning_results,
                                            prob_column = "weighted_mean_DE_normalized_scaled",
                                            n_points = 100,
                                            fish_select = NULL,
                                            time_select = NULL,
                                            min_prob_threshold = 0.001,
                                            max_prob_threshold = 1.0,
                                            seed = NULL,
                                            by_group = TRUE,
                                            crs = NULL) {
  
  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)
  
  # Extract position probabilities
  position_data <- positioning_results$position_probabilities %>%
    as.data.frame()
  
  # Validate prob_column exists
  if (!prob_column %in% names(position_data)) {
    stop("Column '", prob_column, "' not found in position probabilities. ",
         "Available columns: ", paste(names(position_data), collapse = ", "))
  }
  
  cat("Using probability column:", prob_column, "\n")
  cat("Original data:", nrow(position_data), "cells\n")
  
  # Filter by fish if specified
  if (!is.null(fish_select)) {
    # Check if all fish_select values exist
    missing_fish <- setdiff(fish_select, unique(position_data$fish_id))
    if (length(missing_fish) > 0) {
      stop("Fish ID(s) not found: ", paste(missing_fish, collapse = ", "), 
           ". Available fish IDs: ", paste(unique(position_data$fish_id), collapse = ", "))
    }
    position_data <- position_data %>%
      dplyr::filter(fish_id %in% fish_select)
    cat("Filtered to fish:", paste(fish_select, collapse = ", "), 
        "-", nrow(position_data), "cells\n")
  }
  
  # Filter by time if specified
  if (!is.null(time_select)) {
    # Handle different time selection formats
    if (is.numeric(time_select)) {
      # Assume time_period numeric value
      position_data <- position_data %>%
        dplyr::filter(time_period %in% time_select)
    } else if (is.character(time_select)) {
      # Assume time_period_label
      position_data <- position_data %>%
        dplyr::filter(time_period_label %in% time_select)
    } else if (inherits(time_select, "POSIXt")) {
      # Convert to character and match
      time_chars <- as.character(as.Date(time_select))
      position_data <- position_data %>%
        dplyr::filter(time_period_label %in% time_chars)
    }
    cat("Filtered to time periods:", length(unique(position_data$time_period)), 
        "-", nrow(position_data), "cells\n")
  }
  
  # Filter by probability thresholds
  position_data <- position_data %>%
    dplyr::filter(.data[[prob_column]] > min_prob_threshold, 
                  .data[[prob_column]] < max_prob_threshold,
                  !is.na(.data[[prob_column]]))
  
  cat("After probability thresholds (", min_prob_threshold, "< prob <", max_prob_threshold, "):", 
      nrow(position_data), "cells\n")
  
  if (nrow(position_data) == 0) {
    stop("No cells remain after filtering. Try adjusting min_prob_threshold (", min_prob_threshold, 
         ") or max_prob_threshold (", max_prob_threshold, ") or check your selection criteria.")
  }
  
  # Get unique fish-time combinations
  fish_time_groups <- position_data %>%
    dplyr::distinct(fish_id, time_period, time_period_posix, time_period_label) %>%
    dplyr::arrange(fish_id, time_period)
  
  cat("\nProcessing", nrow(fish_time_groups), "fish-time combinations\n")
  
  # Initialize results list
  sampled_list <- list()
  
  # Process each fish-time combination
  for (i in 1:nrow(fish_time_groups)) {
    current_fish <- fish_time_groups$fish_id[i]
    current_time <- fish_time_groups$time_period[i]
    current_time_posix <- fish_time_groups$time_period_posix[i]
    current_time_label <- fish_time_groups$time_period_label[i]
    
    # Get data for this combination
    group_data <- position_data %>%
      dplyr::filter(fish_id == current_fish, 
                    time_period == current_time)
    
    if (nrow(group_data) == 0) next
    
    # Determine number of points to sample for this group
    if (by_group) {
      group_n_points <- n_points
    } else {
      # Distribute points proportionally if sampling total
      total_prob_mass <- sum(group_data[[prob_column]])
      overall_prob_mass <- sum(position_data[[prob_column]])
      group_n_points <- round(n_points * total_prob_mass / overall_prob_mass)
      group_n_points <- max(1, group_n_points)  # At least 1 point per group
    }
    
    # Check if we have enough cells for sampling without replacement
    if (group_n_points > nrow(group_data)) {
      replacement <- TRUE
    } else {
      replacement <- FALSE
    }
    
    # Prepare sampling weights
    sampling_weights <- group_data[[prob_column]]
    
    # Normalize weights to sum to 1 (probability distribution)
    sampling_weights <- sampling_weights / sum(sampling_weights)
    
    # Sample indices based on probabilities
    sampled_indices <- sample(
      1:nrow(group_data), 
      size = group_n_points, 
      replace = replacement,
      prob = sampling_weights
    )
    
    # Extract sampled points
    sampled_data <- group_data[sampled_indices, ] %>%
      dplyr::select(fish_id, time_period, time_period_posix, time_period_label, 
                    x, y, probability = !!rlang::sym(prob_column)) %>%
      dplyr::mutate(
        sample_id = 1:group_n_points,
        group_id = paste(current_fish, current_time, sep = "_")
      )
    
    sampled_list[[i]] <- sampled_data
    
    cat("  Fish", current_fish, "- Time", current_time_label, 
        ": sampled", group_n_points, "points from", nrow(group_data), "cells\n")
  }
  
  # Combine all sampled points
  all_sampled <- dplyr::bind_rows(sampled_list)
  
  # Re-number sample_id across all groups if needed
  if (!by_group) {
    all_sampled$sample_id <- 1:nrow(all_sampled)
  }
  
  # Determine CRS to use
  if (!is.null(crs)) {
    # User specified CRS
    use_crs <- sf::st_crs(crs)
  } else {
    # Use CRS from positioning results
    use_crs <- sf::st_crs(positioning_results$position_probabilities)
  }
  
  # Convert to sf object
  sampled_sf <- sf::st_as_sf(
    all_sampled, 
    coords = c("x", "y"), 
    crs = use_crs
  )
  
  # Add attributes for reference
  attr(sampled_sf, "prob_column") <- prob_column
  attr(sampled_sf, "n_points") <- n_points
  attr(sampled_sf, "min_prob_threshold") <- min_prob_threshold
  attr(sampled_sf, "max_prob_threshold") <- max_prob_threshold
  attr(sampled_sf, "by_group") <- by_group
  
  # Summary
  cat("\n=== SAMPLING SUMMARY ===\n")
  cat("Total points sampled:", nrow(sampled_sf), "\n")
  cat("Probability column:", prob_column, "\n")
  cat("Fish ID(s):", paste(unique(sampled_sf$fish_id), collapse = ", "), "\n")
  cat("Time periods:", length(unique(sampled_sf$time_period)), "\n")
  cat("Fish-time combinations:", length(unique(sampled_sf$group_id)), "\n")
  cat("Sampled probability range:", round(min(sampled_sf$probability), 4), 
      "to", round(max(sampled_sf$probability), 4), "\n")
  cat("Mean probability:", round(mean(sampled_sf$probability), 4), "\n")
  
  return(sampled_sf)
}