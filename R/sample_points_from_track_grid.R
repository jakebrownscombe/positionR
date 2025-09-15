#' Sample points from track data grid cell usage
#'
#' Generates random points based on grid cell usage from track data.
#' Points are sampled with probability proportional to the number of track points
#' in each grid cell, providing a usage-weighted spatial sampling approach.
#'
#' @param track_data Data frame with fish tracks containing columns: fish_id/path_id, 
#'   datetime, x, y. Can be from fish_simulation$tracks or similar track data.
#' @param grid_resolution Numeric. Grid cell size in map units (typically meters).
#'   Default is 100 meters.
#' @param n_points Integer. Number of points to sample per fish-time combination.
#'   Default is 100.
#' @param by_fish Logical. Whether to group by fish_id. Default is TRUE.
#' @param by_time_period Logical. Whether to group by time periods. Default is TRUE.
#' @param time_aggregation Character. How to aggregate time periods: "hour", "day", 
#'   "month", or "none". Default is "day".
#' @param fish_select Integer, character vector, or NULL. Fish ID(s) to sample from.
#'   If NULL, samples from all fish. Default is NULL.
#' @param time_select Character vector, POSIXct vector, or NULL. Time period(s) 
#'   to sample from. If NULL, samples from all time periods. Default is NULL.
#' @param min_count_threshold Numeric. Minimum count threshold.
#'   Only cells with count above this threshold are eligible for sampling.
#'   Default is 1 to exclude empty cells.
#' @param max_count_threshold Numeric. Maximum count threshold.
#'   Only cells with count below this threshold are eligible for sampling.
#'   Default is Inf (no upper limit).
#' @param seed Integer. Random seed for reproducible sampling. Default is NULL.
#' @param by_group Logical. If TRUE, samples n_points for each fish-time combination.
#'   If FALSE, samples n_points total distributed across all combinations.
#'   Default is TRUE.
#' @param crs Coordinate reference system for the output sf object. Can be:
#'   \itemize{
#'     \item NULL (default) - attempts to detect from input data or uses WGS84
#'     \item Numeric EPSG code (e.g., 4326 for WGS84, 32618 for UTM Zone 18N)
#'     \item Character proj4 string
#'     \item An sf/sfc object from which to extract CRS
#'   }
#' @param reference_raster Optional raster object to use for defining grid cells.
#'   If provided, uses actual raster cell boundaries instead of arbitrary grid.
#'
#' @return An sf object containing the sampled points with columns:
#'   \item{fish_id}{Fish identifier}
#'   \item{time_period_label}{Human-readable time period label (if by_time_period = TRUE)}
#'   \item{time_period_posix}{POSIXct datetime for the time period (if available)}
#'   \item{x}{X coordinates (grid cell centers)}
#'   \item{y}{Y coordinates (grid cell centers)}
#'   \item{count}{The count value used for sampling weights}
#'   \item{sample_id}{Sequential sample identifier}
#'   \item{group_id}{Unique identifier for each fish-time combination}
#'   \item{geometry}{sf point geometry}
#'
#' @details
#' This function creates a grid overlay on track data, counts the number of track
#' points in each grid cell, then samples points with probability proportional
#' to usage intensity. The process:
#' \enumerate{
#'   \item Aggregates track data by time periods if requested
#'   \item Creates regular grid cells based on grid_resolution
#'   \item Counts track points in each cell for each fish-time combination
#'   \item Samples points weighted by track point counts
#'   \item Returns usage-weighted spatial sampling
#' }
#'
#' This is useful for:
#' \itemize{
#'   \item Creating habitat models weighted by usage intensity
#'   \item Sampling environmental variables proportionate to space use
#'   \item Generating representative locations for resource selection analysis
#'   \item Monte Carlo analysis of habitat preferences
#'   \item Comparing used vs available habitat
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage-weighted sampling
#' usage_points <- sample_points_from_track_grid(
#'   track_data = fish_simulation$tracks,
#'   grid_resolution = 100,
#'   n_points = 500,
#'   seed = 123
#' )
#'
#' # Sample from specific fish and time periods
#' selected_usage <- sample_points_from_track_grid(
#'   track_data = fish_simulation$tracks,
#'   grid_resolution = 50,  # Finer resolution
#'   n_points = 200,
#'   fish_select = c(1, 2, 3),
#'   time_select = c("2025-07-15", "2025-07-16"),
#'   min_count_threshold = 3  # Only heavily used cells
#' )
#'
#' # Sample without time grouping
#' overall_usage <- sample_points_from_track_grid(
#'   track_data = fish_simulation$tracks,
#'   by_time_period = FALSE,
#'   n_points = 1000
#' )
#'
#' # Use with reference raster for consistent grid
#' raster_usage <- sample_points_from_track_grid(
#'   track_data = fish_simulation$tracks,
#'   reference_raster = depth_raster,
#'   n_points = 300,
#'   crs = 32617
#' )
#'
#' # Plot usage-weighted sampling
#' library(ggplot2)
#' ggplot() +
#'   geom_sf(data = usage_points, aes(color = count, size = count), 
#'           alpha = 0.6) +
#'   scale_color_viridis_c(name = "Track\nPoints", trans = "sqrt") +
#'   scale_size_continuous(name = "Track\nPoints", range = c(0.5, 3), trans = "sqrt") +
#'   facet_wrap(~fish_id + time_period_label) +
#'   theme_minimal() +
#'   labs(title = "Usage-Weighted Spatial Sampling")
#'
#' # Compare usage intensity across fish
#' usage_summary <- usage_points %>%
#'   st_drop_geometry() %>%
#'   group_by(fish_id) %>%
#'   summarise(
#'     mean_usage = mean(count),
#'     max_usage = max(count),
#'     total_samples = n(),
#'     .groups = 'drop'
#'   )
#' }
#'
#' @seealso \code{\link{calculate_space_use}}, \code{\link{sample_points_from_probabilities}}
#'
#' @export
sample_points_from_track_grid <- function(track_data,
                                        grid_resolution = 100,
                                        n_points = 100,
                                        by_fish = TRUE,
                                        by_time_period = TRUE,
                                        time_aggregation = "day",
                                        fish_select = NULL,
                                        time_select = NULL,
                                        min_count_threshold = 1,
                                        max_count_threshold = Inf,
                                        seed = NULL,
                                        by_group = TRUE,
                                        crs = NULL,
                                        reference_raster = NULL) {
  
  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)
  
  cat("=== TRACK GRID SAMPLING ===\n")
  cat("Grid resolution:", grid_resolution, "meters\n")
  cat("Original track data:", nrow(track_data), "points\n")
  
  # Standardize column names
  if ("path_id" %in% names(track_data) && !"fish_id" %in% names(track_data)) {
    track_data$fish_id <- track_data$path_id
  }
  
  # Check required columns
  required_cols <- c("fish_id", "datetime", "x", "y")
  missing_cols <- setdiff(required_cols, names(track_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in track_data: ", paste(missing_cols, collapse = ", "))
  }
  
  # Handle timezone standardization
  if (inherits(track_data$datetime, "POSIXt")) {
    original_tz <- attr(track_data$datetime, "tzone")
    if (is.null(original_tz) || original_tz == "") original_tz <- "UTC"
    track_data$datetime <- as.POSIXct(track_data$datetime, tz = "UTC")
  }
  
  # Add time period aggregation if requested
  if (by_time_period && time_aggregation != "none") {
    if (time_aggregation == "hour") {
      track_data$time_period <- as.POSIXct(format(track_data$datetime, "%Y-%m-%d %H:00:00"), tz = "UTC")
      track_data$time_period_label <- format(track_data$time_period, "%Y-%m-%d %H:00")
    } else if (time_aggregation == "day") {
      track_data$time_period <- as.POSIXct(as.Date(track_data$datetime), tz = "UTC")
      track_data$time_period_label <- as.character(as.Date(track_data$datetime))
    } else if (time_aggregation == "month") {
      track_data$time_period <- as.POSIXct(paste0(format(track_data$datetime, "%Y-%m"), "-01"), tz = "UTC")
      track_data$time_period_label <- format(track_data$time_period, "%Y-%m")
    }
    track_data$time_period_posix <- track_data$time_period
  }
  
  # Filter by fish if specified
  if (!is.null(fish_select)) {
    missing_fish <- setdiff(fish_select, unique(track_data$fish_id))
    if (length(missing_fish) > 0) {
      stop("Fish ID(s) not found: ", paste(missing_fish, collapse = ", "), 
           ". Available fish IDs: ", paste(unique(track_data$fish_id), collapse = ", "))
    }
    track_data <- track_data %>%
      dplyr::filter(fish_id %in% fish_select)
    cat("Filtered to fish:", paste(fish_select, collapse = ", "), 
        "-", nrow(track_data), "points\n")
  }
  
  # Filter by time if specified
  if (!is.null(time_select) && by_time_period) {
    if (is.character(time_select) && "time_period_label" %in% names(track_data)) {
      track_data <- track_data %>%
        dplyr::filter(time_period_label %in% time_select)
    } else if (inherits(time_select, "POSIXt") && "time_period_posix" %in% names(track_data)) {
      time_chars <- as.character(as.Date(time_select))
      track_data <- track_data %>%
        dplyr::filter(as.character(as.Date(time_period_posix)) %in% time_chars)
    }
    cat("Filtered to time periods:", length(unique(track_data$time_period_label)), 
        "-", nrow(track_data), "points\n")
  }
  
  # Create grid cells and count track points
  cat("Creating grid and counting track points...\n")
  
  if (!is.null(reference_raster)) {
    # Use reference raster for grid definition
    if (!requireNamespace("raster", quietly = TRUE)) {
      stop("raster package needed when reference_raster is provided")
    }
    
    # Get cell indices for track points
    track_coords <- cbind(track_data$x, track_data$y)
    cell_indices <- raster::cellFromXY(reference_raster, track_coords)
    
    # Remove points outside raster
    valid_mask <- !is.na(cell_indices)
    track_data <- track_data[valid_mask, ]
    cell_indices <- cell_indices[valid_mask]
    
    if (nrow(track_data) == 0) {
      stop("No track points fall within the reference raster boundaries")
    }
    
    # Get cell center coordinates
    cell_centers <- raster::xyFromCell(reference_raster, cell_indices)
    track_data$grid_x <- cell_centers[,1]
    track_data$grid_y <- cell_centers[,2]
    track_data$cell_id <- cell_indices
    
  } else {
    # Create regular grid
    track_data$grid_x <- floor(track_data$x / grid_resolution) * grid_resolution + grid_resolution/2
    track_data$grid_y <- floor(track_data$y / grid_resolution) * grid_resolution + grid_resolution/2
    track_data$cell_id <- paste(track_data$grid_x, track_data$grid_y, sep = "_")
  }
  
  # Group data for counting
  if (by_fish && by_time_period && "time_period_label" %in% names(track_data)) {
    group_vars <- c("fish_id", "time_period_label", "time_period_posix", "grid_x", "grid_y", "cell_id")
  } else if (by_fish) {
    group_vars <- c("fish_id", "grid_x", "grid_y", "cell_id")
  } else if (by_time_period && "time_period_label" %in% names(track_data)) {
    group_vars <- c("time_period_label", "time_period_posix", "grid_x", "grid_y", "cell_id")
  } else {
    group_vars <- c("grid_x", "grid_y", "cell_id")
  }
  
  # Count track points in each grid cell
  grid_counts <- track_data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(count = dplyr::n(), .groups = 'drop') %>%
    dplyr::rename(x = grid_x, y = grid_y)
  
  cat("Created", nrow(grid_counts), "grid cells with track points\n")
  
  # Filter by count thresholds
  grid_counts <- grid_counts %>%
    dplyr::filter(count >= min_count_threshold, 
                  count <= max_count_threshold)
  
  cat("After count thresholds (", min_count_threshold, "<=count<=", max_count_threshold, "):", 
      nrow(grid_counts), "cells\n")
  
  if (nrow(grid_counts) == 0) {
    stop("No cells remain after filtering. Try adjusting count thresholds.")
  }
  
  # Determine grouping for sampling
  has_fish_id <- "fish_id" %in% names(grid_counts)
  has_time <- any(c("time_period_label", "time_period_posix") %in% names(grid_counts))
  process_by_groups <- by_group && (has_fish_id || has_time)
  
  if (process_by_groups) {
    # Build grouping columns
    sampling_group_cols <- c()
    if (has_fish_id) sampling_group_cols <- c(sampling_group_cols, "fish_id")
    if ("time_period_label" %in% names(grid_counts)) sampling_group_cols <- c(sampling_group_cols, "time_period_label")
    if ("time_period_posix" %in% names(grid_counts)) sampling_group_cols <- c(sampling_group_cols, "time_period_posix")
    
    # Get unique groups
    sampling_groups <- grid_counts %>%
      dplyr::distinct(dplyr::across(dplyr::all_of(sampling_group_cols))) %>%
      dplyr::arrange(dplyr::across(dplyr::all_of(sampling_group_cols)))
    
    cat("\nProcessing", nrow(sampling_groups), "fish-time combinations\n")
    
    # Initialize results list
    sampled_list <- list()
    
    # Process each group
    for (i in 1:nrow(sampling_groups)) {
      current_group <- sampling_groups[i, ]
      
      # Filter grid data for this group
      group_data <- grid_counts
      for (col in sampling_group_cols) {
        if (col %in% names(current_group)) {
          group_data <- group_data %>%
            dplyr::filter(.data[[col]] == current_group[[col]])
        }
      }
      
      if (nrow(group_data) == 0) next
      
      # Determine number of points to sample
      if (by_group) {
        group_n_points <- n_points
      } else {
        # Distribute points proportionally
        total_count_mass <- sum(group_data$count)
        overall_count_mass <- sum(grid_counts$count)
        group_n_points <- round(n_points * total_count_mass / overall_count_mass)
        group_n_points <- max(1, group_n_points)
      }
      
      # Prepare sampling weights (counts as weights)
      sampling_weights <- group_data$count
      
      # Normalize weights to probabilities
      sampling_weights <- sampling_weights / sum(sampling_weights)
      
      # Sample with replacement (reflects usage intensity)
      sampled_indices <- sample(
        1:nrow(group_data), 
        size = group_n_points, 
        replace = TRUE,
        prob = sampling_weights
      )
      
      # Extract sampled points
      sampled_data <- group_data[sampled_indices, ] %>%
        dplyr::mutate(
          sample_id = 1:group_n_points,
          group_id = paste(unlist(current_group), collapse = "_")
        )
      
      sampled_list[[i]] <- sampled_data
      
      # Progress message
      group_desc <- paste(names(current_group), "=", unlist(current_group), collapse = ", ")
      total_usage <- sum(group_data$count)
      cat("  ", group_desc, ": sampled", group_n_points, "points from", nrow(group_data), 
          "cells (total usage:", total_usage, "track points)\n")
    }
    
    # Combine all sampled points
    all_sampled <- dplyr::bind_rows(sampled_list)
    
  } else {
    # Simple sampling without groups
    cat("\nProcessing all data as single group\n")
    
    # Prepare sampling weights
    sampling_weights <- grid_counts$count
    sampling_weights <- sampling_weights / sum(sampling_weights)
    
    # Sample with replacement
    sampled_indices <- sample(
      1:nrow(grid_counts), 
      size = n_points, 
      replace = TRUE,
      prob = sampling_weights
    )
    
    # Extract sampled points
    all_sampled <- grid_counts[sampled_indices, ] %>%
      dplyr::mutate(sample_id = 1:n_points)
    
    total_usage <- sum(grid_counts$count)
    cat("Sampled", n_points, "points from", nrow(grid_counts), 
        "cells (total usage:", total_usage, "track points)\n")
  }
  
  # Re-number sample_id if needed
  if (!by_group && process_by_groups) {
    all_sampled$sample_id <- 1:nrow(all_sampled)
  }
  
  # Determine CRS to use
  if (!is.null(crs)) {
    use_crs <- sf::st_crs(crs)
  } else if (!is.null(reference_raster)) {
    use_crs <- sf::st_crs(reference_raster)
  } else {
    use_crs <- sf::st_crs(4326)
    message("No CRS specified. Using WGS84 (EPSG:4326) as default.")
  }
  
  # Convert to sf object
  sampled_sf <- sf::st_as_sf(
    all_sampled, 
    coords = c("x", "y"), 
    crs = use_crs
  )
  
  # Add attributes for reference
  attr(sampled_sf, "grid_resolution") <- grid_resolution
  attr(sampled_sf, "n_points") <- n_points
  attr(sampled_sf, "min_count_threshold") <- min_count_threshold
  attr(sampled_sf, "max_count_threshold") <- max_count_threshold
  attr(sampled_sf, "by_group") <- by_group
  attr(sampled_sf, "time_aggregation") <- time_aggregation
  
  # Summary
  cat("\n=== SAMPLING SUMMARY ===\n")
  cat("Total points sampled:", nrow(sampled_sf), "\n")
  cat("Grid resolution:", grid_resolution, "meters\n")
  if (has_fish_id) {
    cat("Fish ID(s):", paste(unique(sampled_sf$fish_id), collapse = ", "), "\n")
  }
  if (has_time && "time_period_label" %in% names(sampled_sf)) {
    cat("Time periods:", length(unique(sampled_sf$time_period_label)), "\n")
  }
  if ("group_id" %in% names(sampled_sf)) {
    cat("Groups:", length(unique(sampled_sf$group_id)), "\n")
  }
  cat("Usage range:", min(sampled_sf$count), "to", max(sampled_sf$count), "track points per cell\n")
  cat("Mean usage:", round(mean(sampled_sf$count), 2), "track points per cell\n")
  
  return(sampled_sf)
}