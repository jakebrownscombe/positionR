# Space Use Point Generation Functions
# Generate random points within space use estimates from calculate_space_use

library(dplyr)
library(ggplot2)

# Helper functions for timezone handling (replicated for consistency)
standardize_datetime <- function(datetime_col, target_tz = "UTC") {
  if (is.null(datetime_col) || length(datetime_col) == 0) {
    return(datetime_col)
  }
  
  # Get current timezone
  current_tz <- attr(datetime_col, "tzone")
  
  # If no timezone specified, assume UTC
  if (is.null(current_tz) || current_tz == "") {
    attr(datetime_col, "tzone") <- "UTC"
    current_tz <- "UTC"
  }
  
  # Convert to target timezone
  return(as.POSIXct(datetime_col, tz = target_tz))
}

# Helper to preserve original timezone info
store_timezone_info <- function(datetime_col) {
  tz <- attr(datetime_col, "tzone")
  if (is.null(tz) || tz == "") {
    return("UTC")
  }
  return(tz)
}

# Helper to restore timezone for output
restore_timezone <- function(datetime_col, original_tz) {
  if (is.null(original_tz) || original_tz == "") {
    original_tz <- "UTC"
  }
  return(as.POSIXct(datetime_col, tz = original_tz))
}

#' Generate random points within space use estimates from calculate_space_use
#'
#' @param space_use_results Output list from calculate_space_use function
#' @param track_data Original track data used for calculate_space_use (needed for polygon creation)
#' @param method Space use method to use for point generation ("constrained_convex_hull", "convex_hull", "bounding_box", "95_ellipse", "grid_cell_count")
#' @param n_points Number of random points to generate per group (default: 100)
#' @param fish_ids Which fish IDs to include (default: all available)
#' @param time_periods Which time periods to include (default: all available)
#' @param reference_raster Optional reference raster for grid-based methods
#' @param time_bin_size Time bin size for standardizing track data (default: 3600)
#' @param point_generation_method How to generate points within polygons ("random", "grid", "stratified")
#' @return Data frame with random points including fish_id, time_period, x, y coordinates, and method used
#' @export
generate_random_points_in_space_use <- function(space_use_results,
                                                track_data,
                                                method = "constrained_convex_hull",
                                                n_points = 100,
                                                fish_ids = NULL,
                                                time_periods = NULL,
                                                reference_raster = NULL,
                                                time_bin_size = 3600,
                                                point_generation_method = "random") {
  
  
  # Handle timezone standardization for datetime columns in track_data
  original_timezone <- "UTC"  # Default fallback
  if ("datetime" %in% names(track_data)) {
    original_timezone <- store_timezone_info(track_data$datetime)
    track_data$datetime <- standardize_datetime(track_data$datetime, "UTC")
  }
  
  # Validate inputs
  if (!is.list(space_use_results) || !"space_use_estimates" %in% names(space_use_results)) {
    stop("space_use_results must be output from calculate_space_use function")
  }
  
  if (is.null(track_data)) {
    stop("track_data is required for creating space use polygons")
  }
  
  space_use_df <- space_use_results$space_use_estimates
  
  # Check if method column exists
  method_area_col <- paste0(method, "_area_hectares")
  if (!method_area_col %in% names(space_use_df)) {
    available_methods <- gsub("_area_hectares$", "", names(space_use_df)[grepl("_area_hectares$", names(space_use_df))])
    stop("Method '", method, "' not found in space use results. Available methods: ", 
         paste(available_methods, collapse = ", "))
  }
  
  # Extract time binning information from space_use_results to match time periods
  time_bin_info <- extract_space_use_time_info(space_use_results)
  
  # Standardize track data to match space_use time periods
  track_std <- standardize_track_data_for_points(track_data, time_bin_info)
  
  # Filter space use estimates
  if (!is.null(fish_ids)) {
    space_use_df <- space_use_df[space_use_df$fish_id %in% fish_ids, ]
  }
  
  if (!is.null(time_periods)) {
    if ("time_period_numeric" %in% names(space_use_df)) {
      space_use_df <- space_use_df[space_use_df$time_period_numeric %in% time_periods, ]
    } else if ("time_period" %in% names(space_use_df)) {
      space_use_df <- space_use_df[space_use_df$time_period %in% time_periods, ]
    }
  }
  
  # Remove rows with NA area for the selected method
  space_use_df <- space_use_df[!is.na(space_use_df[[method_area_col]]), ]
  
  if (nrow(space_use_df) == 0) {
    stop("No valid space use estimates found after filtering")
  }
  
  
  # Generate points for each group
  point_list <- list()
  
  for (i in 1:nrow(space_use_df)) {
    group_data <- space_use_df[i, ]
    
    # Get time period for logging
    time_period_display <- if ("time_period_numeric" %in% names(group_data)) {
      group_data$time_period_numeric
    } else if ("time_period" %in% names(group_data)) {
      group_data$time_period
    } else {
      "unknown"
    }
    
    
    # Filter track data for this group - handle mixed column names
    track_filtered <- NULL
    
    # Get the time period value from group_data
    group_time_value <- NULL
    if ("time_period_numeric" %in% names(group_data)) {
      group_time_value <- group_data$time_period_numeric
    } else if ("time_period" %in% names(group_data)) {
      group_time_value <- group_data$time_period
    }
    
    # Filter track data using available columns - improved matching logic
    if (!is.null(group_time_value)) {
      # Space use data has time periods - match accordingly
      # Try multiple column combinations for robust matching
      if ("time_period_numeric" %in% names(track_std)) {
        # Try to match with time_period_numeric first
        track_filtered <- track_std %>%
          dplyr::filter(fish_id == group_data$fish_id & 
                       time_period_numeric == group_time_value)
      } else if ("time_period" %in% names(track_std)) {
        # Fallback to time_period
        track_filtered <- track_std %>%
          dplyr::filter(fish_id == group_data$fish_id & 
                       time_period == group_time_value)
      }
      
      # If no matches found with numeric, try matching by date label if available
      if ((is.null(track_filtered) || nrow(track_filtered) == 0) && 
          "time_period_label" %in% names(group_data) && 
          "time_period_date" %in% names(track_std)) {
        
        track_filtered <- track_std %>%
          dplyr::filter(fish_id == group_data$fish_id & 
                       time_period_date == group_data$time_period_label)
      }
    } else {
      # Space use data has no time periods - filter by fish only (all time combined)
      track_filtered <- track_std %>%
        dplyr::filter(fish_id == group_data$fish_id)
    }
    
    if (is.null(track_filtered) || nrow(track_filtered) == 0) {
      warning("Could not match time periods for fish ", group_data$fish_id, ". Available columns - track_std: ", 
              paste(grep("time_period", names(track_std), value = TRUE), collapse = ", "), 
              ", group_data: ", paste(grep("time_period", names(group_data), value = TRUE), collapse = ", "))
      next
    }
    
    if (nrow(track_filtered) < 3) {
      warning("Insufficient track points for fish ", group_data$fish_id, 
              " time period ", time_period_display, ". Skipping.")
      next
    }
    
    # Generate points based on method
    if (method == "grid_cell_count") {
      points <- generate_points_grid_cells(track_filtered, n_points, reference_raster, point_generation_method)
      
    } else if (method == "constrained_convex_hull") {
      points <- generate_points_constrained_hull(track_filtered, n_points, reference_raster, point_generation_method)
      
    } else if (method %in% c("convex_hull", "mcp")) {
      points <- generate_points_convex_hull(track_filtered, n_points, point_generation_method)
      
    } else if (method == "bounding_box") {
      points <- generate_points_bounding_box_method(track_filtered, n_points, point_generation_method)
      
    } else if (method == "95_ellipse") {
      points <- generate_points_ellipse(track_filtered, n_points, point_generation_method)
      
    } else {
      warning("Unknown method: ", method, ". Using convex hull instead.")
      points <- generate_points_convex_hull(track_filtered, n_points, point_generation_method)
    }
    
    if (!is.null(points) && nrow(points) > 0) {
      # Add group information
      points$fish_id <- group_data$fish_id
      
      # Add time period information
      if ("time_period_numeric" %in% names(group_data)) {
        points$time_period_numeric <- group_data$time_period_numeric
      } else if ("time_period" %in% names(group_data)) {
        points$time_period_numeric <- group_data$time_period
      }
      
      points$method <- method
      points$space_use_area_hectares <- group_data[[method_area_col]]
      points$point_id <- 1:nrow(points)
      
      # Add time date information if available
      if ("time_period_label" %in% names(group_data)) {
        points$time_period_date <- group_data$time_period_label
      } else if ("time_period_group" %in% names(group_data)) {
        points$time_period_date <- group_data$time_period_group
      } else if ("time_period_date" %in% names(group_data)) {
        points$time_period_date <- group_data$time_period_date
      }
      
      point_list[[i]] <- points
    }
  }
  
  # Combine all points
  if (length(point_list) > 0) {
    all_points <- dplyr::bind_rows(point_list)
    
    
    # Display time periods based on available columns
    
    # Restore original timezone in output data if datetime columns exist
    if (!is.null(original_timezone) && original_timezone != "UTC") {
      if ("time_period_date" %in% names(all_points)) {
        # Note: time_period_date is typically a character string, so no direct timezone conversion needed
        # But we preserve the original timezone info for consistency
      }
    }
    
    
    return(all_points)
  } else {
    warning("No points were generated")
    return(data.frame())
  }
}

# Helper functions for different point generation methods

#' Generate points within grid cells occupied by track
generate_points_grid_cells <- function(track_data, n_points, reference_raster = NULL, method = "random") {
  
  x_coords <- track_data$x
  y_coords <- track_data$y
  
  if (!is.null(reference_raster)) {
    # Use actual raster cells
    tryCatch({
      if (!requireNamespace("raster", quietly = TRUE)) {
        stop("raster package needed when reference_raster is provided")
      }
      
      # Find which raster cells contain the track points
      track_coords <- cbind(x_coords, y_coords)
      cell_indices <- raster::cellFromXY(reference_raster, track_coords)
      
      # Remove any NA cell indices (points outside raster)
      valid_cell_indices <- cell_indices[!is.na(cell_indices)]
      unique_cell_indices <- unique(valid_cell_indices)
      
      if (length(unique_cell_indices) == 0) {
        warning("No track points fall within valid raster cells")
        return(NULL)
      }
      
      # Get cell centers and dimensions
      cell_centers <- raster::xyFromCell(reference_raster, unique_cell_indices)
      res_x <- raster::res(reference_raster)[1]
      res_y <- raster::res(reference_raster)[2]
      
      # Distribute points among cells
      points_per_cell <- ceiling(n_points / length(unique_cell_indices))
      
      point_list <- list()
      
      for (i in seq_along(unique_cell_indices)) {
        center_x <- cell_centers[i, 1]
        center_y <- cell_centers[i, 2]
        
        # Generate points within this cell
        if (method == "random") {
          cell_x <- runif(points_per_cell, center_x - res_x/2, center_x + res_x/2)
          cell_y <- runif(points_per_cell, center_y - res_y/2, center_y + res_y/2)
        } else if (method == "grid") {
          # Create a mini-grid within the cell
          n_side <- ceiling(sqrt(points_per_cell))
          x_seq <- seq(center_x - res_x/2, center_x + res_x/2, length.out = n_side)
          y_seq <- seq(center_y - res_y/2, center_y + res_y/2, length.out = n_side)
          grid_coords <- expand.grid(x = x_seq, y = y_seq)
          if (nrow(grid_coords) > points_per_cell) {
            grid_coords <- grid_coords[sample(nrow(grid_coords), points_per_cell), ]
          }
          cell_x <- grid_coords$x
          cell_y <- grid_coords$y
        } else {
          # Default to random
          cell_x <- runif(points_per_cell, center_x - res_x/2, center_x + res_x/2)
          cell_y <- runif(points_per_cell, center_y - res_y/2, center_y + res_y/2)
        }
        
        point_list[[i]] <- data.frame(x = cell_x, y = cell_y)
      }
      
      # Combine and sample to exact n_points
      all_cell_points <- dplyr::bind_rows(point_list)
      if (nrow(all_cell_points) > n_points) {
        all_cell_points <- all_cell_points[sample(nrow(all_cell_points), n_points), ]
      }
      
      return(all_cell_points)
      
    }, error = function(e) {
      warning("Error using reference raster: ", e$message, ". Using grid approach.")
    })
  }
  
  # Fallback: use regular grid
  grid_resolution <- 100  # Default grid size
  x_grid <- floor(x_coords / grid_resolution) * grid_resolution
  y_grid <- floor(y_coords / grid_resolution) * grid_resolution
  
  unique_cells <- unique(data.frame(x_grid = x_grid, y_grid = y_grid))
  
  # Generate points within each grid cell
  point_list <- list()
  points_per_cell <- ceiling(n_points / nrow(unique_cells))
  
  for (i in 1:nrow(unique_cells)) {
    x_min <- unique_cells$x_grid[i]
    x_max <- x_min + grid_resolution
    y_min <- unique_cells$y_grid[i]
    y_max <- y_min + grid_resolution
    
    if (method == "random") {
      cell_x <- runif(points_per_cell, x_min, x_max)
      cell_y <- runif(points_per_cell, y_min, y_max)
    } else {
      # Grid within cell
      n_side <- ceiling(sqrt(points_per_cell))
      x_seq <- seq(x_min, x_max, length.out = n_side)
      y_seq <- seq(y_min, y_max, length.out = n_side)
      grid_coords <- expand.grid(x = x_seq, y = y_seq)
      if (nrow(grid_coords) > points_per_cell) {
        grid_coords <- grid_coords[sample(nrow(grid_coords), points_per_cell), ]
      }
      cell_x <- grid_coords$x
      cell_y <- grid_coords$y
    }
    
    point_list[[i]] <- data.frame(x = cell_x, y = cell_y)
  }
  
  # Combine and sample to exact n_points
  all_points <- dplyr::bind_rows(point_list)
  if (nrow(all_points) > n_points) {
    all_points <- all_points[sample(nrow(all_points), n_points), ]
  }
  
  return(all_points)
}

#' Generate points within constrained convex hull
generate_points_constrained_hull <- function(track_data, n_points, reference_raster = NULL, method = "random") {
  
  x_coords <- track_data$x
  y_coords <- track_data$y
  
  if (length(x_coords) < 3) return(NULL)
  
  if (!is.null(reference_raster)) {
    # Use actual raster cells within convex hull
    tryCatch({
      if (!requireNamespace("raster", quietly = TRUE)) {
        stop("raster package needed when reference_raster is provided")
      }
      
      # Get convex hull
      hull_indices <- grDevices::chull(x_coords, y_coords)
      hull_x <- x_coords[hull_indices]
      hull_y <- y_coords[hull_indices]
      
      # Find raster cells within hull - use same approach as area calculation
      res_x <- raster::res(reference_raster)[1]
      res_y <- raster::res(reference_raster)[2]
      
      # Use same buffer as area calculation to ensure edge cells are included
      hull_extent <- raster::extent(min(hull_x) - res_x, max(hull_x) + res_x,
                                   min(hull_y) - res_y, max(hull_y) + res_y)
      
      all_cells_in_bbox <- raster::cellsFromExtent(reference_raster, hull_extent)
      cell_values <- raster::extract(reference_raster, all_cells_in_bbox)
      valid_cells <- all_cells_in_bbox[!is.na(cell_values)]
      
      if (length(valid_cells) == 0) {
        warning("No valid cells in reference raster within hull extent")
        return(generate_points_convex_hull(track_data, n_points, method))
      }
      
      # Use the same cell intersection test as the area calculation 
      # (this was the issue - it was calling test_cells_in_polygon_full instead of test_cells_in_polygon)
      cells_intersect <- test_cells_in_polygon(reference_raster, valid_cells, hull_x, hull_y)
      hull_cell_indices <- valid_cells[cells_intersect]
      
      if (length(hull_cell_indices) == 0) {
        warning("No raster cells intersect with convex hull")
        return(generate_points_convex_hull(track_data, n_points, method))
      }
      
      # Generate points within these cells
      cell_centers <- raster::xyFromCell(reference_raster, hull_cell_indices)
      points_per_cell <- ceiling(n_points / length(hull_cell_indices))
      
      point_list <- list()
      
      for (i in seq_along(hull_cell_indices)) {
        center_x <- cell_centers[i, 1]
        center_y <- cell_centers[i, 2]
        
        if (method == "random") {
          cell_x <- runif(points_per_cell, center_x - res_x/2, center_x + res_x/2)
          cell_y <- runif(points_per_cell, center_y - res_y/2, center_y + res_y/2)
        } else {
          n_side <- ceiling(sqrt(points_per_cell))
          x_seq <- seq(center_x - res_x/2, center_x + res_x/2, length.out = n_side)
          y_seq <- seq(center_y - res_y/2, center_y + res_y/2, length.out = n_side)
          grid_coords <- expand.grid(x = x_seq, y = y_seq)
          if (nrow(grid_coords) > points_per_cell) {
            grid_coords <- grid_coords[sample(nrow(grid_coords), points_per_cell), ]
          }
          cell_x <- grid_coords$x
          cell_y <- grid_coords$y
        }
        
        point_list[[i]] <- data.frame(x = cell_x, y = cell_y)
      }
      
      all_points <- dplyr::bind_rows(point_list)
      if (nrow(all_points) > n_points) {
        all_points <- all_points[sample(nrow(all_points), n_points), ]
      }
      
      return(all_points)
      
    }, error = function(e) {
      warning("Error using reference raster: ", e$message, ". Using convex hull approach.")
    })
  }
  
  # Fallback to convex hull polygon
  return(generate_points_convex_hull(track_data, n_points, method))
}

#' Generate points within convex hull polygon
generate_points_convex_hull <- function(track_data, n_points, method = "random") {
  
  x_coords <- track_data$x
  y_coords <- track_data$y
  
  if (length(x_coords) < 3) return(NULL)
  
  # Get convex hull
  hull_indices <- grDevices::chull(x_coords, y_coords)
  hull_x <- x_coords[hull_indices]
  hull_y <- y_coords[hull_indices]
  
  # Get bounding box
  x_range <- range(hull_x)
  y_range <- range(hull_y)
  
  if (method == "random") {
    # Generate random points and filter to those inside hull
    # Generate extra points to account for rejection sampling
    n_candidate <- n_points * 3
    
    candidate_x <- runif(n_candidate, x_range[1], x_range[2])
    candidate_y <- runif(n_candidate, y_range[1], y_range[2])
    
    # Test which points are inside the hull
    inside_hull <- point_in_polygon(candidate_x, candidate_y, hull_x, hull_y)
    
    valid_points <- data.frame(
      x = candidate_x[inside_hull],
      y = candidate_y[inside_hull]
    )
    
    # Sample to exact number needed
    if (nrow(valid_points) >= n_points) {
      return(valid_points[sample(nrow(valid_points), n_points), ])
    } else {
      # If not enough points, try again with more candidates
      warning("Not enough points generated in first attempt, trying again")
      n_candidate <- n_points * 10
      candidate_x <- runif(n_candidate, x_range[1], x_range[2])
      candidate_y <- runif(n_candidate, y_range[1], y_range[2])
      inside_hull <- point_in_polygon(candidate_x, candidate_y, hull_x, hull_y)
      
      valid_points <- data.frame(
        x = candidate_x[inside_hull],
        y = candidate_y[inside_hull]
      )
      
      if (nrow(valid_points) >= n_points) {
        return(valid_points[sample(nrow(valid_points), n_points), ])
      } else {
        return(valid_points)  # Return what we have
      }
    }
    
  } else {
    # Grid method: create grid and filter to points inside hull
    n_side <- ceiling(sqrt(n_points * 2))  # Extra points for filtering
    x_seq <- seq(x_range[1], x_range[2], length.out = n_side)
    y_seq <- seq(y_range[1], y_range[2], length.out = n_side)
    
    grid_coords <- expand.grid(x = x_seq, y = y_seq)
    
    # Test which grid points are inside hull
    inside_hull <- point_in_polygon(grid_coords$x, grid_coords$y, hull_x, hull_y)
    
    valid_points <- grid_coords[inside_hull, ]
    
    if (nrow(valid_points) >= n_points) {
      return(valid_points[sample(nrow(valid_points), n_points), ])
    } else {
      return(valid_points)
    }
  }
}

#' Generate points within bounding box
generate_points_bounding_box_method <- function(track_data, n_points, method = "random") {
  
  x_coords <- track_data$x
  y_coords <- track_data$y
  
  x_range <- range(x_coords)
  y_range <- range(y_coords)
  
  if (method == "random") {
    return(data.frame(
      x = runif(n_points, x_range[1], x_range[2]),
      y = runif(n_points, y_range[1], y_range[2])
    ))
  } else {
    # Grid method
    n_side <- ceiling(sqrt(n_points))
    x_seq <- seq(x_range[1], x_range[2], length.out = n_side)
    y_seq <- seq(y_range[1], y_range[2], length.out = n_side)
    
    grid_coords <- expand.grid(x = x_seq, y = y_seq)
    
    if (nrow(grid_coords) > n_points) {
      grid_coords <- grid_coords[sample(nrow(grid_coords), n_points), ]
    }
    
    return(grid_coords)
  }
}

#' Generate points within 95% confidence ellipse
generate_points_ellipse <- function(track_data, n_points, method = "random") {
  
  x_coords <- track_data$x
  y_coords <- track_data$y
  
  if (length(x_coords) < 3) return(NULL)
  
  tryCatch({
    # Calculate ellipse parameters
    center_x <- mean(x_coords)
    center_y <- mean(y_coords)
    
    x_centered <- x_coords - center_x
    y_centered <- y_coords - center_y
    
    cov_matrix <- cov(cbind(x_centered, y_centered))
    eigen_result <- eigen(cov_matrix)
    eigenvalues <- eigen_result$values
    eigenvectors <- eigen_result$vectors
    
    # 95% confidence level
    chi_sq_95 <- qchisq(0.95, df = 2)
    
    # Semi-axes lengths
    a <- sqrt(chi_sq_95 * eigenvalues[1])
    b <- sqrt(chi_sq_95 * eigenvalues[2])
    
    # Rotation angle
    angle <- atan2(eigenvectors[2, 1], eigenvectors[1, 1])
    
    if (method == "random") {
      # Generate random points within ellipse
      # Use rejection sampling within bounding box of ellipse
      
      # Calculate bounding box of rotated ellipse
      cos_angle <- cos(angle)
      sin_angle <- sin(angle)
      
      # Maximum extent in x and y directions
      x_extent <- abs(a * cos_angle) + abs(b * sin_angle)
      y_extent <- abs(a * sin_angle) + abs(b * cos_angle)
      
      x_min <- center_x - x_extent
      x_max <- center_x + x_extent
      y_min <- center_y - y_extent
      y_max <- center_y + y_extent
      
      # Generate candidate points
      n_candidate <- n_points * 3
      candidate_x <- runif(n_candidate, x_min, x_max)
      candidate_y <- runif(n_candidate, y_min, y_max)
      
      # Transform to ellipse coordinate system
      x_translated <- candidate_x - center_x
      y_translated <- candidate_y - center_y
      
      x_rotated <- x_translated * cos_angle + y_translated * sin_angle
      y_rotated <- -x_translated * sin_angle + y_translated * cos_angle
      
      # Test if points are inside ellipse
      ellipse_test <- (x_rotated^2 / a^2) + (y_rotated^2 / b^2) <= 1
      
      valid_points <- data.frame(
        x = candidate_x[ellipse_test],
        y = candidate_y[ellipse_test]
      )
      
      if (nrow(valid_points) >= n_points) {
        return(valid_points[sample(nrow(valid_points), n_points), ])
      } else {
        return(valid_points)
      }
      
    } else {
      # Grid method: create grid and filter to points inside ellipse
      n_side <- ceiling(sqrt(n_points * 2))
      
      # Calculate grid bounds
      cos_angle <- cos(angle)
      sin_angle <- sin(angle)
      x_extent <- abs(a * cos_angle) + abs(b * sin_angle)
      y_extent <- abs(a * sin_angle) + abs(b * cos_angle)
      
      x_seq <- seq(center_x - x_extent, center_x + x_extent, length.out = n_side)
      y_seq <- seq(center_y - y_extent, center_y + y_extent, length.out = n_side)
      
      grid_coords <- expand.grid(x = x_seq, y = y_seq)
      
      # Transform and test
      x_translated <- grid_coords$x - center_x
      y_translated <- grid_coords$y - center_y
      
      x_rotated <- x_translated * cos_angle + y_translated * sin_angle
      y_rotated <- -x_translated * sin_angle + y_translated * cos_angle
      
      ellipse_test <- (x_rotated^2 / a^2) + (y_rotated^2 / b^2) <= 1
      
      valid_points <- grid_coords[ellipse_test, ]
      
      if (nrow(valid_points) >= n_points) {
        return(valid_points[sample(nrow(valid_points), n_points), ])
      } else {
        return(valid_points)
      }
    }
    
  }, error = function(e) {
    warning("Error generating ellipse points: ", e$message)
    return(generate_points_bounding_box_method(track_data, n_points, method))
  })
}

# Utility functions

#' Extract time binning information from space_use_results
extract_space_use_time_info <- function(space_use_results) {
  
  # Initialize default time binning info
  time_bin_info <- list(
    method = "bin_size",           # Default to numeric binning
    bin_size = 3600,              # Default 1 hour
    aggregation_method = "hour",   # Default aggregation
    available_periods = c(),       # Available time periods in space_use_results
    available_labels = c()         # Available time labels
  )
  
  if (!is.null(space_use_results$space_use_estimates)) {
    space_use_df <- space_use_results$space_use_estimates
    
    # Extract available time periods - handle different column names
    if ("time_period_numeric" %in% names(space_use_df)) {
      time_bin_info$available_periods <- unique(space_use_df$time_period_numeric)
    } else if ("time_period" %in% names(space_use_df)) {
      time_bin_info$available_periods <- unique(space_use_df$time_period)
      # Store both the column name and values for better matching
      time_bin_info$space_use_time_col <- "time_period"
    }
    
    # Check if time_period_date or legacy columns exist (indicates POSIX-based grouping)
    if ("time_period_date" %in% names(space_use_df)) {
      time_bin_info$method <- "posix_based"
      time_bin_info$available_labels <- unique(space_use_df$time_period_date)
    } else if ("time_period_label" %in% names(space_use_df)) {
      time_bin_info$method <- "posix_based"
      time_bin_info$available_labels <- unique(space_use_df$time_period_label)
    } else if ("time_period_group" %in% names(space_use_df)) {
      time_bin_info$method <- "posix_based"  
      time_bin_info$available_labels <- unique(space_use_df$time_period_group)
    }
    
    # Try to infer aggregation method from labels (for both time_period_label and time_period_group)
    if (time_bin_info$method == "posix_based") {
      labels <- time_bin_info$available_labels[!is.na(time_bin_info$available_labels)]
      if (length(labels) > 0) {
        sample_label <- labels[1]
        if (grepl("\\d{4}-\\d{2}-\\d{2}$", sample_label)) {
          time_bin_info$aggregation_method <- "day"
        } else if (grepl("\\d{4}-\\d{2}$", sample_label)) {
          time_bin_info$aggregation_method <- "month"
        } else if (grepl("\\d{4}-\\d{2}-\\d{2}\\s+\\d{2}:\\d{2}", sample_label)) {
          time_bin_info$aggregation_method <- "hour"
        }
      }
    }
    
    # Try to infer bin size from numeric time periods if available
    if (time_bin_info$method == "bin_size" && length(time_bin_info$available_periods) > 1) {
      periods <- sort(unique(time_bin_info$available_periods))
      if (is.numeric(periods) && length(periods) >= 2) {
        time_diff <- periods[2] - periods[1]
        if (time_diff > 0) {
          time_bin_info$bin_size <- time_diff
        }
      }
    }
  }
  
  
  return(time_bin_info)
}

#' Alias for test_cells_in_polygon_simple for backward compatibility
test_cells_in_polygon_full <- function(raster_obj, cell_indices, poly_x, poly_y) {
  return(test_cells_in_polygon_simple(raster_obj, cell_indices, poly_x, poly_y))
}

#' Test if raster cells intersect with a polygon using cell centers (original working version)
#' This is the original simple method that was working
test_cells_in_polygon_simple <- function(raster_obj, cell_indices, poly_x, poly_y) {
  
  tryCatch({
    # Get cell centers
    cell_centers <- raster::xyFromCell(raster_obj, cell_indices)
    
    # Test cell centers (simple approach that was working)
    centers_inside <- point_in_polygon(cell_centers[, 1], cell_centers[, 2], poly_x, poly_y)
    
    return(centers_inside)
    
  }, error = function(e) {
    warning("Error testing cell intersection: ", e$message)
    return(rep(FALSE, length(cell_indices)))
  })
}

#' Point-in-polygon test using ray casting algorithm
point_in_polygon <- function(px, py, poly_x, poly_y) {
  
  n_points <- length(px)
  n_poly <- length(poly_x)
  inside <- logical(n_points)
  
  for (i in 1:n_points) {
    x <- px[i]
    y <- py[i]
    
    # Ray casting algorithm
    n_intersections <- 0
    j <- n_poly
    
    for (k in 1:n_poly) {
      if (((poly_y[k] > y) != (poly_y[j] > y)) &&
          (x < (poly_x[j] - poly_x[k]) * (y - poly_y[k]) / (poly_y[j] - poly_y[k]) + poly_x[k])) {
        n_intersections <- n_intersections + 1
      }
      j <- k
    }
    
    inside[i] <- (n_intersections %% 2) == 1
  }
  
  return(inside)
}

#' Standardize track data for point generation to match space_use_results time binning
standardize_track_data_for_points <- function(track_data, time_bin_info) {
  
  track_std <- track_data
  
  # Store original timezone for consistent time binning (before any UTC conversion)
  original_timezone <- "UTC"  # Default fallback
  if ("datetime" %in% names(track_std)) {
    original_timezone <- store_timezone_info(track_std$datetime)
  }
  
  # Handle different fish ID column names
  if ("path_id" %in% names(track_std) && !"fish_id" %in% names(track_std)) {
    track_std$fish_id <- track_std$path_id
  }
  
  # Match track data to space_use_results time bins
  if (!"time_period_numeric" %in% names(track_std) && !"time_period" %in% names(track_std)) {
    
    if (time_bin_info$method == "posix_based" && "datetime" %in% names(track_std)) {
      # Use POSIX datetime binning to match space_use_results
      
      # Create time periods using the same method as calculate_space_use (preserving original timezone)
      if (time_bin_info$aggregation_method == "hour") {
        # Group by date and hour (matching calculate_space_use exactly)
        track_std$time_period_date <- format(track_std$datetime, "%Y-%m-%d %H:00")
        track_std$time_period_numeric <- as.numeric(as.POSIXct(track_std$time_period_date, tz = original_timezone))
      } else if (time_bin_info$aggregation_method == "day") {
        # Group by date only (matching calculate_space_use exactly)
        track_std$time_period_date <- as.character(as.Date(track_std$datetime))
        track_std$time_period_numeric <- as.numeric(as.POSIXct(track_std$time_period_date, tz = original_timezone))
      } else if (time_bin_info$aggregation_method == "month") {
        # Group by year and month (matching calculate_space_use exactly)
        track_std$time_period_date <- format(track_std$datetime, "%Y-%m")
        track_std$time_period_numeric <- as.numeric(as.POSIXct(paste0(track_std$time_period_date, "-01"), tz = original_timezone))
      } else {
        # Fallback for other methods
        track_std$time_period_date <- as.character(track_std$datetime)
        track_std$time_period_numeric <- as.numeric(track_std$datetime)
      }
      
      # Now standardize datetime to UTC for internal processing (after time binning is complete)
      track_std$datetime <- standardize_datetime(track_std$datetime, "UTC")
      
    } else if (time_bin_info$method == "bin_size" && "time_seconds" %in% names(track_std)) {
      # Use numeric bin size method (fallback)
      track_std$time_period_numeric <- floor(track_std$time_seconds / time_bin_info$bin_size) * time_bin_info$bin_size
      
    } else if ("time_bin" %in% names(track_std)) {
      track_std$time_period_numeric <- track_std$time_bin
      
    } else {
      stop("Cannot match time periods: need 'datetime' column for POSIX matching or 'time_seconds' for numeric binning")
    }
  }
  
  # Ensure track time periods match available space_use_results time periods
  if (length(time_bin_info$available_periods) > 0) {
    # Filter to only include time periods that exist in space_use_results
    original_nrow <- nrow(track_std)
    
    # Debug: show what we're trying to match
    
    # Try to match with the appropriate column - multiple approaches
    matched_successfully <- FALSE
    
    # Approach 1: Direct numeric matching
    if ("time_period_numeric" %in% names(track_std)) {
      filtered_track <- track_std %>%
        dplyr::filter(time_period_numeric %in% time_bin_info$available_periods)
      
      if (nrow(filtered_track) > 0) {
        track_std <- filtered_track
        matched_successfully <- TRUE
      }
    } else if ("time_period" %in% names(track_std)) {
      filtered_track <- track_std %>%
        dplyr::filter(time_period %in% time_bin_info$available_periods)
      
      if (nrow(filtered_track) > 0) {
        track_std <- filtered_track
        matched_successfully <- TRUE
      }
    }
    
    # Approach 2: If direct matching failed, try date-based matching
    if (!matched_successfully && "time_period_date" %in% names(track_std) && 
        length(time_bin_info$available_labels) > 0) {
      filtered_track <- track_std %>%
        dplyr::filter(time_period_date %in% time_bin_info$available_labels)
      
      if (nrow(filtered_track) > 0) {
        track_std <- filtered_track
        matched_successfully <- TRUE
      }
    }
    
    # Approach 3: If still no matches, disable filtering (use all data)
    if (!matched_successfully) {
      # Don't filter, keep original track_std
      # But ensure we have the needed time columns for later matching
      if ("time_period_numeric" %in% names(track_std) && !"time_period" %in% names(track_std)) {
        track_std$time_period <- track_std$time_period_numeric
      }
    }
    
    
    # If no matches, this indicates a fundamental mismatch in time binning approach
    if (nrow(track_std) == 0) {
      warning("No matching time periods found between track_data and space_use_results. Time binning methods may be incompatible.")
    }
  }
  
  # Ensure required columns exist
  required_cols <- c("fish_id", "x", "y")
  missing_cols <- setdiff(required_cols, names(track_std))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in track_data: ", paste(missing_cols, collapse = ", "))
  }
  
  # Ensure fish_id exists and maintain its type (character or numeric)
  # No conversion - keep fish_id as-is to support both character and numeric IDs
  
  # Remove duplicate points (same x,y for same fish/time)
  if ("time_period" %in% names(track_std)) {
    track_std <- track_std %>%
      dplyr::distinct(fish_id, time_period, x, y, .keep_all = TRUE)
  } else {
    track_std <- track_std %>%
      dplyr::distinct(fish_id, x, y, .keep_all = TRUE)
  }
  
  return(track_std)
}
