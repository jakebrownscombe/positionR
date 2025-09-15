# Complete Space Use Functions - Fixed Version
# This includes your original calculate_space_use function plus the fixed helper functions

library(dplyr)

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

#' Calculate Space Use Estimates from Track Data
#'
#' @param track_data Data frame with fish tracks (columns: fish_id/path_id, datetime, x, y)
#' @param by_fish Calculate separate estimates for each fish (default: TRUE)
#' @param by_time_period Calculate estimates by time periods (default: FALSE)
#' @param time_aggregation How to aggregate time: "hour", "day", "month", "none" (default: "day")
#' @param methods Vector of methods to use (default: all available)
#' @param grid_resolution Grid resolution for raster-based methods (default: 100m)
#' @param reference_raster Optional raster object to use for grid-based methods
#' @return List with space use estimates and summary statistics
calculate_space_use <- function(track_data,
                                by_fish = TRUE,
                                by_time_period = FALSE,
                                time_aggregation = "day",
                                methods = c("convex_hull", "bounding_box", "mcp", "grid_cell_count", "95_ellipse", "constrained_convex_hull"),
                                grid_resolution = 100,
                                reference_raster = NULL,
                                prob_column = NULL,
                                prob_thresholds = NULL,
                                verbose = TRUE) {

  if (verbose) if(getOption("positionR.verbose", TRUE)) cat("=== SPACE USE ESTIMATION ===\n")

  # Handle probability thresholds if specified
  if (!is.null(prob_column) && !is.null(prob_thresholds)) {
    if (!prob_column %in% names(track_data)) {
      stop("Probability column '", prob_column, "' not found in track_data")
    }
    
    # Convert single threshold to vector
    if (length(prob_thresholds) == 1) {
      prob_thresholds <- c(prob_thresholds)
    }
    
    if (verbose) if(getOption("positionR.verbose", TRUE)) cat("Using probability thresholds:", paste(prob_thresholds, collapse = ", "), "\n")
    
    # Store original data
    track_data_original <- track_data
    
    # Process each threshold
    all_results <- list()
    
    for (threshold in prob_thresholds) {
      if (verbose) if(getOption("positionR.verbose", TRUE)) cat("\nProcessing threshold:", threshold, "\n")
      
      # Filter data by threshold
      track_data_filtered <- track_data_original %>%
        dplyr::filter(.data[[prob_column]] >= threshold)
      
      if (verbose) if(getOption("positionR.verbose", TRUE)) cat("  Points after filtering:", nrow(track_data_filtered), 
          "(", round(nrow(track_data_filtered)/nrow(track_data_original)*100, 1), "% )\n")
      
      # Process with filtered data
      threshold_results <- calculate_space_use_single(
        track_data_filtered,
        by_fish = by_fish,
        by_time_period = by_time_period,
        time_aggregation = time_aggregation,
        methods = methods,
        grid_resolution = grid_resolution,
        reference_raster = reference_raster,
        verbose = verbose
      )
      
      # Add threshold info to results
      threshold_results$threshold <- threshold
      all_results[[paste0("threshold_", threshold)]] <- threshold_results
    }
    
    # If multiple thresholds, return combined results
    if (length(prob_thresholds) > 1) {
      return(list(
        results_by_threshold = all_results,
        thresholds = prob_thresholds,
        prob_column = prob_column
      ))
    } else {
      # Single threshold - return regular results with threshold info
      return(all_results[[1]])
    }
  }

  # Regular processing without probability filtering
  return(calculate_space_use_single(
    track_data,
    by_fish = by_fish,
    by_time_period = by_time_period,
    time_aggregation = time_aggregation,
    methods = methods,
    grid_resolution = grid_resolution,
    reference_raster = reference_raster,
    verbose = verbose
  ))
}

# Internal function for single threshold processing
calculate_space_use_single <- function(track_data,
                                      by_fish = TRUE,
                                      by_time_period = FALSE,
                                      time_aggregation = "day",
                                      methods = c("convex_hull", "bounding_box", "mcp", "grid_cell_count", "95_ellipse", "constrained_convex_hull"),
                                      grid_resolution = 100,
                                      reference_raster = NULL,
                                      verbose = TRUE) {
  
  # Handle timezone standardization for datetime columns
  original_timezone <- "UTC"  # Default fallback
  if ("datetime" %in% names(track_data)) {
    original_timezone <- store_timezone_info(track_data$datetime)
    track_data$datetime <- standardize_datetime(track_data$datetime, "UTC")
    if(getOption("positionR.verbose", TRUE)) cat("Timezone standardized: Input =", original_timezone, ", processing in UTC\n")
  }

  # Standardize track data
  track_std <- standardize_track_data_space(track_data, time_aggregation)

  if(getOption("positionR.verbose", TRUE)) cat("Track data summary:\n")
  if(getOption("positionR.verbose", TRUE)) cat("  Fish IDs:", paste(unique(track_std$fish_id), collapse = ", "), "\n")
  if (by_time_period && "time_period" %in% names(track_std)) {
    if ("time_period_label" %in% names(track_std)) {
      if(getOption("positionR.verbose", TRUE)) cat("  Time periods:", paste(sort(unique(track_std$time_period_label)), collapse = ", "), "\n")
    } else {
      if(getOption("positionR.verbose", TRUE)) cat("  Time periods:", paste(sort(unique(track_std$time_period)), collapse = ", "), "\n")
    }
  }
  if(getOption("positionR.verbose", TRUE)) cat("  Total track points:", nrow(track_std), "\n")
  if(getOption("positionR.verbose", TRUE)) cat("  Methods:", paste(methods, collapse = ", "), "\n")
  if (!is.null(reference_raster)) {
    if(getOption("positionR.verbose", TRUE)) cat("  Using reference raster for grid-based methods\n")
  }

  # Determine grouping variables
  grouping_vars <- c()
  if (by_fish) {
    grouping_vars <- c(grouping_vars, "fish_id")
  }
  if (by_time_period) {
    grouping_vars <- c(grouping_vars, "time_period")
  }

  # Calculate space use for each group
  space_use_results <- list()

  if (length(grouping_vars) > 0) {
    # Group by specified variables
    grouped_data <- track_std %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars))) %>%
      dplyr::group_split()

    group_keys <- track_std %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars))) %>%
      dplyr::group_keys()

    for (i in seq_along(grouped_data)) {
      group_data <- grouped_data[[i]]
      group_key <- group_keys[i, ]

      # Create group identifier
      group_id <- paste(names(group_key), group_key, sep = "=", collapse = "_")
      if(getOption("positionR.verbose", TRUE)) cat("\nProcessing group:", group_id, "with", nrow(group_data), "points\n")

      # Calculate space use metrics for this group
      group_result <- calculate_group_space_use(group_data, methods, grid_resolution, reference_raster)

      # Add group information to summary
      for (col in names(group_key)) {
        group_result$summary[[col]] <- group_key[[col]]
      }
      
      # Also add time_period_label if it exists in the group data
      if ("time_period_label" %in% names(group_data)) {
        group_result$summary$time_period_label <- group_data$time_period_label[1]
      }

      # Add group metadata to spatial data
      for (method_name in names(group_result$spatial_data)) {
        if (!is.null(group_result$spatial_data[[method_name]])) {
          group_result$spatial_data[[method_name]]$group_id <- group_id
          for (col in names(group_key)) {
            group_result$spatial_data[[method_name]][[col]] <- group_key[[col]]
          }
          if ("time_period_label" %in% names(group_data)) {
            group_result$spatial_data[[method_name]]$time_period_label <- group_data$time_period_label[1]
          }
        }
      }

      space_use_results[[group_id]] <- group_result
    }

  } else {
    # Process all data as single group
    if(getOption("positionR.verbose", TRUE)) cat("\nProcessing combined data with", nrow(track_std), "points\n")
    group_result <- calculate_group_space_use(track_std, methods, grid_resolution, reference_raster)
    group_result$summary$group <- "combined"
    
    # Add group metadata to spatial data
    for (method_name in names(group_result$spatial_data)) {
      if (!is.null(group_result$spatial_data[[method_name]])) {
        group_result$spatial_data[[method_name]]$group_id <- "combined"
        group_result$spatial_data[[method_name]]$group <- "combined"
      }
    }
    
    space_use_results[["combined"]] <- group_result
  }

  # Combine summary results
  summary_list <- lapply(space_use_results, function(x) x$summary)
  space_use_df <- dplyr::bind_rows(summary_list, .id = "group_id")
  
  # Combine spatial data by method
  spatial_data_combined <- list()
  for (method in methods) {
    method_data_list <- list()
    for (group_id in names(space_use_results)) {
      if (method %in% names(space_use_results[[group_id]]$spatial_data)) {
        method_data <- space_use_results[[group_id]]$spatial_data[[method]]
        if (!is.null(method_data) && nrow(method_data) > 0) {
          method_data_list[[group_id]] <- method_data
        }
      }
    }
    if (length(method_data_list) > 0) {
      spatial_data_combined[[method]] <- dplyr::bind_rows(method_data_list)
    }
  }
  
  # Also combine hull outline data for constrained convex hull
  if ("constrained_convex_hull" %in% methods) {
    outline_data_list <- list()
    for (group_id in names(space_use_results)) {
      if ("constrained_convex_hull_outline" %in% names(space_use_results[[group_id]]$spatial_data)) {
        outline_data <- space_use_results[[group_id]]$spatial_data[["constrained_convex_hull_outline"]]
        if (!is.null(outline_data) && nrow(outline_data) > 0) {
          outline_data_list[[group_id]] <- outline_data
        }
      }
    }
    if (length(outline_data_list) > 0) {
      spatial_data_combined[["constrained_convex_hull_outline"]] <- dplyr::bind_rows(outline_data_list)
    }
  }

  # Create summary statistics
  summary_stats <- create_space_use_summary(space_use_df, methods, grouping_vars)

  # Restore original timezone in output data if it was processed
  if (!is.null(original_timezone) && original_timezone != "UTC") {
    # Restore timezone in space_use_estimates that contain datetime columns
    if ("time_period_label" %in% names(space_use_df)) {
      # Convert time_period_label back to original timezone if it's a date string
      if(getOption("positionR.verbose", TRUE)) cat("Restoring timezone in output: UTC ->", original_timezone, "\n")
    }
  }

  if(getOption("positionR.verbose", TRUE)) cat("\n=== SPACE USE ESTIMATION COMPLETE ===\n")
  if(getOption("positionR.verbose", TRUE)) cat("Groups processed:", nrow(space_use_df), "\n")

  return(list(
    space_use_estimates = space_use_df,
    spatial_data = spatial_data_combined,
    summary_stats = summary_stats,
    parameters = list(
      methods = methods,
      by_fish = by_fish,
      by_time_period = by_time_period,
      time_aggregation = time_aggregation,
      grid_resolution = grid_resolution,
      grouping_vars = grouping_vars,
      reference_raster = reference_raster,
      original_timezone = original_timezone
    )
  ))
}

#' Standardize track data for space use analysis
standardize_track_data_space <- function(track_data, time_aggregation) {

  track_std <- track_data

  # Handle different fish ID column names
  if ("path_id" %in% names(track_std) && !"fish_id" %in% names(track_std)) {
    track_std$fish_id <- track_std$path_id
  }

  # Handle time period conversion based on datetime column
  if ("datetime" %in% names(track_std) && time_aggregation != "none") {
    # Ensure datetime is POSIXct
    if (!inherits(track_std$datetime, "POSIXct")) {
      track_std$datetime <- as.POSIXct(track_std$datetime)
    }
    
    if (time_aggregation == "hour") {
      # Group by date and hour
      track_std$time_period_label <- format(track_std$datetime, "%Y-%m-%d %H:00")
      track_std$time_period <- as.numeric(as.POSIXct(track_std$time_period_label))
      if(getOption("positionR.verbose", TRUE)) cat("Aggregating time by hour\n")
    } else if (time_aggregation == "day") {
      # Group by date only
      track_std$time_period_label <- as.character(as.Date(track_std$datetime))
      track_std$time_period <- as.numeric(as.POSIXct(track_std$time_period_label))
      if(getOption("positionR.verbose", TRUE)) cat("Aggregating time by day\n")
    } else if (time_aggregation == "month") {
      # Group by year and month
      track_std$time_period_label <- format(track_std$datetime, "%Y-%m")
      track_std$time_period <- as.numeric(as.POSIXct(paste0(track_std$time_period_label, "-01")))
      if(getOption("positionR.verbose", TRUE)) cat("Aggregating time by month\n")
    }
  } else if (!"time_period" %in% names(track_std)) {
    # Fallback for data without datetime column
    if ("time_seconds" %in% names(track_std)) {
      # Default to hourly bins if no datetime
      track_std$time_period <- floor(track_std$time_seconds / 3600) * 3600
      if(getOption("positionR.verbose", TRUE)) cat("No datetime column found. Using hourly bins from time_seconds\n")
    } else if ("time_bin" %in% names(track_std)) {
      track_std$time_period <- track_std$time_bin
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

#' Calculate space use metrics for a single group
calculate_group_space_use <- function(group_data, methods, grid_resolution, reference_raster = NULL) {

  if (nrow(group_data) < 3) {
    warning("Group has fewer than 3 points, skipping complex calculations")
    return(list(
      summary = data.frame(
        n_points = nrow(group_data),
        convex_hull_area_m2 = NA,
        convex_hull_area_hectares = NA,
        bounding_box_area_m2 = NA,
        bounding_box_area_hectares = NA,
        mcp_area_m2 = NA,
        mcp_area_hectares = NA,
        grid_cells_used = NA,
        grid_cell_count_area_m2 = NA,
        grid_cell_count_area_hectares = NA,
        ellipse_95_area_m2 = NA,
        ellipse_95_area_hectares = NA,
        range_x_m = NA,
        range_y_m = NA,
        centroid_x = NA,
        centroid_y = NA
      ),
      spatial_data = list()
    ))
  }

  result <- data.frame(n_points = nrow(group_data))
  spatial_data <- list()

  # Basic metrics
  x_coords <- group_data$x
  y_coords <- group_data$y

  result$range_x_m <- diff(range(x_coords, na.rm = TRUE))
  result$range_y_m <- diff(range(y_coords, na.rm = TRUE))
  result$centroid_x <- mean(x_coords, na.rm = TRUE)
  result$centroid_y <- mean(y_coords, na.rm = TRUE)

  # Method 1: Convex Hull (Minimum Convex Polygon)
  if ("convex_hull" %in% methods || "mcp" %in% methods) {
    hull_area <- calculate_convex_hull_area(x_coords, y_coords)
    result$convex_hull_area_m2 <- hull_area
    result$convex_hull_area_hectares <- hull_area / 10000
    result$mcp_area_m2 <- hull_area  # MCP is same as convex hull
    result$mcp_area_hectares <- hull_area / 10000
    
    # Store spatial data for convex hull
    hull_polygon <- create_space_use_polygon_data(x_coords, y_coords, "convex_hull")
    if (!is.null(hull_polygon)) {
      spatial_data$convex_hull <- hull_polygon
      spatial_data$mcp <- hull_polygon  # Same as convex hull
    }
  }

  # Method 2: Bounding Box (Rectangle)
  if ("bounding_box" %in% methods) {
    bbox_area <- result$range_x_m * result$range_y_m
    result$bounding_box_area_m2 <- bbox_area
    result$bounding_box_area_hectares <- bbox_area / 10000
    
    # Store spatial data for bounding box
    bbox_polygon <- create_space_use_polygon_data(x_coords, y_coords, "bounding_box")
    if (!is.null(bbox_polygon)) {
      spatial_data$bounding_box <- bbox_polygon
    }
  }

  # Method 3: Grid Cell Count (UPDATED)
  if ("grid_cell_count" %in% methods) {
    grid_metrics <- calculate_grid_cell_usage_with_data(x_coords, y_coords, grid_resolution, reference_raster)
    result$grid_cells_used <- grid_metrics$n_cells
    result$grid_cell_count_area_m2 <- grid_metrics$total_area_m2
    result$grid_cell_count_area_hectares <- grid_metrics$total_area_m2 / 10000
    
    # Store spatial data for grid cells
    if (!is.null(grid_metrics$spatial_data)) {
      spatial_data$grid_cell_count <- grid_metrics$spatial_data
    }
  }

  # Method 4: 95% Confidence Ellipse
  if ("95_ellipse" %in% methods) {
    ellipse_area <- calculate_95_ellipse_area(x_coords, y_coords)
    result$ellipse_95_area_m2 <- ellipse_area
    result$ellipse_95_area_hectares <- ellipse_area / 10000
    
    # Store spatial data for ellipse
    ellipse_polygon <- create_space_use_polygon_data(x_coords, y_coords, "95_ellipse")
    if (!is.null(ellipse_polygon)) {
      spatial_data$ellipse_95 <- ellipse_polygon
    }
  }

  # Method 5: Constrained Convex Hull (UPDATED)
  if ("constrained_convex_hull" %in% methods) {
    constrained_metrics <- calculate_constrained_convex_hull_area_with_data(x_coords, y_coords, grid_resolution, reference_raster)
    result$constrained_convex_hull_area_m2 <- constrained_metrics$area
    result$constrained_convex_hull_area_hectares <- constrained_metrics$area / 10000
    
    # Store spatial data for constrained hull
    if (!is.null(constrained_metrics$spatial_data)) {
      spatial_data$constrained_convex_hull <- constrained_metrics$spatial_data
    }
    # Store hull outline for plotting
    if (!is.null(constrained_metrics$hull_outline)) {
      spatial_data$constrained_convex_hull_outline <- constrained_metrics$hull_outline
    }
  }

  return(list(
    summary = result,
    spatial_data = spatial_data
  ))
}

#' Calculate convex hull area
calculate_convex_hull_area <- function(x_coords, y_coords) {

  if (length(x_coords) < 3) return(NA)

  tryCatch({
    # Find convex hull points
    hull_indices <- grDevices::chull(x_coords, y_coords)
    hull_x <- x_coords[hull_indices]
    hull_y <- y_coords[hull_indices]

    # Calculate area using shoelace formula
    n <- length(hull_x)
    area <- 0.5 * abs(sum(hull_x * c(hull_y[-1], hull_y[1]) - c(hull_x[-1], hull_x[1]) * hull_y))

    return(area)

  }, error = function(e) {
    warning("Error calculating convex hull: ", e$message)
    return(NA)
  })
}

#' Calculate 95% confidence ellipse area
calculate_95_ellipse_area <- function(x_coords, y_coords) {

  if (length(x_coords) < 3) return(NA)

  tryCatch({
    # Center the coordinates
    x_centered <- x_coords - mean(x_coords)
    y_centered <- y_coords - mean(y_coords)

    # Calculate covariance matrix
    cov_matrix <- cov(cbind(x_centered, y_centered))

    # Eigenvalues and eigenvectors
    eigen_result <- eigen(cov_matrix)
    eigenvalues <- eigen_result$values

    # 95% confidence level (chi-square with 2 df)
    chi_sq_95 <- qchisq(0.95, df = 2)

    # Semi-axes lengths
    a <- sqrt(chi_sq_95 * eigenvalues[1])  # Semi-major axis
    b <- sqrt(chi_sq_95 * eigenvalues[2])  # Semi-minor axis

    # Ellipse area
    ellipse_area <- pi * a * b

    return(ellipse_area)

  }, error = function(e) {
    warning("Error calculating 95% ellipse: ", e$message)
    return(NA)
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

#' Create summary statistics for space use estimates
create_space_use_summary <- function(space_use_df, methods, grouping_vars) {

  summary_list <- list()

  # Overall summary
  numeric_cols <- sapply(space_use_df, is.numeric)
  numeric_cols <- names(numeric_cols)[numeric_cols]

  summary_list$overall <- space_use_df %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(numeric_cols),
                    list(mean = ~mean(.x, na.rm = TRUE),
                         median = ~median(.x, na.rm = TRUE),
                         min = ~min(.x, na.rm = TRUE),
                         max = ~max(.x, na.rm = TRUE),
                         sd = ~sd(.x, na.rm = TRUE)),
                    .names = "{.col}_{.fn}"),
      .groups = 'drop'
    )

  # By-group summaries if grouping variables exist
  if (length(grouping_vars) > 0) {
    for (var in grouping_vars) {
      if (var %in% names(space_use_df)) {

        # Create summary with only available columns
        summary_data <- space_use_df %>%
          dplyr::group_by(!!dplyr::sym(var)) %>%
          dplyr::summarise(
            n_groups = dplyr::n(),
            mean_n_points = mean(n_points, na.rm = TRUE),
            .groups = 'drop'
          )

        # Add method-specific summaries only if columns exist
        available_area_cols <- intersect(
          c("convex_hull_area_hectares", "bounding_box_area_hectares",
            "grid_cell_count_area_hectares", "ellipse_95_area_hectares",
            "constrained_convex_hull_area_hectares"),
          names(space_use_df)
        )

        if (length(available_area_cols) > 0) {
          method_summaries <- space_use_df %>%
            dplyr::group_by(!!dplyr::sym(var)) %>%
            dplyr::summarise(
              dplyr::across(dplyr::all_of(available_area_cols),
                            ~mean(.x, na.rm = TRUE),
                            .names = "mean_{.col}"),
              .groups = 'drop'
            )

          # Combine summaries
          summary_data <- summary_data %>%
            dplyr::left_join(method_summaries, by = var)
        }

        summary_list[[paste0("by_", var)]] <- summary_data
      }
    }
  }

  return(summary_list)
}

#' Create grid raster data using actual raster cell boundaries
create_grid_raster_data <- function(track_points, reference_raster = NULL) {

  x_coords <- track_points$x
  y_coords <- track_points$y

  if (length(x_coords) == 0) return(NULL)

  if (is.null(reference_raster)) {
    # Fallback to simple grid approach
    grid_resolution <- 100
    x_grid_centers <- floor(x_coords / grid_resolution) * grid_resolution + grid_resolution/2
    y_grid_centers <- floor(y_coords / grid_resolution) * grid_resolution + grid_resolution/2

    unique_cells <- unique(data.frame(
      x = x_grid_centers,
      y = y_grid_centers
    ))

    if(getOption("positionR.verbose", TRUE)) cat("Grid cells used (fallback):", nrow(unique_cells), "\n")
    return(unique_cells)
  }

  tryCatch({
    if (!requireNamespace("raster", quietly = TRUE)) {
      stop("raster package needed when reference_raster is provided")
    }

    if(getOption("positionR.verbose", TRUE)) cat("Finding depth_raster cells that contain track points\n")

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

    # Get the actual cell boundaries from the raster
    cell_data <- get_raster_cells_as_polygons(reference_raster, unique_cell_indices)

    if(getOption("positionR.verbose", TRUE)) cat("Track points occupy", length(unique_cell_indices), "depth_raster cells\n")
    if(getOption("positionR.verbose", TRUE)) cat("Cell resolution:", raster::res(reference_raster), "units\n")

    return(cell_data)

  }, error = function(e) {
    warning("Error finding depth_raster cells: ", e$message)
    return(NULL)
  })
}

#' Create constrained hull by finding which depth_raster cells are within convex hull polygon
create_constrained_grid_raster_data <- function(track_points, reference_raster) {

  x_coords <- track_points$x
  y_coords <- track_points$y

  if (length(x_coords) < 3) return(NULL)

  if (is.null(reference_raster)) {
    warning("reference_raster is required for constrained hull")
    return(NULL)
  }

  tryCatch({
    if (!requireNamespace("raster", quietly = TRUE)) {
      stop("raster package needed when reference_raster is provided")
    }

    if(getOption("positionR.verbose", TRUE)) cat("Finding depth_raster cells within convex hull polygon\n")

    # Get convex hull points
    hull_indices <- grDevices::chull(x_coords, y_coords)
    hull_x <- x_coords[hull_indices]
    hull_y <- y_coords[hull_indices]

    # Get extent of convex hull with small buffer to ensure we don't miss edge cells
    hull_extent <- list(
      xmin = min(hull_x) - raster::res(reference_raster)[1],
      xmax = max(hull_x) + raster::res(reference_raster)[1],
      ymin = min(hull_y) - raster::res(reference_raster)[2],
      ymax = max(hull_y) + raster::res(reference_raster)[2]
    )

    # Create extent object for the hull area
    hull_bbox <- raster::extent(hull_extent$xmin, hull_extent$xmax,
                                hull_extent$ymin, hull_extent$ymax)

    # Get all cell numbers within the bounding box (including NA cells)
    all_cells_in_bbox <- raster::cellsFromExtent(reference_raster, hull_bbox)

    if (length(all_cells_in_bbox) == 0) {
      warning("No raster cells found within hull extent")
      return(NULL)
    }

    # Filter to only non-NA cells (cells with valid depth values)
    cell_values <- raster::extract(reference_raster, all_cells_in_bbox)
    valid_cells <- all_cells_in_bbox[!is.na(cell_values)]

    if(getOption("positionR.verbose", TRUE)) cat("Testing", length(valid_cells), "valid depth raster cells within hull extent\n")

    if (length(valid_cells) == 0) {
      warning("No valid (non-NA) cells in reference raster within hull extent")
      return(NULL)
    }

    # Test each cell to see if it intersects with the convex hull
    # We'll test multiple points within each cell for more accurate results
    cells_in_hull <- test_cells_in_polygon(reference_raster, valid_cells, hull_x, hull_y)

    # Get the cell indices that intersect with the hull
    hull_cell_indices <- valid_cells[cells_in_hull]

    if (length(hull_cell_indices) == 0) {
      warning("No depth raster cells intersect with convex hull")
      return(NULL)
    }

    if(getOption("positionR.verbose", TRUE)) cat("Convex hull intersects with", length(hull_cell_indices), "depth_raster cells\n")

    # Get the actual cell polygons for visualization
    cell_data <- get_raster_cells_as_polygons(reference_raster, hull_cell_indices)

    return(cell_data)

  }, error = function(e) {
    warning("Error finding cells within convex hull: ", e$message)
    return(NULL)
  })
}

#' Get raster cells as proper polygon data for accurate visualization
#' This function creates the actual cell boundaries for each raster cell
get_raster_cells_as_polygons <- function(raster_obj, cell_indices) {

  tryCatch({
    # Get raster resolution
    res_x <- raster::res(raster_obj)[1]
    res_y <- raster::res(raster_obj)[2]

    # Get cell centers
    cell_centers <- raster::xyFromCell(raster_obj, cell_indices)

    # Create polygon data for each cell
    polygon_list <- list()

    for (i in seq_along(cell_indices)) {
      center_x <- cell_centers[i, 1]
      center_y <- cell_centers[i, 2]

      # Calculate exact cell boundaries
      x_min <- center_x - res_x/2
      x_max <- center_x + res_x/2
      y_min <- center_y - res_y/2
      y_max <- center_y + res_y/2

      # Create polygon vertices (rectangle) - fix row names warning
      polygon_list[[i]] <- data.frame(
        cell_id = rep(cell_indices[i], 5),
        x = c(x_min, x_max, x_max, x_min, x_min),  # Close the polygon
        y = c(y_min, y_min, y_max, y_max, y_min),  # Close the polygon
        center_x = rep(center_x, 5),
        center_y = rep(center_y, 5),
        xmin = rep(x_min, 5),
        xmax = rep(x_max, 5),
        ymin = rep(y_min, 5),
        ymax = rep(y_max, 5),
        row.names = NULL
      )
    }

    # Combine all cell polygons
    if (length(polygon_list) > 0) {
      polygons_df <- dplyr::bind_rows(polygon_list)
      if(getOption("positionR.verbose", TRUE)) cat("Created polygon data for", length(unique(polygons_df$cell_id)), "raster cells\n")
      return(polygons_df)
    } else {
      return(NULL)
    }

  }, error = function(e) {
    warning("Error creating cell polygons: ", e$message)
    return(NULL)
  })
}

#' Check if a line segment intersects with a box
#' Uses basic line-box intersection algorithm
line_intersects_box <- function(x1, y1, x2, y2, xmin, ymin, xmax, ymax) {
  # Check if line segment is completely outside box
  if ((x1 < xmin && x2 < xmin) || (x1 > xmax && x2 > xmax) ||
      (y1 < ymin && y2 < ymin) || (y1 > ymax && y2 > ymax)) {
    return(FALSE)
  }
  
  # Check if either endpoint is inside the box
  if ((x1 >= xmin && x1 <= xmax && y1 >= ymin && y1 <= ymax) ||
      (x2 >= xmin && x2 <= xmax && y2 >= ymin && y2 <= ymax)) {
    return(TRUE)
  }
  
  # Check intersection with each box edge
  # Using parametric line equation: point = p1 + t * (p2 - p1)
  dx <- x2 - x1
  dy <- y2 - y1
  
  # Check intersection with vertical edges (x = xmin and x = xmax)
  if (dx != 0) {
    # Left edge (x = xmin)
    t <- (xmin - x1) / dx
    if (t >= 0 && t <= 1) {
      y <- y1 + t * dy
      if (y >= ymin && y <= ymax) return(TRUE)
    }
    # Right edge (x = xmax)
    t <- (xmax - x1) / dx
    if (t >= 0 && t <= 1) {
      y <- y1 + t * dy
      if (y >= ymin && y <= ymax) return(TRUE)
    }
  }
  
  # Check intersection with horizontal edges (y = ymin and y = ymax)
  if (dy != 0) {
    # Bottom edge (y = ymin)
    t <- (ymin - y1) / dy
    if (t >= 0 && t <= 1) {
      x <- x1 + t * dx
      if (x >= xmin && x <= xmax) return(TRUE)
    }
    # Top edge (y = ymax)
    t <- (ymax - y1) / dy
    if (t >= 0 && t <= 1) {
      x <- x1 + t * dx
      if (x >= xmin && x <= xmax) return(TRUE)
    }
  }
  
  return(FALSE)
}

#' Test if raster cells intersect with a polygon using simple corner test
#' Include cell if ANY of its 4 corners are inside the polygon
test_cells_in_polygon <- function(raster_obj, cell_indices, poly_x, poly_y) {

  tryCatch({
    # Get raster resolution
    res_x <- raster::res(raster_obj)[1]
    res_y <- raster::res(raster_obj)[2]

    # Get cell centers
    cell_centers <- raster::xyFromCell(raster_obj, cell_indices)

    # For each cell, test if any corner is inside the polygon
    cells_intersect <- logical(length(cell_indices))

    # Define the 4 corners of each cell relative to center
    corner_offsets <- data.frame(
      x_offset = c(-res_x/2, res_x/2, res_x/2, -res_x/2),  # left, right, right, left
      y_offset = c(-res_y/2, -res_y/2, res_y/2, res_y/2)   # bottom, bottom, top, top
    )

    for (i in seq_along(cell_indices)) {
      center_x <- cell_centers[i, 1]
      center_y <- cell_centers[i, 2]

      # Calculate the 4 corner coordinates
      corner_x <- center_x + corner_offsets$x_offset
      corner_y <- center_y + corner_offsets$y_offset

      # Test each corner individually
      corner_results <- logical(4)
      for (j in 1:4) {
        corner_results[j] <- point_in_polygon(corner_x[j], corner_y[j], poly_x, poly_y)[1]
      }
      
      # Include cell if ANY corner is inside
      cells_intersect[i] <- any(corner_results)
    }

    if(getOption("positionR.verbose", TRUE)) cat("Cell intersection test (corner test): ", sum(cells_intersect), " of ", length(cells_intersect), " cells have corners inside polygon\n")

    return(cells_intersect)

  }, error = function(e) {
    warning("Error testing cell intersection with polygon: ", e$message)
    return(rep(FALSE, length(cell_indices)))
  })
}

#' Updated grid cell usage calculation to work with reference raster
calculate_grid_cell_usage <- function(x_coords, y_coords, grid_resolution, reference_raster = NULL) {

  if (!is.null(reference_raster)) {
    # Use actual raster cells instead of arbitrary grid
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
        return(list(n_cells = 0, total_area_m2 = 0, cell_area_m2 = 0))
      }

      # Calculate actual cell area from raster
      res_x <- raster::res(reference_raster)[1]
      res_y <- raster::res(reference_raster)[2]
      cell_area_m2 <- res_x * res_y

      n_cells <- length(unique_cell_indices)
      total_area_m2 <- n_cells * cell_area_m2

      if(getOption("positionR.verbose", TRUE)) cat("Using reference raster: ", n_cells, " cells, cell area = ", cell_area_m2, " m²\n")

      return(list(
        n_cells = n_cells,
        total_area_m2 = total_area_m2,
        cell_area_m2 = cell_area_m2
      ))

    }, error = function(e) {
      warning("Error using reference raster, falling back to grid approach: ", e$message)
    })
  }

  # Fallback to original grid approach
  x_grid <- floor(x_coords / grid_resolution) * grid_resolution
  y_grid <- floor(y_coords / grid_resolution) * grid_resolution

  # Count unique grid cells
  unique_cells <- unique(paste(x_grid, y_grid, sep = "_"))
  n_cells <- length(unique_cells)

  # Calculate total area
  cell_area_m2 <- grid_resolution^2
  total_area_m2 <- n_cells * cell_area_m2

  return(list(
    n_cells = n_cells,
    total_area_m2 = total_area_m2,
    cell_area_m2 = cell_area_m2
  ))
}

#' Updated constrained convex hull calculation to work with reference raster
calculate_constrained_convex_hull_area <- function(x_coords, y_coords, grid_resolution, reference_raster = NULL) {

  if (length(x_coords) < 3) return(NA)

  if (!is.null(reference_raster)) {
    # Use actual raster cells with improved intersection testing
    tryCatch({
      if (!requireNamespace("raster", quietly = TRUE)) {
        stop("raster package needed when reference_raster is provided")
      }

      # Get convex hull points
      hull_indices <- grDevices::chull(x_coords, y_coords)
      hull_x <- x_coords[hull_indices]
      hull_y <- y_coords[hull_indices]

      # Get extent of convex hull with buffer to ensure edge cells are included
      res_x <- raster::res(reference_raster)[1]
      res_y <- raster::res(reference_raster)[2]

      hull_extent <- list(
        xmin = min(hull_x) - res_x,
        xmax = max(hull_x) + res_x,
        ymin = min(hull_y) - res_y,
        ymax = max(hull_y) + res_y
      )

      # Create extent object for the hull area
      hull_bbox <- raster::extent(hull_extent$xmin, hull_extent$xmax,
                                  hull_extent$ymin, hull_extent$ymax)

      # Get all cell numbers within the bounding box
      all_cells_in_bbox <- raster::cellsFromExtent(reference_raster, hull_bbox)

      if (length(all_cells_in_bbox) == 0) {
        warning("No raster cells found within hull extent")
        return(NA)
      }

      # Filter to only non-NA cells (cells with valid depth values)
      cell_values <- raster::extract(reference_raster, all_cells_in_bbox)
      valid_cells <- all_cells_in_bbox[!is.na(cell_values)]

      if (length(valid_cells) == 0) {
        warning("No valid (non-NA) cells in reference raster within hull extent")
        return(NA)
      }

      # Test each cell to see if it intersects with the convex hull
      cells_intersect <- test_cells_in_polygon(reference_raster, valid_cells, hull_x, hull_y)

      # Get the cell indices that intersect with the hull
      hull_cell_indices <- valid_cells[cells_intersect]

      if (length(hull_cell_indices) == 0) {
        warning("No depth raster cells intersect with convex hull")
        return(NA)
      }

      # Calculate actual cell area from raster
      cell_area_m2 <- res_x * res_y

      n_cells_in_hull <- length(hull_cell_indices)
      constrained_area <- n_cells_in_hull * cell_area_m2

      if(getOption("positionR.verbose", TRUE)) cat("Constrained hull (raster): ", n_cells_in_hull, " cells, area = ", constrained_area, " m²\n")

      return(constrained_area)

    }, error = function(e) {
      warning("Error using reference raster for constrained hull, falling back to grid approach: ", e$message)
    })
  }

  # Fallback to original grid approach
  tryCatch({
    # Get convex hull points
    hull_indices <- grDevices::chull(x_coords, y_coords)
    hull_x <- x_coords[hull_indices]
    hull_y <- y_coords[hull_indices]

    # Find the bounding box of the convex hull
    hull_x_range <- range(hull_x)
    hull_y_range <- range(hull_y)

    # Create a grid covering the entire convex hull area
    x_grid_seq <- seq(from = floor(hull_x_range[1] / grid_resolution) * grid_resolution,
                      to = ceiling(hull_x_range[2] / grid_resolution) * grid_resolution,
                      by = grid_resolution)
    y_grid_seq <- seq(from = floor(hull_y_range[1] / grid_resolution) * grid_resolution,
                      to = ceiling(hull_y_range[2] / grid_resolution) * grid_resolution,
                      by = grid_resolution)

    # Create all possible grid cell combinations
    grid_combinations <- expand.grid(x_grid = x_grid_seq, y_grid = y_grid_seq)

    # Calculate grid cell centers
    cell_centers_x <- grid_combinations$x_grid + grid_resolution/2
    cell_centers_y <- grid_combinations$y_grid + grid_resolution/2

    # Test which grid cell centers are inside the convex hull
    cells_in_hull <- point_in_polygon(cell_centers_x, cell_centers_y, hull_x, hull_y)

    # Calculate area of cells within hull
    n_cells_in_hull <- sum(cells_in_hull)
    cell_area_m2 <- grid_resolution^2
    constrained_area <- n_cells_in_hull * cell_area_m2

    if(getOption("positionR.verbose", TRUE)) cat("Constrained hull (grid): ", n_cells_in_hull, " cells, area = ", constrained_area, " m²\n")

    return(constrained_area)

  }, error = function(e) {
    warning("Error calculating constrained convex hull: ", e$message)
    return(NA)
  })
}

#' Plot space use estimates
plot_space_use <- function(space_use_results,
                           method_compare = NULL,
                           plot_type = "comparison",
                           track_data = NULL,
                           fish_select = NULL,
                           time_select = NULL,
                           method_select = "convex_hull",
                           background_raster = NULL,
                           time_aggregation = "day",
                           point_type = "track") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for plotting")
  }

  if (plot_type == "map") {
    return(plot_space_use_map(space_use_results, track_data, fish_select,
                              time_select, method_select, background_raster,
                              time_aggregation, point_type))
  }

  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' needed for data reshaping")
  }

  space_use_df <- space_use_results$space_use_estimates

  # Auto-detect available methods if not specified
  if (is.null(method_compare)) {
    area_cols <- names(space_use_df)[grepl("_area_hectares$", names(space_use_df))]
    method_compare <- gsub("_area_hectares$", "", area_cols)
    if(getOption("positionR.verbose", TRUE)) cat("Auto-detected methods:", paste(method_compare, collapse = ", "), "\n")
  }

  # Create column names for the specified methods
  area_cols <- paste0(method_compare, "_area_hectares")
  available_cols <- intersect(area_cols, names(space_use_df))

  if (length(available_cols) == 0) {
    stop("No area columns found for specified methods")
  }

  if (plot_type == "comparison") {
    # Reshape data for comparison
    id_cols <- intersect(c("fish_id", "time_period", "group_id"), names(space_use_df))

    long_data <- space_use_df %>%
      dplyr::select(dplyr::all_of(c(id_cols, available_cols))) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(available_cols),
                          names_to = "method",
                          values_to = "area_hectares") %>%
      dplyr::mutate(method = gsub("_area_hectares", "", method)) %>%
      dplyr::filter(!is.na(area_hectares))

    p <- ggplot2::ggplot(long_data, ggplot2::aes(x = method, y = area_hectares)) +
      ggplot2::geom_boxplot(alpha = 0.7) +
      ggplot2::geom_jitter(width = 0.2, alpha = 0.6) +
      ggplot2::scale_y_log10() +
      ggplot2::labs(
        title = "Space Use Estimates by Method",
        x = "Method",
        y = "Area (hectares, log scale)"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  } else if (plot_type == "by_fish" && "fish_id" %in% names(space_use_df)) {
    method_col <- available_cols[1]
    method_name <- gsub("_area_hectares", "", method_col)

    p <- ggplot2::ggplot(space_use_df, ggplot2::aes(x = factor(fish_id), y = !!ggplot2::sym(method_col))) +
      ggplot2::geom_col(alpha = 0.7) +
      ggplot2::labs(
        title = paste("Space Use by Fish (", method_name, ")"),
        x = "Fish ID",
        y = "Area (hectares)"
      ) +
      ggplot2::theme_minimal()
  }

  return(p)
}

#' Updated plot function with proper raster cell visualization
plot_space_use_map <- function(space_use_results, track_data, fish_select,
                               time_select, method_select, background_raster,
                               time_aggregation, point_type = "track") {

  if (is.null(track_data)) {
    stop("track_data is required for map plots")
  }

  space_use_df <- space_use_results$space_use_estimates
  track_std <- standardize_track_data_space(track_data, time_aggregation)

  # Select fish if not specified
  available_fish <- unique(space_use_df$fish_id)
  if (is.null(fish_select)) {
    fish_select <- available_fish[1]
  }

  # Filter data - handle both numeric and date-based time selection
  if (!is.null(time_select)) {
    if (is.character(time_select) && grepl("\\d{4}-\\d{2}-\\d{2}", time_select) && 
        "time_period_label" %in% names(space_use_df)) {
      # Date-based filtering
      space_use_filtered <- space_use_df %>%
        dplyr::filter(fish_id == fish_select & time_period_label == time_select)
      track_filtered <- track_std %>%
        dplyr::filter(fish_id == fish_select & time_period_label == time_select)
      title_suffix <- paste("Fish", fish_select, "- Time Period", time_select)
    } else {
      # Numeric filtering
      space_use_filtered <- space_use_df %>%
        dplyr::filter(fish_id == fish_select & time_period == time_select)
      track_filtered <- track_std %>%
        dplyr::filter(fish_id == fish_select & time_period == time_select)
      title_suffix <- paste("Fish", fish_select, "- Time Period", time_select)
    }
  } else {
    space_use_filtered <- space_use_df %>%
      dplyr::filter(fish_id == fish_select)
    track_filtered <- track_std %>%
      dplyr::filter(fish_id == fish_select)
    title_suffix <- paste("Fish", fish_select, "- All Time Periods")
  }

  if (nrow(track_filtered) < 3) {
    stop("Insufficient track points for mapping")
  }

  # Create base plot
  p <- ggplot2::ggplot()

  # Add background raster if provided
  if (!is.null(background_raster)) {
    if ("RasterLayer" %in% class(background_raster)) {
      if (requireNamespace("raster", quietly = TRUE)) {
        raster_df <- raster::as.data.frame(background_raster, xy = TRUE)
        names(raster_df)[3] <- "values"
        raster_df <- raster_df[!is.na(raster_df$values), ]

        p <- p +
          ggplot2::geom_raster(
            data = raster_df,
            ggplot2::aes(x = x, y = y, alpha = values),
            fill = "lightblue"
          ) +
          ggplot2::scale_alpha_continuous(range = c(0.1, 0.4), guide = "none")
      }
    }
  }

  # Add space use visualization using pre-calculated spatial data
  spatial_data <- space_use_results$spatial_data
  
  if (!is.null(spatial_data) && method_select %in% names(spatial_data)) {
    method_data <- spatial_data[[method_select]]
    
    # Filter spatial data to match the selected fish and time
    if (!is.null(method_data) && nrow(method_data) > 0) {
      # Filter by fish if specified
      if (!is.null(fish_select) && "fish_id" %in% names(method_data)) {
        method_data <- method_data %>% dplyr::filter(fish_id == fish_select)
      }
      
      # Filter by time if specified
      if (!is.null(time_select) && "time_period_label" %in% names(method_data)) {
        method_data <- method_data %>% dplyr::filter(time_period_label == time_select)
      }
      
      if (nrow(method_data) > 0) {
        if (method_select %in% c("grid_cell_count", "constrained_convex_hull")) {
          # Plot grid cells as polygons
          if ("cell_id" %in% names(method_data)) {
            p <- p +
              ggplot2::geom_polygon(
                data = method_data,
                ggplot2::aes(x = x, y = y, group = cell_id),
                fill = "red",
                alpha = if(method_select == "grid_cell_count") 0.4 else 0.3,
                color = "darkred",
                size = if(method_select == "grid_cell_count") 0.3 else 0.2
              )
          } else {
            # Fallback to regular polygon
            p <- p +
              ggplot2::geom_polygon(
                data = method_data,
                ggplot2::aes(x = x, y = y),
                fill = "red",
                alpha = 0.3,
                color = "darkred",
                size = 1
              )
          }
          
          # Add convex hull outline for constrained hull
          if (method_select == "constrained_convex_hull" && 
              "constrained_convex_hull_outline" %in% names(spatial_data)) {
            hull_data <- spatial_data$constrained_convex_hull_outline
            if (!is.null(fish_select) && "fish_id" %in% names(hull_data)) {
              hull_data <- hull_data %>% dplyr::filter(fish_id == fish_select)
            }
            if (!is.null(time_select) && "time_period_label" %in% names(hull_data)) {
              hull_data <- hull_data %>% dplyr::filter(time_period_label == time_select)
            }
            
            if (!is.null(hull_data) && nrow(hull_data) > 0) {
              p <- p +
                ggplot2::geom_polygon(
                  data = hull_data,
                  ggplot2::aes(x = x, y = y),
                  fill = NA,
                  color = "darkred",
                  size = 2,
                  alpha = 0.8
                )
            }
          }
          
        } else {
          # Plot other methods as regular polygons
          p <- p +
            ggplot2::geom_polygon(
              data = method_data,
              ggplot2::aes(x = x, y = y),
              fill = "red",
              alpha = 0.3,
              color = "red",
              size = 1
            )
        }
      } else {
        warning("No spatial data found for selected fish/time combination")
      }
    } else {
      warning("No spatial data available for method: ", method_select)
    }
  } else {
    stop("No spatial data available in space_use_results. Please ensure calculate_space_use was run with appropriate methods.")
  }

  # Add track or points (green) based on point_type
  if (point_type == "track" && nrow(track_filtered) > 1) {
    # Plot as connected track
    p <- p + ggplot2::geom_path(
      data = track_filtered,
      ggplot2::aes(x = x, y = y),
      color = "green", size = 1.5
    )
    p <- p + ggplot2::geom_point(
      data = track_filtered,
      ggplot2::aes(x = x, y = y),
      color = "green", size = 2
    )
  } else if (point_type == "point") {
    # Plot as points only (no connecting lines)
    p <- p + ggplot2::geom_point(
      data = track_filtered,
      ggplot2::aes(x = x, y = y),
      color = "green", size = 2.5, alpha = 0.7
    )
  } else if (point_type == "track") {
    # Single point case for track
    p <- p + ggplot2::geom_point(
      data = track_filtered,
      ggplot2::aes(x = x, y = y),
      color = "green", size = 2
    )
  }

  # Calculate zoom considering both track data and space use spatial data
  current_spatial_data <- NULL
  if (!is.null(spatial_data) && method_select %in% names(spatial_data)) {
    method_data <- spatial_data[[method_select]]
    if (!is.null(method_data) && nrow(method_data) > 0) {
      # Filter spatial data to match the selected fish and time
      if (!is.null(fish_select) && "fish_id" %in% names(method_data)) {
        method_data <- method_data %>% dplyr::filter(fish_id == fish_select)
      }
      if (!is.null(time_select) && "time_period_label" %in% names(method_data)) {
        method_data <- method_data %>% dplyr::filter(time_period_label == time_select)
      }
      if (nrow(method_data) > 0) {
        current_spatial_data <- method_data
      }
    }
  }
  
  zoom_extent <- calculate_space_use_zoom(
    track_points = track_filtered,
    space_use_spatial_data = current_spatial_data,
    method = method_select
  )

  p <- p +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste("Space Use Map -", stringr::str_to_title(gsub("_", " ", method_select))),
      subtitle = paste(title_suffix, "| Green =", ifelse(point_type == "track", "track", "points"), ", Red = space use"),
      x = "X Coordinate (m)",
      y = "Y Coordinate (m)"
    ) +
    ggplot2::coord_fixed(
      xlim = zoom_extent$xlim,
      ylim = zoom_extent$ylim
    )

  return(p)
}

#' Create space use polygon for visualization
create_space_use_polygon <- function(track_points, method) {

  x_coords <- track_points$x
  y_coords <- track_points$y

  if (length(x_coords) < 3) return(NULL)

  tryCatch({
    if (method == "convex_hull" || method == "mcp") {
      # Convex hull
      hull_indices <- grDevices::chull(x_coords, y_coords)
      return(data.frame(x = x_coords[hull_indices], y = y_coords[hull_indices]))

    } else if (method == "bounding_box") {
      # Bounding rectangle
      x_range <- range(x_coords)
      y_range <- range(y_coords)
      return(data.frame(
        x = c(x_range[1], x_range[2], x_range[2], x_range[1], x_range[1]),
        y = c(y_range[1], y_range[1], y_range[2], y_range[2], y_range[1])
      ))

    } else if (method == "95_ellipse") {
      # 95% confidence ellipse
      return(create_ellipse_polygon(x_coords, y_coords))

    } else {
      warning("Unknown method: ", method, ". Using convex hull.")
      hull_indices <- grDevices::chull(x_coords, y_coords)
      return(data.frame(x = x_coords[hull_indices], y = y_coords[hull_indices]))
    }

  }, error = function(e) {
    warning("Error creating polygon for method ", method, ": ", e$message)
    return(NULL)
  })
}

#' Create ellipse polygon for 95% confidence ellipse
create_ellipse_polygon <- function(x_coords, y_coords, n_points = 100) {

  if (length(x_coords) < 3) return(NULL)

  tryCatch({
    # Center coordinates
    center_x <- mean(x_coords)
    center_y <- mean(y_coords)

    # Center the coordinates
    x_centered <- x_coords - center_x
    y_centered <- y_coords - center_y

    # Calculate covariance matrix
    cov_matrix <- cov(cbind(x_centered, y_centered))

    # Eigenvalues and eigenvectors
    eigen_result <- eigen(cov_matrix)
    eigenvalues <- eigen_result$values
    eigenvectors <- eigen_result$vectors

    # 95% confidence level
    chi_sq_95 <- qchisq(0.95, df = 2)

    # Semi-axes lengths
    a <- sqrt(chi_sq_95 * eigenvalues[1])  # Semi-major axis
    b <- sqrt(chi_sq_95 * eigenvalues[2])  # Semi-minor axis

    # Rotation angle
    angle <- atan2(eigenvectors[2, 1], eigenvectors[1, 1])

    # Create ellipse points
    theta <- seq(0, 2*pi, length.out = n_points)
    ellipse_x <- a * cos(theta)
    ellipse_y <- b * sin(theta)

    # Rotate and translate
    cos_angle <- cos(angle)
    sin_angle <- sin(angle)

    rotated_x <- ellipse_x * cos_angle - ellipse_y * sin_angle + center_x
    rotated_y <- ellipse_x * sin_angle + ellipse_y * cos_angle + center_y

    return(data.frame(x = rotated_x, y = rotated_y))

  }, error = function(e) {
    warning("Error creating ellipse polygon: ", e$message)
    return(NULL)
  })
}

#' Create grid cell polygon
create_grid_polygon <- function(x_coords, y_coords, grid_resolution) {

  tryCatch({
    # Find unique grid cells
    x_grid <- floor(x_coords / grid_resolution) * grid_resolution
    y_grid <- floor(y_coords / grid_resolution) * grid_resolution

    unique_cells <- unique(data.frame(x_grid = x_grid, y_grid = y_grid))

    # Create rectangle for each cell
    polygon_list <- list()

    for (i in 1:nrow(unique_cells)) {
      x_min <- unique_cells$x_grid[i]
      x_max <- x_min + grid_resolution
      y_min <- unique_cells$y_grid[i]
      y_max <- y_min + grid_resolution

      polygon_list[[i]] <- data.frame(
        x = c(x_min, x_max, x_max, x_min, x_min),
        y = c(y_min, y_min, y_max, y_max, y_min),
        cell_id = i
      )
    }

    # Combine all cells
    if (length(polygon_list) > 0) {
      return(dplyr::bind_rows(polygon_list))
    } else {
      return(NULL)
    }

  }, error = function(e) {
    warning("Error creating grid polygon: ", e$message)
    return(NULL)
  })
}

#' Calculate zoom extent for map considering both track data and space use area
calculate_space_use_zoom <- function(track_points = NULL, space_use_spatial_data = NULL, method = NULL, buffer_factor = 0.15) {

  all_x <- c()
  all_y <- c()
  
  # Include track points if provided
  if (!is.null(track_points) && nrow(track_points) > 0) {
    all_x <- c(all_x, track_points$x)
    all_y <- c(all_y, track_points$y)
  }
  
  # Include space use spatial data if provided
  if (!is.null(space_use_spatial_data) && nrow(space_use_spatial_data) > 0) {
    all_x <- c(all_x, space_use_spatial_data$x)
    all_y <- c(all_y, space_use_spatial_data$y)
  }
  
  # If no data available, return default extent
  if (length(all_x) == 0 || length(all_y) == 0) {
    warning("No data available for zoom calculation, using default extent")
    return(list(
      xlim = c(0, 1000),
      ylim = c(0, 1000)
    ))
  }

  extent_x <- range(all_x, na.rm = TRUE)
  extent_y <- range(all_y, na.rm = TRUE)

  x_buffer <- max(diff(extent_x) * buffer_factor, 50)
  y_buffer <- max(diff(extent_y) * buffer_factor, 50)

  max_buffer <- max(x_buffer, y_buffer)

  return(list(
    xlim = c(extent_x[1] - max_buffer, extent_x[2] + max_buffer),
    ylim = c(extent_y[1] - max_buffer, extent_y[2] + max_buffer)
  ))
}

#' Create space use polygon data for direct plotting
create_space_use_polygon_data <- function(x_coords, y_coords, method) {
  
  if (length(x_coords) < 3) return(NULL)
  
  tryCatch({
    if (method == "convex_hull" || method == "mcp") {
      # Convex hull
      hull_indices <- grDevices::chull(x_coords, y_coords)
      return(data.frame(
        x = x_coords[hull_indices], 
        y = y_coords[hull_indices],
        method = method
      ))
      
    } else if (method == "bounding_box") {
      # Bounding rectangle
      x_range <- range(x_coords)
      y_range <- range(y_coords)
      return(data.frame(
        x = c(x_range[1], x_range[2], x_range[2], x_range[1], x_range[1]),
        y = c(y_range[1], y_range[1], y_range[2], y_range[2], y_range[1]),
        method = method
      ))
      
    } else if (method == "95_ellipse") {
      # 95% confidence ellipse
      return(create_ellipse_polygon_data(x_coords, y_coords, method))
      
    } else {
      warning("Unknown method: ", method, ". Using convex hull.")
      hull_indices <- grDevices::chull(x_coords, y_coords)
      return(data.frame(
        x = x_coords[hull_indices], 
        y = y_coords[hull_indices],
        method = method
      ))
    }
    
  }, error = function(e) {
    warning("Error creating polygon data for method ", method, ": ", e$message)
    return(NULL)
  })
}

#' Create ellipse polygon data with method label
create_ellipse_polygon_data <- function(x_coords, y_coords, method, n_points = 100) {
  
  if (length(x_coords) < 3) return(NULL)
  
  tryCatch({
    # Center coordinates
    center_x <- mean(x_coords)
    center_y <- mean(y_coords)
    
    # Center the coordinates
    x_centered <- x_coords - center_x
    y_centered <- y_coords - center_y
    
    # Calculate covariance matrix
    cov_matrix <- cov(cbind(x_centered, y_centered))
    
    # Eigenvalues and eigenvectors
    eigen_result <- eigen(cov_matrix)
    eigenvalues <- eigen_result$values
    eigenvectors <- eigen_result$vectors
    
    # 95% confidence level
    chi_sq_95 <- qchisq(0.95, df = 2)
    
    # Semi-axes lengths
    a <- sqrt(chi_sq_95 * eigenvalues[1])  # Semi-major axis
    b <- sqrt(chi_sq_95 * eigenvalues[2])  # Semi-minor axis
    
    # Rotation angle
    angle <- atan2(eigenvectors[2, 1], eigenvectors[1, 1])
    
    # Create ellipse points
    theta <- seq(0, 2*pi, length.out = n_points)
    ellipse_x <- a * cos(theta)
    ellipse_y <- b * sin(theta)
    
    # Rotate and translate
    cos_angle <- cos(angle)
    sin_angle <- sin(angle)
    
    rotated_x <- ellipse_x * cos_angle - ellipse_y * sin_angle + center_x
    rotated_y <- ellipse_x * sin_angle + ellipse_y * cos_angle + center_y
    
    return(data.frame(
      x = rotated_x, 
      y = rotated_y,
      method = method
    ))
    
  }, error = function(e) {
    warning("Error creating ellipse polygon data: ", e$message)
    return(NULL)
  })
}

#' Calculate grid cell usage with spatial data
calculate_grid_cell_usage_with_data <- function(x_coords, y_coords, grid_resolution, reference_raster = NULL) {
  
  if (!is.null(reference_raster)) {
    # Use actual raster cells instead of arbitrary grid
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
        return(list(n_cells = 0, total_area_m2 = 0, cell_area_m2 = 0, spatial_data = NULL))
      }
      
      # Calculate actual cell area from raster
      res_x <- raster::res(reference_raster)[1]
      res_y <- raster::res(reference_raster)[2]
      cell_area_m2 <- res_x * res_y
      
      n_cells <- length(unique_cell_indices)
      total_area_m2 <- n_cells * cell_area_m2
      
      # Create spatial data for plotting
      spatial_data <- get_raster_cells_as_polygons(reference_raster, unique_cell_indices)
      if (!is.null(spatial_data)) {
        spatial_data$method <- "grid_cell_count"
      }
      
      if(getOption("positionR.verbose", TRUE)) cat("Using reference raster: ", n_cells, " cells, cell area = ", cell_area_m2, " m²\n")
      
      return(list(
        n_cells = n_cells,
        total_area_m2 = total_area_m2,
        cell_area_m2 = cell_area_m2,
        spatial_data = spatial_data
      ))
      
    }, error = function(e) {
      warning("Error using reference raster, falling back to grid approach: ", e$message)
    })
  }
  
  # Fallback to original grid approach
  x_grid <- floor(x_coords / grid_resolution) * grid_resolution
  y_grid <- floor(y_coords / grid_resolution) * grid_resolution
  
  # Count unique grid cells
  unique_cells <- unique(paste(x_grid, y_grid, sep = "_"))
  n_cells <- length(unique_cells)
  
  # Calculate total area
  cell_area_m2 <- grid_resolution^2
  total_area_m2 <- n_cells * cell_area_m2
  
  # Create spatial data for fallback grid
  unique_grid_coords <- unique(data.frame(
    x_grid = x_grid + grid_resolution/2,
    y_grid = y_grid + grid_resolution/2
  ))
  
  # Create polygon data for each grid cell
  spatial_data <- create_grid_cell_polygons(unique_grid_coords$x_grid, unique_grid_coords$y_grid, grid_resolution)
  if (!is.null(spatial_data)) {
    spatial_data$method <- "grid_cell_count"
  }
  
  return(list(
    n_cells = n_cells,
    total_area_m2 = total_area_m2,
    cell_area_m2 = cell_area_m2,
    spatial_data = spatial_data
  ))
}

#' Calculate constrained convex hull with spatial data
calculate_constrained_convex_hull_area_with_data <- function(x_coords, y_coords, grid_resolution, reference_raster = NULL) {
  
  if (length(x_coords) < 3) return(list(area = NA, spatial_data = NULL))
  
  if (!is.null(reference_raster)) {
    # Use actual raster cells with improved intersection testing
    tryCatch({
      if (!requireNamespace("raster", quietly = TRUE)) {
        stop("raster package needed when reference_raster is provided")
      }
      
      # Get convex hull points
      hull_indices <- grDevices::chull(x_coords, y_coords)
      hull_x <- x_coords[hull_indices]
      hull_y <- y_coords[hull_indices]
      
      # Get extent of convex hull with buffer to ensure edge cells are included
      res_x <- raster::res(reference_raster)[1]
      res_y <- raster::res(reference_raster)[2]
      
      hull_extent <- list(
        xmin = min(hull_x) - res_x,
        xmax = max(hull_x) + res_x,
        ymin = min(hull_y) - res_y,
        ymax = max(hull_y) + res_y
      )
      
      # Create extent object for the hull area
      hull_bbox <- raster::extent(hull_extent$xmin, hull_extent$xmax,
                                  hull_extent$ymin, hull_extent$ymax)
      
      # Get all cell numbers within the bounding box
      all_cells_in_bbox <- raster::cellsFromExtent(reference_raster, hull_bbox)
      
      if (length(all_cells_in_bbox) == 0) {
        warning("No raster cells found within hull extent")
        return(list(area = NA, spatial_data = NULL))
      }
      
      # Filter to only non-NA cells (cells with valid depth values)
      cell_values <- raster::extract(reference_raster, all_cells_in_bbox)
      valid_cells <- all_cells_in_bbox[!is.na(cell_values)]
      
      if (length(valid_cells) == 0) {
        warning("No valid (non-NA) cells in reference raster within hull extent")
        return(list(area = NA, spatial_data = NULL))
      }
      
      # Test each cell to see if it intersects with the convex hull
      cells_intersect <- test_cells_in_polygon(reference_raster, valid_cells, hull_x, hull_y)
      
      # Get the cell indices that intersect with the hull
      hull_cell_indices <- valid_cells[cells_intersect]
      
      if (length(hull_cell_indices) == 0) {
        warning("No depth raster cells intersect with convex hull")
        return(list(area = NA, spatial_data = NULL))
      }
      
      # Calculate actual cell area from raster
      cell_area_m2 <- res_x * res_y
      
      n_cells_in_hull <- length(hull_cell_indices)
      constrained_area <- n_cells_in_hull * cell_area_m2
      
      # Create spatial data for plotting
      spatial_data <- get_raster_cells_as_polygons(reference_raster, hull_cell_indices)
      if (!is.null(spatial_data)) {
        spatial_data$method <- "constrained_convex_hull"
      }
      
      # Create hull outline data for plotting
      hull_outline <- data.frame(
        x = c(hull_x, hull_x[1]),  # Close the polygon
        y = c(hull_y, hull_y[1]),  # Close the polygon
        method = "convex_hull_outline"
      )
      
      if(getOption("positionR.verbose", TRUE)) cat("Constrained hull (raster): ", n_cells_in_hull, " cells, area = ", constrained_area, " m²\n")
      
      return(list(
        area = constrained_area,
        spatial_data = spatial_data,
        hull_outline = hull_outline
      ))
      
    }, error = function(e) {
      warning("Error using reference raster for constrained hull, falling back to grid approach: ", e$message)
    })
  }
  
  # Fallback to original grid approach
  tryCatch({
    # Get convex hull points
    hull_indices <- grDevices::chull(x_coords, y_coords)
    hull_x <- x_coords[hull_indices]
    hull_y <- y_coords[hull_indices]
    
    # Find the bounding box of the convex hull
    hull_x_range <- range(hull_x)
    hull_y_range <- range(hull_y)
    
    # Create a grid covering the entire convex hull area
    x_grid_seq <- seq(from = floor(hull_x_range[1] / grid_resolution) * grid_resolution,
                      to = ceiling(hull_x_range[2] / grid_resolution) * grid_resolution,
                      by = grid_resolution)
    y_grid_seq <- seq(from = floor(hull_y_range[1] / grid_resolution) * grid_resolution,
                      to = ceiling(hull_y_range[2] / grid_resolution) * grid_resolution,
                      by = grid_resolution)
    
    # Create all possible grid cell combinations
    grid_combinations <- expand.grid(x_grid = x_grid_seq, y_grid = y_grid_seq)
    
    # Calculate grid cell centers
    cell_centers_x <- grid_combinations$x_grid + grid_resolution/2
    cell_centers_y <- grid_combinations$y_grid + grid_resolution/2
    
    # Test which grid cell centers are inside the convex hull
    cells_in_hull <- point_in_polygon(cell_centers_x, cell_centers_y, hull_x, hull_y)
    
    # Calculate area of cells within hull
    n_cells_in_hull <- sum(cells_in_hull)
    cell_area_m2 <- grid_resolution^2
    constrained_area <- n_cells_in_hull * cell_area_m2
    
    # Create spatial data for cells within hull
    hull_cell_centers <- data.frame(
      x_grid = cell_centers_x[cells_in_hull],
      y_grid = cell_centers_y[cells_in_hull]
    )
    
    spatial_data <- create_grid_cell_polygons(hull_cell_centers$x_grid, hull_cell_centers$y_grid, grid_resolution)
    if (!is.null(spatial_data)) {
      spatial_data$method <- "constrained_convex_hull"
    }
    
    # Create hull outline data for plotting
    hull_outline <- data.frame(
      x = c(hull_x, hull_x[1]),  # Close the polygon
      y = c(hull_y, hull_y[1]),  # Close the polygon
      method = "convex_hull_outline"
    )
    
    if(getOption("positionR.verbose", TRUE)) cat("Constrained hull (grid): ", n_cells_in_hull, " cells, area = ", constrained_area, " m²\n")
    
    return(list(
      area = constrained_area,
      spatial_data = spatial_data,
      hull_outline = hull_outline
    ))
    
  }, error = function(e) {
    warning("Error calculating constrained convex hull: ", e$message)
    return(list(area = NA, spatial_data = NULL))
  })
}

#' Create grid cell polygons for plotting
create_grid_cell_polygons <- function(center_x, center_y, grid_resolution) {
  
  if (length(center_x) == 0) return(NULL)
  
  tryCatch({
    polygon_list <- list()
    
    for (i in seq_along(center_x)) {
      x_min <- center_x[i] - grid_resolution/2
      x_max <- center_x[i] + grid_resolution/2
      y_min <- center_y[i] - grid_resolution/2
      y_max <- center_y[i] + grid_resolution/2
      
      polygon_list[[i]] <- data.frame(
        x = c(x_min, x_max, x_max, x_min, x_min),
        y = c(y_min, y_min, y_max, y_max, y_min),
        cell_id = i,
        center_x = rep(center_x[i], 5),
        center_y = rep(center_y[i], 5)
      )
    }
    
    if (length(polygon_list) > 0) {
      return(dplyr::bind_rows(polygon_list))
    } else {
      return(NULL)
    }
    
  }, error = function(e) {
    warning("Error creating grid cell polygons: ", e$message)
    return(NULL)
  })
}
