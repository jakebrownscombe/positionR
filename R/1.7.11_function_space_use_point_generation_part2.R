# Space Use Point Generation Functions - Part 2
# Plotting and summary functions

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

#' Plot space use areas with generated random points
#'
#' @param space_use_points Output from generate_random_points_in_space_use
#' @param track_data Original track data
#' @param space_use_results Original space use results
#' @param fish_select Which fish to plot
#' @param time_select Which time period to plot
#' @param background_raster Optional background raster
#' @param reference_raster Optional reference raster for grid visualization
#' @param show_space_use Whether to show space use area outline/cells
#' @param show_track Whether to show original track
#' @param show_points Whether to show generated points
#' 
#' @return A ggplot2 object showing space use areas with generated points
#' 
#' @export
plot_space_use_points <- function(space_use_points,
                                  track_data,
                                  space_use_results,
                                  fish_select = NULL,
                                  time_select = NULL,
                                  background_raster = NULL,
                                  reference_raster = NULL,
                                  show_space_use = TRUE,
                                  show_track = TRUE,
                                  show_points = TRUE) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for plotting")
  }
  
  
  # Select fish and time if not specified
  if (is.null(fish_select)) {
    if ("fish_id" %in% names(space_use_points)) {
      fish_select <- unique(space_use_points$fish_id)[1]
    } else if ("path_id" %in% names(space_use_points)) {
      fish_select <- unique(space_use_points$path_id)[1]
    } else {
      stop("No fish_id or path_id column found in space_use_points")
    }
  }
  if (is.null(time_select)) {
    if ("time_period_numeric" %in% names(space_use_points)) {
      time_select <- unique(space_use_points$time_period_numeric)[1]
    } else if ("time_period_date" %in% names(space_use_points)) {
      time_select <- unique(space_use_points$time_period_date)[1]
    } else if ("time_period" %in% names(space_use_points)) {
      time_select <- unique(space_use_points$time_period)[1]
    } else {
      # No time aggregation used - use all data for this fish
      time_select <- "no_time_aggregation"
    }
  }
  
  # Store original time_select BEFORE any conversions happen
  original_time_select <- time_select
  
  # Determine fish ID column name and filter points
  fish_col <- if ("fish_id" %in% names(space_use_points)) "fish_id" else "path_id"
  
  # Filter points - simple approach: just use first available time period if none match
  plot_points <- space_use_points %>%
    dplyr::filter(!!rlang::sym(fish_col) == fish_select)
  
  # If no points for this fish, skip
  if (nrow(plot_points) == 0) {
    stop("No points found for fish ", fish_select)
  }
  
  # Skip time filtering if no time aggregation was used
  if (time_select == "no_time_aggregation") {
    # Use all points for this fish, no time filtering needed
  } else {
    # If time_select specified, try to filter, otherwise use first available
    if (!is.null(time_select)) {
    # Try different matching approaches with backward compatibility
    if (is.character(time_select)) {
      # Date-based filtering
      if ("time_period_date" %in% names(plot_points)) {
        filtered_points <- plot_points %>%
          dplyr::filter(!is.na(time_period_date) & as.character(time_period_date) == as.character(time_select))
        if (nrow(filtered_points) > 0) {
          plot_points <- filtered_points
          # Get corresponding numeric value if available
          if ("time_period_numeric" %in% names(plot_points)) {
            time_select <- unique(plot_points$time_period_numeric)[1]
          }
        }
      } else if ("time_period_label" %in% names(plot_points)) {
        filtered_points <- plot_points %>%
          dplyr::filter(!is.na(time_period_label) & as.character(time_period_label) == as.character(time_select))
        if (nrow(filtered_points) > 0) {
          plot_points <- filtered_points
          if ("time_period" %in% names(plot_points)) {
            time_select <- unique(plot_points$time_period)[1]
          }
        }
      }
    } else if (is.numeric(time_select)) {
      # Numeric filtering
      if ("time_period_numeric" %in% names(plot_points)) {
        filtered_points <- plot_points %>%
          dplyr::filter(time_period_numeric == time_select)
        if (nrow(filtered_points) > 0) {
          plot_points <- filtered_points
        }
      } else if ("time_period" %in% names(plot_points)) {
        filtered_points <- plot_points %>%
          dplyr::filter(time_period == time_select)
        if (nrow(filtered_points) > 0) {
          plot_points <- filtered_points
        }
      }
    }
  }
  
  # If still no match, just use the first time period for this fish
  check_condition <- TRUE
  if ("time_period_numeric" %in% names(plot_points)) {
    check_condition <- is.null(time_select) || nrow(plot_points %>% dplyr::filter(time_period_numeric == time_select)) == 0
  } else if ("time_period" %in% names(plot_points)) {
    check_condition <- is.null(time_select) || nrow(plot_points %>% dplyr::filter(time_period == time_select)) == 0
  }
  
  if (check_condition) {
    if ("time_period_numeric" %in% names(plot_points)) {
      available_periods <- unique(plot_points$time_period_numeric)
      time_select <- available_periods[1]
      plot_points <- plot_points %>%
        dplyr::filter(time_period_numeric == time_select)
    } else if ("time_period_date" %in% names(plot_points)) {
      available_periods <- unique(plot_points$time_period_date)
      time_select <- available_periods[1]
      plot_points <- plot_points %>%
        dplyr::filter(time_period_date == time_select)
    } else if ("time_period" %in% names(plot_points)) {
      available_periods <- unique(plot_points$time_period)
      time_select <- available_periods[1]
      plot_points <- plot_points %>%
        dplyr::filter(time_period == time_select)
    }
  }
  
  }  # End of time filtering block
  
  if (nrow(plot_points) == 0) {
    stop("No points found for fish ", fish_select, " and time period ", time_select)
  }
  
  method <- unique(plot_points$method)[1]
  
  # Simplify: Just use the original track data and filter by fish only
  # No complex time period matching needed
  track_fish_col <- if ("fish_id" %in% names(track_data)) "fish_id" else "path_id"
  
  track_filtered <- track_data %>%
    dplyr::filter(!!rlang::sym(track_fish_col) == fish_select)
  
  # Create base plot
  p <- ggplot2::ggplot()
  
  # Add background raster if provided
  if (!is.null(background_raster)) {
    tryCatch({
      # Calculate plot extent
      all_x <- c(track_filtered$x, plot_points$x)
      all_y <- c(track_filtered$y, plot_points$y)
      
      if (length(all_x) > 0) {
        x_range <- range(all_x, na.rm = TRUE)
        y_range <- range(all_y, na.rm = TRUE)
        x_buffer <- diff(x_range) * 0.2
        y_buffer <- diff(y_range) * 0.2
        
        plot_extent <- raster::extent(
          x_range[1] - x_buffer, x_range[2] + x_buffer,
          y_range[1] - y_buffer, y_range[2] + y_buffer
        )
        
        cropped_raster <- raster::crop(background_raster, plot_extent)
        raster_df <- raster::as.data.frame(cropped_raster, xy = TRUE)
        
        if (ncol(raster_df) >= 3) {
          names(raster_df)[3] <- "values"
          raster_df <- raster_df[!is.na(raster_df$values), ]
          
          if (nrow(raster_df) > 0) {
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
    }, error = function(e) {
      warning("Could not add background raster: ", e$message)
    })
  }
  
  # Add space use visualization
  if (show_space_use) {
    if (method == "grid_cell_count" || method == "constrained_convex_hull") {
      # Show actual cells if reference raster available
      if (!is.null(reference_raster)) {
        tryCatch({
          if (method == "grid_cell_count") {
            # Show cells used by track
            track_coords <- cbind(track_filtered$x, track_filtered$y)
            cell_indices <- raster::cellFromXY(reference_raster, track_coords)
            valid_cell_indices <- unique(cell_indices[!is.na(cell_indices)])
          } else {
            # Show cells within convex hull
            if (nrow(track_filtered) < 3) {
              valid_cell_indices <- c()  # Skip hull visualization
            } else {
              hull_indices <- grDevices::chull(track_filtered$x, track_filtered$y)
              hull_x <- track_filtered$x[hull_indices]
              hull_y <- track_filtered$y[hull_indices]
            
              res_x <- raster::res(reference_raster)[1]
              res_y <- raster::res(reference_raster)[2]  
              hull_extent <- raster::extent(min(hull_x) - res_x, max(hull_x) + res_x,
                                           min(hull_y) - res_y, max(hull_y) + res_y)
              
              all_cells <- raster::cellsFromExtent(reference_raster, hull_extent)
              cell_values <- raster::extract(reference_raster, all_cells)
              valid_cells <- all_cells[!is.na(cell_values)]
              
              cells_intersect <- test_cells_in_polygon(reference_raster, valid_cells, hull_x, hull_y)
              valid_cell_indices <- valid_cells[cells_intersect]
            }  # End of hull calculation block
          }
          
          if (length(valid_cell_indices) > 0) {
            # Create cell polygons
            cell_centers <- raster::xyFromCell(reference_raster, valid_cell_indices)
            res_x <- raster::res(reference_raster)[1]
            res_y <- raster::res(reference_raster)[2]
            
            cell_polygons <- data.frame()
            for (i in seq_along(valid_cell_indices)) {
              center_x <- cell_centers[i, 1]
              center_y <- cell_centers[i, 2]
              
              cell_poly <- data.frame(
                cell_id = rep(valid_cell_indices[i], 5),
                x = center_x + c(-res_x/2, res_x/2, res_x/2, -res_x/2, -res_x/2),
                y = center_y + c(-res_y/2, -res_y/2, res_y/2, res_y/2, -res_y/2)
              )
              cell_polygons <- rbind(cell_polygons, cell_poly)
            }
            
            p <- p +
              ggplot2::geom_polygon(
                data = cell_polygons,
                ggplot2::aes(x = x, y = y, group = cell_id),
                fill = "red",
                alpha = 0.3,
                color = "darkred",
                size = 0.2
              )
          }
        }, error = function(e) {
          warning("Could not add cell visualization: ", e$message)
        })
      }
    } else {
      # Show polygon outline for other methods
      if (method %in% c("convex_hull", "mcp")) {
        hull_indices <- grDevices::chull(track_filtered$x, track_filtered$y)
        polygon_data <- data.frame(x = track_filtered$x[hull_indices], 
                                  y = track_filtered$y[hull_indices])
      } else if (method == "bounding_box") {
        x_range <- range(track_filtered$x)
        y_range <- range(track_filtered$y)
        polygon_data <- data.frame(
          x = c(x_range[1], x_range[2], x_range[2], x_range[1], x_range[1]),
          y = c(y_range[1], y_range[1], y_range[2], y_range[2], y_range[1])
        )
      }
      
      if (exists("polygon_data") && nrow(polygon_data) > 0) {
        p <- p +
          ggplot2::geom_polygon(
            data = polygon_data,
            ggplot2::aes(x = x, y = y),
            fill = "red",
            alpha = 0.3,
            color = "darkred",
            size = 1
          )
      }
    }
  }
  
  # Add track
  if (show_track && nrow(track_filtered) > 0) {
    if (nrow(track_filtered) > 1) {
      p <- p +
        ggplot2::geom_path(
          data = track_filtered,
          ggplot2::aes(x = x, y = y),
          color = "green",
          size = 2,
          alpha = 0.8
        )
    }
    
    p <- p +
      ggplot2::geom_point(
        data = track_filtered,
        ggplot2::aes(x = x, y = y),
        color = "darkgreen",
        size = 2.5,
        alpha = 0.9
      )
  }
  
  # Add generated points
  if (show_points) {
    p <- p +
      ggplot2::geom_point(
        data = plot_points,
        ggplot2::aes(x = x, y = y),
        color = "blue",
        alpha = 0.7,
        size = 1.5
      )
  }
  
  # Calculate zoom and styling
  all_x <- c(track_filtered$x, plot_points$x)
  all_y <- c(track_filtered$y, plot_points$y)
  
  if (length(all_x) > 0) {
    x_range <- range(all_x, na.rm = TRUE)
    y_range <- range(all_y, na.rm = TRUE)
    x_buffer <- diff(x_range) * 0.1
    y_buffer <- diff(y_range) * 0.1
    
    p <- p +
      ggplot2::coord_fixed(
        xlim = c(x_range[1] - x_buffer, x_range[2] + x_buffer),
        ylim = c(y_range[1] - y_buffer, y_range[2] + y_buffer)
      )
  }
  
  # Labels and styling
  area_hectares <- round(unique(plot_points$space_use_area_hectares)[1], 2)
  
  # Use time_period_label for display, or format numeric time properly
  display_time <- if("time_period_label" %in% names(plot_points) && !is.na(unique(plot_points$time_period_label)[1])) {
    unique(plot_points$time_period_label)[1]
  } else if("time_period_date" %in% names(plot_points) && !is.na(unique(plot_points$time_period_date)[1])) {
    as.character(unique(plot_points$time_period_date)[1])
  } else if(is.numeric(time_select) && time_select > 1e9) {
    # Looks like a timestamp, convert to readable date
    format(as.POSIXct(time_select, origin = "1970-01-01", tz = "UTC"), "%Y-%m-%d")
  } else if(time_select == "no_time_aggregation") {
    "All Time"
  } else {
    as.character(time_select)
  }
  
  p <- p +
    ggplot2::labs(
      title = paste("Space Use Point Generation:", tools::toTitleCase(gsub("_", " ", method))),
      subtitle = paste0("Fish ", fish_select, " - Time Period ", display_time, 
                       " | Area: ", area_hectares, " ha | Points: ", nrow(plot_points)),
      x = "X Coordinate (m)",
      y = "Y Coordinate (m)",
      caption = "Red = space use area | Green = fish track | Blue = random points"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 12),
      plot.caption = ggplot2::element_text(size = 10, face = "italic")
    )
  
  return(p)
}

#' Summary function for space use points
#'
#' @param space_use_points Output from generate_random_points_in_space_use
#' 
#' @return A data frame with summary statistics for each fish and time period
#' 
#' @export
summarize_space_use_points <- function(space_use_points) {
  
  summary_stats <- space_use_points %>%
    dplyr::group_by(fish_id, time_period, method) %>%
    dplyr::summarise(
      n_points = dplyr::n(),
      space_use_area_hectares = dplyr::first(space_use_area_hectares),
      mean_x = mean(x),
      mean_y = mean(y),
      sd_x = sd(x),
      sd_y = sd(y),
      min_x = min(x),
      max_x = max(x),
      min_y = min(y),
      max_y = max(y),
      range_x = max_x - min_x,
      range_y = max_y - min_y,
      .groups = 'drop'
    )
  
  return(summary_stats)
}

#' Compare multiple space use methods by generating points for each
#'
#' @param space_use_results Output from calculate_space_use function
#' @param track_data Original track data
#' @param methods Vector of methods to compare
#' @param n_points Number of points per method
#' @param fish_ids Which fish to include
#' @param time_periods Which time periods to include
#' @param reference_raster Optional reference raster
#' 
#' @return A data frame containing generated points for each method with comparison metadata
#' 
#' @export
compare_space_use_methods <- function(space_use_results,
                                      track_data,
                                      methods = c("convex_hull", "constrained_convex_hull", "bounding_box"),
                                      n_points = 100,
                                      fish_ids = NULL,
                                      time_periods = NULL,
                                      reference_raster = NULL) {
  
  cat("=== COMPARING SPACE USE METHODS ===\n")
  
  all_points_list <- list()
  
  for (method in methods) {
    cat("\nGenerating points for method:", method, "\n")
    
    tryCatch({
      points <- generate_random_points_in_space_use(
        space_use_results = space_use_results,
        track_data = track_data,
        method = method,
        n_points = n_points,
        fish_ids = fish_ids,
        time_periods = time_periods,
        reference_raster = reference_raster
      )
      
      if (nrow(points) > 0) {
        all_points_list[[method]] <- points
      }
      
    }, error = function(e) {
      warning("Error generating points for method ", method, ": ", e$message)
    })
  }
  
  if (length(all_points_list) > 0) {
    # Combine all points
    combined_points <- dplyr::bind_rows(all_points_list)
    
    cat("\n=== COMPARISON COMPLETE ===\n")
    cat("Methods compared:", paste(names(all_points_list), collapse = ", "), "\n")
    cat("Total points generated:", nrow(combined_points), "\n")
    
    return(combined_points)
  } else {
    warning("No points generated for any method")
    return(data.frame())
  }
}

#' Plot comparison of multiple space use methods
#'
#' @param comparison_points Output from compare_space_use_methods
#' @param track_data Original track data
#' @param fish_select Which fish to plot
#' @param time_select Which time period to plot
#' @param background_raster Optional background raster
#' 
#' @return A ggplot2 object comparing space use methods with generated points
#' 
#' @export
plot_space_use_method_comparison <- function(comparison_points,
                                             track_data,
                                             fish_select = NULL,
                                             time_select = NULL,
                                             background_raster = NULL) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for plotting")
  }
  
  # Select fish and time if not specified
  if (is.null(fish_select)) {
    fish_select <- unique(comparison_points$fish_id)[1]
  }
  if (is.null(time_select)) {
    time_select <- unique(comparison_points$time_period)[1]
  }
  
  # Filter data - handle both numeric and string time_select
  if (is.character(time_select) && "time_period_label" %in% names(comparison_points)) {
    plot_points <- comparison_points %>%
      dplyr::filter(fish_id == fish_select & 
                   (!is.na(time_period_label) & as.character(time_period_label) == as.character(time_select)))
  } else {
    plot_points <- comparison_points %>%
      dplyr::filter(fish_id == fish_select & time_period == time_select)
  }
  
  if (nrow(plot_points) == 0) {
    stop("No points found for fish ", fish_select, " and time period ", time_select)
  }
  
  # Standardize and filter track data to match the same time binning as comparison_points
  time_bin_info <- list(
    method = "posix_based",
    aggregation_method = "day",
    available_periods = unique(comparison_points$time_period),
    available_labels = if("time_period_label" %in% names(comparison_points)) unique(comparison_points$time_period_label) else c()
  )
  
  track_std <- standardize_track_data_for_points(track_data, time_bin_info)
  
  if (is.character(time_select) && "time_period_label" %in% names(track_std)) {
    track_filtered <- track_std %>%
      dplyr::filter(fish_id == fish_select & 
                   (!is.na(time_period_label) & as.character(time_period_label) == as.character(time_select)))
  } else {
    track_filtered <- track_std %>%
      dplyr::filter(fish_id == fish_select & time_period == time_select)
  }
  
  # Create base plot
  p <- ggplot2::ggplot()
  
  # Add background raster if provided
  if (!is.null(background_raster)) {
    tryCatch({
      all_x <- c(track_filtered$x, plot_points$x)
      all_y <- c(track_filtered$y, plot_points$y)
      
      if (length(all_x) > 0) {
        x_range <- range(all_x, na.rm = TRUE)
        y_range <- range(all_y, na.rm = TRUE)
        x_buffer <- diff(x_range) * 0.2
        y_buffer <- diff(y_range) * 0.2
        
        plot_extent <- raster::extent(
          x_range[1] - x_buffer, x_range[2] + x_buffer,
          y_range[1] - y_buffer, y_range[2] + y_buffer
        )
        
        cropped_raster <- raster::crop(background_raster, plot_extent)
        raster_df <- raster::as.data.frame(cropped_raster, xy = TRUE)
        
        if (ncol(raster_df) >= 3) {
          names(raster_df)[3] <- "values"
          raster_df <- raster_df[!is.na(raster_df$values), ]
          
          if (nrow(raster_df) > 0) {
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
    }, error = function(e) {
      warning("Could not add background raster: ", e$message)
    })
  }
  
  # Add track
  if (nrow(track_filtered) > 0) {
    if (nrow(track_filtered) > 1) {
      p <- p +
        ggplot2::geom_path(
          data = track_filtered,
          ggplot2::aes(x = x, y = y),
          color = "black",
          size = 2,
          alpha = 0.8
        )
    }
    
    p <- p +
      ggplot2::geom_point(
        data = track_filtered,
        ggplot2::aes(x = x, y = y),
        color = "black",
        size = 3,
        alpha = 0.9
      )
  }
  
  # Add points colored by method
  p <- p +
    ggplot2::geom_point(
      data = plot_points,
      ggplot2::aes(x = x, y = y, color = method),
      alpha = 0.7,
      size = 1.5
    ) +
    ggplot2::scale_color_viridis_d(name = "Method")
  
  # Calculate zoom and styling
  all_x <- c(track_filtered$x, plot_points$x)
  all_y <- c(track_filtered$y, plot_points$y)
  
  if (length(all_x) > 0) {
    x_range <- range(all_x, na.rm = TRUE)
    y_range <- range(all_y, na.rm = TRUE)
    x_buffer <- diff(x_range) * 0.1
    y_buffer <- diff(y_range) * 0.1
    
    p <- p +
      ggplot2::coord_fixed(
        xlim = c(x_range[1] - x_buffer, x_range[2] + x_buffer),
        ylim = c(y_range[1] - y_buffer, y_range[2] + y_buffer)
      )
  }
  
  # Labels and styling
  methods_used <- paste(unique(plot_points$method), collapse = ", ")
  n_methods <- length(unique(plot_points$method))
  
  # Use time_period_label for display, or format numeric time properly
  display_time <- if("time_period_label" %in% names(plot_points) && !is.na(unique(plot_points$time_period_label)[1])) {
    unique(plot_points$time_period_label)[1]
  } else if("time_period_date" %in% names(plot_points) && !is.na(unique(plot_points$time_period_date)[1])) {
    as.character(unique(plot_points$time_period_date)[1])
  } else if(is.numeric(time_select) && time_select > 1e9) {
    # Looks like a timestamp, convert to readable date
    format(as.POSIXct(time_select, origin = "1970-01-01", tz = "UTC"), "%Y-%m-%d")
  } else if(time_select == "no_time_aggregation") {
    "All Time"
  } else {
    as.character(time_select)
  }
  
  p <- p +
    ggplot2::labs(
      title = "Space Use Method Comparison",
      subtitle = paste0("Fish ", fish_select, " - Time Period ", display_time, 
                       " | Methods: ", methods_used),
      x = "X Coordinate (m)",
      y = "Y Coordinate (m)",
      caption = "Black = fish track | Colors = different space use methods"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 12),
      plot.caption = ggplot2::element_text(size = 10, face = "italic"),
      legend.position = "bottom"
    )
  
  return(p)
}
