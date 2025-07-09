#' Plot Particle Filter Results
#'
#' @description
#' Creates comprehensive visualizations of particle filter positioning results,
#' showing estimated fish tracks, uncertainty bounds, particle distributions,
#' and receiver station locations with detection information.
#'
#' @param pf_results A list object returned by \code{\link{particle_filter_positioning}}
#'   containing position estimates, particle trajectories, and summary statistics.
#' @param fish_select Numeric or character. Fish ID to plot. If NULL (default),
#'   plots the first fish in the results. Use \code{unique(pf_results$position_estimates$fish_id)}
#'   to see available fish IDs.
#' @param show_particles Logical. Whether to display individual particles as points
#'   (default: FALSE). When TRUE, shows a sample of particles with transparency
#'   proportional to their weights.
#' @param show_uncertainty Logical. Whether to display uncertainty bounds as
#'   error bars around position estimates (default: TRUE). Shows 95% confidence
#'   intervals as red cross-shaped error bars.
#' @param depth_raster_df Optional. A data frame or raster object containing
#'   bathymetric/depth data for background visualization. Will be displayed
#'   with improved transparency and color handling.
#' @param station_coords Optional. Data frame or sf object containing receiver
#'   station coordinates for plotting station locations. Same format as used in
#'   \code{\link{particle_filter_positioning}}.
#' @param detection_data Optional. Data frame containing detection events for
#'   enhanced station visualization. Should contain columns: fish_id (or path_id),
#'   station_id, and optionally detected.
#' @param actual_track Optional. Data frame containing the true fish track for
#'   comparison with estimates. Should contain columns: fish_id (or path_id),
#'   x, y, and optionally time information. Plotted in green.
#'
#' @return A ggplot2 object showing:
#'   \itemize{
#'     \item Blue line and points: Estimated fish track
#'     \item Red error bars: 95% confidence intervals (if show_uncertainty = TRUE)
#'     \item Gray points: Individual particles (if show_particles = TRUE)
#'     \item Green line and points: Actual track (if provided)
#'     \item Yellow circles: Detecting stations (sized by detection count)
#'     \item Red circles: Non-detecting stations
#'     \item Light blue background: Bathymetry data (if provided)
#'   }
#'
#' @details
#' The function creates a comprehensive visualization with several components:
#' \itemize{
#'   \item Position estimates shown as blue points connected by lines
#'   \item Station symbols sized by detection count when detection_data provided
#'   \item Automatic zoom to detection area with 25% buffer
#'   \item Bathymetry background with proper transparency handling
#'   \item Uncertainty visualization as confidence interval error bars
#' }
#'
#' @section Automatic Scaling:
#' The plot automatically zooms to show the relevant detection area:
#' \itemize{
#'   \item Primary: Scale to detecting stations (where fish was detected)
#'   \item Fallback: Scale to all stations if no detections
#'   \item Final fallback: Scale to position estimates
#'   \item Adds 25% buffer around the selected area
#' }
#'
#' @examples
#' \dontrun{
#' # Basic plot
#' plot_particle_filter_results(
#'   pf_results = pf_results,
#'   fish_select = 1,
#'   station_coords = receiver_stations
#' )
#'
#' # Full-featured plot
#' plot_particle_filter_results(
#'   pf_results = pf_results,
#'   fish_select = 1,
#'   station_coords = receiver_stations,
#'   detection_data = detection_events,
#'   depth_raster_df = depth_raster,
#'   actual_track = true_tracks,
#'   show_uncertainty = TRUE
#' )
#'
#' # Debug plot with particles
#' plot_particle_filter_results(
#'   pf_results = pf_results,
#'   fish_select = 1,
#'   show_particles = TRUE,
#'   show_uncertainty = FALSE
#' )
#' }
#'
#' @seealso
#' \code{\link{particle_filter_positioning}} for generating the positioning results
#'
#' @export
plot_particle_filter_results <- function(pf_results,
                                         fish_select = NULL,
                                         show_particles = FALSE,
                                         show_uncertainty = TRUE,
                                         depth_raster_df = NULL,
                                         station_coords = NULL,
                                         detection_data = NULL,
                                         actual_track = NULL) {

  # Load required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for this function to work. Please install it.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' needed for this function to work. Please install it.")
  }

  # Validate input structure
  required_components <- c("position_estimates", "all_particles", "summary_stats")
  missing_components <- setdiff(required_components, names(pf_results))
  if (length(missing_components) > 0) {
    stop("pf_results missing required components: ", paste(missing_components, collapse = ", "),
         "\nExpected output from particle_filter_positioning()")
  }

  # Select fish to plot
  available_fish <- unique(pf_results$position_estimates$fish_id)
  if (length(available_fish) == 0) {
    stop("No fish found in position estimates")
  }

  if (is.null(fish_select)) {
    fish_select <- available_fish[1]
    message("No fish specified, plotting fish ID: ", fish_select)
  }

  if (!fish_select %in% available_fish) {
    stop("Fish ID ", fish_select, " not found in results. Available fish: ",
         paste(available_fish, collapse = ", "))
  }

  # Filter data for selected fish
  position_data <- pf_results$position_estimates %>%
    dplyr::filter(fish_id == fish_select)

  particle_data <- pf_results$all_particles %>%
    dplyr::filter(fish_id == fish_select)

  if (nrow(position_data) == 0) {
    stop("No position data found for fish_id: ", fish_select)
  }

  # Create base plot
  p <- ggplot2::ggplot()

  # Add improved depth background if provided
  if (!is.null(depth_raster_df)) {
    # Handle different input types
    if (is.data.frame(depth_raster_df)) {
      raster_df <- depth_raster_df
    } else if ("RasterLayer" %in% class(depth_raster_df)) {
      # Convert raster to dataframe if raster package available
      if (requireNamespace("raster", quietly = TRUE)) {
        raster_df <- raster::as.data.frame(depth_raster_df, xy = TRUE)
      } else {
        warning("Raster package needed to plot raster depth data. Skipping depth background.")
        raster_df <- NULL
      }
    } else {
      warning("depth_raster_df should be a data frame or RasterLayer. Skipping depth background.")
      raster_df <- NULL
    }

    # Add improved raster background if successfully processed
    if (!is.null(raster_df) && nrow(raster_df) > 0) {
      # Remove NA values to avoid artifacts
      raster_df <- raster_df[!is.na(raster_df[[3]]), ]

      # Determine the value column name (more robust detection)
      value_col <- if ("layer" %in% names(raster_df)) {
        "layer"
      } else if ("depth" %in% names(raster_df)) {
        "depth"
      } else if (ncol(raster_df) > 2) {
        names(raster_df)[3]  # Use third column
      } else {
        NULL
      }

      if (!is.null(value_col) && nrow(raster_df) > 0) {
        # Add bathymetry with improved transparency layering
        p <- p +
          # Base layer with depth-based transparency
          ggplot2::geom_raster(
            data = raster_df,
            ggplot2::aes(x = x, y = y, alpha = !!ggplot2::sym(value_col))
          ) +
          ggplot2::scale_alpha_continuous(range = c(0.1, 0.4), guide = "none") +
          # Overlay with consistent light blue color
          ggplot2::geom_raster(
            data = raster_df,
            ggplot2::aes(x = x, y = y),
            fill = "lightblue", alpha = 0.3
          )
      }
    }
  }

  # Add particles if requested
  if (show_particles && nrow(particle_data) > 0) {
    # Sample particles to avoid overplotting (limit to 1000 points)
    n_sample <- min(1000, nrow(particle_data))

    if (nrow(particle_data) > n_sample) {
      particle_sample <- particle_data %>%
        dplyr::slice_sample(n = n_sample)
      message("Showing sample of ", n_sample, " particles (", nrow(particle_data), " total)")
    } else {
      particle_sample <- particle_data
    }

    p <- p +
      ggplot2::geom_point(
        data = particle_sample,
        ggplot2::aes(x = x, y = y, alpha = weight),
        color = "gray", size = 0.5
      ) +
      ggplot2::scale_alpha_continuous(
        guide = "none",
        range = c(0.1, 0.7)
      )
  }

  # Add uncertainty bounds if requested
  if (show_uncertainty && nrow(position_data) > 0) {
    # Check if confidence interval columns exist
    ci_cols <- c("x_ci_lower", "x_ci_upper", "y_ci_lower", "y_ci_upper")
    if (all(ci_cols %in% names(position_data))) {
      p <- p +
        # Horizontal uncertainty bars (x-direction)
        ggplot2::geom_segment(
          data = position_data,
          ggplot2::aes(x = x_ci_lower, xend = x_ci_upper,
                       y = y_est, yend = y_est),
          color = "red", alpha = 0.6, size = 1
        ) +
        # Vertical uncertainty bars (y-direction)
        ggplot2::geom_segment(
          data = position_data,
          ggplot2::aes(x = x_est, xend = x_est,
                       y = y_ci_lower, yend = y_ci_upper),
          color = "red", alpha = 0.6, size = 1
        )
    } else {
      warning("Confidence interval columns not found. Cannot show uncertainty bounds.")
    }
  }

  # Add estimated track
  if (nrow(position_data) > 1) {
    p <- p + ggplot2::geom_path(
      data = position_data,
      ggplot2::aes(x = x_est, y = y_est),
      color = "blue", size = 1.2, alpha = 0.8
    )
  }

  # Add estimated positions
  if (!is.null(detection_data) && !is.null(station_coords)) {
    # If we have detection data, don't use size for uncertainty to avoid scale conflicts
    p <- p + ggplot2::geom_point(
      data = position_data,
      ggplot2::aes(x = x_est, y = y_est),
      color = "blue", alpha = 0.8, size = 3
    )
  } else {
    # Use size for uncertainty only when we don't have detection data
    p <- p + ggplot2::geom_point(
      data = position_data,
      ggplot2::aes(x = x_est, y = y_est, size = position_uncertainty),
      color = "blue", alpha = 0.8
    ) +
      ggplot2::scale_size_continuous(
        name = "Uncertainty\n(meters)",
        range = c(1, 5),
        guide = ggplot2::guide_legend(override.aes = list(alpha = 1))
      )
  }

  # Add actual track if provided
  if (!is.null(actual_track)) {
    # Handle different column naming for fish ID
    actual_track_std <- actual_track
    if ("path_id" %in% names(actual_track) && !"fish_id" %in% names(actual_track)) {
      actual_track_std$fish_id <- actual_track_std$path_id
    }

    # Check if fish exists in actual track
    if ("fish_id" %in% names(actual_track_std)) {
      actual_fish_track <- actual_track_std %>%
        dplyr::filter(fish_id == fish_select)

      if (nrow(actual_fish_track) > 0) {
        # Sort by time if time column exists
        if ("time_bin" %in% names(actual_fish_track)) {
          actual_fish_track <- actual_fish_track %>% dplyr::arrange(time_bin)
        } else if ("time_seconds" %in% names(actual_fish_track)) {
          actual_fish_track <- actual_fish_track %>% dplyr::arrange(time_seconds)
        } else if ("time" %in% names(actual_fish_track)) {
          actual_fish_track <- actual_fish_track %>% dplyr::arrange(time)
        }

        # Add track line
        if (nrow(actual_fish_track) > 1) {
          p <- p + ggplot2::geom_path(
            data = actual_fish_track,
            ggplot2::aes(x = x, y = y),
            color = "green", size = 1.5, alpha = 0.7
          )
        }

        # Add track points
        p <- p + ggplot2::geom_point(
          data = actual_fish_track,
          ggplot2::aes(x = x, y = y),
          color = "green", size = 2, alpha = 0.7
        )
      } else {
        warning("No actual track data found for fish_id: ", fish_select)
      }
    } else {
      warning("actual_track missing fish_id column. Cannot match with selected fish.")
    }
  } else {
    # Try to extract actual track from detection_data if available
    if (!is.null(detection_data)) {
      detection_data_std <- detection_data

      # Handle different column naming conventions
      if ("path_id" %in% names(detection_data) && !"fish_id" %in% names(detection_data)) {
        detection_data_std$fish_id <- detection_data_std$path_id
      }
      if ("time_seconds" %in% names(detection_data) && !"time_bin" %in% names(detection_data)) {
        detection_data_std$time_bin <- detection_data_std$time_seconds
      }

      # Check if detection data has position information (x, y columns)
      if (all(c("fish_id", "x", "y") %in% names(detection_data_std))) {

        # Determine which time column to use
        time_col <- if ("time_bin" %in% names(detection_data_std)) {
          "time_bin"
        } else if ("time_seconds" %in% names(detection_data_std)) {
          "time_seconds"
        } else {
          NULL
        }

        # Extract actual track from detections
        if (!is.null(time_col)) {
          actual_track_from_detections <- detection_data_std %>%
            dplyr::filter(fish_id == fish_select) %>%
            dplyr::select(fish_id, x, y, time = !!dplyr::sym(time_col)) %>%
            dplyr::distinct() %>%
            dplyr::arrange(time)
        } else {
          actual_track_from_detections <- detection_data_std %>%
            dplyr::filter(fish_id == fish_select) %>%
            dplyr::select(fish_id, x, y) %>%
            dplyr::distinct()
        }

        if (nrow(actual_track_from_detections) > 0) {
          # Add track line
          if (nrow(actual_track_from_detections) > 1) {
            p <- p + ggplot2::geom_path(
              data = actual_track_from_detections,
              ggplot2::aes(x = x, y = y),
              color = "green", size = 1.5, alpha = 0.7
            )
          }

          # Add track points
          p <- p + ggplot2::geom_point(
            data = actual_track_from_detections,
            ggplot2::aes(x = x, y = y),
            color = "green", size = 2, alpha = 0.7
          )
        }
      }
    }
  }

  # Add enhanced stations visualization
  if (!is.null(station_coords)) {
    # Handle sf objects
    if ("sf" %in% class(station_coords)) {
      if (requireNamespace("sf", quietly = TRUE)) {
        coords_matrix <- sf::st_coordinates(station_coords)
        station_df <- station_coords %>%
          sf::st_drop_geometry() %>%
          dplyr::mutate(
            station_x = coords_matrix[,1],
            station_y = coords_matrix[,2]
          )
      } else {
        warning("sf package needed to handle sf objects. Skipping station coordinates.")
        station_df <- NULL
      }
    } else {
      station_df <- station_coords
      # Handle different column naming conventions
      if ("x" %in% names(station_df) && !"station_x" %in% names(station_df)) {
        station_df$station_x <- station_df$x
        station_df$station_y <- station_df$y
      }
    }

    # Process station visualization with detection data if available
    if (!is.null(station_df) && "station_x" %in% names(station_df) && "station_y" %in% names(station_df)) {

      if (!is.null(detection_data)) {
        # Enhanced station visualization with detection information
        detection_data_std <- detection_data

        # Handle different column naming conventions for detection data
        if ("path_id" %in% names(detection_data) && !"fish_id" %in% names(detection_data)) {
          detection_data_std$fish_id <- detection_data_std$path_id
        }
        if ("time_seconds" %in% names(detection_data) && !"time_bin" %in% names(detection_data)) {
          detection_data_std$time_bin <- detection_data_std$time_seconds
        }

        # Handle different column naming conventions for station coordinates
        station_df_std <- station_df
        if ("point_id" %in% names(station_df) && !"station_id" %in% names(station_df)) {
          station_df_std$station_id <- station_df_std$point_id
        }

        # Check required columns for detection data (simplified)
        req_det_cols <- c("fish_id", "station_id")
        has_station_id <- "station_id" %in% names(station_df_std)

        if (all(req_det_cols %in% names(detection_data_std)) && has_station_id) {

          # Get detection summary for selected fish (across all times)
          detection_summary <- detection_data_std %>%
            dplyr::filter(fish_id == fish_select) %>%
            dplyr::group_by(station_id) %>%
            dplyr::summarise(
              n_detections = if("detected" %in% names(detection_data_std)) {
                sum(detected, na.rm = TRUE)
              } else {
                dplyr::n()  # Count number of records as proxy for detections
              },
              .groups = 'drop'
            )

          # Merge with station coordinates
          station_with_detections <- station_df_std %>%
            dplyr::left_join(detection_summary, by = "station_id") %>%
            dplyr::mutate(
              n_detections = ifelse(is.na(n_detections), 0, n_detections),
              detecting = n_detections > 0
            )

          # Plot detecting stations (yellow with size based on detections)
          detecting_stations <- station_with_detections %>%
            dplyr::filter(detecting)

          if (nrow(detecting_stations) > 0) {
            p <- p + ggplot2::geom_point(
              data = detecting_stations,
              ggplot2::aes(x = station_x, y = station_y, size = n_detections),
              color = "yellow", fill = NA, shape = 21, stroke = 1.5
            ) +
              ggplot2::scale_size_continuous(
                name = "Detections",
                range = c(2, 8),
                guide = ggplot2::guide_legend(override.aes = list(stroke = 1.5))
              )
          }

          # Plot non-detecting stations (red, smaller)
          non_detecting_stations <- station_with_detections %>%
            dplyr::filter(!detecting)

          if (nrow(non_detecting_stations) > 0) {
            p <- p + ggplot2::geom_point(
              data = non_detecting_stations,
              ggplot2::aes(x = station_x, y = station_y),
              color = "red", size = 2, alpha = 0.7
            )
          }

        } else {
          # Fallback to basic station plotting if detection data incomplete or station_id missing
          missing_cols <- c()
          if (!all(req_det_cols %in% names(detection_data_std))) {
            missing_det_cols <- setdiff(req_det_cols, names(detection_data_std))
            missing_cols <- c(missing_cols, paste("detection_data:", paste(missing_det_cols, collapse = ", ")))
          }
          if (!has_station_id) {
            available_station_cols <- names(station_df_std)[grepl("id|station|point", names(station_df_std), ignore.case = TRUE)]
            missing_cols <- c(missing_cols, paste("station_coords: station_id (available ID columns:", paste(available_station_cols, collapse = ", "), ")"))
          }

          warning("Missing required columns for detection visualization: ", paste(missing_cols, collapse = "; "),
                  ". Using basic station visualization.")
          p <- p + ggplot2::geom_point(
            data = station_df,
            ggplot2::aes(x = station_x, y = station_y),
            shape = 17, size = 3, color = "black"
          )
        }

      } else {
        # Basic station plotting without detection information
        p <- p + ggplot2::geom_point(
          data = station_df,
          ggplot2::aes(x = station_x, y = station_y),
          shape = 17, size = 3, color = "black"
        )
      }

    } else {
      warning("station_coords missing required columns (station_x, station_y or x, y). Skipping stations.")
    }
  }

  # Add formatting and labels
  subtitle_parts <- c("Blue = estimated track")
  if (!is.null(actual_track)) subtitle_parts <- c(subtitle_parts, "Green = actual track")
  if (!is.null(detection_data) && !is.null(station_coords)) {
    subtitle_parts <- c(subtitle_parts, "Yellow = detecting stations (sized by total count), Red = non-detecting")
  } else if (!is.null(station_coords)) {
    subtitle_parts <- c(subtitle_parts, "Black triangles = stations")
  }

  p <- p +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste("Particle Filter Positioning Results - Fish", fish_select),
      subtitle = paste(subtitle_parts, collapse = ", "),
      x = "X Coordinate (m)",
      y = "Y Coordinate (m)"
    )

  # Calculate automatic zoom based on detecting stations only
  if (!is.null(station_coords) && !is.null(detection_data)) {
    # Try to use only detecting stations for scaling if available
    if (!is.null(station_df) && "station_x" %in% names(station_df) && "station_y" %in% names(station_df)) {

      # Check if we have detection information processed
      if (exists("detecting_stations") && nrow(detecting_stations) > 0) {
        # Use only detecting stations for scaling
        x_range <- range(detecting_stations$station_x, na.rm = TRUE)
        y_range <- range(detecting_stations$station_y, na.rm = TRUE)

        # Add buffer (25% of range on each side - increased from 15%)
        x_buffer <- diff(x_range) * 0.25
        y_buffer <- diff(y_range) * 0.25

        # Ensure minimum buffer size (increased from 50m)
        x_buffer <- max(x_buffer, 100)  # minimum 100m buffer
        y_buffer <- max(y_buffer, 100)

        # Set coordinate limits based on detecting stations only
        xlim <- c(x_range[1] - x_buffer, x_range[2] + x_buffer)
        ylim <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)

        p <- p + ggplot2::coord_equal(xlim = xlim, ylim = ylim)
      } else {
        # Fallback to all stations if no detecting stations found
        x_range <- range(station_df$station_x, na.rm = TRUE)
        y_range <- range(station_df$station_y, na.rm = TRUE)

        x_buffer <- diff(x_range) * 0.25
        y_buffer <- diff(y_range) * 0.25
        x_buffer <- max(x_buffer, 100)
        y_buffer <- max(y_buffer, 100)

        xlim <- c(x_range[1] - x_buffer, x_range[2] + x_buffer)
        ylim <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)

        p <- p + ggplot2::coord_equal(xlim = xlim, ylim = ylim)
      }
    } else {
      p <- p + ggplot2::coord_equal()
    }
  } else if (!is.null(station_coords)) {
    # Use all station coordinates for scaling if no detection data
    if (!is.null(station_df) && "station_x" %in% names(station_df) && "station_y" %in% names(station_df)) {
      x_range <- range(station_df$station_x, na.rm = TRUE)
      y_range <- range(station_df$station_y, na.rm = TRUE)

      x_buffer <- diff(x_range) * 0.25
      y_buffer <- diff(y_range) * 0.25
      x_buffer <- max(x_buffer, 100)
      y_buffer <- max(y_buffer, 100)

      xlim <- c(x_range[1] - x_buffer, x_range[2] + x_buffer)
      ylim <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)

      p <- p + ggplot2::coord_equal(xlim = xlim, ylim = ylim)
    } else {
      p <- p + ggplot2::coord_equal()
    }
  } else if (nrow(position_data) > 0) {
    # Fallback to position estimates if no stations provided
    x_cols <- c("x_est")
    y_cols <- c("y_est")

    if (show_uncertainty && all(c("x_ci_lower", "x_ci_upper") %in% names(position_data))) {
      x_cols <- c(x_cols, "x_ci_lower", "x_ci_upper")
      y_cols <- c(y_cols, "y_ci_lower", "y_ci_upper")
    }

    x_range <- range(unlist(position_data[x_cols]), na.rm = TRUE)
    y_range <- range(unlist(position_data[y_cols]), na.rm = TRUE)

    x_buffer <- diff(x_range) * 0.25
    y_buffer <- diff(y_range) * 0.25
    x_buffer <- max(x_buffer, 100)
    y_buffer <- max(y_buffer, 100)

    xlim <- c(x_range[1] - x_buffer, x_range[2] + x_buffer)
    ylim <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)

    p <- p + ggplot2::coord_equal(xlim = xlim, ylim = ylim)
  } else {
    p <- p + ggplot2::coord_equal()
  }

  return(p)
}
