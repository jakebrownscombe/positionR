#' Simulate fish tracks with correlated random walks and detection events
#'
#' Generates realistic fish movement paths using correlated random walks within
#' a raster boundary and simulates acoustic detections based on a detection
#' system's coverage probabilities.
#'
#' @param raster A RasterLayer object defining the study area boundaries.
#'   Fish paths are constrained to non-NA cells.
#' @param detection_system Data frame output from \code{\link{calculate_detection_system}}
#'   containing spatial detection probabilities.
#' @param station_distances Data frame output from \code{\link{calculate_station_distances}}
#'   with individual receiver detection probabilities. Default is NULL.
#' @param n_paths Integer. Number of fish paths to simulate. Default is 1.
#' @param n_steps Integer. Number of steps per path. Default is 100.
#' @param step_length_mean Numeric. Mean step length in map units. Default is 50.
#' @param step_length_sd Numeric. Standard deviation of step length. Default is 20.
#' @param turning_angle_mean Numeric. Mean turning angle in degrees. Default is 0.
#' @param turning_angle_sd Numeric. Standard deviation of turning angle in degrees.
#'   Default is 45.
#' @param time_step Numeric. Time between steps in seconds. Default is 60.
#' @param detection_type Character. Type of detection probability to use:
#'   "cumulative" (any receiver) or "3_plus" (3+ receivers for positioning).
#'   Default is "cumulative".
#' @param start_locations Matrix or data frame with x,y coordinates for starting
#'   locations. If NULL, random start locations are chosen. Default is NULL.
#' @param seed Numeric. Random seed for reproducible results. Default is NULL.
#'
#' @return A list containing:
#'   \item{tracks}{Data frame with path coordinates, step lengths, and bearings}
#'   \item{detections}{Data frame with system-level detection events}
#'   \item{station_detections}{Data frame with individual receiver detections}
#'   \item{parameters}{List of simulation parameters used}
#'
#' @details
#' The simulation process involves:
#' \enumerate{
#'   \item Generating correlated random walks using von Mises distributed turning angles
#'   \item Constraining paths to valid raster cells with boundary reflection
#'   \item Evaluating detection probabilities at each location
#'   \item Simulating detection events for individual receivers and system-level
#' }
#'
#' The correlated random walk uses:
#' \itemize{
#'   \item Normally distributed step lengths (truncated at minimum 5 units)
#'   \item Von Mises distributed turning angles for realistic directional persistence
#'   \item Boundary reflection when paths encounter raster edges or invalid cells
#' }
#'
#' Detection simulation requires both detection_system (for system probabilities)
#' and station_distances (for individual receiver probabilities) to generate
#' realistic detection patterns.
#'
#' @examples
#' \dontrun{
#' # Generate receiver network and detection system
#' stations <- generate_random_points(depth_raster, n_points = 6, seed = 123)
#' distances <- calculate_station_distances(depth_raster, stations, max_distance = 500)
#' de_model <- create_logistic_curve_depth(
#'   min_depth = 2, max_depth = 30,
#'   d50_min_depth = 50, d95_min_depth = 20,
#'   d50_max_depth = 150, d95_max_depth = 60,
#'   plot = FALSE
#' )
#' detection_sys <- calculate_detection_system(distances, stations, de_model$log_model)
#'
#' # Simulate fish tracks
#' fish_sim <- simulate_fish_tracks(
#'   raster = depth_raster,
#'   detection_system = detection_sys$data,
#'   station_distances = distances,
#'   n_paths = 3,
#'   n_steps = 50,
#'   seed = 456
#' )
#'
#' # Analyze results
#' detection_rate <- nrow(fish_sim$detections[fish_sim$detections$detected == 1, ]) /
#'                   nrow(fish_sim$detections)
#'
#' # Plot tracks and detections
#' track_plot <- plot_fish_tracks(fish_sim, depth_raster, stations)
#' }
#'
#' @seealso \code{\link{calculate_detection_system}}, \code{\link{calculate_station_distances}}, \code{\link{plot_fish_tracks}}
#'
#' @export
simulate_fish_tracks <- function(raster,
                                 detection_system,
                                 station_distances = NULL,  # Add station_distances_df parameter
                                 n_paths = 1,
                                 n_steps = 100,
                                 step_length_mean = 50,
                                 step_length_sd = 20,
                                 turning_angle_mean = 0,
                                 turning_angle_sd = 45,
                                 time_step = 60,  # seconds between steps
                                 detection_type = "cumulative",  # "cumulative" or "3_plus"
                                 start_locations = NULL,
                                 seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Validate detection_type
  if (!detection_type %in% c("cumulative", "3_plus")) {
    stop("detection_type must be 'cumulative' or '3_plus'")
  }

  # Convert degrees to radians for calculations
  turning_angle_sd_rad <- turning_angle_sd * pi / 180
  turning_angle_mean_rad <- turning_angle_mean * pi / 180

  # Get raster properties
  raster_extent <- raster::extent(raster)
  valid_cells <- which(!is.na(raster::values(raster)))
  valid_coords <- raster::xyFromCell(raster, valid_cells)

  # Initialize results
  all_tracks <- list()
  all_detections <- list()
  all_station_detections <- list()

  for (path_id in 1:n_paths) {
    cat("Generating path", path_id, "of", n_paths, "\n")

    # Set starting location
    if (is.null(start_locations)) {
      # Random start location within valid raster area
      start_idx <- sample(nrow(valid_coords), 1)
      start_x <- valid_coords[start_idx, 1]
      start_y <- valid_coords[start_idx, 2]
    } else {
      # Use provided start locations
      start_x <- start_locations[path_id, 1]
      start_y <- start_locations[path_id, 2]
    }

    # Initialize path tracking
    track_x <- rep(NA, n_steps + 1)
    track_y <- rep(NA, n_steps + 1)
    track_time <- seq(0, n_steps * time_step, by = time_step)

    # Set starting position and direction
    track_x[1] <- start_x
    track_y[1] <- start_y
    current_bearing <- runif(1, 0, 2 * pi)  # Random initial direction

    # Generate correlated random walk
    for (step in 1:n_steps) {
      # Generate step length (truncated normal to avoid negative values)
      step_length <- stats::rnorm(1, step_length_mean, step_length_sd)
      step_length <- max(step_length, 5)  # Minimum step length

      # Generate turning angle (von Mises distribution for circular data)
      # Suppress circular warnings
      suppressWarnings({
        turning_angle <- circular::rvonmises(1, mu = turning_angle_mean_rad,
                                             kappa = 1 / (turning_angle_sd_rad^2))
      })

      # Update bearing
      current_bearing <- current_bearing + turning_angle

      # Calculate new position
      new_x <- track_x[step] + step_length * cos(current_bearing)
      new_y <- track_y[step] + step_length * sin(current_bearing)

      # Check if new position is within raster bounds and valid
      if (new_x >= raster_extent@xmin && new_x <= raster_extent@xmax &&
          new_y >= raster_extent@ymin && new_y <= raster_extent@ymax) {

        # Check if position has valid raster value
        cell_value <- raster::extract(raster, matrix(c(new_x, new_y), ncol = 2))

        if (!is.na(cell_value)) {
          track_x[step + 1] <- new_x
          track_y[step + 1] <- new_y
        } else {
          # Bounce off invalid areas
          current_bearing <- current_bearing + pi  # Reverse direction
          new_x <- track_x[step] + step_length * cos(current_bearing)
          new_y <- track_y[step] + step_length * sin(current_bearing)
          track_x[step + 1] <- new_x
          track_y[step + 1] <- new_y
        }
      } else {
        # Bounce off raster boundaries
        current_bearing <- current_bearing + pi  # Reverse direction
        new_x <- track_x[step] + step_length * cos(current_bearing)
        new_y <- track_y[step] + step_length * sin(current_bearing)
        track_x[step + 1] <- new_x
        track_y[step + 1] <- new_y
      }
    }

    # Create track dataframe
    track_df <- data.frame(
      path_id = path_id,
      step = 0:n_steps,
      time_seconds = track_time,
      x = track_x,
      y = track_y
    ) %>%
      dplyr::filter(!is.na(x) & !is.na(y)) %>%
      dplyr::mutate(
        step_length = c(0, sqrt(diff(x)^2 + diff(y)^2)),
        bearing = c(0, atan2(diff(y), diff(x)) * 180 / pi)
      )

    all_tracks[[path_id]] <- track_df

    # Generate detections at each track point
    detections_list <- list()
    station_detections_list <- list()

    # Check if detection_system is valid
    if (!is.null(detection_system) && is.data.frame(detection_system) && nrow(detection_system) > 0) {

      for (i in 1:nrow(track_df)) {
        point_x <- track_df$x[i]
        point_y <- track_df$y[i]
        point_time <- track_df$time_seconds[i]

        # Find closest raster cell for this point
        distances_to_cells <- sqrt((detection_system$x - point_x)^2 +
                                     (detection_system$y - point_y)^2)
        closest_cell_idx <- which.min(distances_to_cells)
        closest_cell <- detection_system[closest_cell_idx, ]

        # Generate station-level detections for ALL stations at this location
        system_detected <- 0  # Track if any station detects

        if (!is.null(station_distances)) {
          # Find all station distances for this location
          fish_distances <- station_distances %>%
            dplyr::filter(cell_id == closest_cell$cell_id)

          if (nrow(fish_distances) > 0) {
            # For each station, use the pre-calculated detection probability
            for (j in 1:nrow(fish_distances)) {
              station_row <- fish_distances[j, ]
              station_id <- station_row$station_no
              station_distance <- station_row$cost_distance
              station_detection_prob <- station_row$DE_pred  # Use pre-calculated probability

              # Skip if distance is NA or detection probability is NA
              if (!is.na(station_distance) && !is.na(station_detection_prob)) {

                # Simulate detection for this specific station
                station_detected <- stats::rbinom(1, 1, station_detection_prob)

                if (station_detected == 1) {
                  system_detected <- 1  # Mark that system detected

                  station_detections_list[[length(station_detections_list) + 1]] <- data.frame(
                    path_id = path_id,
                    step = track_df$step[i],
                    time_seconds = point_time,
                    station_id = station_id,
                    x = point_x,
                    y = point_y,
                    distance_to_station = station_distance,
                    detection_prob = station_detection_prob,
                    detected = 1
                  )
                }
              }
            }
          }
        }

        # Record system-level detection based on whether ANY station detected
        if (detection_type == "cumulative") {
          prob_type <- "system"
          # Use the original system probability for recording, but detection status from stations
          if ("cumulative_prob" %in% names(detection_system)) {
            system_prob <- closest_cell$cumulative_prob
          } else {
            system_prob <- NA
          }
        } else if (detection_type == "3_plus") {
          prob_type <- "positioning"
          if ("prob_3_plus" %in% names(detection_system)) {
            system_prob <- closest_cell$prob_3_plus
          } else {
            system_prob <- NA
          }
        }

        # Record detection event (using actual station detection result)
        if (!is.na(system_prob)) {
          detections_list[[length(detections_list) + 1]] <- data.frame(
            path_id = path_id,
            step = track_df$step[i],
            time_seconds = point_time,
            x = point_x,
            y = point_y,
            cell_id = closest_cell$cell_id,
            detection_type = prob_type,
            detection_prob = system_prob,
            detected = system_detected  # Use result from station simulations
          )
        }
      }
    } else {
      cat("No valid detection system provided - generating tracks only\n")
    }

    # Combine detections for this path
    if (length(detections_list) > 0) {
      path_detections <- do.call(rbind, detections_list)
      all_detections[[path_id]] <- path_detections
    }

    if (length(station_detections_list) > 0) {
      path_station_detections <- do.call(rbind, station_detections_list)
      all_station_detections[[path_id]] <- path_station_detections
    }
  }

  # Combine all tracks and detections
  final_tracks <- do.call(rbind, all_tracks)

  if (length(all_detections) > 0) {
    final_detections <- do.call(rbind, all_detections)
  } else {
    final_detections <- data.frame()
  }

  if (length(all_station_detections) > 0) {
    final_station_detections <- do.call(rbind, all_station_detections)
  } else {
    final_station_detections <- data.frame()
  }

  # Return results
  results <- list(
    tracks = final_tracks,
    detections = final_detections,
    station_detections = final_station_detections,
    parameters = list(
      n_paths = n_paths,
      n_steps = n_steps,
      step_length_mean = step_length_mean,
      step_length_sd = step_length_sd,
      turning_angle_mean = turning_angle_mean,
      turning_angle_sd = turning_angle_sd,
      time_step = time_step,
      detection_type = detection_type
    )
  )

  return(results)
}

#' Plot simulated fish tracks with detection events
#'
#' Creates a visualization of simulated fish movement paths overlaid on a raster
#' background, showing detection events and receiver station performance.
#'
#' @param simulation_results List output from \code{\link{simulate_fish_tracks}}
#'   containing tracks, detections, and station_detections.
#' @param raster A RasterLayer object used as the background for plotting.
#' @param receiver_frame An sf object containing receiver station locations.
#'   Default is NULL.
#' @param show_detections Logical. Whether to display detection events on the plot.
#'   Default is TRUE.
#' @param path_alpha Numeric. Transparency level for path lines (0-1). Default is 0.7.
#'
#' @return A ggplot2 object showing fish tracks, detection events, and receiver stations.
#'
#' @details
#' The plot includes:
#' \itemize{
#'   \item Raster background (typically depth or habitat)
#'   \item Fish movement paths colored by path ID
#'   \item Detection events: yellow circles (successful), red X marks (missed)
#'   \item Receiver stations: sized by detection count, colored by activity
#' }
#'
#' Receiver station visualization:
#' \itemize{
#'   \item Green circles indicate stations with detections
#'   \item Red circles indicate stations with no detections
#'   \item Circle size scales with number of detections
#' }
#'
#' The function automatically handles different simulation outputs and adapts
#' the visualization based on available data (tracks only vs. full detection simulation).
#'
#' @examples
#' \dontrun{
#' # Generate and plot fish tracks
#' fish_sim <- simulate_fish_tracks(
#'   raster = depth_raster,
#'   detection_system = detection_sys$data,
#'   station_distances = distances,
#'   n_paths = 2,
#'   n_steps = 100
#' )
#'
#' # Basic plot
#' plot_fish_tracks(fish_sim, depth_raster, stations)
#'
#' # Plot without detection events
#' plot_fish_tracks(fish_sim, depth_raster, stations, show_detections = FALSE)
#'
#' # Customize path transparency
#' plot_fish_tracks(fish_sim, depth_raster, stations, path_alpha = 0.5)
#' }
#'
#' @seealso \code{\link{simulate_fish_tracks}}
#'
#' @export
plot_fish_tracks <- function(simulation_results, raster, receiver_frame = NULL,
                             show_detections = TRUE, path_alpha = 0.7) {
  # Convert raster to dataframe for plotting
  raster_df <- raster::as.data.frame(raster, xy = TRUE)
  # Base plot with raster
  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = raster_df, ggplot2::aes(x = x, y = y, fill = layer)) +
    ggplot2::scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                                 na.value = "transparent", name = "Depth") +
    ggplot2::theme_minimal()
  # Add tracks
  tracks <- simulation_results$tracks
  p <- p + ggplot2::geom_path(data = tracks,
                              ggplot2::aes(x = x, y = y, group = path_id, color = factor(path_id)),
                              alpha = path_alpha, size = 0.8) +
    ggplot2::scale_color_discrete(name = "Path ID")
  # Add detections if requested and available
  if (show_detections && nrow(simulation_results$detections) > 0) {
    detections <- simulation_results$detections
    # Separate detected and missed detections
    detected_events <- detections %>% dplyr::filter(detected == 1)
    missed_events <- detections %>% dplyr::filter(detected == 0)

    # Add missed events FIRST (failed detections) - X marks plotted first so they appear below
    if (nrow(missed_events) > 0) {
      p <- p + ggplot2::geom_point(data = missed_events, ggplot2::aes(x = x, y = y),
                                   color = "red", size = 2, alpha = 0.8,
                                   stroke = 1, shape = 4)  # X shape
    }

    # Add detected events LAST (successful detections) - plotted on top for visibility
    if (nrow(detected_events) > 0) {
      p <- p + ggplot2::geom_point(data = detected_events, ggplot2::aes(x = x, y = y),
                                   color = "black", size = 1.5, alpha = 0.8,
                                   stroke = 1, shape = 21, fill = "yellow")
    }
  }
  # Add receiver locations last so they're visible on top
  if (!is.null(receiver_frame)) {
    # Calculate detection counts per receiver station if station detection data exists
    if (!is.null(simulation_results$station_detections) && nrow(simulation_results$station_detections) > 0) {
      # Count detections per station
      station_counts <- simulation_results$station_detections %>%
        dplyr::group_by(station_id) %>%
        dplyr::summarise(detection_count = dplyr::n(), .groups = 'drop')
      # Get receiver coordinates and merge with detection counts
      receiver_coords <- sf::st_coordinates(receiver_frame)
      receiver_df <- data.frame(
        station_id = receiver_frame$point_id,
        x = receiver_coords[,1],
        y = receiver_coords[,2]
      ) %>%
        dplyr::left_join(station_counts, by = "station_id") %>%
        dplyr::mutate(detection_count = ifelse(is.na(detection_count), 0, detection_count))
      # Plot receivers with size based on detection count
      p <- p + ggplot2::geom_point(data = receiver_df,
                                   ggplot2::aes(x = x, y = y, size = detection_count,
                                                color = ifelse(detection_count == 0, "No Detections", "Detections")),
                                   alpha = 1, stroke = 1, shape = 21, fill = NA) +
        ggplot2::scale_size_continuous(name = "Station\nDetections", range = c(2, 8)) +
        ggplot2::scale_color_manual(values = c("No Detections" = "red", "Detections" = "green"), name="")
    } else {
      # Default receiver plotting if no station detection data
      p <- p + ggplot2::geom_sf(data = receiver_frame, color = "green", size = 2, alpha = 1,
                                stroke = 2, shape = 21, fill = "transparent")
    }
  }
  p <- p + ggplot2::coord_sf() + ggplot2::labs(title = "Fish Tracks and Detections",
                                               subtitle = "Yellow = Detected, Red X = Missed, Green circles = Receivers (size = detections)")
  return(p)
}








#' Calculate detection rate summaries from fish track simulations
#'
#' Computes detection performance metrics including overall detection rates,
#' per-path statistics, and summary statistics across all simulated fish tracks.
#'
#' @param simulation_results List output from \code{\link{simulate_fish_tracks}}
#'   containing tracks, detections, and station_detections data frames.
#'
#' @return A list containing:
#'   \item{overall}{List with total_steps, detected_steps, detection_rate, and detection_percentage}
#'   \item{by_path}{Data frame with detection statistics for each individual path}
#'   \item{summary_stats}{List with mean, median, min, max, and standard deviation of detection rates}
#'
#'   Returns NULL if no detection data is available.
#'
#' @details
#' The function analyzes detection success at each step of simulated fish tracks,
#' providing both aggregate and individual path performance metrics. This is useful
#' for evaluating receiver array effectiveness and optimizing system design.
#'
#' @examples
#' \dontrun{
#' # Simulate fish tracks with detections
#' fish_sim <- simulate_fish_tracks(
#'   raster = depth_raster,
#'   detection_system = system_DE$data,
#'   n_paths = 5,
#'   n_steps = 100
#' )
#'
#' # Calculate detection summaries
#' detection_stats <- calculate_detection_summaries(fish_sim)
#'
#' # Access overall performance
#' overall_rate <- detection_stats$overall$detection_percentage
#'
#' # Access per-path performance
#' path_performance <- detection_stats$by_path
#' }
#'
#' @seealso \code{\link{simulate_fish_tracks}}, \code{\link{print_detection_summary}}, \code{\link{plot_detection_rates}}
#'
#' @export
calculate_detection_summaries <- function(simulation_results) {

  if (is.null(simulation_results$detections) || nrow(simulation_results$detections) == 0) {
    cat("No detection data available in simulation results.\n")
    return(NULL)
  }

  detections <- simulation_results$detections

  # Overall detection rate
  overall_detected <- sum(detections$detected == 1, na.rm = TRUE)
  overall_total <- nrow(detections)
  overall_rate <- overall_detected / overall_total

  # Detection rate by path
  path_summary <- detections %>%
    dplyr::group_by(path_id) %>%
    dplyr::summarise(
      total_steps = dplyr::n(),
      detected_steps = sum(detected == 1, na.rm = TRUE),
      detection_rate = detected_steps / total_steps,
      .groups = 'drop'
    )

  # Detection rate statistics
  detection_stats <- list(
    overall = list(
      total_steps = overall_total,
      detected_steps = overall_detected,
      detection_rate = overall_rate,
      detection_percentage = round(overall_rate * 100, 1)
    ),
    by_path = path_summary,
    summary_stats = list(
      mean_detection_rate = mean(path_summary$detection_rate),
      median_detection_rate = median(path_summary$detection_rate),
      min_detection_rate = min(path_summary$detection_rate),
      max_detection_rate = max(path_summary$detection_rate),
      sd_detection_rate = sd(path_summary$detection_rate)
    )
  )

  return(detection_stats)
}

#' Print detection summary report in readable format
#'
#' Displays a formatted report of detection performance statistics including
#' overall rates, summary statistics, and individual path performance.
#'
#' @param detection_stats List output from \code{\link{calculate_detection_summaries}}
#'   containing detection performance metrics.
#'
#' @return Invisibly returns NULL. Called for side effect of printing summary report.
#'
#' @details
#' The printed report includes:
#' \itemize{
#'   \item Overall detection performance across all paths
#'   \item Summary statistics (mean, median, range, standard deviation)
#'   \item Individual path performance breakdown
#' }
#'
#' @examples
#' \dontrun{
#' # Calculate and print detection summaries
#' detection_stats <- calculate_detection_summaries(fish_simulation)
#' print_detection_summary(detection_stats)
#' }
#'
#' @seealso \code{\link{calculate_detection_summaries}}
#'
#' @export
print_detection_summary <- function(detection_stats) {

  if (is.null(detection_stats)) {
    return(invisible())
  }

  cat("=== DETECTION SUMMARY REPORT ===\n\n")

  # Overall summary
  cat("OVERALL DETECTION PERFORMANCE:\n")
  cat(sprintf("  Total steps simulated: %d\n", detection_stats$overall$total_steps))
  cat(sprintf("  Steps with detections: %d\n", detection_stats$overall$detected_steps))
  cat(sprintf("  Overall detection rate: %.1f%%\n\n", detection_stats$overall$detection_percentage))

  # Summary statistics across paths
  cat("DETECTION RATE STATISTICS ACROSS PATHS:\n")
  cat(sprintf("  Mean detection rate: %.1f%%\n", detection_stats$summary_stats$mean_detection_rate * 100))
  cat(sprintf("  Median detection rate: %.1f%%\n", detection_stats$summary_stats$median_detection_rate * 100))
  cat(sprintf("  Range: %.1f%% - %.1f%%\n",
              detection_stats$summary_stats$min_detection_rate * 100,
              detection_stats$summary_stats$max_detection_rate * 100))
  cat(sprintf("  Standard deviation: %.1f%%\n\n", detection_stats$summary_stats$sd_detection_rate * 100))

  # Individual path performance
  cat("DETECTION RATE BY INDIVIDUAL PATH:\n")
  for (i in 1:nrow(detection_stats$by_path)) {
    path_data <- detection_stats$by_path[i, ]
    cat(sprintf("  Path %d: %d/%d steps detected (%.1f%%)\n",
                path_data$path_id,
                path_data$detected_steps,
                path_data$total_steps,
                path_data$detection_rate * 100))
  }
  cat("\n")
}

#' Create visualization plots of detection rate performance
#'
#' Generates multiple plots to visualize detection performance including
#' per-path detection rates, distribution of rates, and time series analysis.
#'
#' @param detection_stats List output from \code{\link{calculate_detection_summaries}}
#'   containing detection performance metrics.
#' @param simulation_results Optional. List output from \code{\link{simulate_fish_tracks}}.
#'   If provided, creates additional time series plots. Default is NULL.
#'
#' @return A list of ggplot2 objects:
#'   \item{by_path}{Bar chart showing detection rate by individual path}
#'   \item{distribution}{Histogram of detection rates across all paths}
#'   \item{time_series}{Cumulative detection rate over time (if simulation_results provided)}
#'
#'   Returns NULL if no detection statistics are available.
#'
#' @details
#' The function creates three types of visualizations:
#' \enumerate{
#'   \item Bar chart comparing detection rates across individual fish paths
#'   \item Histogram showing the distribution of detection rates
#'   \item Time series plot showing how cumulative detection rates evolve during tracks
#' }
#'
#' All plots include reference lines showing mean detection rates for context.
#'
#' @examples
#' \dontrun{
#' # Calculate detection statistics
#' detection_stats <- calculate_detection_summaries(fish_simulation)
#'
#' # Create all available plots
#' plots <- plot_detection_rates(detection_stats, fish_simulation)
#'
#' # Display individual plots
#' print(plots$by_path)
#' print(plots$distribution)
#' print(plots$time_series)
#'
#' # Create plots without time series
#' basic_plots <- plot_detection_rates(detection_stats)
#' }
#'
#' @seealso \code{\link{calculate_detection_summaries}}, \code{\link{simulate_fish_tracks}}
#'
#' @export
plot_detection_rates <- function(detection_stats, simulation_results = NULL) {

  if (is.null(detection_stats)) {
    cat("No detection statistics to plot.\n")
    return(NULL)
  }

  library(ggplot2)

  # Create plots
  plot_list <- list()

  # 1. Bar plot of detection rates by path
  path_data <- detection_stats$by_path

  p1 <- ggplot2::ggplot(path_data, ggplot2::aes(x = factor(path_id), y = detection_rate)) +
    ggplot2::geom_col(fill = "steelblue", alpha = 0.7) +
    ggplot2::geom_hline(yintercept = detection_stats$summary_stats$mean_detection_rate,
                        color = "red", linetype = "dashed", size = 1) +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    ggplot2::labs(
      title = "Detection Rate by Fish Path",
      subtitle = paste0("Red line = mean detection rate (",
                        round(detection_stats$summary_stats$mean_detection_rate * 100, 1), "%)"),
      x = "Path ID",
      y = "Detection Rate"
    ) +
    ggplot2::theme_minimal()

  plot_list$by_path <- p1

  # 2. Histogram of detection rates
  p2 <- ggplot2::ggplot(path_data, ggplot2::aes(x = detection_rate)) +
    ggplot2::geom_histogram(bins = max(5, nrow(path_data)/2), fill = "lightblue",
                            color = "black", alpha = 0.7) +
    ggplot2::geom_vline(xintercept = detection_stats$summary_stats$mean_detection_rate,
                        color = "red", linetype = "dashed", size = 1) +
    ggplot2::scale_x_continuous(labels = scales::percent) +
    ggplot2::labs(
      title = "Distribution of Detection Rates Across Paths",
      subtitle = paste0("Red line = mean (",
                        round(detection_stats$summary_stats$mean_detection_rate * 100, 1), "%)"),
      x = "Detection Rate",
      y = "Number of Paths"
    ) +
    ggplot2::theme_minimal()

  plot_list$distribution <- p2

  # 3. Time series of detections (if simulation results provided)
  if (!is.null(simulation_results) && !is.null(simulation_results$detections)) {

    detection_data <- simulation_results$detections %>%
      dplyr::group_by(path_id) %>%
      dplyr::arrange(time_seconds) %>%
      dplyr::mutate(
        cumulative_detections = cumsum(detected),
        cumulative_rate = cumulative_detections / dplyr::row_number()
      )

    p3 <- ggplot2::ggplot(detection_data, ggplot2::aes(x = time_seconds, y = cumulative_rate,
                                                       color = factor(path_id))) +
      ggplot2::geom_line(size = 1, alpha = 0.8) +
      ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      ggplot2::scale_color_discrete(name = "Path ID") +
      ggplot2::labs(
        title = "Cumulative Detection Rate Over Time",
        subtitle = "How detection rate evolves during each fish track",
        x = "Time (seconds)",
        y = "Cumulative Detection Rate"
      ) +
      ggplot2::theme_minimal()

    plot_list$time_series <- p3
  }

  # Return plots
  return(plot_list)
}

#' Comprehensive detection performance analysis with summary and plots
#'
#' Performs complete analysis of detection performance from fish track simulations,
#' including statistical summaries and visualization plots. This is a convenience
#' function that combines calculation, printing, and plotting in one call.
#'
#' @param simulation_results List output from \code{\link{simulate_fish_tracks}}
#'   containing tracks, detections, and station_detections data frames.
#' @param create_plots Logical. Whether to generate and display visualization plots.
#'   Default is TRUE.
#'
#' @return A list containing:
#'   \item{statistics}{Complete detection statistics from \code{\link{calculate_detection_summaries}}}
#'   \item{plots}{List of ggplot2 objects from \code{\link{plot_detection_rates}} (if create_plots = TRUE)}
#'
#'   Returns NULL if no detection data is available.
#'
#' @details
#' This function provides a complete workflow for detection analysis:
#' \enumerate{
#'   \item Calculates detection performance statistics
#'   \item Prints formatted summary report to console
#'   \item Creates and displays visualization plots (optional)
#'   \item Returns all results for further analysis
#' }
#'
#' The function is designed for quick analysis and reporting of receiver array
#' performance, making it easy to evaluate different system configurations.
#'
#' @examples
#' \dontrun{
#' # Complete analysis with plots
#' analysis <- analyze_detection_performance(fish_simulation)
#'
#' # Access overall detection rate
#' overall_rate <- analysis$statistics$overall$detection_percentage
#'
#' # Access individual plots
#' path_plot <- analysis$plots$by_path
#' distribution_plot <- analysis$plots$distribution
#'
#' # Analysis without plots (faster)
#' stats_only <- analyze_detection_performance(fish_simulation, create_plots = FALSE)
#' }
#'
#' @seealso \code{\link{simulate_fish_tracks}}, \code{\link{calculate_detection_summaries}},
#'   \code{\link{print_detection_summary}}, \code{\link{plot_detection_rates}}
#'
#' @export
analyze_detection_performance <- function(simulation_results, create_plots = TRUE) {

  cat("Analyzing detection performance...\n\n")

  # Calculate summaries
  detection_stats <- calculate_detection_summaries(simulation_results)

  if (is.null(detection_stats)) {
    return(NULL)
  }

  # Print summary
  print_detection_summary(detection_stats)

  # Create plots if requested
  plots <- NULL
  if (create_plots) {
    cat("Creating visualization plots...\n")
    plots <- plot_detection_rates(detection_stats, simulation_results)

    # Display plots
    if (!is.null(plots$by_path)) print(plots$by_path)
    if (!is.null(plots$distribution)) print(plots$distribution)
    if (!is.null(plots$time_series)) print(plots$time_series)
  }

  # Return everything
  return(list(
    statistics = detection_stats,
    plots = plots
  ))
}
