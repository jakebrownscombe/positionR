#' Calculate fish positions from acoustic telemetry data
#'
#' Estimates fish positions by combining detection and non-detection data from
#' acoustic telemetry arrays. The function aggregates detection events into time
#' bins, models detection efficiency, and calculates weighted position probabilities.
#'
#' @param station_detections A data frame containing detection data with fish tracks
#'   and detection events at receiver stations.
#' @param station_distances_df A data frame containing pre-calculated distances
#'   between receiver stations and spatial grid cells, typically from
#'   \code{\link{calculate_station_distances}}.
#' @param receiver_stations An sf object containing receiver station locations
#'   and metadata, typically from point generation functions.
#' @param de_model A fitted detection efficiency model object (e.g., from
#'   \code{\link{create_logistic_curve_depth}}). The model should accept
#'   dist_m and depth_m as predictors. Default is NULL, which requires
#'   DE_pred column to already exist in station_distances_df.
#' @param bin_size_seconds Numeric. Time bin size in seconds for aggregating
#'   detections. Default is 3600 (1 hour).
#' @param detection_weight Numeric. Weight given to detection events in the
#'   integrated probability calculation (0-1). Default is 0.5.
#' @param non_detection_weight Numeric. Weight given to non-detection events
#'   in the integrated probability calculation (0-1). Default is 0.5.
#' @param max_non_detection_distance Numeric. Maximum distance (in meters) from
#'   detecting stations to consider non-detecting stations. Set to NULL to
#'   include all stations. Default is 2000.
#' @param normalization_method Character. Method for normalizing detection
#'   efficiency values. Options are "min_max", "z_score", or "robust".
#'   Default is "min_max".
#' @param fish_id_col Character. Name of the column containing fish identifiers.
#'   Default is "path_id".
#' @param time_col Character. Name of the column containing time values in seconds.
#'   Default is "time_seconds".
#' @param station_col Character. Name of the column containing station identifiers.
#'   Default is "station_id".
#' @param verbose Logical. Whether to print progress messages. Default is TRUE.
#'
#' @return A list containing:
#'   \item{position_probabilities}{Data frame with integrated position probabilities
#'     for each fish, time period, and spatial cell}
#'   \item{detection_data}{Data frame with processed detection probability data}
#'   \item{non_detection_data}{Data frame with processed non-detection probability data}
#'   \item{station_detections_binned}{Data frame with time-binned detection events
#'     and station coordinates}
#'   \item{station_coordinates}{Data frame with station coordinate information}
#'   \item{summary}{List with summary statistics of the positioning analysis}
#'
#' @details
#' The function implements a multi-step positioning algorithm:
#' \enumerate{
#'   \item Time binning of detection events
#'   \item Aggregation of detection data with spatial distance information
#'   \item Creation of non-detection events for nearby stations
#'   \item Normalization of detection efficiency across receivers
#'   \item Integration of detection and non-detection probabilities
#'   \item Calculation of weighted position estimates
#' }
#'
#' Detection and non-detection weights must sum to 1. The algorithm focuses
#' non-detection analysis on stations within a specified distance of detecting
#' stations to maintain biological realism and computational efficiency.
#'
#' @examples
#' \dontrun{
#' # Basic positioning analysis
#' results <- calculate_fish_positions(
#'   station_detections = fish_tracks$detections,
#'   station_distances_df = distances,
#'   receiver_stations = stations,
#'   bin_size_seconds = 3600
#' )
#'
#' # Custom weighting and filtering
#' results <- calculate_fish_positions(
#'   station_detections = fish_tracks$detections,
#'   station_distances_df = distances,
#'   receiver_stations = stations,
#'   detection_weight = 0.4,
#'   non_detection_weight = 0.6,
#'   max_non_detection_distance = 1500
#' )
#' }
#'
#' @seealso \code{\link{plot_fish_positions}}, \code{\link{analyze_positioning_performance}}
#'
#' @export
calculate_fish_positions <- function(station_detections,
                                     station_distances_df,
                                     receiver_stations,
                                     de_model = NULL,
                                     bin_size_seconds = 3600,
                                     detection_weight = 0.5,
                                     non_detection_weight = 0.5,
                                     max_non_detection_distance = 2000,
                                     normalization_method = "min_max",
                                     fish_id_col = "path_id",
                                     time_col = "time_seconds",
                                     station_col = "station_id",
                                     verbose = TRUE) {

  # Check required packages
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # Validate inputs
  if (abs(detection_weight + non_detection_weight - 1) > 1e-10) {
    stop("detection_weight and non_detection_weight must sum to 1")
  }

  if (!normalization_method %in% c("min_max", "z_score", "robust")) {
    stop("normalization_method must be 'min_max', 'z_score', or 'robust'")
  }

  # Check required columns
  required_cols <- c(fish_id_col, time_col, station_col)
  missing_cols <- setdiff(required_cols, names(station_detections))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns in station_detections:", paste(missing_cols, collapse = ", ")))
  }

  if (verbose) cat("=== CALCULATING FISH POSITIONS ===\n")

  # Step 1: Create time bins
  if (verbose) cat("Step 1: Creating time bins...\n")
  station_detections_binned <- station_detections %>%
    dplyr::mutate(
      time_bin = floor(!!dplyr::sym(time_col) / bin_size_seconds) * bin_size_seconds,
      .after = dplyr::all_of(time_col)
    )

  # Step 2: Aggregate detections for prediction
  if (verbose) cat("Step 2: Aggregating detection data...\n")
  # Step 3: Aggregate detections for prediction
  if (verbose) cat("Step 3: Aggregating detection data...\n")
  detection_probs <- aggregate_detections_for_prediction(
    station_detections = station_detections_binned,
    station_distances_df = station_distances_df,
    fish_id_col = fish_id_col,
    time_col = "time_bin",
    station_col = station_col
  )

  # Step 4: Create non-detections
  if (verbose) cat("Step 4: Creating non-detection data...\n")
  non_detections <- create_non_detections(
    station_detections = station_detections_binned,
    points_regular = receiver_stations,
    max_distance_from_detecting = max_non_detection_distance
  )

  # Step 5: Aggregate non-detections
  if (verbose) cat("Step 5: Aggregating non-detection data...\n")
  non_detection_probs <- aggregate_non_detections(
    non_detections = non_detections,
    station_distances_df = station_distances_df,
    fish_id_col = fish_id_col,
    time_col = "time_bin",
    station_col = station_col
  )

  # Step 6: Normalize detection efficiency by station
  if (verbose) cat("Step 6: Normalizing detection efficiency...\n")
  detection_probs_norm <- normalize_DE_by_station(
    data = detection_probs,
    DE_col = "DE_pred",
    station_col = station_col,
    method = normalization_method
  )

  non_detection_probs_norm <- normalize_DE_by_station(
    data = non_detection_probs,
    DE_col = "DE_pred",
    station_col = station_col,
    method = normalization_method
  )

  # Step 7: Combine detection and non-detection data
  if (verbose) cat("Step 7: Combining detection and non-detection data...\n")
  position_probs_combined <- dplyr::bind_rows(
    detection_probs_norm %>% dplyr::mutate(type = "detection"),
    non_detection_probs_norm %>% dplyr::mutate(type = "non-detection")
  )

  # Step 8: Calculate integrated positioning probabilities
  if (verbose) cat("Step 8: Calculating position probabilities...\n")
  position_probs <- aggregate_probability(
    df = position_probs_combined,
    detection_weight = detection_weight,
    non_detection_weight = non_detection_weight,
    normalize_method = "global"
  )

  # Step 9: Add station coordinates for plotting
  if (verbose) cat("Step 9: Adding station coordinates...\n")

  # Extract station coordinates (handle both sf and data.frame formats)
  if ("sf" %in% class(receiver_stations)) {
    station_coords <- receiver_stations %>%
      sf::st_drop_geometry() %>%
      dplyr::bind_cols(sf::st_coordinates(receiver_stations)) %>%
      dplyr::select(station_id = dplyr::contains("point_id"), station_x = X, station_y = Y)
  } else {
    station_coords <- receiver_stations %>%
      dplyr::select(station_id = dplyr::contains("point_id"), station_x = x, station_y = y)
  }

  # Add station coordinates to binned detections for plotting
  join_by <- stats::setNames("station_id", station_col)
  station_detections_plot <- station_detections_binned %>%
    dplyr::left_join(station_coords, by = join_by)

  # Also add coordinates to detection and non-detection data
  detection_probs_norm <- detection_probs_norm %>%
    dplyr::left_join(station_coords, by = c("station_id" = "station_id"))

  non_detection_probs_norm <- non_detection_probs_norm %>%
    dplyr::left_join(station_coords, by = c("station_id" = "station_id"))

  # Create summary statistics
  summary_stats <- list(
    n_fish = dplyr::n_distinct(position_probs$fish_id),
    n_time_periods = dplyr::n_distinct(position_probs$time_period),
    n_cells = dplyr::n_distinct(position_probs$cell_id),
    n_stations = dplyr::n_distinct(station_coords$station_id),
    time_bin_size = bin_size_seconds,
    detection_weight = detection_weight,
    non_detection_weight = non_detection_weight,
    total_detections = nrow(station_detections_binned),
    total_position_estimates = nrow(position_probs)
  )

  if (verbose) {
    cat("=== POSITIONING COMPLETE ===\n")
    cat("Fish tracked:", summary_stats$n_fish, "\n")
    cat("Time periods:", summary_stats$n_time_periods, "\n")
    cat("Spatial cells:", summary_stats$n_cells, "\n")
    cat("Position estimates generated:", summary_stats$total_position_estimates, "\n")
  }

  # Return comprehensive results
  return(list(
    position_probabilities = position_probs,
    detection_data = detection_probs_norm,
    non_detection_data = non_detection_probs_norm,
    station_detections_binned = station_detections_plot,
    station_coordinates = station_coords,
    summary = summary_stats
  ))
}


#' Plot fish position estimates from acoustic telemetry
#'
#' Creates visualization plots of fish position estimates showing detection
#' probabilities, non-detection probabilities, and integrated positioning results.
#'
#' @param positioning_results A list returned by \code{\link{calculate_fish_positions}}
#'   containing position probabilities and associated data.
#' @param depth_raster_df Optional depth/bathymetry data for background visualization.
#'   Can be either a data frame with 'x', 'y' columns or a raster object which will be
#'   automatically converted to a data frame.
#' @param fish_select Numeric. Fish ID to plot. Default is 1.
#' @param time_select Numeric. Time bin to plot. Default is 0.
#' @param prob_threshold Numeric. Minimum probability threshold for display (0-1).
#'   Cells below this threshold are not plotted. Default is 0.05.
#' @param detection_threshold Numeric. Minimum detection probability threshold (0-1)
#'   for displaying integrated probabilities. Only cells with detection probability
#'   above this threshold will show integrated probability values. This prevents
#'   artificially high integrated probabilities in areas with no acoustic coverage.
#'   Default is 0.05.
#' @param xlim Numeric vector of length 2. X-axis limits for the plot. Default is NULL.
#' @param ylim Numeric vector of length 2. Y-axis limits for the plot. Default is NULL.
#' @param plot_type Character. Type of plot(s) to create. Options are:
#'   \itemize{
#'     \item "detection" - Detection probability only
#'     \item "non_detection" - Non-detection probability only
#'     \item "integrated" - Integrated position probability only
#'     \item "all" - All three plots combined vertically
#'   }
#'   Default is "all".
#' @param return_list Logical. If TRUE, returns a named list of individual ggplot
#'   objects instead of a combined plot. Default is FALSE.
#'
#' @return Depending on \code{plot_type} and \code{return_list}:
#'   \itemize{
#'     \item Single ggplot object (when plot_type is specific type)
#'     \item Combined patchwork plot (when plot_type = "all" and return_list = FALSE)
#'     \item Named list of ggplot objects (when return_list = TRUE)
#'   }
#'
#' @details
#' The function creates up to three types of visualizations:
#' \enumerate{
#'   \item Detection probability plot showing weighted mean detection efficiency
#'   \item Non-detection probability plot showing non-detection patterns
#'   \item Integrated probability plot combining both detection and non-detection data
#' }
#'
#' Plot elements include:
#' \itemize{
#'   \item Background bathymetry (if provided)
#'   \item Probability surfaces as colored rasters
#'   \item Detecting stations as yellow circles (size = number of detections)
#'   \item Non-detecting stations as red circles
#'   \item Actual fish positions as green circles (if available)
#' }
#'
#' @examples
#' \dontrun{
#' # Calculate positions
#' results <- calculate_fish_positions(station_detections, distances, stations)
#'
#' # Create all plots combined
#' plot_fish_positions(results, depth_raster_df = depth_data)
#'
#' # Create only integrated probability plot
#' plot_fish_positions(results, plot_type = "integrated")
#'
#' # Get individual plots as list
#' plots <- plot_fish_positions(results, return_list = TRUE)
#' plots$detection
#' plots$integrated
#'
#' # Custom thresholds
#' plot_fish_positions(
#'   results,
#'   fish_select = 2,
#'   time_select = 3600,
#'   prob_threshold = 0.1,
#'   detection_threshold = 0.02,  # Lower detection threshold
#'   xlim = c(722000, 728000),
#'   ylim = c(4935000, 4940000)
#' )
#' }
#'
#' @seealso \code{\link{calculate_fish_positions}}, \code{\link{analyze_positioning_performance}}
#'
#' @export
plot_fish_positions <- function(positioning_results,
                                depth_raster_df = NULL,
                                fish_select = 1,
                                time_select = 0,
                                prob_threshold = 0.05,
                                detection_threshold = 0.05,
                                xlim = NULL,
                                ylim = NULL,
                                plot_type = "all",
                                return_list = FALSE) {

  # Check required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("raster", quietly = TRUE)) {
    stop("Package 'raster' needed for raster conversion. Please install it.",
         call. = FALSE)
  }
  if (plot_type == "all" && !return_list && !requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' needed for combined plots. Please install it or set return_list = TRUE.",
         call. = FALSE)
  }

  # Validate plot_type
  valid_types <- c("detection", "non_detection", "integrated", "all")
  if (!plot_type %in% valid_types) {
    stop("plot_type must be one of: ", paste(valid_types, collapse = ", "))
  }

  # Extract data from results
  position_probs <- positioning_results$position_probabilities
  station_detections <- positioning_results$station_detections_binned
  detection_data <- positioning_results$detection_data
  non_detection_data <- positioning_results$non_detection_data
  station_coords <- positioning_results$station_coordinates

  # Get the correct fish_id column name
  fish_id_col <- names(station_detections)[grepl("path_id|fish_id", names(station_detections))][1]

  # Filter data for selected fish and time
  position_data <- position_probs %>%
    dplyr::filter(fish_id == fish_select, time_period == time_select)

  # Get detection summary for this fish/time
  detection_summary <- station_detections %>%
    dplyr::filter(!!dplyr::sym(fish_id_col) == fish_select, time_bin == time_select) %>%
    dplyr::group_by(station_id) %>%
    dplyr::summarise(
      n_detections = sum(detected, na.rm = TRUE),
      station_x = dplyr::first(station_x, na_rm = TRUE),
      station_y = dplyr::first(station_y, na_rm = TRUE),
      .groups = 'drop'
    ) %>%
    dplyr::filter(n_detections > 0, !is.na(station_x), !is.na(station_y))

  # Get non-detecting stations for this fish/time
  non_detection_summary <- non_detection_data %>%
    dplyr::filter(fish_id == fish_select, time_period == time_select) %>%
    dplyr::group_by(station_id) %>%
    dplyr::summarise(
      station_x = dplyr::first(station_x, na_rm = TRUE),
      station_y = dplyr::first(station_y, na_rm = TRUE),
      .groups = 'drop'
    ) %>%
    dplyr::filter(!is.na(station_x), !is.na(station_y))

  # Get actual fish position if available
  actual_position <- station_detections %>%
    dplyr::filter(!!dplyr::sym(fish_id_col) == fish_select, time_bin == time_select) %>%
    dplyr::slice(1) %>%
    dplyr::select(x, y)

  # Get fish track for context (show some track around the time period)
  fish_track <- station_detections %>%
    dplyr::filter(!!dplyr::sym(fish_id_col) == fish_select) %>%
    dplyr::select(x, y, time_bin, !!dplyr::sym(fish_id_col)) %>%
    dplyr::distinct() %>%
    dplyr::arrange(time_bin)

  # Helper function to create base plot elements
  create_base_plot <- function() {
    p <- ggplot2::ggplot()

    # Add depth/bathymetry if provided
    if (!is.null(depth_raster_df)) {
      # Convert raster to dataframe if needed
      if (is.data.frame(depth_raster_df)) {
        raster_df <- depth_raster_df
      } else {
        # Convert raster to dataframe for plotting
        raster_df <- raster::as.data.frame(depth_raster_df, xy = TRUE)
      }

      # Remove NA values to avoid blue squares
      raster_df <- raster_df[!is.na(raster_df[[3]]), ]

      # Determine the value column name
      value_col <- if ("layer" %in% names(raster_df)) {
        "layer"
      } else if (ncol(raster_df) > 2) {
        names(raster_df)[3]  # Use third column
      } else {
        NULL
      }

      # Add raster layer with actual depth values but muted colors
      if (!is.null(value_col) && nrow(raster_df) > 0) {
        p <- p + ggplot2::geom_raster(data = raster_df,
                                      ggplot2::aes(x = x, y = y, alpha = !!ggplot2::sym(value_col))) +
          ggplot2::scale_alpha_continuous(range = c(0.1, 0.4), guide = "none") +
          ggplot2::geom_raster(data = raster_df,
                               ggplot2::aes(x = x, y = y),
                               fill = "lightblue")
      }
    }

    return(p)
  }

  # Helper function to add station and fish elements
  add_station_elements <- function(p) {
    # Add fish track first (so it appears behind other elements)
    if (nrow(fish_track) > 1) {
      p <- p +
        ggplot2::geom_path(data = fish_track,
                           ggplot2::aes(x = x, y = y),
                           color = "green", alpha = 0.8, size = 1.2) +
        ggplot2::geom_point(data = fish_track,
                            ggplot2::aes(x = x, y = y),
                            color = "green", size = 1.2, alpha = 0.8)
    }

    # Add detecting stations
    if (nrow(detection_summary) > 0) {
      p <- p +
        ggplot2::geom_point(data = detection_summary,
                            ggplot2::aes(x = station_x, y = station_y, size = n_detections),
                            color = "yellow", shape = 21, stroke = 1.5) +
        ggplot2::scale_size_continuous(name = "Detections", range = c(2, 8))
    }

    # Add non-detecting stations
    if (nrow(non_detection_summary) > 0) {
      p <- p +
        ggplot2::geom_point(data = non_detection_summary,
                            ggplot2::aes(x = station_x, y = station_y),
                            color = "red", size = 2, alpha = 0.7)
    }

    return(p)
  }

  # Helper function to apply final formatting
  apply_formatting <- function(p, title) {
    # Calculate automatic zoom based on stations if xlim/ylim not provided
    if (is.null(xlim) || is.null(ylim)) {
      # Get all station coordinates for this fish/time
      all_stations <- rbind(
        detection_summary %>% dplyr::select(x = station_x, y = station_y),
        non_detection_summary %>% dplyr::select(x = station_x, y = station_y)
      ) %>% dplyr::distinct()

      if (nrow(all_stations) > 0) {
        # Calculate extent with buffer (10% of range on each side)
        x_range <- range(all_stations$x, na.rm = TRUE)
        y_range <- range(all_stations$y, na.rm = TRUE)

        x_buffer <- diff(x_range) * 0.15  # 15% buffer
        y_buffer <- diff(y_range) * 0.15  # 15% buffer

        auto_xlim <- c(x_range[1] - x_buffer, x_range[2] + x_buffer)
        auto_ylim <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)

        # Use automatic limits if not provided by user
        final_xlim <- if (is.null(xlim)) auto_xlim else xlim
        final_ylim <- if (is.null(ylim)) auto_ylim else ylim

        p <- p + ggplot2::coord_sf(xlim = final_xlim, ylim = final_ylim)
      }
    } else {
      # Use user-provided limits
      p <- p + ggplot2::coord_sf(xlim = xlim, ylim = ylim)
    }

    # Apply theme and labels
    p <- p +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = title,
        subtitle = paste("Fish", fish_select, "| Time bin:", time_select, "| Actual track = green"),
        x = "X Coordinate",
        y = "Y Coordinate"
      )

    return(p)
  }

  # Create plots based on plot_type to avoid unnecessary computation
  plots_to_create <- switch(plot_type,
                            "detection" = "detection",
                            "non_detection" = "non_detection",
                            "integrated" = "integrated",
                            "all" = c("detection", "non_detection", "integrated")
  )

  plot_list_result <- list()

  # Plot 1: Detection probability (only if needed)
  if ("detection" %in% plots_to_create) {
    plot_detection <- create_base_plot()

    detection_plot_data <- position_data %>%
      dplyr::filter(!is.na(weighted_mean_DE_normalized),
                    weighted_mean_DE_normalized > prob_threshold)

    if (nrow(detection_plot_data) > 0) {
      plot_detection <- plot_detection +
        ggplot2::geom_raster(data = detection_plot_data,
                             ggplot2::aes(x = x, y = y, fill = weighted_mean_DE_normalized)) +
        ggplot2::scale_fill_viridis_c(option = "magma", name = "Detection\nProbability")
    }

    plot_detection <- add_station_elements(plot_detection)
    plot_detection <- apply_formatting(plot_detection, "Detection Probability")
    plot_list_result$detection <- plot_detection
  }

  # Plot 2: Non-detection probability (only if needed)
  if ("non_detection" %in% plots_to_create) {
    plot_non_detection <- create_base_plot()

    non_detection_plot_data <- position_data %>%
      dplyr::filter(!is.na(non_det_DE_normalized))

    if (nrow(non_detection_plot_data) > 0) {
      plot_non_detection <- plot_non_detection +
        ggplot2::geom_raster(data = non_detection_plot_data,
                             ggplot2::aes(x = x, y = y, fill = non_det_DE_normalized)) +
        ggplot2::scale_fill_viridis_c(name = "Non-Detection\nProbability")
    }

    plot_non_detection <- add_station_elements(plot_non_detection)
    plot_non_detection <- apply_formatting(plot_non_detection, "Non-Detection Probability")
    plot_list_result$non_detection <- plot_non_detection
  }

  # Plot 3: Integrated probability (only if needed)
  if ("integrated" %in% plots_to_create) {
    plot_integrated <- create_base_plot()

    integrated_plot_data <- position_data %>%
      dplyr::filter(!is.na(integrated_prob),
                    !is.na(weighted_mean_DE_normalized),
                    integrated_prob > prob_threshold,
                    weighted_mean_DE_normalized > detection_threshold)  # Only show where detection prob > threshold

    if (nrow(integrated_plot_data) > 0) {
      plot_integrated <- plot_integrated +
        ggplot2::geom_raster(data = integrated_plot_data,
                             ggplot2::aes(x = x, y = y, fill = integrated_prob)) +
        ggplot2::scale_fill_viridis_c(option = "magma", name = "Integrated\nProbability")
    }

    plot_integrated <- add_station_elements(plot_integrated)
    plot_integrated <- apply_formatting(plot_integrated, "Integrated Position Probability")
    plot_list_result$integrated <- plot_integrated
  }

  # Return based on plot_type and return_list parameters
  if (return_list) {
    return(plot_list_result)
  }

  if (plot_type == "detection") {
    return(plot_list_result$detection)
  } else if (plot_type == "non_detection") {
    return(plot_list_result$non_detection)
  } else if (plot_type == "integrated") {
    return(plot_list_result$integrated)
  } else if (plot_type == "all") {
    # Create combined plot using patchwork
    combined_plot <- plot_list_result$detection / plot_list_result$non_detection / plot_list_result$integrated
    return(combined_plot)
  }
}


#' Analyze positioning performance from acoustic telemetry
#'
#' Calculates performance metrics for fish positioning estimates to evaluate
#' the quality and reliability of position estimates.
#'
#' @param positioning_results A list returned by \code{\link{calculate_fish_positions}}
#'   containing position probabilities and summary statistics.
#'
#' @return A list containing:
#'   \item{by_fish_time}{Data frame with performance metrics for each fish-time combination:
#'     \itemize{
#'       \item max_prob - Maximum probability value in the position estimate
#'       \item mean_prob - Mean probability across all cells
#'       \item prob_95_quantile - 95th percentile of probability values
#'       \item n_high_prob_cells - Number of cells with probability > 0.1
#'       \item total_cells - Total number of cells with position estimates
#'     }}
#'   \item{overall}{List with overall performance statistics:
#'     \itemize{
#'       \item mean_max_prob - Average of maximum probabilities across all estimates
#'       \item mean_concentration - Average probability concentration ratio
#'       \item positioning_success_rate - Proportion of estimates with max_prob > 0.1
#'     }}
#'   \item{summary_stats}{Summary statistics from the original positioning analysis}
#'
#' @details
#' Performance metrics help evaluate positioning quality:
#' \itemize{
#'   \item High max_prob values indicate confident position estimates
#'   \item Low concentration ratios suggest widely distributed probability
#'   \item Success rate shows proportion of meaningful position estimates
#' }
#'
#' The 0.1 probability threshold for "high probability cells" and success rate
#' can be adjusted based on the specific application and required confidence levels.
#'
#' @examples
#' \dontrun{
#' # Calculate positions
#' results <- calculate_fish_positions(station_detections, distances, stations)
#'
#' # Analyze performance
#' performance <- analyze_positioning_performance(results)
#'
#' # View overall performance
#' performance$overall
#'
#' # View detailed metrics
#' head(performance$by_fish_time)
#'
#' # Get summary statistics
#' performance$summary_stats
#' }
#'
#' @seealso \code{\link{calculate_fish_positions}}, \code{\link{plot_fish_positions}}
#'
#' @export
analyze_positioning_performance <- function(positioning_results) {

  # Check required packages
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' needed for this function to work. Please install it.",
         call. = FALSE)
  }

  position_probs <- positioning_results$position_probabilities
  summary_stats <- positioning_results$summary

  # Calculate performance metrics by fish and time
  performance_summary <- position_probs %>%
    dplyr::group_by(fish_id, time_period) %>%
    dplyr::summarise(
      max_prob = max(integrated_prob, na.rm = TRUE),
      mean_prob = mean(integrated_prob, na.rm = TRUE),
      prob_95_quantile = stats::quantile(integrated_prob, 0.95, na.rm = TRUE),
      n_high_prob_cells = sum(integrated_prob > 0.1, na.rm = TRUE),
      total_cells = dplyr::n(),
      .groups = 'drop'
    )

  # Overall performance metrics
  overall_performance <- list(
    mean_max_prob = mean(performance_summary$max_prob, na.rm = TRUE),
    mean_concentration = mean(performance_summary$n_high_prob_cells / performance_summary$total_cells, na.rm = TRUE),
    positioning_success_rate = mean(performance_summary$max_prob > 0.1, na.rm = TRUE)
  )

  cat("=== POSITIONING PERFORMANCE ANALYSIS ===\n")
  cat("Mean maximum probability per estimate:", round(overall_performance$mean_max_prob, 3), "\n")
  cat("Mean probability concentration:", round(overall_performance$mean_concentration, 3), "\n")
  cat("Positioning success rate (>0.1 max prob):", round(overall_performance$positioning_success_rate, 3), "\n")

  return(list(
    by_fish_time = performance_summary,
    overall = overall_performance,
    summary_stats = summary_stats
  ))
}


# Helper functions (not exported) ----

# Helper function to create time bins for aggregation
create_time_bins <- function(station_detections,
                             time_col = "time_seconds",
                             bin_size_seconds = 300) {

  station_detections %>%
    dplyr::mutate(
      time_bin = floor(!!dplyr::sym(time_col) / bin_size_seconds) * bin_size_seconds,
      .after = dplyr::all_of(time_col)
    )
}

# Function to aggregate station detections and generate cells for prediction
aggregate_detections_for_prediction <- function(station_detections,
                                                station_distances_df,
                                                fish_id_col = "path_id",
                                                time_col = "time_seconds",
                                                station_col = "station_id") {

  # Validate inputs
  if (!fish_id_col %in% names(station_detections)) {
    stop(paste("Column", fish_id_col, "not found in station_detections"))
  }
  if (!time_col %in% names(station_detections)) {
    stop(paste("Column", time_col, "not found in station_detections"))
  }
  if (!station_col %in% names(station_detections)) {
    stop(paste("Column", station_col, "not found in station_detections"))
  }

  # Rename columns to standard names for processing
  station_detections_renamed <- station_detections %>%
    dplyr::rename(
      fish_id = dplyr::all_of(fish_id_col),
      time_period = dplyr::all_of(time_col),
      station_id = dplyr::all_of(station_col)
    )

  # Aggregate detections by fish, time, and station
  detection_summary <- station_detections_renamed %>%
    dplyr::group_by(fish_id, time_period, station_id) %>%
    dplyr::summarise(
      n_detections = dplyr::n(),
      first_detection_time = min(time_period),
      last_detection_time = max(time_period),
      mean_detection_prob = mean(detection_prob, na.rm = TRUE),
      total_distance = mean(distance_to_station, na.rm = TRUE),
      .groups = 'drop'
    )

  # Merge with station_distances_df to get all possible cells for each combination
  prediction_data <- detection_summary %>%
    dplyr::left_join(
      station_distances_df,
      by = c("station_id" = "station_no"),
      relationship = "many-to-many"
    ) %>%
    dplyr::filter(!is.na(cell_id)) %>%
    dplyr::select(
      fish_id, time_period, station_id, n_detections,
      first_detection_time, last_detection_time, mean_detection_prob,
      total_distance, cell_id, x, y, raster_value, cost_distance,
      straight_distance, tortuosity, DE_pred
    ) %>%
    dplyr::arrange(fish_id, time_period, station_id, cell_id)

  return(prediction_data)
}

# Function to normalize DE values by station
normalize_DE_by_station <- function(data,
                                    DE_col = "DE_pred",
                                    station_col = "station_id",
                                    method = "min_max") {

  # Validate method
  if (!method %in% c("min_max", "z_score", "robust")) {
    stop("method must be 'min_max', 'z_score', or 'robust'")
  }

  # Create normalized column name
  normalized_col <- paste0(DE_col, "_normalized")

  if (method == "min_max") {
    result <- data %>%
      dplyr::group_by(!!dplyr::sym(station_col)) %>%
      dplyr::mutate(
        !!normalized_col := (!!dplyr::sym(DE_col) - min(!!dplyr::sym(DE_col), na.rm = TRUE)) /
          (max(!!dplyr::sym(DE_col), na.rm = TRUE) - min(!!dplyr::sym(DE_col), na.rm = TRUE))
      ) %>%
      dplyr::ungroup()
  } else if (method == "z_score") {
    result <- data %>%
      dplyr::group_by(!!dplyr::sym(station_col)) %>%
      dplyr::mutate(
        !!normalized_col := (!!dplyr::sym(DE_col) - mean(!!dplyr::sym(DE_col), na.rm = TRUE)) /
          stats::sd(!!dplyr::sym(DE_col), na.rm = TRUE)
      ) %>%
      dplyr::ungroup()
  } else if (method == "robust") {
    result <- data %>%
      dplyr::group_by(!!dplyr::sym(station_col)) %>%
      dplyr::mutate(
        !!normalized_col := (!!dplyr::sym(DE_col) - stats::median(!!dplyr::sym(DE_col), na.rm = TRUE)) /
          stats::mad(!!dplyr::sym(DE_col), na.rm = TRUE)
      ) %>%
      dplyr::ungroup()
  }

  return(result)
}

# Create dataset of non-detecting receivers for each fish-time combination
create_non_detections <- function(station_detections, points_regular, max_distance_from_detecting = 2000) {

  # Get all unique fish-time combinations from station_detections
  fish_time_data <- station_detections %>%
    dplyr::select(path_id, step, time_seconds, time_bin, x, y) %>%
    dplyr::distinct()

  # Get all available stations from points_regular
  all_stations <- points_regular %>%
    sf::st_drop_geometry() %>%
    dplyr::select(station_id = point_id, station_x = x, station_y = y, depth = raster_value)

  # Get stations that DID detect for each fish-time
  detecting_stations <- station_detections %>%
    dplyr::select(path_id, time_bin, station_id) %>%
    dplyr::distinct()

  # Get coordinates of detecting stations for each fish-time combination
  detecting_station_coords <- detecting_stations %>%
    dplyr::left_join(all_stations, by = "station_id") %>%
    dplyr::select(path_id, time_bin, detecting_station_id = station_id,
                  detecting_x = station_x, detecting_y = station_y)

  # Filter non-detecting stations based on distance
  if (!is.null(max_distance_from_detecting) && max_distance_from_detecting > 0) {
    candidate_non_detections <- fish_time_data %>%
      dplyr::left_join(detecting_station_coords, by = c("path_id", "time_bin")) %>%
      dplyr::cross_join(all_stations %>%
                          dplyr::rename(candidate_station_id = station_id,
                                        candidate_x = station_x,
                                        candidate_y = station_y,
                                        candidate_depth = depth)) %>%
      dplyr::mutate(
        distance_between_stations = sqrt((candidate_x - detecting_x)^2 + (candidate_y - detecting_y)^2)
      ) %>%
      dplyr::filter(distance_between_stations <= max_distance_from_detecting) %>%
      dplyr::select(path_id, time_bin, station_id = candidate_station_id,
                    station_x = candidate_x, station_y = candidate_y, depth = candidate_depth,
                    step, time_seconds, x, y) %>%
      dplyr::distinct() %>%
      dplyr::mutate(
        distance_to_station = sqrt((x - station_x)^2 + (y - station_y)^2)
      )
  } else {
    candidate_non_detections <- fish_time_data %>%
      dplyr::cross_join(all_stations) %>%
      dplyr::mutate(
        distance_to_station = sqrt((x - station_x)^2 + (y - station_y)^2)
      )
  }

  # Identify non-detecting combinations
  non_detections <- candidate_non_detections %>%
    dplyr::anti_join(detecting_stations, by = c("path_id", "time_bin", "station_id")) %>%
    dplyr::mutate(
      detected = 0,
      detection_prob = NA
    ) %>%
    dplyr::select(path_id, step, time_seconds, time_bin, station_id,
                  x, y, distance_to_station, detection_prob, detected,
                  station_x, station_y, depth) %>%
    dplyr::arrange(path_id, time_bin, station_id)

  return(non_detections)
}

# Function to aggregate non-detections and generate cells for prediction
aggregate_non_detections <- function(non_detections,
                                     station_distances_df,
                                     fish_id_col = "path_id",
                                     time_col = "time_bin",
                                     station_col = "station_id") {

  # Validate inputs
  if (!fish_id_col %in% names(non_detections)) {
    stop(paste("Column", fish_id_col, "not found in non_detections"))
  }
  if (!time_col %in% names(non_detections)) {
    stop(paste("Column", time_col, "not found in non_detections"))
  }
  if (!station_col %in% names(non_detections)) {
    stop(paste("Column", station_col, "not found in non_detections"))
  }

  # Rename columns to standard names for processing
  non_detections_renamed <- non_detections %>%
    dplyr::rename(
      fish_id = dplyr::all_of(fish_id_col),
      time_period = dplyr::all_of(time_col),
      station_id = dplyr::all_of(station_col)
    )

  # Aggregate non-detections by fish, time, and station
  non_detection_summary <- non_detections_renamed %>%
    dplyr::group_by(fish_id, time_period, station_id) %>%
    dplyr::summarise(
      n_detections = 0,
      first_detection_time = NA,
      last_detection_time = NA,
      mean_detection_prob = mean(detection_prob, na.rm = TRUE),
      total_distance = mean(distance_to_station, na.rm = TRUE),
      .groups = 'drop'
    )

  # Merge with station_distances_df
  prediction_data <- non_detection_summary %>%
    dplyr::left_join(
      station_distances_df,
      by = c("station_id" = "station_no"),
      relationship = "many-to-many"
    ) %>%
    dplyr::filter(!is.na(cell_id)) %>%
    dplyr::select(
      fish_id, time_period, station_id, n_detections,
      first_detection_time, last_detection_time, mean_detection_prob,
      total_distance, cell_id, x, y, raster_value, cost_distance,
      straight_distance, tortuosity, DE_pred
    ) %>%
    dplyr::arrange(fish_id, time_period, station_id, cell_id)

  return(prediction_data)
}

# Function to aggregate probabilities
aggregate_probability <- function(df, detection_weight = 0.5, non_detection_weight = 0.5, normalize_method = "global") {

  # Validate weights sum to 1
  if (abs(detection_weight + non_detection_weight - 1) > 1e-10) {
    stop("Weights must sum to 1")
  }

  # Get unique x,y coordinates for each cell_id
  cell_coords <- df %>%
    dplyr::select(cell_id, x, y) %>%
    dplyr::distinct()

  # Separate detection and non-detection data
  detection_data <- df %>%
    dplyr::filter(type == "detection") %>%
    dplyr::group_by(fish_id, time_period, cell_id) %>%
    dplyr::summarise(
      weighted_mean_DE_normalized = stats::weighted.mean(DE_pred_normalized, n_detections),
      detections = dplyr::first(n_detections),
      .groups = "drop"
    )

  non_detection_data <- df %>%
    dplyr::filter(type == "non-detection") %>%
    dplyr::group_by(fish_id, time_period, cell_id) %>%
    dplyr::summarise(
      non_det_DE_normalized = mean(DE_pred_normalized),
      .groups = "drop"
    )

  # Join the two datasets
  combined_data <- dplyr::full_join(
    detection_data,
    non_detection_data,
    by = c("fish_id", "time_period", "cell_id")
  )

  # Apply normalization
  if (normalize_method == "global") {
    combined_data <- combined_data %>%
      dplyr::mutate(
        det_min = min(weighted_mean_DE_normalized, na.rm = TRUE),
        det_max = max(weighted_mean_DE_normalized, na.rm = TRUE),
        det_range = det_max - det_min,
        weighted_mean_DE_normalized_scaled = ifelse(det_range > 0,
                                                    (weighted_mean_DE_normalized - det_min) / det_range,
                                                    weighted_mean_DE_normalized),
        non_det_min = min(non_det_DE_normalized, na.rm = TRUE),
        non_det_max = max(non_det_DE_normalized, na.rm = TRUE),
        non_det_range = non_det_max - non_det_min,
        non_det_DE_normalized_scaled = ifelse(non_det_range > 0,
                                              (non_det_DE_normalized - non_det_min) / non_det_range,
                                              non_det_DE_normalized)
      ) %>%
      dplyr::select(-det_min, -det_max, -det_range, -non_det_min, -non_det_max, -non_det_range)
  } else {
    combined_data <- combined_data %>%
      dplyr::mutate(
        weighted_mean_DE_normalized_scaled = weighted_mean_DE_normalized,
        non_det_DE_normalized_scaled = non_det_DE_normalized
      )
  }

  # Calculate integrated probability
  result <- combined_data %>%
    dplyr::left_join(cell_coords, by = "cell_id") %>%
    dplyr::mutate(
      weighted_mean_DE_normalized_scaled = ifelse(is.na(weighted_mean_DE_normalized_scaled), 0, weighted_mean_DE_normalized_scaled),
      non_det_DE_normalized_scaled = ifelse(is.na(non_det_DE_normalized_scaled), 0, non_det_DE_normalized_scaled),
      integrated_prob = (weighted_mean_DE_normalized_scaled * detection_weight) +
        ((1 - non_det_DE_normalized_scaled) * non_detection_weight)
    ) %>%
    dplyr::select(
      fish_id, time_period, cell_id, x, y, detections,
      weighted_mean_DE_normalized, non_det_DE_normalized, integrated_prob
    )

  return(result)
}
