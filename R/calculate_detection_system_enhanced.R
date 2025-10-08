#' Calculate detection system probabilities with deployment filtering (Enhanced)
#'
#' Enhanced version that handles character station IDs and filters by deployment periods.
#' Computes cumulative detection probabilities and fine-scale positioning
#' probabilities (3+ receivers) across a spatial area for a specific date,
#' considering only receivers that were deployed at that time.
#'
#' @param distance_frame Data frame containing distance calculations from
#'   \code{\link{calculate_station_distances}}. Must include columns: cell_id,
#'   x, y, raster_value, cost_distance, station_no (can be character or numeric).
#' @param receiver_frame An sf object or data frame containing receiver station information.
#'   Must include deployment information and coordinates.
#' @param model A fitted model object (e.g., from \code{\link{create_logistic_curve_depth}})
#'   that can predict detection efficiency. Must accept 'dist_m' and 'depth_m'
#'   as predictors.
#' @param target_date Date or POSIXct. The date for which to calculate detection probabilities.
#'   Only receivers deployed on this date will be included. If NULL, uses all receivers.
#' @param output_type Character. Type of probabilities to calculate:
#'   \itemize{
#'     \item "cumulative" - Probability of detection by at least one receiver
#'     \item "3_plus" - Probability of detection by 3+ receivers (for positioning)
#'     \item "both" - Calculate both types (default)
#'   }
#' @param plots Logical. Whether to create and display visualization plots.
#'   Default is TRUE.
#' @param include_barriers Logical. If TRUE, sets detection efficiency to zero
#'   for receiver-cell pairs where the path crosses barriers (land/obstacles).
#'   Requires the 'crosses_barrier' column in distance_frame (added by updated
#'   calculate_station_distances function). Default is FALSE.
#' @param station_col Character. Column name in receiver_frame containing station IDs.
#'   Default is "station_id". Must match the station identifiers in distance_frame$station_no.
#' @param deploy_col Character. Column name in receiver_frame containing deployment dates.
#'   Default is "deploy_datetime_UTC".
#' @param recover_col Character. Column name in receiver_frame containing recovery dates.
#'   Default is "recover_datetime_UTC".
#'
#' @return If plots = TRUE, returns a list containing:
#'   \item{data}{Data frame with spatial coordinates and calculated probabilities}
#'   \item{plots}{List of individual ggplot objects}
#'   \item{combined_plot}{Combined plot using patchwork (if both types calculated)}
#'   \item{active_stations}{Character vector of station IDs active on target_date}
#'   \item{n_active_stations}{Number of active stations}
#'
#'   If plots = FALSE, returns only the data frame with calculated probabilities.
#'
#' @details
#' This enhanced version:
#' \itemize{
#'   \item Supports both numeric and character station IDs
#'   \item Filters receivers by deployment periods for a specific date
#'   \item Handles flexible column naming for station IDs and deployment dates
#'   \item Provides deployment statistics in the output
#'   \item Maintains compatibility with existing workflows
#' }
#'
#' The function performs the following steps:
#' \enumerate{
#'   \item Identifies receivers active on the target date
#'   \item Filters distance_frame to include only active receivers
#'   \item Predicts detection efficiency for each cell-receiver pair
#'   \item (Optional) Masks detections where barriers are present (if include_barriers = TRUE)
#'   \item Calculates cumulative and/or 3+ receiver probabilities
#'   \item Creates spatial visualizations showing active receiver coverage
#' }
#'
#' When \code{include_barriers = TRUE}, detection efficiency is set to zero for
#' any receiver-cell pair where the straight-line path crosses land or other
#' barriers (NA cells in the raster). This prevents unrealistic detections through
#' obstacles and provides more accurate detection coverage predictions.
#'
#' @examples
#' \dontrun{
#' # With character station IDs and deployment filtering
#' target_date <- as.Date("2023-07-15")
#'
#' system_DE <- calculate_detection_system(
#'   distance_frame = station_distances,  # From calculate_station_distances
#'   receiver_frame = stoney_rx_deploy,   # With deployment information
#'   model = logistic_DE$log_model,
#'   target_date = target_date,
#'   output_type = "both",
#'   station_col = "station_id",
#'   deploy_col = "deploy_datetime_UTC",
#'   recover_col = "recover_datetime_UTC"
#' )
#'
#' # Check which stations were active
#' print(system_DE$active_stations)
#' print(paste("Active stations:", system_DE$n_active_stations))
#'
#' # Use all stations (no date filtering)
#' system_DE_all <- calculate_detection_system(
#'   distance_frame = station_distances,
#'   receiver_frame = stoney_rx_deploy,
#'   model = logistic_DE$log_model,
#'   target_date = NULL,  # No filtering
#'   output_type = "cumulative"
#' )
#' }
#'
#' @seealso \code{\link{calculate_station_distances}}, \code{\link{create_logistic_curve_depth}}
#'
#' @export
calculate_detection_system <- function(distance_frame,
                                       receiver_frame,
                                       model,
                                       target_date = NULL,
                                       output_type = "both",
                                       plots = TRUE,
                                       include_barriers = FALSE,
                                       station_col = "station_id",
                                       deploy_col = "deploy_datetime_UTC",
                                       recover_col = "recover_datetime_UTC") {

  # Validate inputs
  if (!output_type %in% c("cumulative", "3_plus", "both")) {
    stop("output_type must be 'cumulative', '3_plus', or 'both'")
  }

  # Check required columns in distance_frame
  required_distance_cols <- c("cell_id", "x", "y", "raster_value", "cost_distance", "station_no")
  missing_distance_cols <- setdiff(required_distance_cols, names(distance_frame))
  if (length(missing_distance_cols) > 0) {
    stop("Missing columns in distance_frame: ", paste(missing_distance_cols, collapse = ", "))
  }

  # Check required columns in receiver_frame
  if (!station_col %in% names(receiver_frame)) {
    stop("Column '", station_col, "' not found in receiver_frame. ",
         "Available columns: ", paste(names(receiver_frame), collapse = ", "))
  }

  # Check for crosses_barrier column if barrier masking requested
  if (include_barriers) {
    if (!"crosses_barrier" %in% names(distance_frame)) {
      stop("include_barriers = TRUE requires 'crosses_barrier' column in distance_frame.\n",
           "This column is added by the updated calculate_station_distances() function.\n",
           "Please recalculate distances with the barrier-aware version.")
    }
    cat("Barrier masking enabled - will set DE = 0 where paths cross barriers\n")
  }

  # Determine what to calculate
  calc_cumulative <- output_type %in% c("cumulative", "both")
  calc_3_plus <- output_type %in% c("3_plus", "both")

  # Handle deployment filtering if target_date is provided
  active_stations <- NULL
  n_active_stations <- 0

  if (!is.null(target_date)) {
    # Check for deployment columns
    if (!deploy_col %in% names(receiver_frame) || !recover_col %in% names(receiver_frame)) {
      stop("Deployment columns '", deploy_col, "' and '", recover_col, 
           "' not found in receiver_frame for date filtering. ",
           "Set target_date = NULL to use all stations.")
    }

    # Convert target_date to Date if needed
    if (inherits(target_date, "POSIXct")) {
      target_date <- as.Date(target_date)
    } else if (!inherits(target_date, "Date")) {
      target_date <- as.Date(target_date)
    }

    cat("Filtering receivers active on", as.character(target_date), "...\n")

    # Convert deployment dates to Date for comparison
    deploy_dates <- as.Date(receiver_frame[[deploy_col]])
    recover_dates <- as.Date(receiver_frame[[recover_col]])

    # Find active stations
    active_mask <- deploy_dates <= target_date & recover_dates >= target_date
    active_stations <- receiver_frame[[station_col]][active_mask & !is.na(active_mask)]
    n_active_stations <- length(active_stations)

    cat("Found", n_active_stations, "active stations on", as.character(target_date), "\n")

    if (n_active_stations == 0) {
      stop("No stations were active on ", as.character(target_date), 
           ". Check deployment dates or choose a different date.")
    }

    # Filter distance_frame to include only active stations
    original_rows <- nrow(distance_frame)
    distance_frame <- distance_frame %>%
      dplyr::filter(station_no %in% active_stations)

    filtered_rows <- nrow(distance_frame)
    cat("Filtered distance_frame from", original_rows, "to", filtered_rows, "rows\n")

    # Filter receiver_frame for plotting
    receiver_frame_filtered <- receiver_frame %>%
      dplyr::filter(!!rlang::sym(station_col) %in% active_stations)

  } else {
    cat("Using all", length(unique(distance_frame$station_no)), "stations (no date filtering)\n")
    active_stations <- unique(distance_frame$station_no)
    n_active_stations <- length(active_stations)
    receiver_frame_filtered <- receiver_frame
  }

  # Predict detection probabilities using the specified model
  cat("Calculating detection probabilities...\n")

  # Check if DE_pred already exists
  if (!"DE_pred" %in% names(distance_frame)) {
    distance_frame$DE_pred <- stats::predict(model,
                                             newdata = distance_frame %>%
                                               dplyr::rename(dist_m = cost_distance) %>%
                                               dplyr::mutate(depth_m = abs(raster_value)),
                                             type = "response")
  } else {
    cat("Using existing DE_pred column\n")
  }

  # Apply barrier masking if requested
  if (include_barriers) {
    barrier_mask <- distance_frame$crosses_barrier
    n_total <- nrow(distance_frame)
    n_masked <- sum(barrier_mask, na.rm = TRUE)

    # Set DE to zero where barriers detected
    distance_frame$DE_pred[barrier_mask & !is.na(barrier_mask)] <- 0

    cat(sprintf("  Masked %d of %d station-cell pairs (%.1f%%) due to barriers\n",
                n_masked, n_total, 100 * n_masked / n_total))
  }

  # Calculate system detection probabilities
  cat("Calculating system probabilities...\n")

  # Base summary
  system_summary <- distance_frame %>%
    dplyr::group_by(cell_id) %>%
    dplyr::summarise(
      x = dplyr::first(x),
      y = dplyr::first(y),
      raster_value = dplyr::first(raster_value),
      n_stations = dplyr::n(),
      .groups = 'drop'
    )

  # Add cumulative probability if requested
  if (calc_cumulative) {
    cumulative_data <- distance_frame %>%
      dplyr::group_by(cell_id) %>%
      dplyr::summarise(
        cumulative_prob = 1 - prod(1 - DE_pred, na.rm = TRUE),
        .groups = 'drop'
      )

    system_summary <- system_summary %>%
      dplyr::left_join(cumulative_data, by = "cell_id")
  }

  # Add 3+ probability if requested
  if (calc_3_plus) {
    prob_3_plus_data <- distance_frame %>%
      dplyr::group_by(cell_id) %>%
      dplyr::summarise(
        prob_3_plus = calculate_prob_3_plus(DE_pred),
        .groups = 'drop'
      )

    system_summary <- system_summary %>%
      dplyr::left_join(prob_3_plus_data, by = "cell_id")
  }

  # Extract coordinates for plotting from receiver_frame_filtered
  if ("sf" %in% class(receiver_frame_filtered)) {
    # Extract coordinates from sf object
    coords_matrix <- sf::st_coordinates(receiver_frame_filtered)
    plot_coords <- data.frame(
      x = coords_matrix[, 1],
      y = coords_matrix[, 2]
    )
  } else {
    # Use existing coordinate columns
    coord_cols <- c("x", "y")
    if (!all(coord_cols %in% names(receiver_frame_filtered))) {
      # Try alternative names
      if (all(c("station_x", "station_y") %in% names(receiver_frame_filtered))) {
        coord_cols <- c("station_x", "station_y")
      } else {
        warning("Could not find coordinate columns for plotting. Using x and y from distance_frame.")
        plot_coords <- distance_frame %>%
          dplyr::group_by(station_no) %>%
          dplyr::summarise(x = mean(x), y = mean(y), .groups = 'drop')
      }
    }
    
    if (exists("coord_cols")) {
      plot_coords <- receiver_frame_filtered[, coord_cols]
      names(plot_coords) <- c("x", "y")
    }
  }

  # Create plots if requested
  plot_list <- list()

  if (plots) {
    cat("Creating plots...\n")

    # Create title suffix based on filtering
    title_suffix <- if (!is.null(target_date)) {
      paste0(" (", n_active_stations, " stations active on ", as.character(target_date), ")")
    } else {
      paste0(" (", n_active_stations, " stations)")
    }

    if (calc_cumulative) {
      p1 <- ggplot2::ggplot(system_summary, ggplot2::aes(x, y, fill = cumulative_prob)) +
        ggplot2::geom_raster() +
        ggplot2::scale_fill_viridis_c(option = "magma", name = "Probability") +
        ggplot2::geom_point(data = plot_coords, ggplot2::aes(x, y),
                            col = "green", fill = NA, pch = 21, size = 2, inherit.aes = FALSE) +
        ggplot2::ggtitle(paste0("Single Detection Probability", title_suffix)) +
        ggplot2::coord_sf() +
        ggplot2::theme_minimal()

      plot_list[["cumulative"]] <- p1
    }

    if (calc_3_plus) {
      p2 <- ggplot2::ggplot(system_summary, ggplot2::aes(x, y, fill = prob_3_plus)) +
        ggplot2::geom_raster() +
        ggplot2::scale_fill_viridis_c(option = "magma", name = "Probability") +
        ggplot2::geom_point(data = plot_coords, ggplot2::aes(x, y),
                            col = "green", fill = NA, pch = 21, size = 2, inherit.aes = FALSE) +
        ggplot2::ggtitle(paste0("Fine Scale Positioning Probability (3+ Receivers)", title_suffix)) +
        ggplot2::coord_sf() +
        ggplot2::theme_minimal()

      plot_list[["prob_3_plus"]] <- p2
    }

    # Combine plots if both calculated
    if (length(plot_list) == 2) {
      if (requireNamespace("patchwork", quietly = TRUE)) {
        combined_plot <- plot_list$cumulative / plot_list$prob_3_plus
      } else {
        # Fallback if patchwork not available
        combined_plot <- plot_list$cumulative
        warning("Package 'patchwork' not available. Showing only cumulative plot.")
      }
    } else if (length(plot_list) == 1) {
      combined_plot <- plot_list[[1]]
    } else {
      combined_plot <- NULL
    }

    # Display the plot
    if (!is.null(combined_plot)) {
      print(combined_plot)
    }

    # Return comprehensive results with plots
    results <- list(
      data = system_summary,
      plots = plot_list,
      combined_plot = combined_plot,
      active_stations = active_stations,
      n_active_stations = n_active_stations,
      target_date = target_date
    )

  } else {
    # Return results without plots but with station info
    results <- list(
      data = system_summary,
      active_stations = active_stations,
      n_active_stations = n_active_stations,
      target_date = target_date
    )
  }

  # Summary statistics
  cat("\n=== Detection System Summary ===\n")
  cat("Target date:", if(is.null(target_date)) "All dates" else as.character(target_date), "\n")
  cat("Active stations:", n_active_stations, "\n")
  if (calc_cumulative) {
    cat("Mean cumulative detection probability:", 
        round(mean(system_summary$cumulative_prob, na.rm = TRUE), 3), "\n")
  }
  if (calc_3_plus) {
    cat("Mean 3+ receiver probability:", 
        round(mean(system_summary$prob_3_plus, na.rm = TRUE), 3), "\n")
    cat("Cells with >0.5 positioning probability:", 
        sum(system_summary$prob_3_plus > 0.5, na.rm = TRUE), 
        "/", nrow(system_summary), "\n")
  }

  return(results)
}

#' Helper function for 3+ receiver probability calculation
#' @keywords internal
calculate_prob_3_plus <- function(probs) {
  n <- length(probs)
  if (n < 3) return(0)  # Can't have 3+ detections with fewer than 3 receivers

  # Calculate using binomial probability
  # P(X >= 3) = 1 - P(X <= 2) = 1 - [P(X=0) + P(X=1) + P(X=2)]

  prob_0 <- prod(1 - probs, na.rm = TRUE)  # Probability of 0 detections

  prob_1 <- 0  # Probability of exactly 1 detection
  for (i in 1:n) {
    if (!is.na(probs[i])) {
      prob_1 <- prob_1 + probs[i] * prod(1 - probs[-i], na.rm = TRUE)
    }
  }

  prob_2 <- 0  # Probability of exactly 2 detections
  if (n >= 2) {
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        if (!is.na(probs[i]) && !is.na(probs[j])) {
          remaining_probs <- probs[-c(i,j)]
          prob_2 <- prob_2 + probs[i] * probs[j] * prod(1 - remaining_probs, na.rm = TRUE)
        }
      }
    }
  }

  prob_3_plus <- 1 - (prob_0 + prob_1 + prob_2)
  return(max(0, prob_3_plus))  # Ensure non-negative
}