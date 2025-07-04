#' Calculate probability of detecting signal on 3 or more receivers
#'
#' Helper function that computes the probability of a signal being detected
#' by at least 3 receivers using binomial probability calculations.
#'
#' @param probs Numeric vector of detection probabilities for individual receivers.
#'
#' @return Numeric value representing the probability of detection on 3+ receivers.
#'   Returns 0 if fewer than 3 receivers are available.
#'
#' @details
#' Calculates P(X >= 3) = 1 - [P(X=0) + P(X=1) + P(X=2)] where X is the number
#' of receivers detecting the signal. Uses exact binomial calculations rather
#' than approximations.
#'
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

#' Calculate detection system probabilities and create visualization
#'
#' Computes cumulative detection probabilities and fine-scale positioning
#' probabilities (3+ receivers) across a spatial area using a detection
#' efficiency model. Creates heatmap visualizations of detection coverage.
#'
#' @param distance_frame Data frame containing distance calculations from
#'   \code{\link{calculate_station_distances}}. Must include columns: cell_id,
#'   x, y, raster_value, cost_distance.
#' @param receiver_frame An sf object containing receiver station locations
#'   with x, y coordinates for plotting.
#' @param model A fitted model object (e.g., from \code{\link{create_logistic_curve_depth}})
#'   that can predict detection efficiency. Must accept 'dist_m' and 'depth_m'
#'   as predictors.
#' @param output_type Character. Type of probabilities to calculate:
#'   \itemize{
#'     \item "cumulative" - Probability of detection by at least one receiver
#'     \item "3_plus" - Probability of detection by 3+ receivers (for positioning)
#'     \item "both" - Calculate both types (default)
#'   }
#' @param plots Logical. Whether to create and display visualization plots.
#'   Default is TRUE.
#'
#' @return If plots = TRUE, returns a list containing:
#'   \item{data}{Data frame with spatial coordinates and calculated probabilities}
#'   \item{plots}{List of individual ggplot objects}
#'   \item{combined_plot}{Combined plot using patchwork (if both types calculated)}
#'
#'   If plots = FALSE, returns only the data frame with calculated probabilities.
#'
#' @details
#' The function performs the following calculations:
#' \enumerate{
#'   \item Predicts detection efficiency for each cell-receiver pair using the model
#'   \item Calculates cumulative detection probability: 1 - prod(1 - individual_probs)
#'   \item Calculates 3+ receiver probability using exact binomial calculations
#'   \item Creates spatial heatmaps showing detection coverage
#' }
#'
#' The cumulative probability represents the likelihood of detecting a signal
#' anywhere in the system, while the 3+ receiver probability indicates areas
#' where fine-scale positioning is possible.
#'
#' The model predictions use 'cost_distance' as 'dist_m' and absolute
#' 'raster_value' as 'depth_m' to match expected model inputs.
#'
#' @examples
#' \dontrun{
#' # Generate receiver stations and calculate distances
#' stations <- generate_random_points(depth_raster, n_points = 6, seed = 123)
#' distances <- calculate_station_distances(depth_raster, stations, max_distance = 500)
#'
#' # Create detection efficiency model
#' de_model <- create_logistic_curve_depth(
#'   min_depth = 2, max_depth = 30,
#'   d50_min_depth = 50, d95_min_depth = 20,
#'   d50_max_depth = 150, d95_max_depth = 60,
#'   plot = FALSE
#' )
#'
#' # Calculate system detection probabilities
#' detection_system <- calculate_detection_system(
#'   distance_frame = distances,
#'   receiver_frame = stations,
#'   model = de_model$log_model,
#'   output_type = "both"
#' )
#'
#' # Extract results
#' detection_data <- detection_system$data
#' positioning_plot <- detection_system$plots$prob_3_plus
#'
#' # Calculate only positioning probabilities without plots
#' positioning_only <- calculate_detection_system(
#'   distance_frame = distances,
#'   receiver_frame = stations,
#'   model = de_model$log_model,
#'   output_type = "3_plus",
#'   plots = FALSE
#' )
#' }
#'
#' @seealso \code{\link{calculate_station_distances}}, \code{\link{create_logistic_curve_depth}}
#'
#' @export
calculate_detection_system <- function(distance_frame,
                                       receiver_frame,
                                       model,
                                       output_type = "both",  # "cumulative", "3_plus", or "both"
                                       plots = TRUE) {

  # Validate output_type
  if (!output_type %in% c("cumulative", "3_plus", "both")) {
    stop("output_type must be 'cumulative', '3_plus', or 'both'")
  }

  # Determine what to calculate
  calc_cumulative <- output_type %in% c("cumulative", "both")
  calc_3_plus <- output_type %in% c("3_plus", "both")

  # Predict detection probabilities using the specified model
  cat("Calculating detection probabilities...\n")
  distance_frame$DE_pred <- stats::predict(model,
                                           newdata = distance_frame %>%
                                             dplyr::rename(dist_m = cost_distance) %>%
                                             dplyr::mutate(depth_m = abs(raster_value)),
                                           type = "response")

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

  # Create plots if requested
  plot_list <- list()

  if (plots) {
    cat("Creating plots...\n")

    if (calc_cumulative) {
      p1 <- ggplot2::ggplot(system_summary, ggplot2::aes(x, y, fill = cumulative_prob)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_viridis_c(option = "magma", name = "Probability") +
        ggplot2::geom_point(data = receiver_frame, ggplot2::aes(x, y),
                            col = "green", fill = NA, pch = 21, inherit.aes = FALSE) +
        ggplot2::ggtitle("Single Detection Probability") +
        ggplot2::coord_sf() +
        ggplot2::theme_minimal()

      plot_list[["cumulative"]] <- p1
    }

    if (calc_3_plus) {
      p2 <- ggplot2::ggplot(system_summary, ggplot2::aes(x, y, fill = prob_3_plus)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_viridis_c(option = "magma", name = "Probability") +
        ggplot2::geom_point(data = receiver_frame, ggplot2::aes(x, y),
                            col = "green", fill = NA, pch = 21, inherit.aes = FALSE) +
        ggplot2::ggtitle("Fine Scale Positioning Probability (3+ Receivers)") +
        ggplot2::coord_sf() +
        ggplot2::theme_minimal()

      plot_list[["prob_3_plus"]] <- p2
    }

    # Combine plots
    if (length(plot_list) == 2) {
      combined_plot <- plot_list$cumulative / plot_list$prob_3_plus
    } else if (length(plot_list) == 1) {
      combined_plot <- plot_list[[1]]
    } else {
      combined_plot <- NULL
    }

    # Display the plot
    if (!is.null(combined_plot)) {
      print(combined_plot)
    }

    # Return results with plots
    results <- list(
      data = system_summary,
      plots = plot_list,
      combined_plot = combined_plot
    )

  } else {
    # Return results without plots
    results <- system_summary
  }

  return(results)
}
