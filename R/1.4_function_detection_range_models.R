#'generate detection range models

#' Create logistic curve for detection efficiency with depth extrapolation
#'
#' Generates a logistic detection efficiency model that varies with both distance
#' and depth. The function interpolates d50 and d95 parameters linearly between
#' minimum and maximum depths, creating a 3D detection efficiency surface.
#'
#' @param min_depth Numeric. Minimum depth in meters for the model.
#' @param max_depth Numeric. Maximum depth in meters for the model.
#' @param d50_min_depth Numeric. Distance (m) at which detection efficiency = 50%
#'   at minimum depth.
#' @param d95_min_depth Numeric. Distance (m) at which detection efficiency = 95%
#'   at minimum depth.
#' @param d50_max_depth Numeric. Distance (m) at which detection efficiency = 50%
#'   at maximum depth.
#' @param d95_max_depth Numeric. Distance (m) at which detection efficiency = 95%
#'   at maximum depth.
#' @param dist_step Numeric. Distance step size for grid predictions. Default is 5.
#' @param depth_step Numeric. Depth step size for grid predictions. Default is 1.
#' @param color_option Character. Viridis color palette option for plots.
#'   Default is "magma".
#' @param plot Logical. Whether to display visualization plots. Default is TRUE.
#' @param return_model Logical. Whether to fit and return a logistic regression model.
#'   Default is TRUE.
#' @param return_object Logical. Whether to return the complete function object.
#'   Default is TRUE.
#'
#' @return If return_object = TRUE, returns a list containing:
#'   \item{func}{The logistic function for calculating detection efficiency}
#'   \item{predict_grid}{Function for generating prediction grids}
#'   \item{predictions}{Data frame with distance, depth, and DE predictions}
#'   \item{params}{List of input parameters}
#'   \item{log_model}{Fitted logistic regression model (if return_model = TRUE)}
#'   \item{line_plot}{ggplot2 line plot object (if plot = TRUE)}
#'   \item{heatmap_plot}{ggplot2 heatmap object (if plot = TRUE)}
#'   \item{combined_plot}{Combined plot using patchwork (if plot = TRUE)}
#'
#'   If return_object = FALSE, returns invisibly or displays plots only.
#'
#' @details
#' The function creates a logistic detection efficiency model where:
#' \itemize{
#'   \item d50 and d95 parameters vary linearly with depth
#'   \item Detection efficiency = 1 / (1 + exp(a + b * distance))
#'   \item b = 2.944 / (d95 - d50)
#'   \item a = -b * d50
#' }
#'
#' The model assumes that detection efficiency decreases with distance and that
#' the detection range changes with depth in a predictable manner.
#'
#' @examples
#' \dontrun{
#' # Create a basic detection efficiency model
#' de_model <- create_logistic_curve_depth(
#'   min_depth = 2,
#'   max_depth = 30,
#'   d50_min_depth = 50,
#'   d95_min_depth = 20,
#'   d50_max_depth = 150,
#'   d95_max_depth = 60
#' )
#'
#' # Use the model to predict detection efficiency
#' efficiency <- de_model$func(dist_m = 100, depth_m = 15)
#'
#' # Generate custom predictions
#' custom_grid <- de_model$predict_grid(
#'   dist_range = seq(0, 200, by = 10),
#'   depth_range = seq(5, 25, by = 2)
#' )
#'
#' # Create model without plots
#' de_model_quiet <- create_logistic_curve_depth(
#'   min_depth = 2, max_depth = 30,
#'   d50_min_depth = 50, d95_min_depth = 20,
#'   d50_max_depth = 150, d95_max_depth = 60,
#'   plot = FALSE
#' )
#' }
#'
#' @seealso \code{\link{generate_random_points}}, \code{\link{generate_spaced_points}}
#'
#' @export
create_logistic_curve_depth <- function(min_depth, max_depth,
                                        d50_min_depth, d95_min_depth,
                                        d50_max_depth, d95_max_depth,
                                        dist_step = 5,             # Distance step size for predictions
                                        depth_step = 1,            # Depth step size for predictions
                                        color_option = "magma",    # Viridis color option
                                        plot = TRUE,               # Display plots by default
                                        return_model = TRUE,       # Option to return logistic model fit
                                        return_object = TRUE) {   # Whether to return the function object

  # Function to interpolate d50 and d95 based on depth
  interpolate_d50_d95 <- function(depth) {
    # Ensure depth is within bounds
    depth <- pmax(min_depth, pmin(max_depth, depth))

    # Calculate proportion of the way from min to max depth
    depth_prop <- (depth - min_depth) / (max_depth - min_depth)

    # Linearly interpolate d50 and d95
    d50 <- d50_min_depth + depth_prop * (d50_max_depth - d50_min_depth)
    d95 <- d95_min_depth + depth_prop * (d95_max_depth - d95_min_depth)

    return(list(d50 = d50, d95 = d95))
  }

  # Main logistic function
  logistic_function <- function(dist_m, depth_m) {
    # Get d50 and d95 for the given depth
    params <- interpolate_d50_d95(depth_m)
    d50 <- params$d50
    d95 <- params$d95

    # Calculate logistic parameters
    b <- 2.944 / (d95 - d50)
    a <- -b * d50

    # Calculate DE
    DE <- 1 / (1 + exp(a + b * dist_m))
    return(DE)
  }

  # Grid prediction function
  predict_grid <- function(dist_range = NULL, depth_range = NULL) {
    # Set default ranges if not provided
    if (is.null(dist_range)) {
      dist_range <- seq(0, d95_max_depth, by = dist_step)
    }
    if (is.null(depth_range)) {
      depth_range <- seq(min_depth, max_depth, by = depth_step)
    }

    grid <- expand.grid(dist_m = dist_range, depth_m = depth_range)
    grid$DE <- mapply(logistic_function, grid$dist_m, grid$depth_m)
    return(grid)
  }

  # Generate predictions for plotting/modeling
  dist_range <- seq(0, d95_max_depth, by = dist_step)
  depth_range <- seq(min_depth, max_depth, by = depth_step)
  predictions <- predict_grid(dist_range, depth_range)

  message("Created detection efficiency model with the following parameters:")
  message(sprintf("  Depth range: %.1f to %.1f m", min_depth, max_depth))
  message(sprintf("  d50 at min depth (%.1f m): %.1f m", min_depth, d50_min_depth))
  message(sprintf("  d95 at min depth (%.1f m): %.1f m", min_depth, d95_min_depth))
  message(sprintf("  d50 at max depth (%.1f m): %.1f m", max_depth, d50_max_depth))
  message(sprintf("  d95 at max depth (%.1f m): %.1f m", max_depth, d95_max_depth))

  # Generate logistic model if requested
  log_model <- NULL
  if (return_model) {
    message("\nFitting logistic regression model...")
    # Fit logistic regression model to the generated data
    log_model <- stats::glm(DE ~ dist_m * depth_m,
                            family = stats::binomial(link = "logit"),
                            data = predictions)
  }

  # Generate plots if requested
  if (plot) {
    # Check if required packages are available
    required_packages <- c("ggplot2", "patchwork", "scales")
    missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

    if (length(missing_packages) > 0) {
      message("The following packages are required for plotting but not installed: ",
              paste(missing_packages, collapse = ", "),
              ". Please install them.")
      if (return_object) {
        message("\nReturning model object without plots.")
        return(list(
          func = logistic_function,
          predict_grid = predict_grid,
          predictions = predictions,
          log_model = log_model,
          params = list(
            min_depth = min_depth,
            max_depth = max_depth,
            d50_min_depth = d50_min_depth,
            d95_min_depth = d95_min_depth,
            d50_max_depth = d50_max_depth,
            d95_max_depth = d95_max_depth
          )
        ))
      } else {
        return(invisible(NULL))
      }
    }

    message("\nGenerating visualization plots...")

    # Calculate points for d50 and d95 lines for annotation
    d50_points <- data.frame()
    d95_points <- data.frame()

    for (d in depth_range) {
      params <- interpolate_d50_d95(d)
      d50_points <- rbind(d50_points, data.frame(dist_m = params$d50, depth_m = d))
      d95_points <- rbind(d95_points, data.frame(dist_m = params$d95, depth_m = d))
    }

    # Select annotation points
    d50_anno <- d50_points[ceiling(nrow(d50_points)/2), , drop = FALSE]
    d95_anno <- d95_points[ceiling(nrow(d95_points)/2), , drop = FALSE]

    # Add DE value to annotation data frames
    d50_anno$DE <- 0.5
    d95_anno$DE <- 0.05

    # Create plots - suppress warnings for the entire plotting section
    suppressWarnings({
      # Line plot
      line_plot <- ggplot2::ggplot(predictions, ggplot2::aes(x = dist_m, y = DE, group = depth_m, color = factor(depth_m))) +
        ggplot2::geom_line(size = 1) +
        ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.5) +
        ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", alpha = 0.5) +
        ggplot2::labs(title = "Detection Efficiency Curves by Depth",
                      x = "Distance (m)",
                      y = "Detection Efficiency",
                      color = "Depth (m)") +
        ggplot2::theme_minimal() +
        ggplot2::scale_y_continuous(labels = scales::percent) +
        ggplot2::scale_color_viridis_d(direction = -1)

      # Heatmap with contour lines
      heatmap_plot <- ggplot2::ggplot(predictions, ggplot2::aes(x = dist_m, y = depth_m, fill = DE)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_viridis_c(option = color_option, labels = scales::percent) +
        ggplot2::labs(title = "Detection Efficiency Heatmap",
                      x = "Distance (m)",
                      y = "Depth (m)",
                      fill = "DE") +
        ggplot2::scale_y_reverse() +
        ggplot2::theme_minimal() +
        # Add contour lines at DE = 0.5 and DE = 0.05
        ggplot2::geom_contour(ggplot2::aes(z = DE), breaks = c(0.05, 0.5), color = "white", linewidth = 0.8) +
        # Add annotations for the DE = 0.5 contour line
        ggplot2::geom_text(data = d50_anno,
                           ggplot2::aes(x = dist_m + 15, y = depth_m),
                           label = "DE = 50%",
                           color = "white",
                           fontface = "bold",
                           size = 3.5,
                           hjust = 0,
                           inherit.aes = FALSE) +
        # Add annotations for the DE = 0.05 contour line
        ggplot2::geom_text(data = d95_anno,
                           ggplot2::aes(x = dist_m + 15, y = depth_m),
                           label = "DE = 5%",
                           color = "white",
                           fontface = "bold",
                           size = 3.5,
                           hjust = 0,
                           inherit.aes = FALSE)

      # Combine plots
      combined_plot <- line_plot / heatmap_plot

      # Display the combined plot
      message("Displaying detection efficiency visualization...")
      print(combined_plot)

      # Add model summary if requested
      if (return_model) {
        message("\nLogistic Model Summary:")
        print(summary(log_model))
      }

      # Store plots in function environment for possible return
      assign("line_plot", line_plot, environment())
      assign("heatmap_plot", heatmap_plot, environment())
      assign("combined_plot", combined_plot, environment())
    })
  }

  # Return results if requested
  if (return_object) {
    message("\nReturning complete model object with functions and data.")
    result <- list(
      func = logistic_function,
      predict_grid = predict_grid,
      predictions = predictions,
      params = list(
        min_depth = min_depth,
        max_depth = max_depth,
        d50_min_depth = d50_min_depth,
        d95_min_depth = d95_min_depth,
        d50_max_depth = d50_max_depth,
        d95_max_depth = d95_max_depth
      )
    )

    if (plot) {
      message("  - Including visualization plots")
      result$line_plot <- line_plot
      result$heatmap_plot <- heatmap_plot
      result$combined_plot <- combined_plot
    }

    if (return_model) {
      message("  - Including fitted logistic regression model")
      result$log_model <- log_model
    }

    return(result)
  } else {
    if (!plot) {
      message("\nNo output displayed or returned. Set plot=TRUE to view visualizations or return_object=TRUE to return the model object.")
    } else {
      message("\nVisualization displayed. Set return_object=TRUE if you want to save the model object for further use.")
    }
    # Return invisibly if no object requested
    return(invisible(NULL))
  }
}
