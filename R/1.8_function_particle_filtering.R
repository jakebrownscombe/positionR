#' Particle Filter Fish Positioning
#'
#' @description
#' Estimates fish positions using particle filtering with correlated random walk
#' movement models and detection efficiency models. Provides biologically realistic
#' positioning that maintains temporal continuity and accounts for detection uncertainty.
#'
#' @param detection_data A data frame containing detection events with the following columns:
#'   \describe{
#'     \item{fish_id or path_id}{Unique identifier for each fish}
#'     \item{time_bin or time_seconds}{Time of detection (numeric, in seconds or time bins)}
#'     \item{station_id}{Identifier for the detecting station/receiver}
#'     \item{x}{X coordinate of the fish at detection time}
#'     \item{y}{Y coordinate of the fish at detection time}
#'   }
#' @param station_coords A data frame or sf object containing receiver station information:
#'   \describe{
#'     \item{station_id or point_id}{Unique identifier for each station}
#'     \item{station_x, station_y or x, y}{X and Y coordinates of the station}
#'     \item{raster_value or depth_m}{Depth at station location (if raster_value is negative, will be converted to positive)}
#'   }
#' @param de_model A detection efficiency model object (typically from logistic regression)
#'   that can predict detection probability based on distance and depth. Should accept a data frame
#'   with columns 'dist_m' and 'depth_m'.
#' @param movement_params A named list containing correlated random walk parameters:
#'   \describe{
#'     \item{step_length_mean}{Mean step length in meters (default: 50)}
#'     \item{step_length_sd}{Standard deviation of step length in meters (default: 20)}
#'     \item{turning_angle_mean}{Mean turning angle in degrees (default: 0)}
#'     \item{turning_angle_sd}{Standard deviation of turning angle in degrees (default: 45)}
#'   }
#' @param boundary_raster Optional. A RasterLayer object defining valid movement areas.
#'   If provided, particles will be constrained to non-NA cells using boundary
#'   reflection (same logic as simulate_fish_tracks). Fish cannot move onto land
#'   or outside the raster extent.
#' @param n_particles Integer. Number of particles to use in the filter (default: 1000).
#'   More particles provide better accuracy but increase computational cost.
#' @param resample_threshold Numeric. Effective sample size threshold for triggering
#'   resampling (default: n_particles/2). Lower values resample more frequently.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing of particles.
#'   If NULL (default), automatically detects available cores and uses a conservative
#'   number (max 4, leaving 1 core free). Set to 1 to disable parallel processing.
#' @param parallel Logical. Whether to enable parallel processing of particles
#'   (default: TRUE). Set to FALSE for debugging or compatibility issues.
#' @param verbose Logical. Whether to print progress messages and status updates
#'   (default: TRUE).
#'
#' @return A list containing three elements:
#'   \describe{
#'     \item{position_estimates}{Data frame with estimated positions and uncertainty metrics:
#'       \itemize{
#'         \item fish_id: Fish identifier
#'         \item time_step: Time point
#'         \item x_est, y_est: Weighted mean estimated position
#'         \item x_var, y_var: Position variance estimates
#'         \item position_uncertainty: Combined position uncertainty (meters)
#'         \item x_ci_lower, x_ci_upper, y_ci_lower, y_ci_upper: 95% confidence intervals
#'         \item n_particles_effective: Effective sample size
#'         \item max_weight: Maximum particle weight
#'       }}
#'     \item{all_particles}{Data frame containing all particle trajectories with columns:
#'       particle_id, fish_id, x, y, heading, weight, time_step}
#'     \item{summary_stats}{List with overall performance metrics and parameters}
#'   }
#'
#' @details
#' The particle filter implements a sequential Monte Carlo method for fish positioning
#' that combines:
#' \itemize{
#'   \item \strong{Movement Model}: Correlated random walk (CRW) with biologically
#'     realistic parameters for step lengths and turning angles
#'   \item \strong{Detection Model}: Uses provided detection efficiency model with
#'     station-specific depths for realistic detection probabilities
#'   \item \strong{Uncertainty Quantification}: Provides position estimates with
#'     confidence intervals and uncertainty metrics
#'   \item \strong{Boundary Constraints}: Optional raster-based movement boundaries
#'     with reflection at edges and invalid cells
#'   \item \strong{Parallel Processing}: Automatic parallelization of particle
#'     likelihood calculations for improved performance
#' }
#'
#' The algorithm consists of three main steps repeated at each time point:
#' \enumerate{
#'   \item \strong{Predict}: Move particles forward using CRW movement model
#'   \item \strong{Update}: Weight particles based on detection likelihood
#'   \item \strong{Resample}: Resample particles when effective sample size drops
#'     below threshold
#' }
#'
#' @section Movement Model:
#' The correlated random walk model simulates realistic fish movement by:
#' \itemize{
#'   \item Drawing step lengths from a truncated normal distribution
#'   \item Drawing turning angles from a von Mises distribution (if circular package
#'     available) or normal distribution
#'   \item Scaling movement by time intervals between detections
#'   \item Maintaining heading persistence between time steps
#'   \item Applying boundary constraints if raster provided (particles bounce off
#'     raster edges and invalid/NA cells)
#' }
#'
#' @section Detection Model:
#' Detection probability is calculated for each particle at each receiver station
#' based on distance and station-specific depth:
#' \itemize{
#'   \item Uses actual depths from station_coords (raster_value or depth_m)
#'   \item Detection probability decreases with distance from receiver
#'   \item Stations that detected the fish contribute positive likelihood
#'   \item Stations that didn't detect contribute (1 - detection_probability)
#'   \item Fallback to 10m depth if no depth information available
#' }
#'
#' @section Parallel Processing:
#' For better performance with large particle counts:
#' \itemize{
#'   \item Parallel processing is automatically enabled on multi-core systems
#'   \item Cross-platform compatibility (Windows, Mac, Linux)
#'   \item Graceful fallback to sequential processing if parallel fails
#'   \item Most effective with 1000+ particles
#'   \item Conservative core detection (max 4 cores) for package compatibility
#' }
#'
#' @section Boundary Handling:
#' When boundary_raster is provided, particle movement is constrained using:
#' \itemize{
#'   \item Raster extent boundaries: particles bounce off edges
#'   \item Invalid cell boundaries: particles bounce off NA cells (land areas)
#'   \item Direction reversal: same "bounce" logic as simulate_fish_tracks
#'   \item Realistic habitat constraints: keeps fish in suitable areas
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with automatic parallel processing
#' movement_params <- list(
#'   step_length_mean = 50,
#'   step_length_sd = 20,
#'   turning_angle_mean = 0,
#'   turning_angle_sd = 45
#' )
#'
#' pf_results <- particle_filter_positioning(
#'   detection_data = fish_detections,
#'   station_coords = receiver_stations,
#'   de_model = logistic_DE$log_model,
#'   movement_params = movement_params,
#'   n_particles = 1000
#' )
#'
#' # High-performance usage with more cores
#' pf_results_fast <- particle_filter_positioning(
#'   detection_data = fish_detections,
#'   station_coords = receiver_stations,
#'   de_model = logistic_DE$log_model,
#'   movement_params = movement_params,
#'   n_particles = 5000,    # More particles
#'   n_cores = 8           # More cores for high-end systems
#' )
#'
#' # Use boundary constraints to keep fish in water
#' pf_results_bounded <- particle_filter_positioning(
#'   detection_data = fish_detections,
#'   station_coords = receiver_stations,
#'   de_model = logistic_DE$log_model,
#'   movement_params = movement_params,
#'   boundary_raster = depth_raster,  # Keep fish within raster bounds
#'   n_particles = 1000
#' )
#'
#' # Disable parallel processing for debugging
#' pf_results_debug <- particle_filter_positioning(
#'   detection_data = fish_detections,
#'   station_coords = receiver_stations,
#'   de_model = logistic_DE$log_model,
#'   movement_params = movement_params,
#'   parallel = FALSE
#' )
#'
#' # Extract results
#' positions <- pf_results$position_estimates
#' summary(pf_results$summary_stats)
#' }
#'
#' @references
#' \itemize{
#'   \item Doucet, A., & Johansen, A. M. (2009). A tutorial on particle filtering and
#'     smoothing: Fifteen years later. Handbook of nonlinear filtering, 12(656-704), 3.
#'   \item Codling, E. A., Plank, M. J., & Benhamou, S. (2008). Random walk models in
#'     biology. Journal of the Royal Society Interface, 5(25), 813-834.
#' }
#'
#' @seealso
#' \code{\link{plot_particle_filter_results}} for visualization,
#' \code{\link{simulate_fish_tracks}} for generating synthetic movement data
#'
#' @export
particle_filter_positioning <- function(detection_data,
                                        station_coords,
                                        de_model,
                                        movement_params,
                                        boundary_raster = NULL,
                                        n_particles = 1000,
                                        resample_threshold = NULL,
                                        n_cores = NULL,
                                        parallel = TRUE,
                                        verbose = TRUE) {

  # Load required packages
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' needed for this function to work.")
  }
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' needed for this function to work.")
  }

  if (verbose) cat("=== PARTICLE FILTER POSITIONING ===\n")

  # Set up parallel processing
  if (is.null(n_cores)) {
    available_cores <- parallel::detectCores(logical = FALSE)
    if (is.na(available_cores) || available_cores < 2) {
      n_cores <- 1
    } else {
      n_cores <- max(1, min(4, available_cores - 1))  # Conservative: max 4 cores, leave 1 free
    }
  }

  # Determine if parallel processing should be used
  use_parallel <- parallel && n_cores > 1

  if (verbose && use_parallel) {
    cat("Using", n_cores, "cores for parallel particle processing\n")
  } else if (verbose) {
    cat("Using sequential processing\n")
  }

  # Set default resample threshold
  if (is.null(resample_threshold)) {
    resample_threshold <- n_particles / 2
  }

  # Validate inputs and standardize column names
  detection_data_std <- detection_data

  # Handle different column naming conventions
  if ("path_id" %in% names(detection_data) && !"fish_id" %in% names(detection_data)) {
    detection_data_std$fish_id <- detection_data_std$path_id
  }
  if ("time_seconds" %in% names(detection_data) && !"time_bin" %in% names(detection_data)) {
    detection_data_std$time_bin <- detection_data_std$time_seconds
  }

  required_detection_cols <- c("fish_id", "time_bin", "station_id", "x", "y")
  missing_detection_cols <- setdiff(required_detection_cols, names(detection_data_std))
  if (length(missing_detection_cols) > 0) {
    stop("Missing columns in detection_data: ", paste(missing_detection_cols, collapse = ", "),
         "\nNote: path_id can substitute for fish_id, time_seconds can substitute for time_bin")
  }

  # Standardize station coordinates column names
  station_coords_std <- station_coords

  # Handle sf objects (extract coordinates)
  if ("sf" %in% class(station_coords)) {
    coords_matrix <- sf::st_coordinates(station_coords)
    station_coords_std <- station_coords %>%
      sf::st_drop_geometry() %>%
      dplyr::mutate(
        station_x = coords_matrix[,1],
        station_y = coords_matrix[,2]
      )
  }

  # Handle different column naming conventions
  if ("point_id" %in% names(station_coords_std) && !"station_id" %in% names(station_coords_std)) {
    station_coords_std$station_id <- station_coords_std$point_id
  }
  if ("x" %in% names(station_coords_std) && !"station_x" %in% names(station_coords_std)) {
    station_coords_std$station_x <- station_coords_std$x
  }
  if ("y" %in% names(station_coords_std) && !"station_y" %in% names(station_coords_std)) {
    station_coords_std$station_y <- station_coords_std$y
  }
  # Handle depth values from raster_value (convert to positive depth)
  if ("raster_value" %in% names(station_coords_std) && !"depth_m" %in% names(station_coords_std)) {
    station_coords_std$depth_m <- abs(station_coords_std$raster_value)
  }

  required_station_cols <- c("station_id", "station_x", "station_y")
  missing_station_cols <- setdiff(required_station_cols, names(station_coords_std))
  if (length(missing_station_cols) > 0) {
    stop("Missing columns in station_coords: ", paste(missing_station_cols, collapse = ", "),
         "\nNote: point_id can substitute for station_id, x/y can substitute for station_x/station_y")
  }

  # Check for depth information (warn if missing but don't stop)
  if (!"depth_m" %in% names(station_coords_std)) {
    warning("No depth information found in station_coords (raster_value or depth_m). Using default 10m depth for all stations.")
    station_coords_std$depth_m <- 10
  }

  required_movement_params <- c("step_length_mean", "step_length_sd", "turning_angle_mean", "turning_angle_sd")
  missing_movement_params <- setdiff(required_movement_params, names(movement_params))
  if (length(missing_movement_params) > 0) {
    stop("Missing movement parameters: ", paste(missing_movement_params, collapse = ", "))
  }

  # =========================================================================
  # HELPER FUNCTIONS
  # =========================================================================

  # Platform-safe parallel processing function
  safe_parallel_apply <- function(X, FUN, ...) {
    if (!use_parallel || length(X) == 1) {
      return(lapply(X, FUN, ...))
    }

    tryCatch({
      if (.Platform$OS.type == "windows") {
        # Windows: use parLapply with cluster
        cl <- parallel::makeCluster(n_cores)
        on.exit(parallel::stopCluster(cl), add = TRUE)

        # Export necessary objects to cluster workers
        parallel::clusterExport(cl, c("de_model", "movement_params", "station_coords_std",
                                      "boundary_raster", "calculate_particle_likelihood"),
                                envir = environment())

        # Load required packages on workers
        parallel::clusterEvalQ(cl, {
          if (!requireNamespace("dplyr", quietly = TRUE)) {
            stop("Package 'dplyr' needed on worker nodes.")
          }
        })

        result <- parallel::parLapply(cl, X, FUN, ...)
      } else {
        # Unix/Mac: use mclapply with better error handling
        result <- parallel::mclapply(X, FUN, ..., mc.cores = n_cores, mc.preschedule = FALSE)

        # Check for errors in mclapply results
        error_indices <- which(sapply(result, function(x) inherits(x, "try-error")))
        if (length(error_indices) > 0) {
          if (verbose) {
            cat("Parallel processing errors detected, falling back to sequential for failed particles\n")
          }
          # Reprocess failed particles sequentially
          for (i in error_indices) {
            result[[i]] <- FUN(X[[i]], ...)
          }
        }
      }
      return(result)
    }, error = function(e) {
      if (verbose) warning("Parallel processing failed, falling back to sequential: ", e$message)
      return(lapply(X, FUN, ...))
    })
  }

  # Calculate likelihood for a single particle
  calculate_particle_likelihood <- function(particle_index, particles_df, current_detections, all_stations) {
    tryCatch({
      fish_x <- particles_df$x[particle_index]
      fish_y <- particles_df$y[particle_index]

      likelihood <- 1
      detecting_stations <- unique(current_detections$station_id)

      # For each station, calculate likelihood
      for (station_id in all_stations$station_id) {
        station_info <- all_stations[all_stations$station_id == station_id, ]

        # Calculate distance from fish to station
        distance <- sqrt((fish_x - station_info$station_x)^2 + (fish_y - station_info$station_y)^2)

        # Get station depth (use actual depth from raster_value if available)
        station_depth <- station_info$depth_m

        # Prepare data for DE model using actual station depth
        newdata <- data.frame(
          dist_m = distance,
          depth_m = station_depth
        )

        # Predict detection probability using station-specific depth
        detection_prob <- tryCatch({
          prob <- stats::predict(de_model, newdata, type = "response")
          pmax(pmin(prob, 0.999), 0.001)  # Bound between 0.1% and 99.9%
        }, error = function(e) {
          # Fallback: simple exponential decay if model fails
          pmax(exp(-distance / 200), 0.001)
        })

        # Check if station detected the fish
        was_detected <- station_id %in% detecting_stations

        # Update likelihood
        if (was_detected) {
          likelihood <- likelihood * detection_prob
        } else {
          likelihood <- likelihood * (1 - detection_prob)
        }
      }

      # Ensure likelihood is a valid numeric value
      if (!is.numeric(likelihood) || is.na(likelihood) || !is.finite(likelihood)) {
        likelihood <- 1e-10  # Small positive value
      }

      return(likelihood)
    }, error = function(e) {
      # Return small positive value if calculation fails
      return(1e-10)
    })
  }

  # Simulate single CRW step with boundary handling (CONSISTENT WITH simulate_fish_tracks)
  simulate_crw_step <- function(current_x, current_y, current_heading, time_diff_hours) {
    # Convert time to seconds for consistency with simulate_fish_tracks
    time_diff_seconds <- time_diff_hours * 3600

    # Generate step length (truncated normal - SAME as simulate_fish_tracks)
    step_length <- stats::rnorm(1, movement_params$step_length_mean, movement_params$step_length_sd)
    step_length <- max(step_length, 5)  # Minimum step length (same as original)

    # Generate turning angle using von Mises distribution (SAME as simulate_fish_tracks)
    # Convert turning angle SD from degrees to radians
    turning_angle_sd_rad <- movement_params$turning_angle_sd * pi / 180
    turning_angle_mean_rad <- movement_params$turning_angle_mean * pi / 180

    # Use von Mises distribution with same kappa calculation
    suppressWarnings({
      if (!requireNamespace("circular", quietly = TRUE)) {
        # Fallback to normal if circular package not available
        turning_angle <- rnorm(1, mean = turning_angle_mean_rad, sd = turning_angle_sd_rad)
      } else {
        turning_angle <- circular::rvonmises(1, mu = turning_angle_mean_rad,
                                             kappa = 1 / (turning_angle_sd_rad^2))
      }
    })

    # Update heading (SAME as simulate_fish_tracks)
    new_heading <- current_heading + turning_angle

    # Calculate new position
    # Scale step length by time (assuming step_length is per time_step from simulate_fish_tracks)
    # Default time_step in simulate_fish_tracks is 60 seconds
    distance <- step_length * (time_diff_seconds / 60)  # Scale by time
    new_x <- current_x + distance * cos(new_heading)
    new_y <- current_y + distance * sin(new_heading)

    # Apply boundary handling if raster is provided
    if (!is.null(boundary_raster)) {
      # Get raster extent
      raster_extent <- raster::extent(boundary_raster)

      # Check if new position is within raster bounds and valid
      if (new_x >= raster_extent@xmin && new_x <= raster_extent@xmax &&
          new_y >= raster_extent@ymin && new_y <= raster_extent@ymax) {

        # Check if position has valid raster value
        cell_value <- raster::extract(boundary_raster, matrix(c(new_x, new_y), ncol = 2))

        if (is.na(cell_value)) {
          # Bounce off invalid areas (like land)
          new_heading <- new_heading + pi  # Reverse direction
          new_x <- current_x + distance * cos(new_heading)
          new_y <- current_y + distance * sin(new_heading)
        }
      } else {
        # Bounce off raster boundaries
        new_heading <- new_heading + pi  # Reverse direction
        new_x <- current_x + distance * cos(new_heading)
        new_y <- current_y + distance * sin(new_heading)
      }
    }

    return(list(x = new_x, y = new_y, heading = new_heading))
  }

  # Calculate detection probability given distance (DEPRECATED - now handled in calculate_detection_likelihood)
  predict_detection_probability <- function(distance) {
    # This function is kept for backwards compatibility but is no longer used
    # Detection probabilities are now calculated with station-specific depths
    # in calculate_detection_likelihood function

    # Prepare data for DE model (using default depth)
    newdata <- data.frame(
      dist_m = distance,
      depth_m = 10  # Default depth - only used as fallback
    )

    # Predict detection probability
    tryCatch({
      prob <- stats::predict(de_model, newdata, type = "response")
      return(pmax(pmin(prob, 0.999), 0.001))  # Bound between 0.1% and 99.9%
    }, error = function(e) {
      # Fallback: simple exponential decay if model fails
      return(pmax(exp(-distance / 200), 0.001))
    })
  }

  # Calculate likelihood of observations given fish position
  calculate_detection_likelihood <- function(fish_x, fish_y, observed_detections, all_stations) {
    likelihood <- 1

    # Get list of stations that detected the fish
    detecting_stations <- unique(observed_detections$station_id)

    # For each station, calculate likelihood
    for (station_id in all_stations$station_id) {
      station_info <- all_stations[all_stations$station_id == station_id, ]

      # Calculate distance from fish to station
      distance <- sqrt((fish_x - station_info$station_x)^2 + (fish_y - station_info$station_y)^2)

      # Get station depth (use actual depth from raster_value if available)
      station_depth <- station_info$depth_m

      # Prepare data for DE model using actual station depth
      newdata <- data.frame(
        dist_m = distance,
        depth_m = station_depth
      )

      # Predict detection probability using station-specific depth
      tryCatch({
        detection_prob <- stats::predict(de_model, newdata, type = "response")
        detection_prob <- pmax(pmin(detection_prob, 0.999), 0.001)  # Bound between 0.1% and 99.9%
      }, error = function(e) {
        # Fallback: simple exponential decay if model fails
        detection_prob <- pmax(exp(-distance / 200), 0.001)
      })

      # Check if station detected the fish
      was_detected <- station_id %in% detecting_stations

      # Update likelihood
      if (was_detected) {
        likelihood <- likelihood * detection_prob
      } else {
        likelihood <- likelihood * (1 - detection_prob)
      }
    }

    return(likelihood)
  }

  # Weighted quantile function
  weighted_quantile <- function(x, weights, probs) {
    if (length(x) == 0) return(NA)

    # Sort by x values
    ord <- order(x)
    x_sorted <- x[ord]
    weights_sorted <- weights[ord]

    # Calculate cumulative weights
    cumsum_weights <- cumsum(weights_sorted) / sum(weights_sorted)

    # Interpolate to find quantiles
    result <- approx(cumsum_weights, x_sorted, xout = probs, rule = 2)$y
    return(result)
  }

  # =========================================================================
  # MAIN PARTICLE FILTER ALGORITHM
  # =========================================================================

  # Get unique fish and time steps
  unique_fish <- sort(unique(detection_data_std$fish_id))
  all_results <- vector("list", length(unique_fish))
  names(all_results) <- paste0("fish_", unique_fish)

  # Process each fish separately
  for (fish_idx in seq_along(unique_fish)) {
    current_fish <- unique_fish[fish_idx]

    if (verbose) cat("Processing fish", current_fish, "...\n")

    # Get detection data for this fish
    fish_detections <- detection_data_std %>%
      dplyr::filter(fish_id == current_fish) %>%
      dplyr::arrange(time_bin)

    # Get unique time steps
    time_steps <- sort(unique(fish_detections$time_bin))
    n_times <- length(time_steps)

    if (n_times < 2) {
      if (verbose) cat("  Warning: Fish", current_fish, "has < 2 time steps. Skipping.\n")
      next
    }

    # Initialize particles at first detection location
    first_detection <- fish_detections %>% dplyr::filter(time_bin == time_steps[1]) %>% dplyr::slice(1)

    particles <- data.frame(
      particle_id = 1:n_particles,
      fish_id = current_fish,
      x = rnorm(n_particles, first_detection$x, 50),  # Initial uncertainty: 50m
      y = rnorm(n_particles, first_detection$y, 50),
      heading = runif(n_particles, 0, 2 * pi),        # Random initial headings
      weight = rep(1 / n_particles, n_particles),
      time_step = time_steps[1]
    )

    # Store all particles
    all_particles_fish <- particles

    # Main particle filter loop with progress bar
    if (verbose) {
      pb <- txtProgressBar(min = 1, max = n_times, style = 3,
                           title = paste("Fish", current_fish))
      cat("  Progress: ")
    }

    for (t in 2:n_times) {
      current_time <- time_steps[t]
      previous_time <- time_steps[t - 1]
      time_diff_hours <- as.numeric(current_time - previous_time) / 3600  # Convert to hours

      # Update progress bar
      if (verbose) {
        setTxtProgressBar(pb, t)
      }

      # PREDICT: Move all particles forward using CRW
      for (p in 1:n_particles) {
        new_state <- simulate_crw_step(
          particles$x[p],
          particles$y[p],
          particles$heading[p],
          time_diff_hours
        )

        particles$x[p] <- new_state$x
        particles$y[p] <- new_state$y
        particles$heading[p] <- new_state$heading
      }

      # UPDATE: Weight particles based on detection data (PARALLEL VERSION)
      current_detections <- fish_detections %>% dplyr::filter(time_bin == current_time)

      # Calculate likelihoods for all particles in parallel
      particle_indices <- 1:n_particles
      particle_likelihoods <- safe_parallel_apply(particle_indices, function(p) {
        calculate_particle_likelihood(p, particles, current_detections, station_coords_std)
      })

      # Check for and handle any errors in parallel results
      if (any(sapply(particle_likelihoods, function(x) !is.numeric(x) || is.na(x) || length(x) != 1))) {
        if (verbose) cat("  Warning: Some parallel likelihood calculations failed, using sequential fallback\n")
        # Fall back to sequential processing for this time step
        particle_likelihoods <- lapply(particle_indices, function(p) {
          calculate_particle_likelihood(p, particles, current_detections, station_coords_std)
        })
      }

      # Convert to numeric vector and ensure all values are valid
      likelihood_values <- as.numeric(unlist(particle_likelihoods))
      if (any(is.na(likelihood_values) | !is.finite(likelihood_values))) {
        if (verbose) cat("  Warning: Invalid likelihood values detected, replacing with small positive values\n")
        likelihood_values[is.na(likelihood_values) | !is.finite(likelihood_values)] <- 1e-10
      }

      # Update particle weights
      particles$weight <- particles$weight * likelihood_values

      # Normalize weights
      total_weight <- sum(particles$weight)
      if (total_weight > 0) {
        particles$weight <- particles$weight / total_weight
      } else {
        # If all weights are 0, reset to uniform
        particles$weight <- rep(1 / n_particles, n_particles)
      }

      # RESAMPLE: Check effective sample size
      eff_sample_size <- 1 / sum(particles$weight^2)

      if (eff_sample_size < resample_threshold) {
        # Resample particles
        resample_indices <- sample(1:n_particles, n_particles, replace = TRUE, prob = particles$weight)
        particles <- particles[resample_indices, ]
        particles$particle_id <- 1:n_particles
        particles$weight <- rep(1 / n_particles, n_particles)
      }

      # Store particles for this time step
      particles$time_step <- current_time
      all_particles_fish <- rbind(all_particles_fish, particles)
    }

    # Close progress bar
    if (verbose) {
      close(pb)
      cat("\n")
    }

    # Calculate position estimates for this fish
    position_estimates_fish <- all_particles_fish %>%
      dplyr::group_by(time_step) %>%
      dplyr::summarise(
        fish_id = dplyr::first(fish_id),

        # Weighted mean position
        x_est = sum(x * weight) / sum(weight),
        y_est = sum(y * weight) / sum(weight),

        # Uncertainty measures
        x_var = sum(weight * (x - x_est)^2) / sum(weight),
        y_var = sum(weight * (y - y_est)^2) / sum(weight),
        position_uncertainty = sqrt(x_var + y_var),

        # Confidence intervals (95%)
        x_ci_lower = weighted_quantile(x, weight, 0.025),
        x_ci_upper = weighted_quantile(x, weight, 0.975),
        y_ci_lower = weighted_quantile(y, weight, 0.025),
        y_ci_upper = weighted_quantile(y, weight, 0.975),

        # Additional metrics
        n_particles_effective = 1 / sum(weight^2),
        max_weight = max(weight),

        .groups = 'drop'
      )

    # Store results for this fish
    all_results[[fish_idx]] <- list(
      position_estimates = position_estimates_fish,
      all_particles = all_particles_fish
    )
  }

  # Combine results across all fish
  combined_position_estimates <- dplyr::bind_rows(lapply(all_results, function(x) x$position_estimates))
  combined_all_particles <- dplyr::bind_rows(lapply(all_results, function(x) x$all_particles))

  # Calculate summary statistics
  summary_stats <- list(
    n_fish = length(unique_fish),
    n_particles = n_particles,
    total_time_steps = nrow(combined_position_estimates),
    mean_position_uncertainty = mean(combined_position_estimates$position_uncertainty, na.rm = TRUE),
    movement_params = movement_params,
    resample_threshold = resample_threshold
  )

  if (verbose) {
    cat("=== PARTICLE FILTER COMPLETE ===\n")
    cat("Fish processed:", summary_stats$n_fish, "\n")
    cat("Total position estimates:", summary_stats$total_time_steps, "\n")
    cat("Mean position uncertainty:", round(summary_stats$mean_position_uncertainty, 1), "meters\n")
  }

  # Return comprehensive results
  return(list(
    position_estimates = combined_position_estimates,
    all_particles = combined_all_particles,
    summary_stats = summary_stats
  ))
}
