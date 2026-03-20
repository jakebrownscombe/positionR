#' Create Detection Efficiency Lookup Grid
#'
#' Pre-computes a dense grid of detection efficiency values from a DE model,
#' enabling fast array-indexed lookups instead of repeated \code{stats::predict()} calls.
#'
#' @param de_model A fitted model object (e.g., logistic regression) that predicts
#'   detection probability from \code{dist_m} and \code{depth_m}.
#' @param max_distance Numeric. Maximum distance in metres for the lookup grid.
#' @param min_depth Numeric. Minimum depth in metres (default 0).
#' @param max_depth Numeric. Maximum depth in metres (default 50).
#' @param n_dist Integer. Number of distance bins (default 2000).
#' @param n_depth Integer. Number of depth bins (default 100).
#'
#' @return A list with components:
#'   \item{grid}{Matrix of detection probabilities (n_dist x n_depth)}
#'   \item{dist_step}{Distance bin width}
#'   \item{depth_step}{Depth bin width}
#'   \item{min_depth}{Minimum depth value}
#'   \item{max_distance}{Maximum distance}
#'   \item{n_dist}{Number of distance bins}
#'   \item{n_depth}{Number of depth bins}
#'
#' @export
create_de_lookup <- function(de_model, max_distance, min_depth = 0, max_depth = 50,
                              n_dist = 2000, n_depth = 100) {
  dist_step <- max_distance / n_dist
  depth_step <- (max_depth - min_depth) / n_depth

  # Create grid vectors

  dist_vec <- seq(0, max_distance, length.out = n_dist)
  depth_vec <- seq(min_depth, max_depth, length.out = n_depth)

  # Expand to full grid for a single predict() call
  grid_df <- expand.grid(dist_m = dist_vec, depth_m = depth_vec)
  probs <- stats::predict(de_model, newdata = grid_df, type = "response")
  probs <- pmax(pmin(probs, 0.999), 0.001)

  # Reshape to matrix
  grid_mat <- matrix(probs, nrow = n_dist, ncol = n_depth)

  list(
    grid = grid_mat,
    dist_step = dist_step,
    depth_step = depth_step,
    min_depth = min_depth,
    max_distance = max_distance,
    n_dist = n_dist,
    n_depth = n_depth
  )
}

#' Fast DE Lookup from Pre-computed Grid
#'
#' Vectorized lookup of detection efficiency values using array indexing.
#'
#' @param distances Numeric vector or matrix of distances (metres).
#' @param depths Numeric vector of station depths (recycled across particles).
#' @param de_lookup List from \code{create_de_lookup()}.
#'
#' @return Numeric vector or matrix of detection probabilities (same shape as distances).
#' @keywords internal
de_lookup_fast <- function(distances, depths, de_lookup) {
  # Convert distances and depths to grid indices
  di <- pmin(pmax(round(distances / de_lookup$dist_step) + 1L, 1L), de_lookup$n_dist)
  zi <- pmin(pmax(round((depths - de_lookup$min_depth) / de_lookup$depth_step) + 1L, 1L), de_lookup$n_depth)

  # Matrix indexing for fast lookup
  if (is.matrix(distances)) {
    # distances is n_particles x n_stations, depths is n_stations (recycled)
    result <- matrix(0, nrow = nrow(distances), ncol = ncol(distances))
    for (s in seq_len(ncol(distances))) {
      result[, s] <- de_lookup$grid[cbind(di[, s], zi[s])]
    }
    return(result)
  } else {
    return(de_lookup$grid[cbind(di, zi)])
  }
}


#' Particle Filter Positioning for Acoustic Telemetry
#'
#' Estimates fish positions using a Sequential Monte Carlo (particle filter) approach.
#' Particles represent possible fish locations and are propagated forward using a
#' Correlated Random Walk, then weighted by detection likelihood at each time step.
#'
#' @param detection_data Data frame of detections with columns for fish ID, time,
#'   station ID, and detection status. Should include both detecting and non-detecting
#'   stations per time step (e.g., from \code{prepare_detection_data_for_wade()}).
#' @param station_info Data frame or sf object with station coordinates and depths.
#'   Must include station_id, x/y coordinates, and depth (raster_value or depth_m).
#' @param de_model A fitted detection efficiency model (e.g., from
#'   \code{create_logistic_curve_depth()}).
#' @param raster A RasterLayer defining the study area. NA cells are barriers.
#' @param n_particles Integer. Number of particles per fish (default 1000).
#' @param step_length_mean Numeric. Mean step length in metres per time_step seconds (default 50).
#' @param step_length_sd Numeric. SD of step length (default 30).
#' @param turning_angle_sd Numeric. SD of turning angle in degrees (default 45).
#' @param time_step Numeric. Reference time step in seconds for movement scaling (default 180).
#' @param max_distance Numeric. Maximum DE lookup distance in metres (default 30000).
#' @param fish_id_col Character. Name of fish ID column (default "path_id").
#' @param time_col Character. Name of time column (default "datetime").
#' @param station_col Character. Name of station ID column (default "station_id").
#' @param ess_threshold Numeric. Effective sample size threshold for resampling,
#'   as fraction of n_particles (default 0.5).
#' @param return_particles Logical. Whether to return full particle histories (default FALSE).
#'   Can be memory-intensive for large datasets.
#' @param verbose Logical. Print progress messages (default TRUE).
#'
#' @return A list with components:
#'   \item{positions}{Data frame with weighted mean position estimates per fish/time}
#'   \item{particles}{(Optional) Data frame with all particle positions and weights}
#'
#' @details
#' The algorithm proceeds in three steps at each time point:
#' \enumerate{
#'   \item \strong{Predict}: Move all particles forward using a Correlated Random Walk.
#'     Step lengths scale with time between observations.
#'   \item \strong{Update}: Weight particles by detection likelihood. For each particle,
#'     detection probability is looked up from a pre-computed grid (not model prediction),
#'     and the likelihood is the product across all stations of P(detected|position) or
#'     P(not detected|position).
#'   \item \strong{Resample}: When the effective sample size drops below the threshold,
#'     particles are resampled with replacement proportional to their weights using
#'     systematic resampling.
#' }
#'
#' Detection efficiency lookups use a pre-computed grid (via \code{create_de_lookup()})
#' instead of calling \code{stats::predict()}, providing orders-of-magnitude speedup.
#'
#' @examples
#' \dontrun{
#' results <- particle_filter_positioning(
#'   detection_data = wade_data$station_detections,
#'   station_info = receiver_stations,
#'   de_model = logistic_DE$log_model,
#'   raster = depth_raster,
#'   n_particles = 1000,
#'   fish_id_col = "path_id",
#'   time_col = "datetime",
#'   station_col = "station_id"
#' )
#' }
#'
#' @export
particle_filter_positioning <- function(detection_data,
                                         station_info,
                                         de_model,
                                         raster,
                                         n_particles = 1000,
                                         step_length_mean = 50,
                                         step_length_sd = 30,
                                         turning_angle_sd = 45,
                                         time_step = 180,
                                         max_distance = 30000,
                                         fish_id_col = "path_id",
                                         time_col = "datetime",
                                         station_col = "station_id",
                                         ess_threshold = 0.5,
                                         return_particles = FALSE,
                                         verbose = TRUE) {

  # --- Setup ---
  dt <- data.table::as.data.table(detection_data)

  # Standardize column names
  if (fish_id_col != "fish_id" && fish_id_col %in% names(dt)) {
    data.table::setnames(dt, fish_id_col, "fish_id")
  }
  if (time_col != "time" && time_col %in% names(dt)) {
    data.table::setnames(dt, time_col, "time")
  }
  if (station_col != "station_id" && station_col %in% names(dt)) {
    data.table::setnames(dt, station_col, "station_id")
  }

  # Extract station coordinates and depths
  if ("sf" %in% class(station_info)) {
    coords <- sf::st_coordinates(station_info)
    stn <- data.table::data.table(
      station_id = if ("station_id" %in% names(station_info)) station_info$station_id
                   else station_info$point_id,
      station_x = coords[, 1],
      station_y = coords[, 2],
      depth_m = abs(if ("raster_value" %in% names(station_info)) station_info$raster_value
                    else if ("depth_m" %in% names(station_info)) station_info$depth_m
                    else rep(10, nrow(station_info)))
    )
  } else {
    stn <- data.table::as.data.table(station_info)
    if (!"station_x" %in% names(stn) && "x" %in% names(stn)) {
      data.table::setnames(stn, c("x", "y"), c("station_x", "station_y"))
    }
    if (!"depth_m" %in% names(stn) && "raster_value" %in% names(stn)) {
      stn[, depth_m := abs(raster_value)]
    }
  }

  # Pre-compute DE lookup grid
  depth_range <- range(stn$depth_m, na.rm = TRUE)
  de_lookup <- create_de_lookup(
    de_model = de_model,
    max_distance = max_distance,
    min_depth = max(0, depth_range[1] - 5),
    max_depth = depth_range[2] + 5
  )

  if (verbose) cat("DE lookup grid pre-computed (",
                    de_lookup$n_dist, "x", de_lookup$n_depth, "bins)\n")

  # Raster extent for boundary checking
  raster_ext <- raster::extent(raster)

  # Get unique fish

  fish_ids <- unique(dt$fish_id)
  n_fish <- length(fish_ids)

  # Movement parameters
  turning_angle_sd_rad <- turning_angle_sd * pi / 180
  kappa <- 1 / (turning_angle_sd_rad^2)

  # --- Helper: vectorized validity check ---
  is_valid <- function(x, y) {
    in_extent <- x >= raster_ext@xmin & x <= raster_ext@xmax &
                 y >= raster_ext@ymin & y <= raster_ext@ymax
    valid <- rep(FALSE, length(x))
    valid[in_extent] <- !is.na(raster::extract(raster, cbind(x[in_extent], y[in_extent])))
    valid
  }

  # --- Helper: vectorized CRW move ---
  move_particles <- function(px, py, headings, time_diff_sec) {
    n <- length(px)
    time_scale <- time_diff_sec / time_step

    # Generate step lengths and turning angles (vectorized)
    steps <- pmax(stats::rnorm(n, step_length_mean, step_length_sd), 5) * time_scale
    turns <- as.numeric(circular::rvonmises(n, mu = circular::circular(0), kappa = kappa))
    new_headings <- headings + turns

    # Propose new positions
    new_x <- px + steps * cos(new_headings)
    new_y <- py + steps * sin(new_headings)

    # Check validity
    valid <- is_valid(new_x, new_y)

    # For invalid particles: try reflection (reverse direction, half step)
    if (any(!valid)) {
      inv <- which(!valid)
      for (attempt in 1:3) {
        ref_headings <- new_headings[inv] + pi + stats::runif(length(inv), -0.5, 0.5)
        ref_steps <- steps[inv] * (0.5 ^ attempt)
        new_x[inv] <- px[inv] + ref_steps * cos(ref_headings)
        new_y[inv] <- py[inv] + ref_steps * sin(ref_headings)
        still_invalid <- !is_valid(new_x[inv], new_y[inv])
        if (!any(still_invalid)) break
        inv <- inv[still_invalid]
      }
      # Last resort: stay in place
      still_bad <- which(!is_valid(new_x, new_y))
      if (length(still_bad) > 0) {
        new_x[still_bad] <- px[still_bad]
        new_y[still_bad] <- py[still_bad]
        new_headings[still_bad] <- headings[still_bad]
      }
    }

    list(x = new_x, y = new_y, heading = new_headings)
  }

  # --- Helper: vectorized likelihood ---
  calculate_likelihoods <- function(px, py, detecting_ids, stn_df, de_lookup) {
    n <- length(px)
    n_stn <- nrow(stn_df)

    # Distance matrix: n_particles x n_stations
    dx <- outer(px, stn_df$station_x, "-")
    dy <- outer(py, stn_df$station_y, "-")
    dist_mat <- sqrt(dx^2 + dy^2)

    # DE lookup (vectorized)
    de_mat <- de_lookup_fast(dist_mat, stn_df$depth_m, de_lookup)

    # Log-likelihood via matrix operations
    detected_mask <- stn_df$station_id %in% detecting_ids

    log_lik <- numeric(n)
    if (any(detected_mask)) {
      log_lik <- log_lik + rowSums(log(de_mat[, detected_mask, drop = FALSE]))
    }
    if (any(!detected_mask)) {
      log_lik <- log_lik + rowSums(log(1 - de_mat[, !detected_mask, drop = FALSE]))
    }

    # Log-sum-exp for numerical stability
    log_lik_shifted <- log_lik - max(log_lik)
    weights <- exp(log_lik_shifted)
    weights[!is.finite(weights)] <- 1e-300
    weights
  }

  # --- Helper: systematic resampling (more efficient than multinomial) ---
  systematic_resample <- function(weights, n) {
    cumw <- cumsum(weights / sum(weights))
    u <- (stats::runif(1) + seq(0, n - 1)) / n
    indices <- findInterval(u, cumw) + 1L
    pmin(indices, length(weights))
  }

  # --- Process each fish ---
  all_positions <- list()
  all_particles <- if (return_particles) list() else NULL

  for (fi in seq_along(fish_ids)) {
    current_fish <- fish_ids[fi]
    if (verbose && interactive()) {
      cat("\rProcessing fish", fi, "of", n_fish, "(", current_fish, ")    ")
      flush.console()
    }

    # Get time steps for this fish — only those with at least one detection
    fish_dt <- dt[fish_id == current_fish]

    # Find the first time step with a detection to initialize particles
    det_times <- sort(unique(fish_dt[detected == 1]$time))
    all_times <- sort(unique(fish_dt$time))

    if (length(det_times) < 2) {
      if (verbose) cat(" - skipping (< 2 detection events)\n")
      next
    }

    # Use all time steps from first detection onward
    start_idx <- which(all_times == det_times[1])
    time_steps <- all_times[start_idx:length(all_times)]
    n_times <- length(time_steps)

    # Initialize particles at first detecting station(s)
    first_dets <- fish_dt[time == det_times[1] & detected == 1]
    init_stations <- stn[station_id %in% first_dets$station_id]
    init_x <- mean(init_stations$station_x)
    init_y <- mean(init_stations$station_y)

    # Scatter particles around init point
    px <- init_x + stats::rnorm(n_particles, 0, 200)
    py <- init_y + stats::rnorm(n_particles, 0, 200)
    headings <- stats::runif(n_particles, 0, 2 * pi)
    weights <- rep(1 / n_particles, n_particles)

    # Ensure valid starting positions
    valid <- is_valid(px, py)
    if (any(!valid)) {
      px[!valid] <- init_x
      py[!valid] <- init_y
    }

    # Storage for position estimates
    fish_positions <- data.table::data.table(
      fish_id = rep(current_fish, n_times),
      time = time_steps,
      x_mean = numeric(n_times),
      y_mean = numeric(n_times),
      x_sd = numeric(n_times),
      y_sd = numeric(n_times),
      ess = numeric(n_times),
      n_detecting = integer(n_times),
      resampled = logical(n_times)
    )

    fish_particle_list <- if (return_particles) vector("list", n_times) else NULL

    # First time step — just weight by likelihood (no movement yet)
    dets_t <- fish_dt[time == time_steps[1] & detected == 1]
    detecting_ids <- unique(dets_t$station_id)
    weights <- calculate_likelihoods(px, py, detecting_ids, stn, de_lookup)
    weights <- weights / sum(weights)

    xm <- sum(px * weights)
    ym <- sum(py * weights)
    data.table::set(fish_positions, 1L, "x_mean", xm)
    data.table::set(fish_positions, 1L, "y_mean", ym)
    data.table::set(fish_positions, 1L, "x_sd", sqrt(sum(weights * (px - xm)^2)))
    data.table::set(fish_positions, 1L, "y_sd", sqrt(sum(weights * (py - ym)^2)))
    data.table::set(fish_positions, 1L, "ess", 1 / sum(weights^2))
    data.table::set(fish_positions, 1L, "n_detecting", length(detecting_ids))
    data.table::set(fish_positions, 1L, "resampled", FALSE)

    if (return_particles) {
      fish_particle_list[[1]] <- data.table::data.table(
        fish_id = current_fish, time = time_steps[1],
        particle = seq_len(n_particles), x = px, y = py, weight = weights
      )
    }

    # --- Main filter loop ---
    for (t in 2:n_times) {
      # Time difference
      if (inherits(time_steps[t], "POSIXct")) {
        time_diff_sec <- as.numeric(difftime(time_steps[t], time_steps[t - 1], units = "secs"))
      } else {
        time_diff_sec <- as.numeric(time_steps[t] - time_steps[t - 1])
      }
      time_diff_sec <- max(time_diff_sec, 1)  # Avoid zero

      # PREDICT: move particles
      moved <- move_particles(px, py, headings, time_diff_sec)
      px <- moved$x
      py <- moved$y
      headings <- moved$heading

      # UPDATE: calculate likelihoods
      dets_t <- fish_dt[time == time_steps[t] & detected == 1]
      detecting_ids <- unique(dets_t$station_id)
      lik <- calculate_likelihoods(px, py, detecting_ids, stn, de_lookup)

      # Update weights
      weights <- weights * lik
      w_sum <- sum(weights)
      if (w_sum > 0 && is.finite(w_sum)) {
        weights <- weights / w_sum
      } else {
        weights <- rep(1 / n_particles, n_particles)
      }

      # Effective sample size
      ess <- 1 / sum(weights^2)
      did_resample <- FALSE

      # RESAMPLE if ESS too low
      if (ess < n_particles * ess_threshold) {
        idx <- systematic_resample(weights, n_particles)
        px <- px[idx]
        py <- py[idx]
        headings <- headings[idx]
        weights <- rep(1 / n_particles, n_particles)
        ess <- n_particles
        did_resample <- TRUE
      }

      # Store weighted position estimate
      xm <- sum(px * weights)
      ym <- sum(py * weights)
      ti <- as.integer(t)
      data.table::set(fish_positions, ti, "x_mean", xm)
      data.table::set(fish_positions, ti, "y_mean", ym)
      data.table::set(fish_positions, ti, "x_sd", sqrt(sum(weights * (px - xm)^2)))
      data.table::set(fish_positions, ti, "y_sd", sqrt(sum(weights * (py - ym)^2)))
      data.table::set(fish_positions, ti, "ess", ess)
      data.table::set(fish_positions, ti, "n_detecting", length(detecting_ids))
      data.table::set(fish_positions, ti, "resampled", did_resample)

      if (return_particles) {
        fish_particle_list[[t]] <- data.table::data.table(
          fish_id = current_fish, time = time_steps[t],
          particle = seq_len(n_particles), x = px, y = py, weight = weights
        )
      }
    }

    all_positions[[fi]] <- fish_positions
    if (return_particles) {
      all_particles[[fi]] <- data.table::rbindlist(fish_particle_list)
    }
  }

  if (verbose && interactive()) cat("\n")
  if (verbose) cat("Particle filter complete for", n_fish, "fish.\n")

  # Combine results
  positions <- data.table::rbindlist(all_positions)

  result <- list(
    positions = as.data.frame(positions)
  )
  if (return_particles) {
    result$particles <- as.data.frame(data.table::rbindlist(all_particles))
  }

  return(result)
}
