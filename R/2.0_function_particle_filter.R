# Suppress R CMD check notes for data.table NSE variables
utils::globalVariables(c("fish_id", "detected", "time", "station_id",
                         "depth_m", "raster_value", "total_dets",
                         "has_detections"))

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
  dist_vec <- seq(0, max_distance, length.out = n_dist)
  depth_vec <- seq(min_depth, max_depth, length.out = n_depth)
  grid_df <- expand.grid(dist_m = dist_vec, depth_m = depth_vec)
  probs <- stats::predict(de_model, newdata = grid_df, type = "response")
  probs <- pmax(pmin(probs, 0.999), 0.001)
  grid_mat <- matrix(probs, nrow = n_dist, ncol = n_depth)
  list(grid = grid_mat, dist_step = dist_step, depth_step = depth_step,
       min_depth = min_depth, max_distance = max_distance,
       n_dist = n_dist, n_depth = n_depth)
}

#' Fast DE Lookup from Pre-computed Grid
#'
#' @param distances Numeric vector or matrix of distances (metres).
#' @param depths Numeric vector of station depths (recycled across particles).
#' @param de_lookup List from \code{create_de_lookup()}.
#' @return Numeric vector or matrix of detection probabilities.
#' @keywords internal
de_lookup_fast <- function(distances, depths, de_lookup) {
  di <- pmin(pmax(round(distances / de_lookup$dist_step) + 1L, 1L), de_lookup$n_dist)
  zi <- pmin(pmax(round((depths - de_lookup$min_depth) / de_lookup$depth_step) + 1L, 1L), de_lookup$n_depth)
  if (is.matrix(distances)) {
    result <- matrix(0, nrow = nrow(distances), ncol = ncol(distances))
    for (s in seq_len(ncol(distances))) {
      result[, s] <- de_lookup$grid[cbind(di[, s], zi[s])]
    }
    return(result)
  } else {
    return(de_lookup$grid[cbind(di, zi)])
  }
}


# ============================================================================
# Internal helpers (package-level, not exported)
# ============================================================================

#' @keywords internal
pf_is_valid <- function(x, y, raster, raster_ext) {
  in_extent <- x >= raster_ext@xmin & x <= raster_ext@xmax &
               y >= raster_ext@ymin & y <= raster_ext@ymax
  valid <- rep(FALSE, length(x))
  valid[in_extent] <- !is.na(raster::extract(raster, cbind(x[in_extent], y[in_extent])))
  valid
}

#' @keywords internal
pf_move_particles <- function(px, py, headings, time_diff_sec,
                               step_length_mean, step_length_sd, kappa,
                               time_step, raster, raster_ext) {
  n <- length(px)
  time_scale <- time_diff_sec / time_step
  steps <- pmax(stats::rnorm(n, step_length_mean, step_length_sd), 5) * time_scale
  turns <- as.numeric(circular::rvonmises(n, mu = circular::circular(0), kappa = kappa))
  new_headings <- headings + turns
  new_x <- px + steps * cos(new_headings)
  new_y <- py + steps * sin(new_headings)

  valid <- pf_is_valid(new_x, new_y, raster, raster_ext)
  if (any(!valid)) {
    inv <- which(!valid)
    for (attempt in 1:3) {
      ref_headings <- new_headings[inv] + pi + stats::runif(length(inv), -0.5, 0.5)
      ref_steps <- steps[inv] * (0.5 ^ attempt)
      new_x[inv] <- px[inv] + ref_steps * cos(ref_headings)
      new_y[inv] <- py[inv] + ref_steps * sin(ref_headings)
      still_invalid <- !pf_is_valid(new_x[inv], new_y[inv], raster, raster_ext)
      if (!any(still_invalid)) break
      inv <- inv[still_invalid]
    }
    still_bad <- which(!pf_is_valid(new_x, new_y, raster, raster_ext))
    if (length(still_bad) > 0) {
      new_x[still_bad] <- px[still_bad]
      new_y[still_bad] <- py[still_bad]
      new_headings[still_bad] <- headings[still_bad]
    }
  }
  list(x = new_x, y = new_y, heading = new_headings)
}

#' @keywords internal
pf_calculate_likelihoods <- function(px, py, detecting_ids, stn_df, de_lookup) {
  n <- length(px)
  dx <- outer(px, stn_df$station_x, "-")
  dy <- outer(py, stn_df$station_y, "-")
  dist_mat <- sqrt(dx^2 + dy^2)
  de_mat <- de_lookup_fast(dist_mat, stn_df$depth_m, de_lookup)
  detected_mask <- stn_df$station_id %in% detecting_ids
  log_lik <- numeric(n)
  if (any(detected_mask)) {
    log_lik <- log_lik + rowSums(log(de_mat[, detected_mask, drop = FALSE]))
  }
  if (any(!detected_mask)) {
    log_lik <- log_lik + rowSums(log(1 - de_mat[, !detected_mask, drop = FALSE]))
  }
  log_lik_shifted <- log_lik - max(log_lik)
  weights <- exp(log_lik_shifted)
  weights[!is.finite(weights)] <- 1e-300
  weights
}

#' @keywords internal
pf_systematic_resample <- function(weights, n) {
  cumw <- cumsum(weights / sum(weights))
  u <- (stats::runif(1) + seq(0, n - 1)) / n
  indices <- findInterval(u, cumw) + 1L
  pmin(indices, length(weights))
}

#' Prepare station data for particle filter
#' @keywords internal
pf_prepare_stations <- function(station_info) {
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
  stn
}

#' Standardize detection data column names
#' @keywords internal
pf_prepare_detections <- function(detection_data, fish_id_col, time_col, station_col) {
  dt <- data.table::as.data.table(detection_data)
  if (fish_id_col != "fish_id" && fish_id_col %in% names(dt)) {
    data.table::setnames(dt, fish_id_col, "fish_id")
  }
  if (time_col != "time" && time_col %in% names(dt)) {
    data.table::setnames(dt, time_col, "time")
  }
  if (station_col != "station_id" && station_col %in% names(dt)) {
    data.table::setnames(dt, station_col, "station_id")
  }
  dt
}


# ============================================================================
# Main particle filter
# ============================================================================

#' Particle Filter Positioning for Acoustic Telemetry
#'
#' Estimates fish positions using a Sequential Monte Carlo (particle filter) approach.
#' Particles represent possible fish locations and are propagated forward using a
#' Correlated Random Walk, then weighted by detection likelihood at each time step.
#'
#' @param detection_data Data frame of detections with columns for fish ID, time,
#'   station ID, and detection status.
#' @param station_info Data frame or sf object with station coordinates and depths.
#' @param de_model A fitted detection efficiency model.
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
#' @param ess_threshold Numeric. ESS threshold for resampling as fraction of n_particles (default 0.5).
#' @param return_particles Logical. Return full particle histories (default FALSE).
#' @param verbose Logical. Print progress (default TRUE).
#'
#' @return A list with components:
#'   \item{positions}{Data frame with weighted mean position estimates per fish/time}
#'   \item{particles}{(Optional) Data frame with all particle positions and weights}
#'
#' @export
particle_filter_positioning <- function(detection_data, station_info, de_model, raster,
                                         n_particles = 1000,
                                         step_length_mean = 50, step_length_sd = 30,
                                         turning_angle_sd = 45, time_step = 180,
                                         max_distance = 30000,
                                         fish_id_col = "path_id", time_col = "datetime",
                                         station_col = "station_id",
                                         ess_threshold = 0.5,
                                         return_particles = FALSE, verbose = TRUE) {

  dt <- pf_prepare_detections(detection_data, fish_id_col, time_col, station_col)
  stn <- pf_prepare_stations(station_info)

  depth_range <- range(stn$depth_m, na.rm = TRUE)
  de_lookup <- create_de_lookup(de_model, max_distance,
                                 min_depth = max(0, depth_range[1] - 5),
                                 max_depth = depth_range[2] + 5)
  if (verbose) cat("DE lookup grid pre-computed (",
                    de_lookup$n_dist, "x", de_lookup$n_depth, "bins)\n")

  raster_ext <- raster::extent(raster)
  fish_ids <- unique(dt$fish_id)
  n_fish <- length(fish_ids)
  turning_angle_sd_rad <- turning_angle_sd * pi / 180
  kappa <- 1 / (turning_angle_sd_rad^2)

  all_positions <- list()
  all_particles <- if (return_particles) list() else NULL

  for (fi in seq_along(fish_ids)) {
    current_fish <- fish_ids[fi]
    if (verbose && interactive()) {
      cat("\rProcessing fish", fi, "of", n_fish, "(", current_fish, ")    ")
      utils::flush.console()
    }

    fish_dt <- dt[fish_id == current_fish]
    det_times <- sort(unique(fish_dt[detected == 1]$time))
    all_times <- sort(unique(fish_dt$time))

    if (length(det_times) < 2) {
      if (verbose) cat(" - skipping (< 2 detection events)\n")
      next
    }

    start_idx <- which(all_times == det_times[1])
    time_steps_vec <- all_times[start_idx:length(all_times)]
    n_times <- length(time_steps_vec)

    first_dets <- fish_dt[time == det_times[1] & detected == 1]
    init_stations <- stn[station_id %in% first_dets$station_id]
    init_x <- mean(init_stations$station_x)
    init_y <- mean(init_stations$station_y)

    px <- init_x + stats::rnorm(n_particles, 0, 200)
    py <- init_y + stats::rnorm(n_particles, 0, 200)
    headings <- stats::runif(n_particles, 0, 2 * pi)
    weights <- rep(1 / n_particles, n_particles)

    valid <- pf_is_valid(px, py, raster, raster_ext)
    if (any(!valid)) { px[!valid] <- init_x; py[!valid] <- init_y }

    fish_positions <- data.table::data.table(
      fish_id = rep(current_fish, n_times), time = time_steps_vec,
      x_mean = numeric(n_times), y_mean = numeric(n_times),
      x_sd = numeric(n_times), y_sd = numeric(n_times),
      ess = numeric(n_times), n_detecting = integer(n_times),
      resampled = logical(n_times)
    )
    fish_particle_list <- if (return_particles) vector("list", n_times) else NULL

    # First time step
    dets_t <- fish_dt[time == time_steps_vec[1] & detected == 1]
    detecting_ids <- unique(dets_t$station_id)
    weights <- pf_calculate_likelihoods(px, py, detecting_ids, stn, de_lookup)
    weights <- weights / sum(weights)
    xm <- sum(px * weights); ym <- sum(py * weights)
    data.table::set(fish_positions, 1L, "x_mean", xm)
    data.table::set(fish_positions, 1L, "y_mean", ym)
    data.table::set(fish_positions, 1L, "x_sd", sqrt(sum(weights * (px - xm)^2)))
    data.table::set(fish_positions, 1L, "y_sd", sqrt(sum(weights * (py - ym)^2)))
    data.table::set(fish_positions, 1L, "ess", 1 / sum(weights^2))
    data.table::set(fish_positions, 1L, "n_detecting", length(detecting_ids))
    data.table::set(fish_positions, 1L, "resampled", FALSE)
    if (return_particles) {
      fish_particle_list[[1]] <- data.table::data.table(
        fish_id = current_fish, time = time_steps_vec[1],
        particle = seq_len(n_particles), x = px, y = py, weight = weights)
    }

    # Main loop
    for (t in 2:n_times) {
      if (inherits(time_steps_vec[t], "POSIXct")) {
        time_diff_sec <- as.numeric(difftime(time_steps_vec[t], time_steps_vec[t-1], units = "secs"))
      } else {
        time_diff_sec <- as.numeric(time_steps_vec[t] - time_steps_vec[t-1])
      }
      time_diff_sec <- max(time_diff_sec, 1)

      moved <- pf_move_particles(px, py, headings, time_diff_sec,
                                  step_length_mean, step_length_sd, kappa,
                                  time_step, raster, raster_ext)
      px <- moved$x; py <- moved$y; headings <- moved$heading

      dets_t <- fish_dt[time == time_steps_vec[t] & detected == 1]
      detecting_ids <- unique(dets_t$station_id)
      lik <- pf_calculate_likelihoods(px, py, detecting_ids, stn, de_lookup)

      weights <- weights * lik
      w_sum <- sum(weights)
      if (w_sum > 0 && is.finite(w_sum)) { weights <- weights / w_sum
      } else { weights <- rep(1 / n_particles, n_particles) }

      ess <- 1 / sum(weights^2)
      did_resample <- FALSE
      if (ess < n_particles * ess_threshold) {
        idx <- pf_systematic_resample(weights, n_particles)
        px <- px[idx]; py <- py[idx]; headings <- headings[idx]
        weights <- rep(1 / n_particles, n_particles)
        ess <- n_particles; did_resample <- TRUE
      }

      xm <- sum(px * weights); ym <- sum(py * weights)
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
          fish_id = current_fish, time = time_steps_vec[t],
          particle = seq_len(n_particles), x = px, y = py, weight = weights)
      }
    }
    all_positions[[fi]] <- fish_positions
    if (return_particles) all_particles[[fi]] <- data.table::rbindlist(fish_particle_list)
  }

  if (verbose && interactive()) cat("\n")
  if (verbose) cat("Particle filter complete for", n_fish, "fish.\n")

  positions <- data.table::rbindlist(all_positions)
  result <- list(positions = as.data.frame(positions))
  if (return_particles) result$particles <- as.data.frame(data.table::rbindlist(all_particles))
  return(result)
}


# ============================================================================
# Two-filter smoother
# ============================================================================

#' Smooth Particle Filter Positions (Two-Filter Smoother)
#'
#' Refines particle filter position estimates by running a backward filter and
#' combining forward and backward particle weights. At each time step, the
#' smoothed distribution incorporates both past and future observations,
#' producing more accurate position estimates than the forward filter alone.
#'
#' @param pf_results Output from \code{particle_filter_positioning()} with
#'   \code{return_particles = TRUE}.
#' @param detection_data Original detection data (same as passed to the filter).
#' @param station_info Station coordinates and depths.
#' @param de_model Detection efficiency model.
#' @param raster Study area RasterLayer.
#' @param step_length_mean,step_length_sd,turning_angle_sd,time_step Movement
#'   parameters (should match the forward filter).
#' @param max_distance Maximum DE lookup distance (default 30000).
#' @param n_backward_particles Number of backward particles per fish. Defaults
#'   to the same as the forward filter.
#' @param fish_id_col,time_col,station_col Column name mappings.
#' @param ess_threshold ESS threshold for resampling (default 0.5).
#' @param return_particles Return smoothed particle clouds (default FALSE).
#' @param verbose Print progress (default TRUE).
#'
#' @return A list with components:
#'   \item{positions}{Data frame with smoothed position estimates}
#'   \item{particles}{(Optional) Smoothed particle clouds}
#'
#' @details
#' The two-filter smoother works by:
#' \enumerate{
#'   \item Using the forward particle clouds from the initial filter
#'   \item Running a backward filter from the last time step to the first
#'   \item At each time step, combining forward and backward log-likelihoods
#'     and re-normalizing to produce the smoothed marginal distribution
#' }
#'
#' @export
pf_smooth <- function(pf_results, detection_data, station_info, de_model, raster,
                       step_length_mean = 50, step_length_sd = 30,
                       turning_angle_sd = 45, time_step = 180,
                       max_distance = 30000, n_backward_particles = NULL,
                       fish_id_col = "path_id", time_col = "datetime",
                       station_col = "station_id", ess_threshold = 0.5,
                       return_particles = FALSE, verbose = TRUE) {

  if (is.null(pf_results$particles)) {
    stop("pf_results must include particles (run particle_filter_positioning with return_particles = TRUE)")
  }

  dt <- pf_prepare_detections(detection_data, fish_id_col, time_col, station_col)
  stn <- pf_prepare_stations(station_info)

  depth_range <- range(stn$depth_m, na.rm = TRUE)
  de_lookup <- create_de_lookup(de_model, max_distance,
                                 min_depth = max(0, depth_range[1] - 5),
                                 max_depth = depth_range[2] + 5)

  raster_ext <- raster::extent(raster)
  turning_angle_sd_rad <- turning_angle_sd * pi / 180
  kappa <- 1 / (turning_angle_sd_rad^2)

  fwd_particles <- data.table::as.data.table(pf_results$particles)
  fish_ids <- unique(fwd_particles$fish_id)
  n_fish <- length(fish_ids)

  if (is.null(n_backward_particles)) {
    n_backward_particles <- max(fwd_particles[fish_id == fish_ids[1] & time == fwd_particles$time[1], .N])
  }
  n_bp <- n_backward_particles

  all_positions <- list()
  all_particles_out <- if (return_particles) list() else NULL

  for (fi in seq_along(fish_ids)) {
    current_fish <- fish_ids[fi]
    if (verbose && interactive()) {
      cat("\rSmoothing fish", fi, "of", n_fish, "(", current_fish, ")    ")
      utils::flush.console()
    }

    fish_fwd <- fwd_particles[fish_id == current_fish]
    time_steps_vec <- sort(unique(fish_fwd$time))
    n_times <- length(time_steps_vec)
    if (n_times < 2) next

    fish_dt <- dt[fish_id == current_fish]

    # --- Backward filter ---
    # Initialize at last time step from forward particles
    last_fwd <- fish_fwd[time == time_steps_vec[n_times]]
    bpx <- last_fwd$x
    bpy <- last_fwd$y
    if (length(bpx) != n_bp) {
      # Resample to match n_backward_particles
      idx <- sample(length(bpx), n_bp, replace = TRUE, prob = last_fwd$weight)
      bpx <- bpx[idx]; bpy <- bpy[idx]
    }
    b_headings <- stats::runif(n_bp, 0, 2 * pi)
    b_weights <- rep(1 / n_bp, n_bp)

    # Store backward log-weights per time step
    backward_log_weights <- vector("list", n_times)
    backward_x <- vector("list", n_times)
    backward_y <- vector("list", n_times)

    # Last time step — weight by likelihood
    dets_t <- fish_dt[time == time_steps_vec[n_times] & detected == 1]
    detecting_ids <- unique(dets_t$station_id)
    b_lik <- pf_calculate_likelihoods(bpx, bpy, detecting_ids, stn, de_lookup)
    b_weights <- b_lik / sum(b_lik)
    backward_log_weights[[n_times]] <- log(pmax(b_weights, 1e-300))
    backward_x[[n_times]] <- bpx
    backward_y[[n_times]] <- bpy

    # Backward loop
    for (t in (n_times - 1):1) {
      if (inherits(time_steps_vec[t], "POSIXct")) {
        time_diff_sec <- as.numeric(difftime(time_steps_vec[t + 1], time_steps_vec[t], units = "secs"))
      } else {
        time_diff_sec <- as.numeric(time_steps_vec[t + 1] - time_steps_vec[t])
      }
      time_diff_sec <- max(time_diff_sec, 1)

      # Move backward (same CRW, different headings naturally explore backward)
      moved <- pf_move_particles(bpx, bpy, b_headings, time_diff_sec,
                                  step_length_mean, step_length_sd, kappa,
                                  time_step, raster, raster_ext)
      bpx <- moved$x; bpy <- moved$y; b_headings <- moved$heading

      # Weight by likelihood at this time step
      dets_t <- fish_dt[time == time_steps_vec[t] & detected == 1]
      detecting_ids <- unique(dets_t$station_id)
      b_lik <- pf_calculate_likelihoods(bpx, bpy, detecting_ids, stn, de_lookup)

      b_weights <- b_weights * b_lik
      w_sum <- sum(b_weights)
      if (w_sum > 0 && is.finite(w_sum)) { b_weights <- b_weights / w_sum
      } else { b_weights <- rep(1 / n_bp, n_bp) }

      # Resample if needed
      ess <- 1 / sum(b_weights^2)
      if (ess < n_bp * ess_threshold) {
        idx <- pf_systematic_resample(b_weights, n_bp)
        bpx <- bpx[idx]; bpy <- bpy[idx]; b_headings <- b_headings[idx]
        b_weights <- rep(1 / n_bp, n_bp)
      }

      backward_log_weights[[t]] <- log(pmax(b_weights, 1e-300))
      backward_x[[t]] <- bpx
      backward_y[[t]] <- bpy
    }

    # --- Combine forward and backward ---
    smoothed_positions <- data.table::data.table(
      fish_id = rep(current_fish, n_times), time = time_steps_vec,
      x_mean = numeric(n_times), y_mean = numeric(n_times),
      x_sd = numeric(n_times), y_sd = numeric(n_times),
      ess = numeric(n_times), n_detecting = integer(n_times),
      resampled = logical(n_times)
    )
    smoothed_particle_list <- if (return_particles) vector("list", n_times) else NULL

    for (t in seq_len(n_times)) {
      # Forward particles and weights at time t
      fwd_t <- fish_fwd[time == time_steps_vec[t]]
      fwd_x <- fwd_t$x; fwd_y <- fwd_t$y; fwd_w <- fwd_t$weight

      # Backward particles and weights at time t
      bwd_x <- backward_x[[t]]; bwd_y <- backward_y[[t]]
      bwd_log_w <- backward_log_weights[[t]]

      # Combine: pool forward and backward particles
      all_x <- c(fwd_x, bwd_x)
      all_y <- c(fwd_y, bwd_y)
      # Combined log-weight = forward log-weight + backward log-weight (for matched)
      # Since forward and backward have different particles, combine by pooling
      fwd_log_w <- log(pmax(fwd_w, 1e-300))
      all_log_w <- c(fwd_log_w, bwd_log_w)
      all_log_w <- all_log_w - max(all_log_w)
      all_w <- exp(all_log_w)
      all_w <- all_w / sum(all_w)

      xm <- sum(all_x * all_w)
      ym <- sum(all_y * all_w)

      dets_t <- fish_dt[time == time_steps_vec[t] & detected == 1]

      ti <- as.integer(t)
      data.table::set(smoothed_positions, ti, "x_mean", xm)
      data.table::set(smoothed_positions, ti, "y_mean", ym)
      data.table::set(smoothed_positions, ti, "x_sd", sqrt(sum(all_w * (all_x - xm)^2)))
      data.table::set(smoothed_positions, ti, "y_sd", sqrt(sum(all_w * (all_y - ym)^2)))
      data.table::set(smoothed_positions, ti, "ess", 1 / sum(all_w^2))
      data.table::set(smoothed_positions, ti, "n_detecting", nrow(dets_t))
      data.table::set(smoothed_positions, ti, "resampled", FALSE)

      if (return_particles) {
        smoothed_particle_list[[t]] <- data.table::data.table(
          fish_id = current_fish, time = time_steps_vec[t],
          particle = seq_along(all_x), x = all_x, y = all_y, weight = all_w)
      }
    }

    all_positions[[fi]] <- smoothed_positions
    if (return_particles) all_particles_out[[fi]] <- data.table::rbindlist(smoothed_particle_list)
  }

  if (verbose && interactive()) cat("\n")
  if (verbose) cat("Smoothing complete for", n_fish, "fish.\n")

  positions <- data.table::rbindlist(all_positions)
  result <- list(positions = as.data.frame(positions))
  if (return_particles) result$particles <- as.data.frame(data.table::rbindlist(all_particles_out))
  return(result)
}


# ============================================================================
# Utilization distribution
# ============================================================================

#' Compute Utilization Distribution from Particle Filter Results
#'
#' Converts weighted particle clouds into probability-of-use rasters and
#' extracts home range contours at specified levels.
#'
#' @param pf_results Output from \code{particle_filter_positioning()} or
#'   \code{pf_smooth()} with \code{return_particles = TRUE}.
#' @param raster Study area RasterLayer used as the spatial grid template.
#' @param contour_levels Numeric vector of contour levels for home range
#'   extraction (default \code{c(0.5, 0.95)}).
#' @param by_time Logical. Compute UD per time step (TRUE) or aggregate across
#'   all time steps per fish (FALSE, default).
#' @param verbose Logical. Print progress (default TRUE).
#'
#' @return A list with components:
#'   \item{home_ranges}{Data frame with area estimates per fish/contour level}
#'   \item{ud_rasters}{Named list of raster objects per fish (probability-of-use)}
#'
#' @export
pf_utilization_distribution <- function(pf_results, raster,
                                         contour_levels = c(0.5, 0.95),
                                         by_time = FALSE,
                                         verbose = TRUE) {

  if (is.null(pf_results$particles)) {
    stop("pf_results must include particles (use return_particles = TRUE)")
  }

  particles <- data.table::as.data.table(pf_results$particles)
  fish_ids <- unique(particles$fish_id)

  # Get cell resolution for area calculation
  cell_area <- raster::res(raster)[1] * raster::res(raster)[2]

  all_hr <- list()
  ud_rasters <- list()

  for (fi in seq_along(fish_ids)) {
    current_fish <- fish_ids[fi]
    if (verbose) cat("Computing UD for fish", current_fish, "\n")

    if (by_time) {
      fish_p <- particles[fish_id == current_fish]
      time_vals <- unique(fish_p$time)
      for (tv in time_vals) {
        tp <- fish_p[time == tv]
        ud_result <- pf_compute_single_ud(tp$x, tp$y, tp$weight, raster,
                                            cell_area, contour_levels)
        for (cl in contour_levels) {
          all_hr[[length(all_hr) + 1]] <- data.frame(
            fish_id = current_fish, time = tv,
            contour_level = cl,
            area_m2 = ud_result$areas[[as.character(cl)]])
        }
      }
    } else {
      fish_p <- particles[fish_id == current_fish]
      ud_result <- pf_compute_single_ud(fish_p$x, fish_p$y, fish_p$weight,
                                          raster, cell_area, contour_levels)
      ud_rasters[[as.character(current_fish)]] <- ud_result$ud_raster
      for (cl in contour_levels) {
        all_hr[[length(all_hr) + 1]] <- data.frame(
          fish_id = current_fish, time = NA,
          contour_level = cl,
          area_m2 = ud_result$areas[[as.character(cl)]])
      }
    }
  }

  home_ranges <- do.call(rbind, all_hr)
  home_ranges$area_ha <- home_ranges$area_m2 / 10000

  list(home_ranges = home_ranges, ud_rasters = ud_rasters)
}

#' Compute single UD from weighted particles
#' @keywords internal
pf_compute_single_ud <- function(x, y, weights, raster, cell_area, contour_levels) {
  # Assign particles to raster cells
  cell_ids <- raster::cellFromXY(raster, cbind(x, y))
  valid <- !is.na(cell_ids)
  cell_ids <- cell_ids[valid]
  w <- weights[valid]

  # Aggregate weights per cell
  dt <- data.table::data.table(cell = cell_ids, w = w)
  cell_sums <- dt[, .(prob = sum(w)), by = cell]

  # Normalize to sum to 1
  cell_sums[, prob := prob / sum(prob)]

  # Create UD raster
  ud_raster <- raster
  raster::values(ud_raster) <- NA
  ud_raster[cell_sums$cell] <- cell_sums$prob

  # Extract home range areas at each contour level
  areas <- list()
  # Sort cells by probability (descending) for isopleth extraction
  cell_sums <- cell_sums[order(-prob)]
  cell_sums[, cum_prob := cumsum(prob)]

  for (cl in contour_levels) {
    # Home range = smallest area containing cl proportion of total probability
    n_cells <- sum(cell_sums$cum_prob <= cl) + 1
    n_cells <- min(n_cells, nrow(cell_sums))
    areas[[as.character(cl)]] <- n_cells * cell_area
  }

  list(ud_raster = ud_raster, areas = areas)
}
