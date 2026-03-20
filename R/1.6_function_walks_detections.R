# Helper functions for timezone handling
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

# Helper to preserve original timezone info
store_timezone_info <- function(datetime_col) {
  tz <- attr(datetime_col, "tzone")
  if (is.null(tz) || tz == "") {
    return("UTC")
  }
  return(tz)
}

# Helper to restore timezone for output
restore_timezone <- function(datetime_col, original_tz) {
  if (is.null(original_tz) || original_tz == "") {
    original_tz <- "UTC"
  }
  return(as.POSIXct(datetime_col, tz = original_tz))
}

# Helper function for flexible date parsing
parse_date_flexible <- function(date_strings) {
  # Common date formats to try
  formats <- c("%Y-%m-%d", "%m/%d/%Y", "%d/%m/%Y", "%Y/%m/%d", 
               "%m-%d-%Y", "%d-%m-%Y", "%Y.%m.%d", "%m.%d.%Y")
  
  parsed_dates <- as.Date(rep(NA, length(date_strings)))
  
  for (format in formats) {
    missing_indices <- is.na(parsed_dates)
    if (!any(missing_indices)) break
    
    try({
      parsed_dates[missing_indices] <- as.Date(date_strings[missing_indices], format = format)
    }, silent = TRUE)
  }
  
  if (any(is.na(parsed_dates))) {
    warning(paste("Could not parse", sum(is.na(parsed_dates)), "dates out of", length(date_strings)))
  }
  
  return(parsed_dates)
}

# Helper function to calculate day length (photoperiod)
calculate_daylight_hours <- function(date, latitude = 45) {
  # Convert date to day of year
  day_of_year <- as.numeric(format(as.Date(date), "%j"))
  
  # Solar declination angle
  declination <- 23.45 * sin((360 * (284 + day_of_year) / 365) * pi / 180)
  
  # Hour angle for sunrise/sunset
  lat_rad <- latitude * pi / 180
  decl_rad <- declination * pi / 180
  
  hour_angle <- acos(-tan(lat_rad) * tan(decl_rad))
  
  # Daylight hours
  daylight_hours <- 2 * hour_angle * 12 / pi
  
  return(daylight_hours)
}

# Helper function to check if daylight is increasing (spring cue)
is_daylight_increasing <- function(date, latitude = 45, window_days = 7) {
  current_daylight <- calculate_daylight_hours(date, latitude)
  past_daylight <- calculate_daylight_hours(date - window_days, latitude)
  
  return(current_daylight > past_daylight)
}

# Helper function to calculate spawning probability based on temperature and photoperiod
calculate_spawning_probability <- function(temperature, date, species_params, latitude = 45) {
  # Check photoperiod condition
  if (!is_daylight_increasing(date, latitude)) {
    return(list(probability = 0, phase = "none"))
  }
  
  # Extract spawning temperature parameters
  min_temp <- species_params$spawning_min_temp
  opt_min <- species_params$spawning_optimal_temp_min
  opt_max <- species_params$spawning_optimal_temp_max
  max_temp <- species_params$spawning_max_temp
  
  if (temperature < min_temp || temperature > max_temp) {
    return(list(probability = 0, phase = "none"))
  }
  
  # Determine spawning phase and probability
  if (temperature >= min_temp && temperature < opt_min) {
    # Pre-spawning search phase
    prob <- (temperature - min_temp) / (opt_min - min_temp) * 0.7  # Max 0.7 probability in search phase
    return(list(probability = prob, phase = "search"))
  } else if (temperature >= opt_min && temperature <= opt_max) {
    # Peak spawning phase
    return(list(probability = 1.0, phase = "spawning"))
  } else if (temperature > opt_max && temperature <= max_temp) {
    # Post-spawning phase - return to normal behavior
    return(list(probability = 0, phase = "none"))
  }
  
  return(list(probability = 0, phase = "none"))
}

# Helper function to get depth bias for spawning behavior
get_depth_bias <- function(current_depth, spawning_phase, species_params) {
  if (spawning_phase == "none") {
    return(1.0)  # No bias
  }
  
  # Extract depth preferences
  preferred_min <- species_params$spawning_depth_preferred_min
  preferred_max <- species_params$spawning_depth_preferred_max
  acceptable_min <- species_params$spawning_depth_min
  acceptable_max <- species_params$spawning_depth_max
  
  # Check if in preferred range
  if (current_depth >= preferred_min && current_depth <= preferred_max) {
    return(2.0)  # Strong bias to stay in preferred depths
  }
  
  # Check if in acceptable range
  if (current_depth >= acceptable_min && current_depth <= acceptable_max) {
    return(1.5)  # Moderate bias to stay in acceptable depths
  }
  
  # Outside acceptable range - bias toward spawning depths
  return(0.3)  # Strong bias to move toward spawning areas
}

#' Simulate fish tracks with correlated random walks and detection events (ROBUST)
#'
#' Automatically detects column names and optimizes for speed.
#'
#' @param raster A RasterLayer object defining the study area boundaries.
#' @param station_distances Data frame with receiver detection probabilities.
#' @param n_paths Integer. Number of fish paths to simulate. Default is 1.
#' @param n_steps Integer. Number of steps per path. Default is 100.
#' @param step_length_mean Numeric. Mean step length in map units. Default is 50.
#' @param step_length_sd Numeric. Standard deviation of step length. Default is 20.
#' @param turning_angle_mean Numeric. Mean turning angle in degrees. Default is 0.
#' @param turning_angle_sd Numeric. Standard deviation of turning angle in degrees. Default is 45.
#' @param time_step Numeric. Time between steps in seconds. Default is 60.
#' @param start_locations Matrix or data frame with x,y coordinates for starting locations. Default is NULL.
#' @param start_time POSIXct object for simulation start time. Default is July 1, 2025 08:00:00 EST.
#' @param seed Numeric. Random seed for reproducible results. Default is NULL.
#' @param station_info Data frame with station deployment information including start/end dates. Default is NULL.
#' @param temporal_info Data frame with daily environmental conditions. Default is NULL.
#' @param de_model Model object for predicting detection efficiency based on temporal conditions. Default is NULL.
#' @param spawning_behavior Logical. Whether to enable spawning behavior modifications. Default is FALSE.
#' @param depth_state_bias Logical. Whether to enable depth-dependent behavioural
#'   state transition bias. When TRUE, fish in shallow water are more likely to
#'   transition to search state, and fish in deep water are more likely to transition
#'   to rest state. Depth thresholds are defined in species_behavioral_params.
#'   Default is FALSE.
#' @param species Character. Species name for preset movement parameters ("Walleye", "Smallmouth Bass", "Muskellunge"). Default is NULL.
#' @param fish_size_cm Numeric. Fish length in centimeters for size-scaled movement parameters. Default is NULL.
#' @param species_params Data frame with custom species movement parameters. Default is NULL.
#' @param temperature_data Data frame with daily temperature data. Default is NULL.
#' @param behavioral_states Logical. Whether to use 3-state behavioral model (cruise/search/rest). Default is TRUE when species specified.
#' @param goal_locations Matrix or data frame with x,y goal coordinates, one row per path. When provided,
#'   fish navigate toward the goal using a biased correlated random walk and stop upon arrival.
#'   Requires \code{start_locations} to be specified. Default is NULL.
#' @param goal_bias Numeric between 0 and 1. Strength of directional bias toward goal.
#'   0 = pure CRW (no bias), 1 = beeline toward goal. Default is 0.5.
#' @param goal_tolerance Numeric. Distance (map units) at which fish is considered "arrived" at goal.
#'   Defaults to \code{step_length_mean} if NULL.
#' @param include_barriers Logical. Whether to apply barrier masking to prevent detections through land obstacles. Default is FALSE.
#'   When TRUE, detection efficiency is set to 0 for any receiver where the line-of-sight crosses a barrier.
#'   Requires \code{crosses_barrier} column in \code{station_distances} (generated by \code{calculate_station_distances()} with \code{barrier_raster} parameter).
#'
#' @return A list containing tracks and station_detections data frames.
#'
#' @details
#' When \code{include_barriers = TRUE}, the function masks detections at receivers where the direct path
#' from the fish location crosses a land barrier. This prevents physically impossible detections through
#' islands, shorelines, or other obstacles. The barrier information must be pre-computed in the
#' \code{station_distances} data frame using \code{calculate_station_distances()} with a barrier raster.
#'
#' @export
simulate_fish_tracks <- function(raster, station_distances, n_paths = 1, n_steps = 100,
                                 step_length_mean = NULL, step_length_sd = NULL,
                                 turning_angle_mean = NULL, turning_angle_sd = NULL,
                                 time_step = 60, start_locations = NULL,
                                 start_time = as.POSIXct("2025-07-01 08:00:00", tz = "America/Toronto"),
                                 seed = NULL, station_info = NULL, temporal_info = NULL, de_model = NULL,
                                 species = NULL, fish_size_cm = NULL, species_params = NULL,
                                 temperature_data = NULL, behavioral_states = NULL, spawning_behavior = FALSE,
                                 depth_state_bias = FALSE,
                                 goal_locations = NULL, goal_bias = 0.5, goal_tolerance = NULL,
                                 include_barriers = FALSE) {

  if (!is.null(seed)) set.seed(seed)
  
  # Handle species-specific movement parameters and behavioral states
  movement_params_resolved <- FALSE
  use_behavioral_states <- FALSE
  species_row <- NULL
  
  if (!is.null(species)) {
    # Load species behavioral parameters
    if (is.null(species_params)) {
      data("species_behavioral_params", envir = environment())
      species_params <- species_behavioral_params
    }
    
    # Find species in parameters
    species_row <- species_params[species_params$species == species, ]
    
    if (nrow(species_row) == 0) {
      stop(paste("Species", species, "not found in species_params. Available species:", 
                paste(unique(species_params$species), collapse = ", ")))
    }
    
    # Calculate size-adjusted parameters
    if (is.null(fish_size_cm)) {
      fish_size_cm <- species_row$typical_size_cm
      cat("Using typical size for", species, ":", fish_size_cm, "cm\n")
    }
    
    # Determine if using behavioral states
    if (is.null(behavioral_states)) {
      use_behavioral_states <- TRUE  # Default to behavioral states when species specified
    } else {
      use_behavioral_states <- behavioral_states
    }
    
    if (use_behavioral_states) {
      cat("Using 3-state behavioral model for", species, "(size:", fish_size_cm, "cm)\n")
      cat("  States: cruise, search, rest\n")
      cat("  Temperature-dependent state switching enabled\n")
      if (depth_state_bias) cat("  Depth-dependent state bias enabled\n")
      cat("  Movement style:", species_row$movement_description, "\n")
    } else {
      # Use simple speed-based parameters (backward compatibility)
      speed_mean_base <- (species_row$cruise_speed_mean_ms + species_row$search_speed_mean_ms) / 2
      speed_sd_base <- (species_row$cruise_speed_sd_ms + species_row$search_speed_sd_ms) / 2
      
      # Calculate size-adjusted speed
      size_speed_mean <- speed_mean_base + (species_row$speed_mean_size_scalar * fish_size_cm)
      size_speed_sd <- speed_sd_base + (species_row$speed_sd_size_scalar * fish_size_cm)
      
      # Convert to step lengths (fallback for simple mode)
      if (is.null(step_length_mean)) {
        step_length_mean <- size_speed_mean * time_step
      }
      if (is.null(step_length_sd)) {
        step_length_sd <- size_speed_sd * time_step
      }
      if (is.null(turning_angle_mean)) {
        turning_angle_mean <- 0  # Default
      }
      if (is.null(turning_angle_sd)) {
        turning_angle_sd <- (species_row$cruise_turn_sd + species_row$search_turn_sd) / 2
      }
      
      cat("Using", species, "simple movement parameters (size:", fish_size_cm, "cm):\n")
      cat("  Step length mean:", round(step_length_mean, 1), "m\n")
      cat("  Step length SD:", round(step_length_sd, 1), "m\n")
      cat("  Turning angle SD:", turning_angle_sd, "degrees\n")
    }
    
    movement_params_resolved <- TRUE
  }

  # Validate depth_state_bias columns
  if (depth_state_bias) {
    required_cols <- c("search_depth_preferred_max", "search_depth_max",
                       "rest_depth_min", "rest_depth_preferred_min", "depth_bias_strength")
    missing <- setdiff(required_cols, names(species_row))
    if (length(missing) > 0) {
      stop("depth_state_bias requires columns in species_params: ", paste(missing, collapse = ", "))
    }
  }

  # Set default values if not specified and no species used
  if (!movement_params_resolved) {
    if (is.null(step_length_mean)) step_length_mean <- 50
    if (is.null(step_length_sd)) step_length_sd <- 20
    if (is.null(turning_angle_mean)) turning_angle_mean <- 0
    if (is.null(turning_angle_sd)) turning_angle_sd <- 45
  }

  # Validate inputs
  if (is.null(station_distances)) {
    stop("station_distances is required for receiver-specific detection simulation")
  }
  
  # Check if using temporal DE prediction or pre-computed DE values
  use_temporal_de <- !is.null(station_info) && !is.null(temporal_info) && !is.null(de_model)

  if (!use_temporal_de && !"DE_pred" %in% names(station_distances)) {
    stop("station_distances must contain 'DE_pred' column with detection probabilities, or provide station_info, temporal_info, and de_model for temporal DE prediction")
  }

  # Validate barrier data if barrier masking is requested
  if (include_barriers) {
    if (!"crosses_barrier" %in% names(station_distances)) {
      stop("include_barriers = TRUE requires 'crosses_barrier' column in station_distances.\n",
           "Run calculate_station_distances() with barrier_raster parameter to generate this column.")
    }
    cat("Barrier masking enabled: DE will be set to 0 where crosses_barrier = TRUE\n")
  }

  # Validate goal-directed movement parameters
  if (!is.null(goal_locations)) {
    if (is.null(start_locations)) {
      stop("goal_locations requires start_locations to be specified — goal-directed walks need a defined start")
    }
    goal_locations <- as.matrix(goal_locations)
    if (ncol(goal_locations) < 2) {
      stop("goal_locations must have at least 2 columns (x, y)")
    }
    if (nrow(goal_locations) != n_paths) {
      stop("goal_locations must have one row per path (", n_paths, " rows expected, got ", nrow(goal_locations), ")")
    }
    if (goal_bias < 0 || goal_bias > 1) {
      stop("goal_bias must be between 0 and 1")
    }
    if (is.null(goal_tolerance)) {
      goal_tolerance <- step_length_mean
    }
    cat("Goal-directed movement enabled (bias:", goal_bias, ", tolerance:", goal_tolerance, "m)\n")
  }

  # Load temporal data if using temporal DE prediction
  if (use_temporal_de) {
    cat("Using temporal DE prediction with provided model...\n")
    
    # Load data files if paths provided as character strings
    if (is.character(station_info)) {
      station_info <- utils::read.csv(station_info, stringsAsFactors = FALSE)
    }
    if (is.character(temporal_info)) {
      temporal_info <- utils::read.csv(temporal_info, stringsAsFactors = FALSE)
    }
    
    # Helper function to parse dates in multiple formats
    parse_date_flexible <- function(date_vector) {
      # Common date formats to try
      formats <- c(
        "%Y-%m-%d",     # 2025-01-15 (ISO standard)
        "%m/%d/%Y",     # 1/15/2025 (US Excel format)
        "%d/%m/%Y",     # 15/1/2025 (European Excel format)
        "%m-%d-%Y",     # 1-15-2025
        "%d-%m-%Y",     # 15-1-2025
        "%Y/%m/%d",     # 2025/1/15
        "%m/%d/%y",     # 1/15/25 (2-digit year)
        "%d/%m/%y"      # 15/1/25 (2-digit year)
      )
      
      # Try each format
      for (fmt in formats) {
        result <- try(as.Date(as.character(date_vector), format = fmt), silent = TRUE)
        if (!inherits(result, "try-error") && !any(is.na(result))) {
          return(result)
        }
      }
      
      # If all specific formats fail, try automatic parsing
      result <- try(as.Date(as.character(date_vector)), silent = TRUE)
      if (!inherits(result, "try-error")) {
        return(result)
      }
      
      # If everything fails, warn and return NAs
      warning(paste("Could not parse dates:", paste(unique(date_vector[1:min(5, length(date_vector))]), collapse = ", "), 
                   "... Please use YYYY-MM-DD format or common Excel formats."))
      return(as.Date(rep(NA, length(date_vector))))
    }
    
    # Convert date columns to Date format with flexible parsing
    station_info$start_date <- parse_date_flexible(station_info$start_date)
    station_info$end_date <- parse_date_flexible(station_info$end_date)
    temporal_info$date <- parse_date_flexible(temporal_info$date)
  }

  # Validate start_time
  if (!inherits(start_time, "POSIXct")) {
    stop("start_time must be a POSIXct object")
  }
  
  # Store original timezone and standardize to UTC for internal processing
  original_timezone <- store_timezone_info(start_time)
  start_time_utc <- standardize_datetime(start_time, "UTC")

  # Load temperature data if using behavioral states
  if (use_behavioral_states && !is.null(species_row)) {
    if (is.null(temperature_data)) {
      data("daily_temperature", envir = environment())
      temperature_data <- daily_temperature
    }
    
    # Load data files if paths provided as character strings
    if (is.character(temperature_data)) {
      temperature_data <- utils::read.csv(temperature_data, stringsAsFactors = FALSE)
    }
    
    # Convert date column to Date format with flexible parsing
    if (!inherits(temperature_data$date, "Date")) {
      temperature_data$date <- parse_date_flexible(temperature_data$date)
    }
    
    cat("Temperature data loaded:", nrow(temperature_data), "days\n")
    cat("Temperature range:", round(min(temperature_data$water_temp_c), 1), "to", 
        round(max(temperature_data$water_temp_c), 1), "\u00B0C\n")
  }
  
  # Load required library for fast nearest neighbor search
  if (!requireNamespace("FNN", quietly = TRUE)) {
    stop("Package 'FNN' is required for fast spatial lookup. Install with: install.packages('FNN')")
  }

  cat("Preprocessing spatial lookup table...\n")

  # Prepare daily DE predictions if using temporal model
  daily_de_lookup <- NULL
  if (use_temporal_de) {
    cat("Pre-computing daily DE predictions...\n")
    
    # Get simulation date range
    simulation_end_time <- start_time + (n_steps * n_paths * time_step)
    unique_dates <- seq(from = as.Date(start_time), to = as.Date(simulation_end_time), by = "day")
    
    daily_de_lookup <- list()
    
    for (date in unique_dates) {
      date_obj <- as.Date(date, origin = "1970-01-01")
      
      # Get active stations for this date
      active_stations <- station_info %>%
        dplyr::filter(start_date <= date_obj & end_date >= date_obj)
      
      # Get environmental conditions for this date
      daily_conditions <- temporal_info %>%
        dplyr::filter(date == date_obj)
      
      if (nrow(active_stations) > 0 && nrow(daily_conditions) > 0) {
        # Create prediction data for all active station-cell combinations
        prediction_data <- station_distances %>%
          dplyr::filter(station_no %in% active_stations$station_id) %>%
          dplyr::mutate(
            dist_m = cost_distance,
            depth_m = abs(raster_value)
          ) %>%
          dplyr::cross_join(daily_conditions %>% dplyr::select(-date))
        
        # Check if model uses these variables (excluding response variable)
        model_formula <- stats::formula(de_model)
        # Get predictor variables only (right side of formula)
        rhs_formula <- model_formula[[3]]  # Right-hand side of formula
        model_vars <- all.vars(rhs_formula)  # Only predictor variables
        available_vars <- names(prediction_data)
        missing_vars <- setdiff(model_vars, available_vars)
        
        if (length(missing_vars) > 0) {
          warning(paste("Missing variables for DE model:", paste(missing_vars, collapse = ", ")))
        }
        
        # Predict DE for this day
        prediction_data$DE_pred <- tryCatch({
          stats::predict(de_model, newdata = prediction_data, type = "response")
        }, error = function(e) {
          warning(paste("DE prediction failed for date", date_obj, ":", e$message))
          rep(0.5, nrow(prediction_data))  # Default fallback
        })
        
        # Store in lookup
        daily_de_lookup[[as.character(date_obj)]] <- prediction_data %>%
          dplyr::select(cell_id, station_no, DE_pred)
      }
    }
    
    cat("Completed daily DE predictions for", length(daily_de_lookup), "days\n")
  }

  # Create spatial lookup table
  if (use_temporal_de) {
    # Select columns including barriers if available
    if (include_barriers) {
      station_lookup <- station_distances %>%
        dplyr::select(cell_id, x, y, station_no, cost_distance, crosses_barrier) %>%
        dplyr::arrange(cell_id, station_no)
    } else {
      station_lookup <- station_distances %>%
        dplyr::select(cell_id, x, y, station_no, cost_distance) %>%
        dplyr::arrange(cell_id, station_no)
    }
  } else {
    # Select columns including barriers if available
    if (include_barriers) {
      station_lookup <- station_distances %>%
        dplyr::select(cell_id, x, y, station_no, cost_distance, DE_pred, crosses_barrier) %>%
        dplyr::mutate(DE_pred = ifelse(crosses_barrier, 0, DE_pred)) %>%
        dplyr::arrange(cell_id, station_no)
    } else {
      station_lookup <- station_distances %>%
        dplyr::select(cell_id, x, y, station_no, cost_distance, DE_pred) %>%
        dplyr::arrange(cell_id, station_no)
    }
  }

  # Create unique cell coordinates for fast lookup
  unique_cells <- station_distances %>%
    dplyr::select(cell_id, x, y) %>%
    dplyr::distinct()

  lookup_coords <- as.matrix(unique_cells[, c("x", "y")])

  # Convert parameters
  turning_angle_sd_rad <- turning_angle_sd * pi / 180
  turning_angle_mean_rad <- turning_angle_mean * pi / 180

  # Get raster properties
  raster_extent <- raster::extent(raster)
  valid_cells <- which(!is.na(raster::values(raster)))
  valid_coords <- raster::xyFromCell(raster, valid_cells)

  # Initialize behavioral state tracking if using behavioral states
  if (use_behavioral_states && !is.null(species_row)) {
    # Initialize state for all paths (start in cruise state)
    current_states <- rep("cruise", n_paths)
    state_history <- list()  # Track state changes for analysis
  }
  
  # Initialize results
  all_tracks <- list()
  all_station_detections <- list()

  for (path_id in 1:n_paths) {
    cat("Generating path", path_id, "of", n_paths, "\n")

    # Set starting location
    if (is.null(start_locations)) {
      start_idx <- sample(nrow(valid_coords), 1)
      start_x <- valid_coords[start_idx, 1]
      start_y <- valid_coords[start_idx, 2]
    } else {
      start_x <- start_locations[path_id, 1]
      start_y <- start_locations[path_id, 2]
    }

    # Initialize path tracking
    track_x <- rep(NA, n_steps + 1)
    track_y <- rep(NA, n_steps + 1)

    # Create time sequences for both seconds (for backward compatibility) and POSIX
    track_time_seconds <- seq(0, n_steps * time_step, by = time_step)
    track_time_posix <- start_time_utc + track_time_seconds

    # Set starting position and direction
    track_x[1] <- start_x
    track_y[1] <- start_y

    # Set initial bearing — toward goal if goal-directed, else random
    if (!is.null(goal_locations)) {
      goal_x <- goal_locations[path_id, 1]
      goal_y <- goal_locations[path_id, 2]
      current_bearing <- atan2(goal_y - start_y, goal_x - start_x)
    } else {
      current_bearing <- runif(1, 0, 2 * pi)
    }

    # Initialize path-specific state tracking
    if (use_behavioral_states && !is.null(species_row)) {
      path_state_history <- rep("cruise", n_steps + 1)  # Track state for each step
    }
    
    # Initialize spawning phase tracking
    if (spawning_behavior && !is.null(species_row)) {
      path_spawning_phase <- rep("none", n_steps + 1)  # Track spawning phase for each step
    }
    
    # Generate correlated random walk with behavioral states
    for (step in 1:n_steps) {

      # Check if fish has arrived at goal
      if (!is.null(goal_locations)) {
        dist_to_goal <- sqrt((goal_x - track_x[step])^2 + (goal_y - track_y[step])^2)
        if (dist_to_goal <= goal_tolerance) {
          # Fish arrived — fill remaining steps with goal position
          track_x[(step + 1):(n_steps + 1)] <- goal_x
          track_y[(step + 1):(n_steps + 1)] <- goal_y
          cat("  Path", path_id, "arrived at goal on step", step, "\n")
          break
        }
      }

      # Handle behavioral state transitions and movement
      if (use_behavioral_states && !is.null(species_row)) {
        # Get current date for temperature lookup
        current_time <- track_time_posix[step + 1]
        current_date <- as.Date(current_time)
        
        # Get temperature for this date
        temp_row <- temperature_data[temperature_data$date == current_date, ]
        if (nrow(temp_row) > 0) {
          current_temp <- temp_row$water_temp_c[1]
        } else {
          current_temp <- species_row$optimal_temp_c  # Default if no temperature data
        }
        
        # Calculate spawning behavior if enabled
        spawning_probability <- 0
        spawning_phase <- "none"
        depth_bias <- 1.0
        
        if (depth_state_bias) {
          current_depth <- abs(raster::extract(raster, cbind(track_x[step], track_y[step])))
          if (is.na(current_depth)) current_depth <- 5
        }

        if (spawning_behavior) {
          # Get current depth (simplified - use raster value at current position)
          if (!depth_state_bias) {
            current_depth <- abs(raster::extract(raster, cbind(track_x[step], track_y[step])))
            if (is.na(current_depth)) current_depth <- 5  # Default depth if NA
          }
          
          # Calculate spawning probability and phase
          spawn_result <- calculate_spawning_probability(current_temp, current_date, species_row)
          spawning_probability <- spawn_result$probability
          spawning_phase <- spawn_result$phase
          
          # Calculate depth bias for movement
          depth_bias <- get_depth_bias(current_depth, spawning_phase, species_row)
          
          # Store spawning phase in history
          path_spawning_phase[step + 1] <- spawning_phase
        }
        
        # Calculate temperature-dependent activity factor
        if (current_temp < species_row$low_temp_threshold) {
          activity_multiplier <- species_row$cold_activity_factor
        } else if (current_temp > species_row$high_temp_threshold) {
          activity_multiplier <- species_row$hot_activity_factor
        } else {
          # Linear interpolation around optimal temperature
          temp_diff <- abs(current_temp - species_row$optimal_temp_c)
          activity_multiplier <- max(0.1, 1 - (temp_diff / species_row$temp_tolerance_c))
        }
        
        # State transitions with spawning behavior overlay
        current_state <- current_states[path_id]
        transition_prob <- runif(1)
        
        # Modify transition probabilities based on spawning behavior
        search_multiplier <- 1.0
        rest_multiplier <- 1.0
        
        if (spawning_behavior && spawning_probability > 0) {
          if (spawning_phase == "search") {
            # Pre-spawning: increase search behavior
            search_multiplier <- 1.0 + spawning_probability * 2.0  # Up to 3x search probability
          } else if (spawning_phase == "spawning") {
            # Peak spawning: increase rest/territory holding behavior
            rest_multiplier <- 1.0 + spawning_probability * 1.5  # Up to 2.5x rest probability
          }
        }

        if (depth_state_bias) {
          # Shallow water: boost search transitions
          if (current_depth <= species_row$search_depth_max) {
            if (current_depth <= species_row$search_depth_preferred_max) {
              search_multiplier <- search_multiplier * species_row$depth_bias_strength
            } else {
              frac <- (species_row$search_depth_max - current_depth) /
                      (species_row$search_depth_max - species_row$search_depth_preferred_max)
              search_multiplier <- search_multiplier * (1 + frac * (species_row$depth_bias_strength - 1))
            }
          }

          # Deep water: boost rest transitions
          if (current_depth >= species_row$rest_depth_min) {
            if (current_depth >= species_row$rest_depth_preferred_min) {
              rest_multiplier <- rest_multiplier * species_row$depth_bias_strength
            } else {
              frac <- (current_depth - species_row$rest_depth_min) /
                      (species_row$rest_depth_preferred_min - species_row$rest_depth_min)
              rest_multiplier <- rest_multiplier * (1 + frac * (species_row$depth_bias_strength - 1))
            }
          }
        }

        if (current_state == "cruise") {
          search_prob <- species_row$base_cruise_to_search * activity_multiplier * search_multiplier
          rest_prob <- species_row$base_cruise_to_rest * (2 - activity_multiplier) * rest_multiplier
          
          if (transition_prob < search_prob) {
            current_states[path_id] <- "search"
          } else if (transition_prob < (search_prob + rest_prob)) {
            current_states[path_id] <- "rest"
          }
        } else if (current_state == "search") {
          cruise_prob <- species_row$base_search_to_cruise * activity_multiplier
          rest_prob <- species_row$base_search_to_rest * (2 - activity_multiplier) * rest_multiplier
          
          if (transition_prob < cruise_prob) {
            current_states[path_id] <- "cruise"
          } else if (transition_prob < (cruise_prob + rest_prob)) {
            current_states[path_id] <- "rest"
          }
        } else { # rest state
          cruise_prob <- species_row$base_rest_to_cruise * activity_multiplier
          search_prob <- species_row$base_rest_to_search * activity_multiplier * search_multiplier
          
          if (transition_prob < cruise_prob) {
            current_states[path_id] <- "cruise"
          } else if (transition_prob < (cruise_prob + search_prob)) {
            current_states[path_id] <- "search"
          }
        }
        
        # Store current state in path history
        path_state_history[step + 1] <- current_states[path_id]
        
        # Get state-specific movement parameters
        current_state <- current_states[path_id]
        if (current_state == "cruise") {
          speed_mean <- species_row$cruise_speed_mean_ms + (species_row$speed_mean_size_scalar * fish_size_cm)
          speed_sd <- species_row$cruise_speed_sd_ms + (species_row$speed_sd_size_scalar * fish_size_cm)
          turn_sd <- species_row$cruise_turn_sd
        } else if (current_state == "search") {
          speed_mean <- species_row$search_speed_mean_ms + (species_row$speed_mean_size_scalar * fish_size_cm)
          speed_sd <- species_row$search_speed_sd_ms + (species_row$speed_sd_size_scalar * fish_size_cm)
          turn_sd <- species_row$search_turn_sd
        } else { # rest state
          speed_mean <- species_row$rest_speed_mean_ms + (species_row$speed_mean_size_scalar * fish_size_cm)
          speed_sd <- species_row$rest_speed_sd_ms + (species_row$speed_sd_size_scalar * fish_size_cm)
          turn_sd <- species_row$rest_turn_sd
        }
        
        # Apply temperature effect on speed
        temp_adjusted_speed_mean <- speed_mean * (1 + species_row$temp_speed_scalar * (activity_multiplier - 0.5))
        
        # Generate movement parameters for this step
        swimming_speed <- max(stats::rnorm(1, temp_adjusted_speed_mean, speed_sd), 0.01)
        
        # Modify movement based on spawning phase
        if (spawning_behavior && spawning_phase == "spawning" && depth_bias >= 1.5) {
          # Fish is in suitable spawning habitat - greatly reduce movement
          swimming_speed <- swimming_speed * 0.1  # 90% reduction in speed
        } else if (spawning_behavior && spawning_phase == "search") {
          # Searching for spawning habitat - maintain or increase movement
          swimming_speed <- swimming_speed * 1.2  # 20% increase in search speed
        }
        
        step_length <- swimming_speed * time_step
        
        # Generate turning angle with state-specific variability and spawning bias
        turn_sd_rad <- turn_sd * pi / 180
        
        # Apply depth bias to turning angle if spawning behavior is active
        if (spawning_behavior && spawning_phase == "spawning" && depth_bias >= 1.5) {
          # Fish is in spawning habitat - reduce turning to stay in area
          turn_sd_rad <- turn_sd_rad * 0.3  # Much less random turning
        } else if (spawning_behavior && spawning_probability > 0 && depth_bias < 1.0) {
          # Fish is outside preferred spawning depth - bias movement toward preferred depths
          # This is a simplified implementation - in reality would need gradient calculation
          turn_sd_rad <- turn_sd_rad * depth_bias  # Reduce randomness when seeking depths
        }
        
        suppressWarnings({
          turning_angle <- circular::rvonmises(1, mu = 0, kappa = 1 / (turn_sd_rad^2))
        })
        
      } else {
        # Traditional fixed-parameter approach
        step_length <- max(stats::rnorm(1, step_length_mean, step_length_sd), 5)
        suppressWarnings({
          turning_angle <- circular::rvonmises(1, mu = turning_angle_mean_rad,
                                               kappa = 1 / (turning_angle_sd_rad^2))
        })
      }

      current_bearing <- current_bearing + turning_angle

      # Apply goal-directed bias via weighted circular mean
      if (!is.null(goal_locations)) {
        goal_bearing <- atan2(goal_y - track_y[step], goal_x - track_x[step])
        crw_weight <- 1 - goal_bias
        goal_weight <- goal_bias
        mean_x <- crw_weight * cos(current_bearing) + goal_weight * cos(goal_bearing)
        mean_y <- crw_weight * sin(current_bearing) + goal_weight * sin(goal_bearing)
        current_bearing <- atan2(mean_y, mean_x)
      }

      new_x <- track_x[step] + step_length * cos(current_bearing)
      new_y <- track_y[step] + step_length * sin(current_bearing)

      # Boundary checking
      if (new_x >= raster_extent@xmin && new_x <= raster_extent@xmax &&
          new_y >= raster_extent@ymin && new_y <= raster_extent@ymax) {
        cell_value <- raster::extract(raster, matrix(c(new_x, new_y), ncol = 2))
        if (!is.na(cell_value)) {
          track_x[step + 1] <- new_x
          track_y[step + 1] <- new_y
        } else {
          current_bearing <- current_bearing + pi
          new_x <- track_x[step] + step_length * cos(current_bearing)
          new_y <- track_y[step] + step_length * sin(current_bearing)
          track_x[step + 1] <- new_x
          track_y[step + 1] <- new_y
        }
      } else {
        current_bearing <- current_bearing + pi
        new_x <- track_x[step] + step_length * cos(current_bearing)
        new_y <- track_y[step] + step_length * sin(current_bearing)
        track_x[step + 1] <- new_x
        track_y[step + 1] <- new_y
      }
    }

    # Create track dataframe with both time formats
    track_df <- data.frame(
      path_id = path_id,
      step = 0:n_steps,
      time_seconds = track_time_seconds,
      datetime = restore_timezone(track_time_posix, original_timezone),
      x = track_x,
      y = track_y
    ) %>%
      dplyr::filter(!is.na(x) & !is.na(y)) %>%
      dplyr::mutate(
        step_length = c(0, sqrt(diff(x)^2 + diff(y)^2)),
        bearing = c(0, atan2(diff(y), diff(x)) * 180 / pi)
      )
    
    # Add behavioral state information if using behavioral states
    if (use_behavioral_states && !is.null(species_row)) {
      # Add state information (trim to match track_df length after filtering)
      if (nrow(track_df) <= length(path_state_history)) {
        track_df$behavioral_state <- path_state_history[1:nrow(track_df)]
      } else {
        track_df$behavioral_state <- rep("cruise", nrow(track_df))  # Fallback
      }
    }
    
    # Add spawning phase information if using spawning behavior
    if (spawning_behavior && !is.null(species_row)) {
      if (exists("path_spawning_phase")) {
        # Add spawning phase information (trim to match track_df length after filtering)
        if (nrow(track_df) <= length(path_spawning_phase)) {
          track_df$spawning_phase <- path_spawning_phase[1:nrow(track_df)]
        } else {
          track_df$spawning_phase <- rep("none", nrow(track_df))  # Fallback
        }
      } else {
        # If path_spawning_phase doesn't exist, create default
        track_df$spawning_phase <- rep("none", nrow(track_df))
      }
    }
    
    # Add temperature information if using behavioral states with species
    if (use_behavioral_states && !is.null(species_row) && !is.null(temperature_data)) {
      track_df$water_temp_c <- NA
      for (i in 1:nrow(track_df)) {
        current_date <- as.Date(track_df$datetime[i])
        temp_row <- temperature_data[temperature_data$date == current_date, ]
        if (nrow(temp_row) > 0) {
          track_df$water_temp_c[i] <- temp_row$water_temp_c[1]
        }
      }
    }

    all_tracks[[path_id]] <- track_df

    # Batch process all points for this path
    if (nrow(track_df) > 0) {

      # Get all track coordinates
      track_coords <- as.matrix(track_df[, c("x", "y")])

      # Fast nearest neighbor lookup for ALL points at once
      closest_indices <- FNN::get.knnx(lookup_coords, track_coords, k = 1)$nn.index

      # Add cell IDs to track data
      track_df$cell_id <- unique_cells$cell_id[closest_indices]

      # Join with station data and get DE predictions
      if (use_temporal_de) {
        # Use temporal DE prediction
        detection_data_list <- list()
        
        for (i in 1:nrow(track_df)) {
          current_date <- as.character(as.Date(track_df$datetime[i]))
          
          if (current_date %in% names(daily_de_lookup)) {
            # Get DE predictions for this date
            daily_de <- daily_de_lookup[[current_date]]
            
            # Join track point with station data and daily DE
            point_data <- track_df[i, ] %>%
              dplyr::left_join(station_lookup %>% dplyr::select(-x, -y),
                               by = "cell_id", relationship = "many-to-many") %>%
              dplyr::left_join(daily_de, by = c("cell_id", "station_no"), relationship = "many-to-many") %>%
              dplyr::filter(!is.na(cost_distance) & !is.na(DE_pred))

            # Apply barrier masking if enabled
            if (include_barriers && "crosses_barrier" %in% names(point_data)) {
              point_data <- point_data %>%
                dplyr::mutate(DE_pred = ifelse(crosses_barrier, 0, DE_pred))
            }

            detection_data_list[[i]] <- point_data
          }
        }
        
        detection_data <- do.call(rbind, detection_data_list) %>%
          dplyr::rename(station_id = station_no,
                        distance_to_station = cost_distance,
                        detection_prob = DE_pred)
        
      } else {
        # Use pre-computed DE values
        detection_data <- track_df %>%
          dplyr::left_join(station_lookup %>% dplyr::select(-x, -y),
                           by = "cell_id", relationship = "many-to-many") %>%
          dplyr::filter(!is.na(cost_distance) & !is.na(DE_pred)) %>%
          dplyr::rename(station_id = station_no,
                        distance_to_station = cost_distance,
                        detection_prob = DE_pred)
      }

      # Vectorized detection simulation
      if (nrow(detection_data) > 0) {
        detection_data$detected <- stats::rbinom(nrow(detection_data), 1, detection_data$detection_prob)

        # Select final columns and aggregate by station-time-fish combination
        path_detections <- detection_data %>%
          dplyr::group_by(path_id, step, time_seconds, datetime, station_id) %>%
          dplyr::summarise(
            x = dplyr::first(x),  # Use first position for this time step
            y = dplyr::first(y),
            distance_to_station = mean(distance_to_station),  # Average distance
            detection_prob = mean(detection_prob),  # Average detection probability
            detected = max(detected),  # 1 if any detection occurred, 0 otherwise
            .groups = 'drop'
          )

        all_station_detections[[path_id]] <- path_detections
      }
    }
  }

  # Combine results
  final_tracks <- do.call(rbind, all_tracks)

  if (length(all_station_detections) > 0) {
    final_station_detections <- do.call(rbind, all_station_detections)
  } else {
    final_station_detections <- data.frame()
  }

  # Return results
  results <- list(
    tracks = final_tracks,
    station_detections = final_station_detections,
    parameters = list(
      n_paths = n_paths,
      n_steps = n_steps,
      step_length_mean = step_length_mean,
      step_length_sd = step_length_sd,
      turning_angle_mean = turning_angle_mean,
      turning_angle_sd = turning_angle_sd,
      time_step = time_step,
      start_time = start_time,
      species = species,
      fish_size_cm = fish_size_cm,
      behavioral_states = use_behavioral_states,
      temperature_dependent = use_behavioral_states && !is.null(temperature_data),
      spawning_behavior = spawning_behavior,
      depth_state_bias = depth_state_bias,
      goal_directed = !is.null(goal_locations),
      goal_bias = if (!is.null(goal_locations)) goal_bias else NULL,
      goal_tolerance = if (!is.null(goal_locations)) goal_tolerance else NULL,
      timezone = original_timezone
    )
  )

  cat("Simulation complete!\n")
  cat("Start time:", format(start_time, "%Y-%m-%d %H:%M:%S %Z"), "\n")
  cat("End time:", format(max(final_tracks$datetime), "%Y-%m-%d %H:%M:%S %Z"), "\n")

  return(results)
}


#' Plot simulated fish tracks with detection events
#'
#' Creates a visualization of simulated fish movement paths overlaid on a raster
#' background, showing detection events and receiver station performance.
#'
#' @param simulation_results List output from \code{\link{simulate_fish_tracks}}
#'   containing tracks and station_detections.
#' @param raster A RasterLayer object used as the background for plotting.
#' @param receiver_frame An sf object containing receiver station locations.
#'   Default is NULL.
#' @param show_detections Logical. Whether to display detection events on the plot.
#'   Default is TRUE.
#' @param path_alpha Numeric. Transparency level for path lines (0-1). Default is 0.7.
#' @param color_by Character. What to color the tracks by: "path_id" (default),
#'   "behavioral_state" (if behavioral states used), or "spawning_phase" (if spawning behavior used).
#'   Default is "path_id".
#' @param path_color Character. Optional single colour for all paths (e.g., "yellow").
#'   When specified, overrides color_by for path colouring. Default is NULL.
#' @param point_size Numeric. Size of points when plotting states/phases. Default is 1.5.
#' @param sample_rate Numeric. Fraction of points to plot (1 = all points, 0.1 = 10% of points). 
#'   Useful for long tracks to reduce overplotting. Default is 1.
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
                             show_detections = TRUE, path_alpha = 0.7, color_by = "path_id",
                             path_color = NULL, point_size = 1.5, sample_rate = 1) {
  # Convert raster to dataframe for plotting
  raster_df <- raster::as.data.frame(raster, xy = TRUE)

  # Base plot with raster
  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = raster_df, ggplot2::aes(x = x, y = y, fill = layer)) +
    ggplot2::scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                                 na.value = "transparent", name = "Depth") +
    ggplot2::theme_minimal()

  # Add tracks with appropriate coloring
  tracks <- simulation_results$tracks
  
  # Apply sampling if requested
  if (sample_rate < 1 && sample_rate > 0) {
    # Sample points for state/phase visualization
    n_points <- nrow(tracks)
    sample_indices <- sort(sample(1:n_points, size = round(n_points * sample_rate)))
    tracks_sampled <- tracks[sample_indices, ]
  } else {
    tracks_sampled <- tracks
  }
  
  # Determine coloring based on user choice and available data
  if (color_by == "behavioral_state" && "behavioral_state" %in% names(tracks)) {
    # Color by behavioral state
    p <- p + ggplot2::geom_path(data = tracks,
                                ggplot2::aes(x = x, y = y, group = path_id),
                                alpha = path_alpha * 0.5, size = 0.3, color = "grey70") +
             ggplot2::geom_point(data = tracks_sampled,
                                ggplot2::aes(x = x, y = y, color = behavioral_state),
                                size = point_size, alpha = path_alpha) +
             ggplot2::scale_color_manual(name = "Behavioral State",
                                       values = c("cruise" = "yellow", 
                                                 "search" = "orange", 
                                                 "rest" = "red"),
                                       breaks = c("cruise", "search", "rest"))
  } else if (color_by == "spawning_phase" && "spawning_phase" %in% names(tracks)) {
    # Color by spawning phase
    p <- p + ggplot2::geom_path(data = tracks,
                                ggplot2::aes(x = x, y = y, group = path_id),
                                alpha = path_alpha * 0.5, size = 0.3, color = "grey70") +
             ggplot2::geom_point(data = tracks_sampled,
                                ggplot2::aes(x = x, y = y, color = spawning_phase),
                                size = point_size, alpha = path_alpha) +
             ggplot2::scale_color_manual(name = "Spawning Phase",
                                       values = c("none" = "yellow", 
                                                 "search" = "orange", 
                                                 "spawning" = "green"),
                                       breaks = c("none", "search", "spawning"))
  } else if (!is.null(path_color)) {
    # Single color for all paths
    p <- p + ggplot2::geom_path(data = tracks,
                                ggplot2::aes(x = x, y = y, group = path_id),
                                color = path_color, alpha = path_alpha, size = 0.8)
  } else {
    # Default: color by path ID
    p <- p + ggplot2::geom_path(data = tracks,
                                ggplot2::aes(x = x, y = y, group = path_id, color = factor(path_id)),
                                alpha = path_alpha, size = 0.8) +
             ggplot2::scale_color_discrete(name = "Path ID")
  }

  # Add detections if requested and available
  if (show_detections && !is.null(simulation_results$station_detections) &&
      nrow(simulation_results$station_detections) > 0) {

    # Create system-level detections from station detections
    # A location is "detected" if ANY receiver detected it
    system_detections <- simulation_results$station_detections %>%
      dplyr::group_by(path_id, step, time_seconds, x, y) %>%
      dplyr::summarise(detected = max(detected), .groups = 'drop')

    # Separate detected and missed detections
    detected_events <- system_detections %>% dplyr::filter(detected == 1)
    missed_events <- system_detections %>% dplyr::filter(detected == 0)

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
        dplyr::filter(detected == 1) %>%
        dplyr::group_by(station_id) %>%
        dplyr::summarise(detection_count = dplyr::n(), .groups = 'drop')

      # Get receiver coordinates and merge with detection counts
      receiver_coords <- sf::st_coordinates(receiver_frame)
      receiver_df <- data.frame(
        station_id = receiver_frame$station_id,
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

  p <- p + ggplot2::coord_sf() +
    ggplot2::labs(title = "Fish Tracks and Detections")

  return(p)
}

#' Calculate detection rate summaries from fish track simulations
#'
#' Computes detection performance metrics including overall detection rates,
#' per-path statistics, and summary statistics across all simulated fish tracks.
#' Works with receiver-specific detection data to create system-level summaries.
#'
#' @param simulation_results List output from \code{\link{simulate_fish_tracks}}
#'   containing tracks and station_detections data frames.
#'
#' @return A list containing:
#'   \item{overall}{List with total_steps, detected_steps, detection_rate, and detection_percentage}
#'   \item{by_path}{Data frame with detection statistics for each individual path}
#'   \item{by_station}{Data frame with detection statistics for each receiver}
#'   \item{summary_stats}{List with mean, median, min, max, and standard deviation of detection rates}
#'
#'   Returns NULL if no detection data is available.
#'
#' @details
#' The function analyzes detection success at each step of simulated fish tracks,
#' providing both aggregate and individual path performance metrics. System-level
#' detections are derived from station-level data (a location is "detected" if
#' ANY receiver detected the fish at that location).
#'
#' @examples
#' \dontrun{
#' # Simulate fish tracks with detections
#' fish_sim <- simulate_fish_tracks(
#'   raster = depth_raster,
#'   station_distances = distances,
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
#'
#' # Access per-station performance
#' station_performance <- detection_stats$by_station
#' }
#'
#' @seealso \code{\link{simulate_fish_tracks}}, \code{\link{print_detection_summary}}, \code{\link{plot_detection_rates}}
#'
#' @export
calculate_detection_summaries <- function(simulation_results) {

  if (is.null(simulation_results$station_detections) || nrow(simulation_results$station_detections) == 0) {
    cat("No detection data available in simulation results.\n")
    return(NULL)
  }

  station_detections <- simulation_results$station_detections

  # Create system-level detections (detected if ANY receiver detected)
  system_detections <- station_detections %>%
    dplyr::group_by(path_id, step, time_seconds, x, y) %>%
    dplyr::summarise(detected = max(detected), .groups = 'drop')

  # Overall detection rate
  overall_detected <- sum(system_detections$detected == 1, na.rm = TRUE)
  overall_total <- nrow(system_detections)
  overall_rate <- overall_detected / overall_total

  # Detection rate by path
  path_summary <- system_detections %>%
    dplyr::group_by(path_id) %>%
    dplyr::summarise(
      total_steps = dplyr::n(),
      detected_steps = sum(detected == 1, na.rm = TRUE),
      detection_rate = detected_steps / total_steps,
      .groups = 'drop'
    )

  # Detection rate by station
  station_summary <- station_detections %>%
    dplyr::group_by(station_id) %>%
    dplyr::summarise(
      total_opportunities = dplyr::n(),
      detections = sum(detected == 1, na.rm = TRUE),
      detection_rate = detections / total_opportunities,
      mean_detection_prob = mean(detection_prob, na.rm = TRUE),
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
    by_station = station_summary,
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
#' overall rates, summary statistics, individual path performance, and
#' receiver-specific performance.
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
#'   \item Individual receiver performance breakdown
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

  cat("\n")
}

#' Create visualization plots of detection rate performance
#'
#' Generates multiple plots to visualize detection performance including
#' per-path detection rates, per-station rates, distribution of rates,
#' and time series analysis.
#'
#' @param detection_stats List output from \code{\link{calculate_detection_summaries}}
#'   containing detection performance metrics.
#' @param simulation_results Optional. List output from \code{\link{simulate_fish_tracks}}.
#'   If provided, creates additional time series plots. Default is NULL.
#'
#' @return A list of ggplot2 objects:
#'   \item{by_path}{Bar chart showing detection rate by individual path}
#'   \item{by_station}{Bar chart showing detection rate by individual receiver}
#'   \item{distribution}{Histogram of detection rates across all paths}
#'   \item{time_series}{Cumulative detection rate over time (if simulation_results provided)}
#'
#'   Returns NULL if no detection statistics are available.
#'
#' @details
#' The function creates four types of visualizations:
#' \enumerate{
#'   \item Bar chart comparing detection rates across individual fish paths
#'   \item Bar chart comparing detection rates across receivers
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
#' print(plots$by_station)
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

  # 2. Bar plot of detection rates by station
  station_data <- detection_stats$by_station

  p2 <- ggplot2::ggplot(station_data, ggplot2::aes(x = factor(station_id), y = detection_rate)) +
    ggplot2::geom_col(fill = "lightgreen", alpha = 0.7) +
    ggplot2::geom_point(ggplot2::aes(y = mean_detection_prob), color = "red", size = 2) +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, max(c(station_data$detection_rate, station_data$mean_detection_prob), na.rm = TRUE) * 1.1)) +
    ggplot2::labs(
      title = "Detection Rate by Receiver Station",
      subtitle = "Bars = actual detection rate, Red dots = expected detection probability",
      x = "Station ID",
      y = "Detection Rate"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  plot_list$by_station <- p2

  # 3. Histogram of detection rates across paths
  p3 <- ggplot2::ggplot(path_data, ggplot2::aes(x = detection_rate)) +
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

  plot_list$distribution <- p3

  # 4. Time series of detections (if simulation results provided)
  if (!is.null(simulation_results) && !is.null(simulation_results$station_detections)) {

    # Create system-level detections for time series
    system_detections <- simulation_results$station_detections %>%
      dplyr::group_by(path_id, step, time_seconds, x, y) %>%
      dplyr::summarise(detected = max(detected), .groups = 'drop')

    detection_data <- system_detections %>%
      dplyr::group_by(path_id) %>%
      dplyr::arrange(time_seconds) %>%
      dplyr::mutate(
        cumulative_detections = cumsum(detected),
        cumulative_rate = cumulative_detections / dplyr::row_number()
      )

    p4 <- ggplot2::ggplot(detection_data, ggplot2::aes(x = time_seconds, y = cumulative_rate,
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

    plot_list$time_series <- p4
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
#'   containing tracks and station_detections data frames.
#' @param create_plots Logical. Whether to generate visualization plots.
#'   Default is TRUE.
#' @param display_plots Logical. Whether to display the generated plots.
#'   Default is TRUE. Set to FALSE to create plots without displaying them.
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
#' station_plot <- analysis$plots$by_station
#' distribution_plot <- analysis$plots$distribution
#'
#' # Analysis without plots (faster)
#' stats_only <- analyze_detection_performance(fish_simulation, create_plots = FALSE)
#' 
#' # Create plots but don't display them (useful in scripts/vignettes)
#' analysis_quiet <- analyze_detection_performance(fish_simulation, 
#'                                                 create_plots = TRUE, 
#'                                                 display_plots = FALSE)
#' }
#'
#' @seealso \code{\link{simulate_fish_tracks}}, \code{\link{calculate_detection_summaries}},
#'   \code{\link{print_detection_summary}}, \code{\link{plot_detection_rates}}
#'
#' @export
analyze_detection_performance <- function(simulation_results, create_plots = TRUE, display_plots = TRUE) {

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

    # Display plots if requested
    if (display_plots) {
      if (!is.null(plots$by_path)) print(plots$by_path)
      if (!is.null(plots$by_station)) print(plots$by_station)
      if (!is.null(plots$distribution)) print(plots$distribution)
      if (!is.null(plots$time_series)) print(plots$time_series)
    }
  }

  # Return everything
  return(list(
    statistics = detection_stats,
    plots = plots
  ))
}

#' Plot behavioral states over time with temperature coloring
#'
#' Creates a time series visualization showing behavioral state transitions
#' for each simulated fish path, with background temperature information.
#' Only works with simulation results that include behavioral state data.
#'
#' @param simulation_results List output from \code{\link{simulate_fish_tracks}}
#'   containing tracks with behavioral_state and water_temp_c columns.
#' @param show_temperature Logical. Whether to show temperature as background color.
#'   Default is TRUE.
#' @param facet_by_path Logical. Whether to create separate panels for each fish path.
#'   Default is TRUE for multiple paths, FALSE for single path.
#' @param time_units Character. Time units for x-axis: "seconds", "minutes", "hours", or "days".
#'   Default is "hours".
#'
#' @return A ggplot2 object showing behavioral states over time.
#'   Returns NULL if behavioral state data is not available.
#'
#' @details
#' The plot shows:
#' \itemize{
#'   \item Behavioral states on y-axis (rest/search/cruise)
#'   \item Gray lines connecting state transitions over time
#'   \item Points colored by temperature (if show_temperature = TRUE) or behavioral state
#'   \item Time progression on x-axis in specified units
#' }
#'
#' When show_temperature = FALSE, state colors are:
#' \itemize{
#'   \item Cruise: Blue (active directed movement)
#'   \item Search: Orange (active foraging movement)  
#'   \item Rest: Red (minimal movement)
#' }
#'
#' When show_temperature = TRUE, point coloring uses the plasma viridis color scale (purple=cold, yellow=warm).
#'
#' @examples
#' \dontrun{
#' # Simulate fish with behavioral states
#' fish_sim <- simulate_fish_tracks(
#'   raster = depth_raster,
#'   station_distances = distances,
#'   species = "Walleye",
#'   fish_size_cm = 45,
#'   n_paths = 3,
#'   n_steps = 200
#' )
#'
#' # Plot behavioral states with temperature
#' plot_behavioral_states(fish_sim)
#'
#' # Plot without temperature background
#' plot_behavioral_states(fish_sim, show_temperature = FALSE)
#'
#' # Plot as single panel with time in minutes
#' plot_behavioral_states(fish_sim, facet_by_path = FALSE, time_units = "minutes")
#' }
#'
#' @seealso \code{\link{simulate_fish_tracks}}
#'
#' @export
plot_behavioral_states <- function(simulation_results, show_temperature = TRUE, 
                                  facet_by_path = NULL, time_units = "hours") {
  
  # Check if behavioral state data is available
  if (is.null(simulation_results$tracks) || 
      !"behavioral_state" %in% names(simulation_results$tracks)) {
    cat("No behavioral state data available in simulation results.\n")
    cat("Make sure to run simulate_fish_tracks() with a species and behavioral_states = TRUE.\n")
    return(NULL)
  }
  
  tracks <- simulation_results$tracks
  
  # Check if temperature data is available
  has_temperature <- "water_temp_c" %in% names(tracks) && 
                    !all(is.na(tracks$water_temp_c))
  
  if (show_temperature && !has_temperature) {
    warning("Temperature data not available, plotting without temperature background.")
    show_temperature <- FALSE
  }
  
  # Set default faceting based on number of paths
  if (is.null(facet_by_path)) {
    facet_by_path <- length(unique(tracks$path_id)) > 1
  }
  
  # Convert time units
  time_divisor <- switch(time_units,
                        "seconds" = 1,
                        "minutes" = 60,
                        "hours" = 3600,
                        "days" = 86400,
                        3600)  # Default to hours
  
  tracks$time_converted <- tracks$time_seconds / time_divisor
  
  # Create numeric version of behavioral state for plotting
  tracks$state_numeric <- as.numeric(factor(tracks$behavioral_state, 
                                           levels = c("rest", "search", "cruise")))
  
  # Define state colors
  state_colors <- c("rest" = "#d62728",      # Red
                   "search" = "#ff7f0e",     # Orange  
                   "cruise" = "#1f77b4")     # Blue
  
  # Start building the plot
  p <- ggplot2::ggplot(tracks, ggplot2::aes(x = time_converted))
  
  # Add behavioral state lines connecting states
  p <- p + ggplot2::geom_line(ggplot2::aes(y = state_numeric, group = path_id),
                             color = "gray60", size = 0.8, alpha = 0.7)
  
  # Add behavioral state points colored by temperature or state
  if (show_temperature) {
    p <- p + ggplot2::geom_point(ggplot2::aes(y = state_numeric, color = water_temp_c),
                                 size = 2.5, alpha = 0.9) +
      ggplot2::scale_color_viridis_c(option = "plasma", name = "Water\nTemp (\u00B0C)",
                                    guide = ggplot2::guide_colorbar(order = 1))
  } else {
    p <- p + ggplot2::geom_point(ggplot2::aes(y = state_numeric, color = behavioral_state),
                                 size = 2.5, alpha = 0.9) +
      ggplot2::scale_color_manual(values = state_colors, name = "Behavioral\nState",
                                 guide = ggplot2::guide_legend(order = 1))
  }
  
  # Set y-axis
  p <- p + ggplot2::scale_y_continuous(breaks = 1:3, 
                                      labels = c("Rest", "Search", "Cruise"),
                                      limits = c(0.5, 3.5))
  
  # Add faceting if requested
  if (facet_by_path) {
    p <- p + ggplot2::facet_wrap(~ paste("Fish Path", path_id), 
                                scales = "free_x", ncol = 1)
  }
  
  # Add labels and theme
  time_label <- paste0("Time (", time_units, ")")
  
  p <- p + ggplot2::labs(
    title = "Behavioral States Over Time",
    subtitle = if (show_temperature) "Point color shows water temperature, lines show state transitions" else "Point color shows behavioral state, lines show state transitions",
    x = time_label,
    y = "Behavioral State"
  ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
  
  # Add path legend if not faceted and multiple paths
  if (!facet_by_path && length(unique(tracks$path_id)) > 1) {
    p <- p + ggplot2::aes(linetype = factor(path_id)) +
      ggplot2::scale_linetype_discrete(name = "Path ID",
                                      guide = ggplot2::guide_legend(order = 3))
  }
  
  return(p)
}

#' Analyze behavioral state patterns by temperature bins
#'
#' Creates visualizations showing time spent in each behavioral state and 
#' state transition probabilities across different temperature ranges.
#' Only works with simulation results that include behavioral state data.
#'
#' @param simulation_results List output from \code{\link{simulate_fish_tracks}}
#'   containing tracks with behavioral_state and water_temp_c columns.
#' @param n_temp_bins Integer. Number of temperature bins to create. Default is 5.
#' @param create_plots Logical. Whether to generate and return plots. Default is TRUE.
#' @param time_units Character. Time units for analysis: "seconds", "minutes", "hours", or "days".
#'   Default is "hours".
#' @param by_fish_id Logical. Whether to create separate plots for each fish path. Default is FALSE.
#'
#' @return A list containing:
#'   \item{data}{Data frame with behavioral statistics by temperature bin}
#'   \item{plots}{List of ggplot2 objects (if create_plots = TRUE):}
#'   \itemize{
#'     \item time_in_state: Stacked bar chart showing time in each state by temperature
#'     \item transition_probs: Heatmap showing transition probabilities by temperature
#'   }
#'   Returns NULL if behavioral state data is not available.
#'
#' @details
#' The function analyzes behavioral patterns across temperature ranges by:
#' \enumerate{
#'   \item Binning temperature data into equal-width ranges
#'   \item Calculating time spent in each behavioral state per temperature bin
#'   \item Computing state transition probabilities within each temperature bin
#'   \item Creating visualizations of these patterns
#' }
#'
#' @examples
#' \dontrun{
#' # Analyze behavioral patterns by temperature
#' temp_analysis <- analyze_behavioral_temperature(fish_sim)
#'
#' # View the plots
#' print(temp_analysis$plots$time_in_state)
#' print(temp_analysis$plots$transition_probs)
#'
#' # Access the data
#' behavioral_data <- temp_analysis$data
#'
#' # Custom temperature bins
#' analysis_10bins <- analyze_behavioral_temperature(fish_sim, n_temp_bins = 10)
#' 
#' # Individual fish analysis
#' individual_analysis <- analyze_behavioral_temperature(fish_sim, by_fish_id = TRUE)
#' }
#'
#' @seealso \code{\link{simulate_fish_tracks}}, \code{\link{plot_behavioral_states}}
#'
#' @export
analyze_behavioral_temperature <- function(simulation_results, n_temp_bins = 5, 
                                          create_plots = TRUE, time_units = "hours",
                                          by_fish_id = FALSE) {
  
  # Check if behavioral state data is available
  if (is.null(simulation_results$tracks) || 
      !"behavioral_state" %in% names(simulation_results$tracks)) {
    cat("No behavioral state data available in simulation results.\n")
    cat("Make sure to run simulate_fish_tracks() with a species and behavioral_states = TRUE.\n")
    return(NULL)
  }
  
  tracks <- simulation_results$tracks
  
  # Check if temperature data is available
  if (!"water_temp_c" %in% names(tracks) || all(is.na(tracks$water_temp_c))) {
    cat("No temperature data available in simulation results.\n")
    return(NULL)
  }
  
  # Convert time units
  time_divisor <- switch(time_units,
                        "seconds" = 1,
                        "minutes" = 60,
                        "hours" = 3600,
                        "days" = 86400,
                        3600)  # Default to hours
  
  time_step_converted <- simulation_results$parameters$time_step / time_divisor
  
  # Create temperature bins
  temp_range <- range(tracks$water_temp_c, na.rm = TRUE)
  
  # Check if temperature range is too small for binning
  if (diff(temp_range) < 0.1) {
    # If temperature is essentially constant, create a single bin
    tracks$temp_bin <- factor(paste0(round(temp_range[1], 1), "\u00B0C"), 
                             levels = paste0(round(temp_range[1], 1), "\u00B0C"))
    tracks$temp_bin_mid <- 1
    cat("Temperature range very small (", round(diff(temp_range), 2), "\u00B0C). Using single temperature bin.\n")
  } else {
    # Create proper temperature bins
    temp_breaks <- seq(temp_range[1], temp_range[2], length.out = n_temp_bins + 1)
    tracks$temp_bin <- cut(tracks$water_temp_c, breaks = temp_breaks, 
                          include.lowest = TRUE, dig.lab = 1)
    tracks$temp_bin_mid <- as.numeric(tracks$temp_bin)
  }
  
  # Calculate time in each state by temperature bin
  if (by_fish_id) {
    # Keep fish ID separate for individual analysis
    time_in_state <- tracks %>%
      dplyr::group_by(path_id, temp_bin, behavioral_state) %>%
      dplyr::summarise(
        n_steps = dplyr::n(),
        time_spent = n_steps * time_step_converted,
        mean_temp = mean(water_temp_c, na.rm = TRUE),
        .groups = 'drop'
      )
    
    # Calculate proportions within each temperature bin and fish
    time_proportions <- time_in_state %>%
      dplyr::group_by(path_id, temp_bin) %>%
      dplyr::mutate(
        total_time_bin = sum(time_spent),
        proportion = time_spent / total_time_bin
      ) %>%
      dplyr::ungroup()
      
  } else {
    # Aggregate across all fish
    time_in_state <- tracks %>%
      dplyr::group_by(path_id, temp_bin, behavioral_state) %>%
      dplyr::summarise(
        n_steps = dplyr::n(),
        time_spent = n_steps * time_step_converted,
        mean_temp = mean(water_temp_c, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      dplyr::group_by(temp_bin, behavioral_state) %>%
      dplyr::summarise(
        total_time = sum(time_spent),
        mean_temp = mean(mean_temp),
        n_observations = sum(n_steps),
        .groups = 'drop'
      )
    
    # Calculate proportions within each temperature bin
    time_proportions <- time_in_state %>%
      dplyr::group_by(temp_bin) %>%
      dplyr::mutate(
        total_time_bin = sum(total_time),
        proportion = total_time / total_time_bin
      ) %>%
      dplyr::ungroup()
  }
  
  # Calculate state transitions by temperature bin
  transition_data <- tracks %>%
    dplyr::arrange(path_id, step) %>%
    dplyr::group_by(path_id) %>%
    dplyr::mutate(
      next_state = dplyr::lead(behavioral_state),
      transition = paste(behavioral_state, "->", next_state)
    ) %>%
    dplyr::filter(!is.na(next_state)) %>%
    dplyr::ungroup()
  
  # Calculate transition probabilities by temperature bin
  if (by_fish_id) {
    # Keep fish ID separate for individual analysis
    transition_probs <- transition_data %>%
      dplyr::group_by(path_id, temp_bin, behavioral_state, next_state) %>%
      dplyr::summarise(
        n_transitions = dplyr::n(),
        mean_temp = mean(water_temp_c, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      dplyr::group_by(path_id, temp_bin, behavioral_state) %>%
      dplyr::mutate(
        total_from_state = sum(n_transitions),
        transition_prob = n_transitions / total_from_state
      ) %>%
      dplyr::ungroup() %>%
      dplyr::filter(behavioral_state != next_state)  # Remove self-transitions
  } else {
    # Aggregate across all fish
    transition_probs <- transition_data %>%
      dplyr::group_by(temp_bin, behavioral_state, next_state) %>%
      dplyr::summarise(
        n_transitions = dplyr::n(),
        mean_temp = mean(water_temp_c, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      dplyr::group_by(temp_bin, behavioral_state) %>%
      dplyr::mutate(
        total_from_state = sum(n_transitions),
        transition_prob = n_transitions / total_from_state
      ) %>%
      dplyr::ungroup() %>%
      dplyr::filter(behavioral_state != next_state)  # Remove self-transitions
  }
  
  # Combine results
  analysis_data <- list(
    time_in_state = time_proportions,
    transition_probs = transition_probs,
    temp_bins = levels(tracks$temp_bin),
    time_units = time_units
  )
  
  # Create plots if requested
  plots <- NULL
  if (create_plots) {
    
    # Plot 1: Time in each state by temperature
    if (by_fish_id) {
      p1 <- ggplot2::ggplot(time_proportions, ggplot2::aes(x = mean_temp, y = proportion, 
                                                           fill = behavioral_state)) +
        ggplot2::geom_col(position = "stack", alpha = 0.8) +
        ggplot2::scale_fill_manual(values = c("rest" = "#d62728", "search" = "#ff7f0e", 
                                             "cruise" = "#1f77b4"), name = "Behavioral\nState") +
        ggplot2::scale_y_continuous(labels = scales::percent) +
        ggplot2::facet_wrap(~ paste("Fish", path_id), scales = "free_x") +
        ggplot2::labs(
          title = "Time Spent in Each Behavioral State by Temperature (Individual Fish)",
          x = "Temperature (\u00B0C)",
          y = "Proportion of Time",
          subtitle = paste("Separate panels for each fish across", n_temp_bins, "temperature bins")
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "right", 
                      strip.text = ggplot2::element_text(face = "bold"))
    } else {
      p1 <- ggplot2::ggplot(time_proportions, ggplot2::aes(x = mean_temp, y = proportion, 
                                                           fill = behavioral_state)) +
        ggplot2::geom_col(position = "stack", alpha = 0.8) +
        ggplot2::scale_fill_manual(values = c("rest" = "#d62728", "search" = "#ff7f0e", 
                                             "cruise" = "#1f77b4"), name = "Behavioral\nState") +
        ggplot2::scale_y_continuous(labels = scales::percent) +
        ggplot2::labs(
          title = "Time Spent in Each Behavioral State by Temperature",
          x = "Temperature (\u00B0C)",
          y = "Proportion of Time",
          subtitle = paste("Stacked bars show relative time allocation across", n_temp_bins, "temperature bins")
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "right")
    }
    
    # Plot 2: Transition probabilities by temperature
    if (nrow(transition_probs) > 0) {
      # Create complete transition matrix for better visualization
      if (by_fish_id) {
        complete_transitions <- transition_probs %>%
          dplyr::select(path_id, temp_bin, mean_temp, behavioral_state, next_state, transition_prob) %>%
          tidyr::complete(path_id, temp_bin, behavioral_state, next_state, 
                         fill = list(transition_prob = 0)) %>%
          dplyr::group_by(path_id, temp_bin) %>%
          dplyr::mutate(mean_temp = mean(mean_temp, na.rm = TRUE)) %>%
          dplyr::ungroup() %>%
          dplyr::filter(behavioral_state != next_state)
        
        p2 <- ggplot2::ggplot(complete_transitions, 
                             ggplot2::aes(x = mean_temp, y = paste(behavioral_state, "→", next_state),
                                         fill = transition_prob)) +
          ggplot2::geom_tile(color = "white", size = 0.5) +
          ggplot2::scale_fill_viridis_c(option = "plasma", name = "Transition\nProbability",
                                       labels = scales::percent) +
          ggplot2::facet_wrap(~ paste("Fish", path_id), scales = "free_x") +
          ggplot2::labs(
            title = "State Transition Probabilities by Temperature (Individual Fish)",
            x = "Temperature (\u00B0C)",
            y = "State Transition",
            subtitle = "Probability of switching between behavioral states for each fish"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 10),
            legend.position = "right",
            strip.text = ggplot2::element_text(face = "bold")
          )
      } else {
        complete_transitions <- transition_probs %>%
          dplyr::select(temp_bin, mean_temp, behavioral_state, next_state, transition_prob) %>%
          tidyr::complete(temp_bin, behavioral_state, next_state, 
                         fill = list(transition_prob = 0)) %>%
          dplyr::group_by(temp_bin) %>%
          dplyr::mutate(mean_temp = mean(mean_temp, na.rm = TRUE)) %>%
          dplyr::ungroup() %>%
          dplyr::filter(behavioral_state != next_state)
        
        p2 <- ggplot2::ggplot(complete_transitions, 
                             ggplot2::aes(x = mean_temp, y = paste(behavioral_state, "→", next_state),
                                         fill = transition_prob)) +
          ggplot2::geom_tile(color = "white", size = 0.5) +
          ggplot2::scale_fill_viridis_c(option = "plasma", name = "Transition\nProbability",
                                       labels = scales::percent) +
          ggplot2::labs(
            title = "State Transition Probabilities by Temperature",
            x = "Temperature (\u00B0C)",
            y = "State Transition",
            subtitle = "Probability of switching between behavioral states"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 10),
            legend.position = "right"
          )
      }
    } else {
      p2 <- ggplot2::ggplot() + 
        ggplot2::annotate("text", x = 0.5, y = 0.5, 
                         label = "No state transitions detected", size = 5) +
        ggplot2::theme_void()
    }
    
    plots <- list(
      time_in_state = p1,
      transition_probs = p2
    )
  }
  
  return(list(
    data = analysis_data,
    plots = plots
  ))
}
