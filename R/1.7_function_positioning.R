# Helper functions for timezone handling (replicated from walks_detections for consistency)
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
#'   and metadata, typically from point generation functions. This parameter
#'   is optional when station_info is provided, as receiver_stations will be
#'   automatically created from station_info coordinates.
#' @param de_model A fitted detection efficiency model object (e.g., from
#'   \code{\link{create_logistic_curve_depth}}). The model should accept
#'   dist_m and depth_m as predictors. Default is NULL, which requires
#'   DE_pred column to already exist in station_distances_df.
#' @param time_aggregation Character. Method for time aggregation. Options are:
#'   \itemize{
#'     \item "seconds" - Numeric seconds with bin_size_seconds (default)
#'     \item "hour" - Hourly aggregation using POSIX datetime
#'     \item "day" - Daily aggregation using POSIX datetime  
#'     \item "month" - Monthly aggregation using POSIX datetime
#'   }
#'   Default is "seconds".
#' @param bin_size_seconds Numeric. Time bin size in seconds for aggregating
#'   detections when time_aggregation = "seconds". Default is 3600 (1 hour).
#' @param detection_weight Numeric. Weight given to detection events in the
#'   integrated probability calculation (0-1). Default is 0.5.
#' @param non_detection_weight Numeric. Weight given to non-detection events
#'   in the integrated probability calculation (0-1). Default is 0.5.
#' @param integration_method Character. Method for integrating detection and
#'   non-detection evidence. Options are:
#'   \itemize{
#'     \item "subtractive" (default): Detection field is the base; non-detection
#'       evidence subtracts from it: \code{det - (nondet * non_detection_weight)},
#'       clamped to 0. Produces tight, detection-anchored position estimates.
#'     \item "multiplicative": Detection field scaled down by non-detection
#'       evidence: \code{det * (1 - nondet * non_detection_weight)}. Smoother
#'       penalty than subtractive; stays non-negative naturally.
#'     \item "additive": Original WADE formula: weighted sum of detection and
#'       inverted non-detection probabilities. Can inflate spatial footprint
#'       beyond detection zones.
#'   }
#'   For "subtractive" and "multiplicative", \code{detection_weight} is ignored
#'   (detection is always the base); only \code{non_detection_weight} controls
#'   the strength of non-detection penalty.
#' @param max_non_detection_distance Numeric. Maximum distance (in meters) from
#'   detecting stations to consider non-detecting stations. Set to NULL to
#'   include all stations. Default is 2000.
#' @param normalization_method Character. Method for normalizing detection
#'   efficiency values. Options are "min_max", "z_score", or "robust".
#'   Default is "min_max".
#' @param fish_id_col Character. Name of the column containing fish identifiers.
#'   Default is "path_id".
#' @param time_col Character. Name of the column containing time values. Can be
#'   numeric seconds (for time_aggregation = "seconds") or POSIX datetime
#'   (for time_aggregation = "hour"/"day"/"month"). Default is "time_seconds".
#' @param station_col Character. Name of the column containing station identifiers.
#'   Default is "station_id".
#' @param station_info Station information as data frame, CSV file path, or sf object.
#'   Must contain station identifier (station_id or point_id) and coordinates (x,y).
#'   If sf object (e.g., from generate_exact_regular_points), coordinates are
#'   automatically extracted. Optional temporal columns: start_date, end_date for
#'   deployment windows. When provided, receiver_stations is auto-created. Default is NULL.
#' @param temporal_info Data frame with daily environmental conditions for 
#'   temporal DE prediction. Must contain date column and environmental predictors 
#'   used by de_model. Only used when station_info contains temporal columns
#'   and de_model is provided. Default is NULL.
#' @param crs Character. Coordinate reference system for station_info coordinates.
#'   If NULL (default), attempts to auto-detect from station_distances_df coordinate
#'   ranges or falls back to "EPSG:32617". Use same CRS as your depth_raster.
#' @param include_barriers Logical. Whether to apply barrier masking to prevent position
#'   estimates through land obstacles. Default is FALSE. When TRUE, detection efficiency
#'   is set to 0 for any cell-receiver pair where the line-of-sight crosses a barrier.
#'   Requires \code{crosses_barrier} column in \code{station_distances_df} (generated by
#'   \code{calculate_station_distances()} with \code{barrier_raster} parameter).
#' @param verbose Logical. Whether to print progress messages. Default is TRUE.
#'
#' @return A list containing:
#'   \item{position_probabilities}{Data frame with integrated position probabilities
#'     for each fish, time period, and spatial cell. The \code{integrated_prob}
#'     column is rescaled to [0, 1] per fish/time period.}
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
#'   \item Time binning of detection events (standardized to time_period)
#'   \item Aggregation of detection data with spatial distance information
#'   \item Creation of non-detection events for nearby stations
#'   \item Normalization of detection efficiency across receivers
#'   \item Integration of detection and non-detection probabilities
#'   \item Calculation of weighted position estimates
#' }
#'
#' For additive integration, detection and non-detection weights must sum to 1.
#' For subtractive and multiplicative integration, only non_detection_weight
#' controls the penalty strength (0-1); detection_weight is ignored. The algorithm
#' focuses non-detection analysis on stations within a specified distance of
#' detecting stations to maintain biological realism and computational efficiency.
#'
#' Barrier Masking:
#' When \code{include_barriers = TRUE}, the function masks position estimates at cells
#' where the direct path to any receiver crosses a land barrier. This prevents physically
#' impossible position estimates through islands, shorelines, or other obstacles. The
#' barrier information must be pre-computed in the \code{station_distances_df} data frame
#' using \code{calculate_station_distances()} with a barrier raster. Barrier masking works
#' with both static DE mode (pre-computed \code{DE_pred}) and temporal DE mode (on-the-fly
#' DE calculation using \code{de_model}).
#'
#' Time Aggregation Options:
#' \itemize{
#'   \item "seconds": Uses numeric time with bin_size_seconds (backward compatible)
#'   \item "hour": Groups POSIX datetime to hourly bins (e.g., 2024-07-01 15:00:00)
#'   \item "day": Groups POSIX datetime to daily bins (e.g., 2024-07-01)
#'   \item "month": Groups POSIX datetime to monthly bins (e.g., 2024-07-01)
#' }
#'
#' @examples
#' \dontrun{
#' # Basic positioning analysis (backward compatible)
#' results <- calculate_fish_positions(
#'   station_detections = fish_tracks$detections,
#'   station_distances_df = distances,
#'   receiver_stations = stations,
#'   bin_size_seconds = 3600
#' )
#'
#' # POSIX datetime aggregation by day
#' results <- calculate_fish_positions(
#'   station_detections = fish_tracks$detections,
#'   station_distances_df = distances,
#'   receiver_stations = stations,
#'   time_col = "datetime",
#'   time_aggregation = "day"
#' )
#'
#' # POSIX datetime aggregation by hour
#' results <- calculate_fish_positions(
#'   station_detections = fish_tracks$detections,
#'   station_distances_df = distances,
#'   receiver_stations = stations,
#'   time_col = "datetime",
#'   time_aggregation = "hour"
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
                                     time_aggregation = "seconds",
                                     bin_size_seconds = 3600,
                                     detection_weight = 0.5,
                                     non_detection_weight = 0.5,
                                     integration_method = "subtractive",
                                     max_non_detection_distance = 2000,
                                     weighting_method = "information_theoretic",
                                     percentile_cutoff = 0.95,
                                     temporal_grouping = "day",
                                     dampening_factor = 1,
                                     normalization_method = "min_max",
                                     fish_id_col = "path_id",
                                     time_col = "time_seconds",
                                     station_col = "station_id",
                                     station_info = NULL,
                                     temporal_info = NULL,
                                     crs = NULL,
                                     include_barriers = FALSE,
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
  
  # ===== CHARACTER ID SUPPORT =====
  # Handle character station_id and fish_id by converting to integers internally
  # This ensures compatibility with field data while maintaining backward compatibility
  
  station_id_mapping <- NULL
  fish_id_mapping <- NULL
  character_ids_detected <- FALSE
  
  # Check if station_detections has character station IDs
  if (is.character(station_detections[[station_col]])) {
    if (verbose) cat("Converting character station IDs to integers for internal processing...\n")
    character_ids_detected <- TRUE
    
    # Create station mapping
    unique_stations <- unique(station_detections[[station_col]])
    station_id_mapping <- data.frame(
      station_id_char = unique_stations,
      station_id_int = seq_along(unique_stations)
    )
    
    # Convert station_detections
    station_detections <- station_detections %>%
      dplyr::left_join(station_id_mapping, by = setNames("station_id_char", station_col)) %>%
      dplyr::mutate(!!station_col := station_id_int) %>%
      dplyr::select(-station_id_int)
  }
  
  # Check if station_detections has character fish IDs
  if (is.character(station_detections[[fish_id_col]])) {
    if (verbose) cat("Converting character fish IDs to integers for internal processing...\n")
    character_ids_detected <- TRUE
    
    # Create fish mapping
    unique_fish <- unique(station_detections[[fish_id_col]])
    fish_id_mapping <- data.frame(
      fish_id_char = unique_fish,
      fish_id_int = seq_along(unique_fish)
    )
    
    # Convert station_detections
    station_detections <- station_detections %>%
      dplyr::left_join(fish_id_mapping, by = setNames("fish_id_char", fish_id_col)) %>%
      dplyr::mutate(!!fish_id_col := fish_id_int) %>%
      dplyr::select(-fish_id_int)
  }
  
  # Convert station_distances_df if station mapping exists
  if (!is.null(station_id_mapping)) {
    # Check the station column name in station_distances_df (usually "station_no")
    station_dist_col <- if ("station_no" %in% names(station_distances_df)) "station_no" else "station_id"
    
    if (is.character(station_distances_df[[station_dist_col]])) {
      station_distances_df <- station_distances_df %>%
        dplyr::left_join(station_id_mapping, by = setNames("station_id_char", station_dist_col)) %>%
        dplyr::filter(!is.na(station_id_int)) %>%
        dplyr::select(-!!station_dist_col) %>%
        dplyr::rename(!!station_dist_col := station_id_int)
    }
  }
  
  # Convert station_info if provided and station mapping exists
  if (!is.null(station_info) && !is.null(station_id_mapping)) {
    if (is.character(station_info[[station_col]])) {
      station_info <- station_info %>%
        dplyr::left_join(station_id_mapping, by = setNames("station_id_char", station_col)) %>%
        dplyr::filter(!is.na(station_id_int)) %>%
        dplyr::select(-!!station_col) %>%
        dplyr::rename(!!station_col := station_id_int)
    }
  }
  
  # Convert receiver_stations if provided (for backward compatibility)
  if (!missing(receiver_stations) && !is.null(receiver_stations) && !is.null(station_id_mapping)) {
    if (station_col %in% names(receiver_stations)) {
      if ("sf" %in% class(receiver_stations)) {
        # Handle sf object
        receiver_coords <- sf::st_coordinates(receiver_stations)
        receiver_df <- sf::st_drop_geometry(receiver_stations)
        
        if (is.character(receiver_df[[station_col]])) {
          receiver_df <- receiver_df %>%
            dplyr::left_join(station_id_mapping, by = setNames("station_id_char", station_col)) %>%
            dplyr::filter(!is.na(station_id_int)) %>%
            dplyr::select(-!!station_col) %>%
            dplyr::rename(!!station_col := station_id_int)
          
          # Recreate sf object if needed
          receiver_stations <- sf::st_as_sf(
            cbind(receiver_df, receiver_coords),
            coords = c("X", "Y"), 
            crs = sf::st_crs(receiver_stations)
          )
        }
      } else {
        # Handle regular data frame
        if (is.character(receiver_stations[[station_col]])) {
          receiver_stations <- receiver_stations %>%
            dplyr::left_join(station_id_mapping, by = setNames("station_id_char", station_col)) %>%
            dplyr::filter(!is.na(station_id_int)) %>%
            dplyr::select(-!!station_col) %>%
            dplyr::rename(!!station_col := station_id_int)
        }
      }
    }
  }
  
  if (character_ids_detected && verbose) {
    cat("Character ID conversion completed. Processing with integer IDs internally.\n")
  }
  if (time_aggregation %in% c("hour", "day", "month") && !requireNamespace("lubridate", quietly = TRUE)) {
    stop("Package 'lubridate' needed for POSIX time aggregation. Please install it.",
         call. = FALSE)
  }

  # Validate inputs
  if (!integration_method %in% c("subtractive", "multiplicative", "additive")) {
    stop("integration_method must be 'subtractive', 'multiplicative', or 'additive'")
  }

  if (integration_method == "additive") {
    if (abs(detection_weight + non_detection_weight - 1) > 1e-10) {
      stop("For additive integration, detection_weight and non_detection_weight must sum to 1")
    }
  } else {
    if (non_detection_weight < 0 || non_detection_weight > 1) {
      stop("non_detection_weight must be between 0 and 1")
    }
    if (!missing(detection_weight) && detection_weight != 0.5) {
      warning("detection_weight is ignored for '", integration_method,
              "' integration (detection field is always the base). ",
              "Only non_detection_weight controls non-detection strength.")
    }
  }

  if (!weighting_method %in% c("information_theoretic", "normalize_stations", "raw")) {
    stop("weighting_method must be 'information_theoretic', 'normalize_stations', or 'raw'")
  }

  if (!normalization_method %in% c("min_max", "z_score", "robust")) {
    stop("normalization_method must be 'min_max', 'z_score', or 'robust'")
  }

  if (!time_aggregation %in% c("seconds", "hour", "day", "month")) {
    stop("time_aggregation must be one of: 'seconds', 'hour', 'day', 'month'")
  }
  
  if (percentile_cutoff <= 0 || percentile_cutoff > 1) {
    stop("percentile_cutoff must be between 0 and 1")
  }

  # Validate barrier data if barrier masking is requested
  if (include_barriers) {
    if (!"crosses_barrier" %in% names(station_distances_df)) {
      stop("include_barriers = TRUE requires 'crosses_barrier' column in station_distances_df.\n",
           "Run calculate_station_distances() with barrier_raster parameter to generate this column.")
    }
    if (verbose) cat("Barrier masking enabled: DE will be set to 0 where crosses_barrier = TRUE\n")
  }

  # Check required columns
  required_cols <- c(fish_id_col, time_col, station_col)
  missing_cols <- setdiff(required_cols, names(station_detections))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns in station_detections:", paste(missing_cols, collapse = ", ")))
  }
  
  # Determine mode based on available parameters
  use_temporal_de <- !is.null(station_info) && !is.null(temporal_info) && !is.null(de_model)
  use_station_info <- !is.null(station_info)
  
  if (use_station_info) {
    # Handle different input types for station_info
    if (is.character(station_info)) {
      # Load CSV file
      station_info <- utils::read.csv(station_info, stringsAsFactors = FALSE)
    } else if (inherits(station_info, "sf")) {
      # Handle sf object (from generate_exact_regular_points)
      if (verbose) cat("Converting sf object to station_info data frame...\n")
      
      # Extract coordinates and convert to data frame
      coords <- sf::st_coordinates(station_info)
      station_info_df <- sf::st_drop_geometry(station_info)
      station_info_df$x <- coords[,1]
      station_info_df$y <- coords[,2]
      station_info <- station_info_df
    }
    
    # Validate station_info - station identifier is required
    # Accept either station_id or point_id as the identifier column
    has_station_id <- "station_id" %in% names(station_info)
    has_point_id <- "point_id" %in% names(station_info)
    
    if (!has_station_id && !has_point_id) {
      stop("station_info must contain either 'station_id' or 'point_id' column")
    }
    
    # Standardize to station_id column
    if (has_point_id && !has_station_id) {
      station_info$station_id <- station_info$point_id
      if (verbose) cat("Using 'point_id' column as 'station_id'\n")
    }
    
    # Check for x,y coordinates (should exist after sf conversion or be provided)
    if (!"x" %in% names(station_info) || !"y" %in% names(station_info)) {
      stop("station_info must contain 'x' and 'y' coordinate columns")
    }
    
    # Check if temporal columns exist
    has_temporal_cols <- all(c("start_date", "end_date") %in% names(station_info))
    
    if (use_temporal_de) {
      if (verbose) cat("Using temporal DE prediction with station deployment windows...\n")
      
      # Validate temporal columns are present
      if (!has_temporal_cols) {
        stop("station_info must contain 'start_date' and 'end_date' columns for temporal DE prediction")
      }
      
      # Validate temporal_info  
      if (!"date" %in% names(temporal_info)) {
        stop("temporal_info must contain 'date' column")
      }
      
      if (is.character(temporal_info)) {
        temporal_info <- utils::read.csv(temporal_info, stringsAsFactors = FALSE)
      }
      
    } else {
      if (verbose) {
        if (has_temporal_cols) {
          cat("Using station_info with deployment windows but static DE predictions...\n")
        } else {
          cat("Using station_info with static station locations and DE predictions...\n")
        }
      }
      
      # Check for DE_pred column when not using temporal DE
      if (is.null(de_model) && !"DE_pred" %in% names(station_distances_df)) {
        stop("station_distances_df must contain 'DE_pred' column when de_model is NULL and temporal_info not provided")
      }
    }
    
  } else {
    if (verbose) cat("Using receiver_stations parameter with static DE predictions...\n")
    
    # Check for receiver_stations in legacy mode
    if (is.null(receiver_stations)) {
      stop("Either station_info or receiver_stations must be provided")
    }
    
    # Check for DE_pred column in legacy mode
    if (is.null(de_model) && !"DE_pred" %in% names(station_distances_df)) {
      stop("station_distances_df must contain 'DE_pred' column when de_model is NULL and no temporal parameters provided")
    }
  }

  # Validate time data type matches aggregation method
  time_data <- station_detections[[time_col]]
  original_timezone <- NULL
  
  if (time_aggregation == "seconds") {
    if (!is.numeric(time_data)) {
      stop("time_aggregation 'seconds' requires numeric time data in column '", time_col, "'. ",
           "Current column type: ", class(time_data)[1])
    }
  } else {
    if (!lubridate::is.POSIXct(time_data) && !lubridate::is.POSIXlt(time_data)) {
      stop("time_aggregation '", time_aggregation, "' requires POSIX datetime in column '", time_col, "'. ",
           "Current column type: ", class(time_data)[1])
    }
    
    # Store original timezone and standardize to UTC for internal processing
    original_timezone <- store_timezone_info(time_data)
    station_detections[[time_col]] <- standardize_datetime(time_data, "UTC")
    if (verbose) {
      cat("Timezone info: Input timezone =", original_timezone, ", processing in UTC\n")
    }
  }
  
  # Create receiver_stations from station_info if using station_info mode
  if (use_station_info) {
    if (verbose) cat("Creating receiver_stations from station_info coordinates...\n")
    
    # Determine CRS for station_info coordinates
    target_crs <- crs
    
    if (is.null(target_crs)) {
      # Try to auto-detect from station_distances_df coordinate ranges
      if ("x" %in% names(station_distances_df) && "y" %in% names(station_distances_df)) {
        x_range <- range(station_distances_df$x, na.rm = TRUE)
        y_range <- range(station_distances_df$y, na.rm = TRUE)
        
        # Check if coordinates are in UTM range (heuristic)
        if (x_range[1] > 100000 && y_range[1] > 1000000) {
          target_crs <- "EPSG:32617"  # UTM Zone 17N NAD83
          if (verbose) cat("Auto-detected CRS: EPSG:32617 (UTM Zone 17N) based on coordinate ranges\n")
        }
      }
      
      # Default fallback
      if (is.null(target_crs)) {
        target_crs <- "EPSG:32617"
        if (verbose) cat("Using default CRS: EPSG:32617\n")
      }
    } else {
      if (verbose) cat("Using specified CRS:", target_crs, "\n")
    }
    
    receiver_stations <- station_info %>%
      sf::st_as_sf(coords = c("x", "y"), crs = target_crs) %>%
      dplyr::mutate(point_id = station_id) %>%
      dplyr::select(point_id, dplyr::everything())
    
    if (verbose) cat("Created receiver_stations with", nrow(receiver_stations), "stations from station_info\n")
  }
  
  # Setup temporal DE prediction if enabled
  daily_de_lookup <- NULL
  if (use_temporal_de) {
    if (verbose) cat("Setting up temporal DE prediction...\n")
    
    # Helper function to parse dates in multiple formats (from simulate_fish_tracks)
    parse_date_flexible <- function(date_vector) {
      # Return early if already Date class
      if (inherits(date_vector, "Date")) {
        return(date_vector)
      }
      
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
    
    # Get positioning date range from detection data
    if (time_aggregation == "seconds") {
      # Convert numeric seconds to dates for temporal DE (assume origin)
      time_range <- range(station_detections[[time_col]], na.rm = TRUE)
      # For now, assume a reasonable origin - this could be parameterized
      origin_date <- as.Date("2025-01-01")  # Could be made parameter
      positioning_start_date <- origin_date + (time_range[1] / 86400)
      positioning_end_date <- origin_date + (time_range[2] / 86400)
    } else {
      # Extract dates from POSIX datetime
      time_range <- range(station_detections[[time_col]], na.rm = TRUE)
      positioning_start_date <- as.Date(time_range[1])
      positioning_end_date <- as.Date(time_range[2])
    }
    
    unique_dates <- seq(from = positioning_start_date, to = positioning_end_date, by = "day")
    if (verbose) cat("Pre-computing daily DE predictions for", length(unique_dates), "days...\n")
    
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
        prediction_data <- station_distances_df %>%
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
    
    if (verbose) cat("Completed daily DE predictions for", length(daily_de_lookup), "days\n")
  }

  if (verbose) cat("=== CALCULATING FISH POSITIONS ===\n")

  # Step 1: Standardize column names and create time periods
  if (verbose) cat("Step 1: Standardizing columns and creating time periods...\n")
  
  if (time_aggregation == "seconds") {
    # EXISTING APPROACH: numeric seconds with bin_size_seconds
    station_detections_processed <- station_detections %>%
      dplyr::rename(
        fish_id = !!fish_id_col,
        station_id = !!station_col
      ) %>%
      dplyr::mutate(
        time_period = floor(!!dplyr::sym(time_col) / bin_size_seconds) * bin_size_seconds,
        time_period_posix = NA,
        time_period_label = paste0("Bin_", time_period),
        time_aggregation_method = "seconds",
        .after = dplyr::all_of(time_col)
      )
  } else {
    # NEW POSIX APPROACH: detect and aggregate
    station_detections_processed <- station_detections %>%
      dplyr::rename(
        fish_id = !!fish_id_col,
        station_id = !!station_col
      ) %>%
      dplyr::mutate(
        # Create time periods based on aggregation method
        time_period_posix = case_when(
          time_aggregation == "hour" ~ lubridate::floor_date(!!dplyr::sym(time_col), "hour"),
          time_aggregation == "day" ~ lubridate::floor_date(!!dplyr::sym(time_col), "day"),
          time_aggregation == "month" ~ lubridate::floor_date(!!dplyr::sym(time_col), "month"),
          TRUE ~ lubridate::floor_date(!!dplyr::sym(time_col), "day")  # Default to day
        ),
        
        # Convert to numeric for compatibility with existing workflow
        time_period = as.numeric(time_period_posix),
        
        # Create human-readable labels
        time_period_label = case_when(
          time_aggregation == "hour" ~ format(time_period_posix, "%Y-%m-%d %H:00"),
          time_aggregation == "day" ~ format(time_period_posix, "%Y-%m-%d"),
          time_aggregation == "month" ~ format(time_period_posix, "%Y-%m"),
          TRUE ~ format(time_period_posix, "%Y-%m-%d")  # Default to day format
        ),
        
        time_aggregation_method = time_aggregation,
        .after = dplyr::all_of(time_col)
      )
  }

  # Create time metadata lookup table to preserve POSIX info
  time_metadata <- station_detections_processed %>%
    dplyr::select(time_period, time_period_posix, time_period_label, time_aggregation_method) %>%
    dplyr::distinct()

  # Step 2: Aggregate detections for prediction
  if (verbose) cat("Step 2: Aggregating detection data...\n")
  
  # Filter detections to only include deployed stations if using temporal deployment info
  detections_for_processing <- station_detections_processed %>% dplyr::filter(detected == 1)
  
  if (use_station_info && all(c("start_date", "end_date") %in% names(station_info))) {
    # Filter detections by station deployment status on each day
    detections_for_processing <- detections_for_processing %>%
      dplyr::mutate(detection_date = if (time_aggregation == "seconds") {
        as.Date("2025-01-01") + (time_period / 86400)  # Convert seconds to date
      } else {
        as.Date(time_period_posix)
      }) %>%
      dplyr::left_join(
        station_info %>% 
          dplyr::select(station_id, start_date, end_date),
        by = "station_id"
      ) %>%
      dplyr::filter(
        detection_date >= start_date & detection_date <= end_date
      ) %>%
      dplyr::select(-detection_date, -start_date, -end_date)
    
    if (verbose) cat("Filtered detections to", nrow(detections_for_processing), "from deployed stations\n")
  }
  
  detection_probs <- aggregate_detections_for_prediction(
    station_detections = detections_for_processing,
    station_distances_df = station_distances_df,
    daily_de_lookup = daily_de_lookup,
    use_temporal_de = use_temporal_de,
    time_aggregation = time_aggregation,
    include_barriers = include_barriers
  )

  # Step 3: Create non-detections using anti_join approach
  if (verbose) cat("Step 3: Creating non-detection data...\n")
  non_detections <- create_non_detections(
    station_detections = detections_for_processing,  # Use filtered detections
    points_regular = receiver_stations,
    max_distance_from_detecting = max_non_detection_distance,
    station_info = if (use_station_info && all(c("start_date", "end_date") %in% names(station_info))) station_info else NULL,
    use_temporal_filtering = use_station_info && all(c("start_date", "end_date") %in% names(station_info)),
    time_aggregation = time_aggregation
  )

  # Step 4: Aggregate non-detections
  if (verbose) cat("Step 4: Aggregating non-detection data...\n")
  non_detection_probs <- aggregate_non_detections(
    non_detections = non_detections,
    station_distances_df = station_distances_df,
    daily_de_lookup = daily_de_lookup,
    use_temporal_de = use_temporal_de,
    time_aggregation = time_aggregation,
    include_barriers = include_barriers
  )

  # Step 5: Apply weighting method
  if (verbose) cat("Step 5: Applying", weighting_method, "weighting...\n")
  
  if (weighting_method == "information_theoretic") {
    # Calculate robust station characteristics
    station_characteristics <- calculate_station_de_characteristics(
      detection_probs, 
      percentile_cutoff = percentile_cutoff,
      temporal_grouping = temporal_grouping,
      time_aggregation = time_aggregation
    )
    
    # Apply information-theoretic weighting to detections
    detection_probs_norm <- apply_information_weighting(
      detection_probs, 
      station_characteristics,
      weighting_type = "detection",
      dampening_factor = dampening_factor
    )
    
    # Apply raw DE weighting to non-detections (no information adjustment)
    non_detection_probs_norm <- apply_information_weighting(
      non_detection_probs,
      station_characteristics, 
      weighting_type = "non_detection",
      dampening_factor = 1  # No dampening for non-detections
    )
    
    # Print diagnostic information if verbose
    if (verbose) {
      cat("\nStation DE Characteristics Summary:\n")
      cat("  Number of stations:", length(unique(station_characteristics$station_id)), "\n")
      cat("  Mean effective DE:", round(mean(station_characteristics$station_effective_DE, na.rm = TRUE), 3), "\n")
      cat("  DE range:", round(range(station_characteristics$station_effective_DE, na.rm = TRUE), 3), "\n")
      
      # Calculate expected multiplier range
      de_range <- range(station_characteristics$station_effective_DE, na.rm = TRUE)
      cat("  Expected multiplier range:", round(1/de_range[2], 2), "-", round(1/de_range[1], 2), "x\n")
    }
    
  } else if (weighting_method == "normalize_stations") {
    # Keep current approach for backwards compatibility
    # But improve non-detection handling
    detection_probs_norm <- normalize_DE_by_station(
      data = detection_probs,
      DE_col = "DE_pred", 
      station_col = "station_id",
      method = normalization_method
    )
    
    # IMPROVED: Use raw DE for non-detections instead of normalizing
    non_detection_probs_norm <- non_detection_probs %>%
      dplyr::mutate(
        DE_pred_normalized = DE_pred  # Raw DE, no station normalization
      )
    
    if (verbose) cat("  Using station normalization for detections, raw DE for non-detections\n")
    
  } else if (weighting_method == "raw") {
    # No weighting adjustment - use raw values
    detection_probs_norm <- detection_probs %>%
      dplyr::mutate(
        DE_pred_normalized = DE_pred
      )
    
    non_detection_probs_norm <- non_detection_probs %>%
      dplyr::mutate(
        DE_pred_normalized = DE_pred
      )
    
    if (verbose) cat("  Using raw DE values without weighting adjustments\n")
  }

  # Step 6: Combine detection and non-detection data
  if (verbose) cat("Step 6: Combining detection and non-detection data...\n")
  position_probs_combined <- dplyr::bind_rows(
    detection_probs_norm %>% dplyr::mutate(type = "detection"),
    non_detection_probs_norm %>% dplyr::mutate(type = "non-detection")
  )

  # Step 7: Calculate integrated positioning probabilities
  if (verbose) cat("Step 7: Calculating position probabilities...\n")
  position_probs <- aggregate_probability(
    df = position_probs_combined,
    detection_weight = detection_weight,
    non_detection_weight = non_detection_weight,
    integration_method = integration_method,
    normalize_method = if(weighting_method == "information_theoretic") "none" else "global"
  )

  # Step 8: Add time metadata back to position probabilities
  if (verbose) cat("Step 8: Adding time metadata and station coordinates...\n")
  position_probs <- position_probs %>%
    dplyr::left_join(time_metadata, by = "time_period")
  
  # Restore original fish IDs if character conversion was used
  # Skip restoration for field data that already preserved character IDs
  if (!is.null(fish_id_mapping) && "fish_id_int" %in% names(fish_id_mapping) && 
      "fish_id_char" %in% names(fish_id_mapping)) {
    position_probs <- position_probs %>%
      dplyr::left_join(fish_id_mapping, by = setNames("fish_id_int", "fish_id")) %>%
      dplyr::select(-fish_id) %>%
      dplyr::rename(fish_id = fish_id_char)
  }

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

  # Add station coordinates to processed detections for plotting
  station_detections_plot <- station_detections_processed %>%
    dplyr::left_join(station_coords, by = "station_id")

  # Also add coordinates to detection and non-detection data
  detection_probs_norm <- detection_probs_norm %>%
    dplyr::left_join(station_coords, by = "station_id") %>%
    dplyr::left_join(time_metadata, by = "time_period")

  non_detection_probs_norm <- non_detection_probs_norm %>%
    dplyr::left_join(station_coords, by = "station_id") %>%
    dplyr::left_join(time_metadata, by = "time_period")
    
  # Restore original fish IDs in all result components if character conversion was used
  # Skip restoration for field data that already preserved character IDs
  if (!is.null(fish_id_mapping) && "fish_id_int" %in% names(fish_id_mapping) && 
      "fish_id_char" %in% names(fish_id_mapping)) {
    station_detections_plot <- station_detections_plot %>%
      dplyr::left_join(fish_id_mapping, by = setNames("fish_id_int", "fish_id")) %>%
      dplyr::select(-fish_id) %>%
      dplyr::rename(fish_id = fish_id_char)
    
    detection_probs_norm <- detection_probs_norm %>%
      dplyr::left_join(fish_id_mapping, by = setNames("fish_id_int", "fish_id")) %>%
      dplyr::select(-fish_id) %>%
      dplyr::rename(fish_id = fish_id_char)
    
    non_detection_probs_norm <- non_detection_probs_norm %>%
      dplyr::left_join(fish_id_mapping, by = setNames("fish_id_int", "fish_id")) %>%
      dplyr::select(-fish_id) %>%
      dplyr::rename(fish_id = fish_id_char)
  }

  # Create enhanced summary statistics
  summary_stats <- list(
    n_fish = dplyr::n_distinct(position_probs$fish_id),
    n_time_periods = dplyr::n_distinct(position_probs$time_period),
    n_cells = dplyr::n_distinct(position_probs$cell_id),
    n_stations = dplyr::n_distinct(station_coords$station_id),
    time_aggregation = time_aggregation,
    time_bin_size = if(time_aggregation == "seconds") bin_size_seconds else NA,
    detection_weight = detection_weight,
    non_detection_weight = non_detection_weight,
    integration_method = integration_method,
    total_detections = nrow(station_detections_processed),
    total_position_estimates = nrow(position_probs)
  )

  if (verbose) {
    cat("=== POSITIONING COMPLETE ===\n")
    cat("Fish tracked:", summary_stats$n_fish, "\n")
    cat("Time periods:", summary_stats$n_time_periods, "\n")
    cat("Time aggregation:", summary_stats$time_aggregation, "\n")
    if (time_aggregation == "seconds") {
      cat("Time bin size:", summary_stats$time_bin_size, "seconds\n")
    }
    cat("Spatial cells:", summary_stats$n_cells, "\n")
    cat("Position estimates generated:", summary_stats$total_position_estimates, "\n")
  }

  # Restore original timezone in output data
  if (!is.null(original_timezone) && time_aggregation != "seconds") {
    # Restore timezone in datasets that contain datetime columns
    if ("datetime" %in% names(position_probs)) {
      position_probs$datetime <- restore_timezone(position_probs$datetime, original_timezone)
    }
    if ("datetime" %in% names(detection_probs_norm)) {
      detection_probs_norm$datetime <- restore_timezone(detection_probs_norm$datetime, original_timezone)
    }
    if ("datetime" %in% names(non_detection_probs_norm)) {
      non_detection_probs_norm$datetime <- restore_timezone(non_detection_probs_norm$datetime, original_timezone)
    }
    if ("datetime" %in% names(station_detections_plot)) {
      station_detections_plot$datetime <- restore_timezone(station_detections_plot$datetime, original_timezone)
    }
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
#' @param track_data Optional complete fish track data for showing full track context.
#'   Should be a data frame with columns for fish_id (or path_id), time_seconds, x, y.
#'   If not provided, only detection-based track segments will be shown.
#' @param fish_select Numeric. Fish ID to plot. Default is 1.
#' @param time_select Time period to plot. Can be:
#'   \itemize{
#'     \item NULL (default) - Plot all time periods (no filtering)
#'     \item Numeric - Use as time_period value (e.g., 0, 3600, 7200)
#'     \item Character - Interpreted based on format:
#'       \itemize{
#'         \item "YYYY-MM-DD" - Daily time period (e.g., "2024-07-01")
#'         \item "YYYY-MM-DD HH:00" - Hourly time period (e.g., "2024-07-01 15:00")
#'         \item Other formats - Treated as time_period_label
#'       }
#'     \item POSIXct/POSIXlt - Exact datetime match
#'   }
#' @param actual_track_size Numeric. Size of the actual track points (green). Default is 1.2.
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
#' The unified time_select parameter intelligently handles different input types:
#' \itemize{
#'   \item NULL displays all time periods (no filtering)
#'   \item Numeric values match time_period values
#'   \item Date strings ("YYYY-MM-DD") match daily periods
#'   \item Datetime strings ("YYYY-MM-DD HH:00") match hourly periods
#'   \item Other strings match time_period_label values
#'   \item POSIXct objects match exact datetime values
#' }
#'
#' @examples
#' \dontrun{
#' # Calculate positions
#' results <- calculate_fish_positions(station_detections, distances, stations)
#'
#' # Create all plots for all time periods (default)
#' plot_fish_positions(results, depth_raster_df = depth_data)
#'
#' # Create only integrated probability plot for all time
#' plot_fish_positions(results, plot_type = "integrated")
#'
#' # Select specific time using numeric value
#' plot_fish_positions(results, time_select = 3600)
#'
#' # Select time using date string
#' plot_fish_positions(results, time_select = "2024-07-01")
#'
#' # Select time using datetime string
#' plot_fish_positions(results, time_select = "2024-07-01 15:00")
#'
#' # Select time using POSIXct object
#' plot_fish_positions(results, time_select = as.POSIXct("2024-07-01 15:00:00"))
#'
#' # Customize actual track point size
#' plot_fish_positions(results, actual_track_size = 2.5)
#'
#' # Get individual plots as list
#' plot_fish_positions(results, return_list = TRUE)
#' plots$detection
#' plots$integrated
#' }
#'
#' @seealso \code{\link{calculate_fish_positions}}, \code{\link{analyze_positioning_performance}}
#'
#' @export
plot_fish_positions <- function(positioning_results,
                                depth_raster_df = NULL,
                                track_data = NULL,
                                fish_select = 1,
                                time_select = NULL,
                                actual_track_size = 1.2,
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
  if (!is.null(time_select) && 
      (inherits(time_select, "POSIXt") || 
       (is.character(time_select) && grepl("^\\d{4}-\\d{2}-\\d{2}", time_select))) && 
      !requireNamespace("lubridate", quietly = TRUE)) {
    stop("Package 'lubridate' needed for date/time selection. Please install it.",
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

  # SMART TIME SELECTION - Handle different input types
  time_filter_value <- NULL
  apply_time_filter <- !is.null(time_select)
  
  if (!is.null(time_select)) {
    # Determine input type and process accordingly
    if (is.numeric(time_select)) {
      # NUMERIC TIME PERIOD VALUE
      available_times <- position_probs %>%
        dplyr::filter(fish_id == fish_select) %>%
        dplyr::pull(time_period) %>%
        unique() %>%
        sort()
      
      if (!time_select %in% available_times) {
        stop("Time period ", time_select, " not found for fish_id ", fish_select, 
             ". Available times: ", paste(available_times, collapse = ", "))
      }
      time_filter_value <- time_select
      
    } else if (inherits(time_select, "POSIXt")) {
      # POSIX DATETIME OBJECT
      target_posix <- lubridate::as_datetime(time_select)
      
      if ("time_period_posix" %in% names(position_probs)) {
        available_posix <- position_probs %>%
          dplyr::filter(!is.na(time_period_posix)) %>%
          dplyr::pull(time_period_posix) %>%
          unique() %>%
          sort()
        
        matching_times <- position_probs %>%
          dplyr::filter(!is.na(time_period_posix)) %>%
          dplyr::filter(time_period_posix == target_posix) %>%
          dplyr::pull(time_period) %>%
          unique()
        
        if (length(matching_times) == 0) {
          stop("No exact match found for POSIX time: ", target_posix, 
               "\nAvailable POSIX times: ", paste(available_posix, collapse = ", "))
        }
        time_filter_value <- matching_times[1]
      } else {
        stop("POSIX time selection requires position_probabilities to have time_period_posix column")
      }
      
    } else if (is.character(time_select)) {
      # CHARACTER INPUT - Check format
      
      # Check if it's a date/datetime format
      if (grepl("^\\d{4}-\\d{2}-\\d{2}$", time_select) || 
          grepl("^\\d{4}-\\d{2}-\\d{2} \\d{2}:00$", time_select)) {
        
        # Try to use as POSIX time first
        if ("time_period_posix" %in% names(position_probs)) {
          target_posix <- lubridate::as_datetime(time_select)
          
          matching_times <- position_probs %>%
            dplyr::filter(!is.na(time_period_posix)) %>%
            dplyr::filter(time_period_posix == target_posix) %>%
            dplyr::pull(time_period) %>%
            unique()
          
          if (length(matching_times) > 0) {
            time_filter_value <- matching_times[1]
          }
        }
        
        # If no POSIX match, try as label
        if (is.null(time_filter_value) && "time_period_label" %in% names(position_probs)) {
          matching_times <- position_probs %>%
            dplyr::filter(!is.na(time_period_label)) %>%
            dplyr::filter(time_period_label == time_select) %>%
            dplyr::pull(time_period) %>%
            unique()
          
          if (length(matching_times) > 0) {
            time_filter_value <- matching_times[1]
          }
        }
        
        if (is.null(time_filter_value)) {
          available_options <- c()
          if ("time_period_posix" %in% names(position_probs)) {
            posix_times <- position_probs %>%
              dplyr::filter(!is.na(time_period_posix)) %>%
              dplyr::pull(time_period_posix) %>%
              unique() %>%
              sort()
            available_options <- c(available_options, paste("POSIX times:", paste(posix_times, collapse = ", ")))
          }
          if ("time_period_label" %in% names(position_probs)) {
            labels <- position_probs %>%
              dplyr::filter(!is.na(time_period_label)) %>%
              dplyr::pull(time_period_label) %>%
              unique() %>%
              sort()
            available_options <- c(available_options, paste("Labels:", paste(labels, collapse = ", ")))
          }
          stop("No match found for time: ", time_select, 
               "\nAvailable options:\n", paste(available_options, collapse = "\n"))
        }
        
      } else {
        # Treat as time_period_label
        if ("time_period_label" %in% names(position_probs)) {
          available_labels <- position_probs %>%
            dplyr::filter(!is.na(time_period_label)) %>%
            dplyr::pull(time_period_label) %>%
            unique() %>%
            sort()
          
          matching_times <- position_probs %>%
            dplyr::filter(!is.na(time_period_label)) %>%
            dplyr::filter(time_period_label == time_select) %>%
            dplyr::pull(time_period) %>%
            unique()
          
          if (length(matching_times) == 0) {
            stop("No exact match found for label: ", time_select,
                 "\nAvailable labels: ", paste(available_labels, collapse = ", "))
          }
          time_filter_value <- matching_times[1]
        } else {
          stop("Label time selection requires position_probabilities to have time_period_label column")
        }
      }
    } else {
      stop("time_select must be NULL, numeric, character, or POSIXct")
    }
  }

  # Filter data for selected fish (and optionally time)
  if (apply_time_filter) {
    position_data <- position_probs %>%
      dplyr::filter(fish_id == fish_select, time_period == time_filter_value)
    
    # Get detection summary for this fish/time
    # Use actual n_detections from station_detections if available
    if ("n_detections" %in% names(station_detections)) {
      detection_summary <- station_detections %>%
        dplyr::filter(fish_id == fish_select, time_period == time_filter_value) %>%
        dplyr::group_by(station_id) %>%
        dplyr::summarise(n_detections = sum(n_detections, na.rm = TRUE), .groups = 'drop')
    } else {
      detection_summary <- station_detections %>%
        dplyr::filter(fish_id == fish_select, time_period == time_filter_value) %>%
        dplyr::group_by(station_id) %>%
        dplyr::summarise(n_detections = sum(detected, na.rm = TRUE), .groups = 'drop')
    }
    
    # Add station coordinates
    detection_summary <- detection_summary %>%
      dplyr::left_join(
        detection_data %>% 
        dplyr::filter(fish_id == fish_select, time_period == time_filter_value) %>%
        dplyr::distinct(station_id, .keep_all = TRUE) %>%
        dplyr::select(station_id, station_x, station_y),
        by = "station_id"
      ) %>%
      dplyr::filter(n_detections > 0, !is.na(station_x), !is.na(station_y))
    
    # Get non-detecting stations for this fish/time  
    non_detection_summary <- non_detection_data %>%
      dplyr::filter(fish_id == fish_select, time_period == time_filter_value) %>%
      dplyr::distinct(station_id, .keep_all = TRUE) %>%
      dplyr::select(station_id, station_x, station_y) %>%
      dplyr::filter(!is.na(station_x), !is.na(station_y))
    
    # Get actual fish position if available
    actual_position <- detection_data %>%
      dplyr::filter(fish_id == fish_select, time_period == time_filter_value) %>%
      dplyr::slice(1) %>%
      dplyr::select(x, y)
      
  } else {
    # No time filter - get all data for the fish
    position_data <- position_probs %>%
      dplyr::filter(fish_id == fish_select)
    
    # Get all detections for this fish
    # Use actual n_detections from station_detections if available
    if ("n_detections" %in% names(station_detections)) {
      detection_summary <- station_detections %>%
        dplyr::filter(fish_id == fish_select) %>%
        dplyr::group_by(station_id) %>%
        dplyr::summarise(n_detections = sum(n_detections, na.rm = TRUE), .groups = 'drop')
    } else {
      detection_summary <- station_detections %>%
        dplyr::filter(fish_id == fish_select) %>%
        dplyr::group_by(station_id) %>%
        dplyr::summarise(n_detections = sum(detected, na.rm = TRUE), .groups = 'drop')
    }
    
    # Add station coordinates
    detection_summary <- detection_summary %>%
      dplyr::left_join(
        detection_data %>% 
        dplyr::filter(fish_id == fish_select) %>%
        dplyr::distinct(station_id, .keep_all = TRUE) %>%
        dplyr::select(station_id, station_x, station_y),
        by = "station_id"
      ) %>%
      dplyr::filter(n_detections > 0, !is.na(station_x), !is.na(station_y))
    
    # Get stations that NEVER detected this fish across all time periods
    # First, get all unique stations from the analysis
    all_stations_in_analysis <- dplyr::bind_rows(
      detection_data %>%
        dplyr::filter(fish_id == fish_select) %>%
        dplyr::select(station_id, station_x, station_y) %>%
        dplyr::distinct(),
      non_detection_data %>%
        dplyr::filter(fish_id == fish_select) %>%
        dplyr::select(station_id, station_x, station_y) %>%
        dplyr::distinct()
    ) %>%
      dplyr::distinct(station_id, .keep_all = TRUE)
    
    # Identify stations that never detected (not in detection_summary)
    non_detection_summary <- all_stations_in_analysis %>%
      dplyr::anti_join(detection_summary, by = "station_id") %>%
      dplyr::filter(!is.na(station_x), !is.na(station_y))
    
    # Get all fish positions
    actual_position <- detection_data %>%
      dplyr::filter(fish_id == fish_select) %>%
      dplyr::select(x, y) %>%
      dplyr::distinct()
  }

  # Get fish track for context
  if (!is.null(track_data)) {
    # Use complete track data if provided, filtered to the specific time period
    fish_id_col <- if ("fish_id" %in% names(track_data)) "fish_id" else "path_id"
    
    fish_track <- track_data %>%
      dplyr::filter(!!dplyr::sym(fish_id_col) == fish_select) %>%
      dplyr::select(x, y, !!dplyr::sym(fish_id_col)) %>%
      dplyr::distinct()
  } else {
    # No track data available for field data - create empty track to avoid plotting
    # thousands of probability grid cells as track points
    fish_track <- data.frame(x = numeric(0), y = numeric(0), fish_id = integer(0))
  }

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
                            color = "green", size = actual_track_size, alpha = 0.8)
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
    # Calculate automatic zoom based on fish activity if xlim/ylim not provided  
    # Priority order: track → receiver array coverage (for field data)
    if (is.null(xlim) || is.null(ylim)) {
      # Get fish track coordinates for this fish
      track_coords <- NULL

      # Try to get coordinates from track_data if provided
      if (!is.null(track_data)) {
        fish_id_col <- if ("fish_id" %in% names(track_data)) "fish_id" else "path_id"
        track_coords <- track_data %>%
          dplyr::filter(!!dplyr::sym(fish_id_col) == fish_select) %>%
          dplyr::select(x, y) %>%
          dplyr::filter(!is.na(x), !is.na(y))
      }

      # If no track_data, use fish position from actual_position
      if (is.null(track_coords) || nrow(track_coords) == 0) {
        if (nrow(actual_position) > 0 && !is.na(actual_position$x) && !is.na(actual_position$y)) {
          # For field data, don't use single position estimates for zoom - use receivers instead
          if (!is.null(track_data)) {
            track_coords <- actual_position %>% dplyr::select(x, y)
          }
        }
      }

      # Calculate zoom with buffer - focus on fish activity area only
      all_x_coords <- c()
      all_y_coords <- c()
      
      # Priority 1: Track coordinates if available (most reliable)
      if (!is.null(track_coords) && nrow(track_coords) > 0) {
        all_x_coords <- c(all_x_coords, track_coords$x)
        all_y_coords <- c(all_y_coords, track_coords$y)
      }
      
      # Priority 2: For field data without track, use receiver array from selected time period
      if (length(all_x_coords) == 0) {
        # For field data, get stations active in the selected time period for the specific fish
        # This is much more appropriate than using all stations
        if (is.list(positioning_results) && "station_detections_binned" %in% names(positioning_results)) {
          
          # Filter to selected fish and time period
          time_period_stations <- positioning_results$station_detections_binned %>%
            dplyr::filter(fish_id == fish_select, datetime == time_filter_value) %>%
            dplyr::select(station_id, station_x.x, station_y.x, detected, n_detections) %>%
            dplyr::filter(!is.na(station_x.x), !is.na(station_y.x))
          
          if (nrow(time_period_stations) > 0) {
            # Prioritize stations with detections for tighter zoom
            detected_stations <- time_period_stations %>%
              dplyr::filter(detected == 1, n_detections > 0)
            
            if (nrow(detected_stations) > 0) {
              # Use only stations with detections for primary zoom
              all_x_coords <- c(all_x_coords, detected_stations$station_x.x)
              all_y_coords <- c(all_y_coords, detected_stations$station_y.x)
              
              # Add a few nearby non-detection stations for context (if available)
              if (nrow(detected_stations) < 10) {  # Only if we have few detection stations
                non_detected_stations <- time_period_stations %>%
                  dplyr::filter(detected == 0) %>%
                  dplyr::slice_head(n = 5)  # Add up to 5 non-detection stations for context
                all_x_coords <- c(all_x_coords, non_detected_stations$station_x.x)
                all_y_coords <- c(all_y_coords, non_detected_stations$station_y.x)
              }
            } else {
              # Fallback: use all stations if no detections (shouldn't happen in normal use)
              all_x_coords <- c(all_x_coords, time_period_stations$station_x.x)
              all_y_coords <- c(all_y_coords, time_period_stations$station_y.x)
            }
          }
        }
        
        # Fallback to original logic if field data approach doesn't work
        if (length(all_x_coords) == 0) {
          # Get coordinates from detection_data (has full station info)
          if (!is.null(detection_data) && nrow(detection_data) > 0) {
            station_coords <- detection_data %>%
              dplyr::distinct(station_id, .keep_all = TRUE) %>%
              dplyr::select(station_x, station_y) %>%
              dplyr::filter(!is.na(station_x), !is.na(station_y))
            
            if (nrow(station_coords) > 0) {
              all_x_coords <- c(all_x_coords, station_coords$station_x)
              all_y_coords <- c(all_y_coords, station_coords$station_y)
            }
          }
          
          # Fallback: get from non_detection_data if detection_data didn't work
          if (length(all_x_coords) == 0 && !is.null(non_detection_data) && nrow(non_detection_data) > 0) {
            station_coords <- non_detection_data %>%
              dplyr::distinct(station_id, .keep_all = TRUE) %>%
              dplyr::select(station_x, station_y) %>%
              dplyr::filter(!is.na(station_x), !is.na(station_y))
            
            if (nrow(station_coords) > 0) {
              all_x_coords <- c(all_x_coords, station_coords$station_x)
              all_y_coords <- c(all_y_coords, station_coords$station_y)
            }
          }
        }
      }
      
      # Calculate zoom if we have any coordinates
      if (length(all_x_coords) > 0 && length(all_y_coords) > 0) {
        x_range <- range(all_x_coords, na.rm = TRUE)
        y_range <- range(all_y_coords, na.rm = TRUE)

        # Add healthy buffer around receiver array (15% for field data context)
        x_buffer <- diff(x_range) * 0.15
        y_buffer <- diff(y_range) * 0.15

        # Ensure minimum buffer size (1000m) for good context
        x_buffer <- pmax(x_buffer, 1000)
        y_buffer <- pmax(y_buffer, 1000)

        auto_xlim <- c(x_range[1] - x_buffer, x_range[2] + x_buffer)
        auto_ylim <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)

        # Use automatic limits if not provided by user
        final_xlim <- if (is.null(xlim)) auto_xlim else xlim
        final_ylim <- if (is.null(ylim)) auto_ylim else ylim

        # Use coord_equal for proper aspect ratio with explicit limits
        p <- p + ggplot2::coord_equal(xlim = final_xlim, ylim = final_ylim)
      } else {
        # Fallback: use coord_equal with user limits if provided, otherwise default
        if (!is.null(xlim) && !is.null(ylim)) {
          p <- p + ggplot2::coord_equal(xlim = xlim, ylim = ylim)
        } else {
          p <- p + ggplot2::coord_equal()
        }
      }
    } else {
      # Use user-provided limits with coord_equal
      p <- p + ggplot2::coord_equal(xlim = xlim, ylim = ylim)
    }

    # Create informative subtitle with enhanced time information
    time_info <- if (!apply_time_filter) {
      "All time periods"
    } else if ("time_period_label" %in% names(position_data) && nrow(position_data) > 0) {
      # Use the time_period_label from the filtered data
      unique_label <- unique(position_data$time_period_label)
      if (length(unique_label) == 1 && !is.na(unique_label)) {
        unique_label
      } else {
        paste("Time period:", time_filter_value)
      }
    } else {
      paste("Time period:", time_filter_value)
    }

    # Create subtitle with conditional track reference
    subtitle_parts <- c(paste("Fish", fish_select), time_info)
    # Only mention track if there's actual track data being plotted
    if (!is.null(track_data) && nrow(fish_track) > 1) {
      subtitle_parts <- c(subtitle_parts, "Actual track = green")
    }
    
    # Apply theme and labels
    p <- p +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = title,
        subtitle = paste(subtitle_parts, collapse = " | "),
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
                             ggplot2::aes(x = x, y = y, fill = weighted_mean_DE_normalized_scaled)) +
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
                             ggplot2::aes(x = x, y = y, fill = non_det_DE_normalized_scaled)) +
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
    combined_plot <- plot_list_result$detection | plot_list_result$non_detection | plot_list_result$integrated
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

# Function to aggregate station detections and generate cells for prediction
aggregate_detections_for_prediction <- function(station_detections,
                                                station_distances_df,
                                                daily_de_lookup = NULL,
                                                use_temporal_de = FALSE,
                                                time_aggregation = "seconds",
                                                include_barriers = FALSE) {

  # Aggregate detections by fish, time, and station
  detection_summary <- station_detections %>%
    dplyr::group_by(fish_id, time_period, station_id) %>%
    dplyr::summarise(
      n_detections = ifelse("n_detections" %in% names(station_detections),
                           sum(n_detections, na.rm = TRUE),
                           dplyr::n()),  # Use n_detections column if available, otherwise count rows
      first_detection_time = min(time_period),
      last_detection_time = max(time_period),
      mean_detection_prob = mean(detection_prob, na.rm = TRUE),
      total_distance = mean(distance_to_station, na.rm = TRUE),
      .groups = 'drop'
    )

  # Use temporal DE prediction if enabled, otherwise use static DE from station_distances_df
  if (use_temporal_de && !is.null(daily_de_lookup)) {
    # Get unique dates from detection data
    if (time_aggregation == "seconds") {
      # Convert seconds to dates (assuming origin date as in main function)
      origin_date <- as.Date("2025-01-01")
      detection_dates <- detection_summary %>%
        dplyr::mutate(detection_date = origin_date + (time_period / 86400)) %>%
        dplyr::select(fish_id, time_period, detection_date) %>%
        dplyr::distinct()
    } else {
      # Extract dates directly from POSIX time_period
      detection_dates <- detection_summary %>%
        dplyr::mutate(detection_date = as.Date(as.POSIXct(time_period, origin = "1970-01-01"))) %>%
        dplyr::select(fish_id, time_period, detection_date) %>%
        dplyr::distinct()
    }
    
    # Merge with temporal DE predictions
    prediction_data_list <- list()
    
    for (i in 1:nrow(detection_dates)) {
      current_date <- detection_dates$detection_date[i]
      current_fish <- detection_dates$fish_id[i]
      current_time <- detection_dates$time_period[i]
      
      # Get DE predictions for this date
      date_key <- as.character(current_date)
      
      if (date_key %in% names(daily_de_lookup)) {
        daily_de <- daily_de_lookup[[date_key]]
        
        # Get detection summary for this fish/time
        current_detections <- detection_summary %>%
          dplyr::filter(fish_id == current_fish, time_period == current_time)
        
        # Merge with temporal DE predictions
        # Select columns including barriers if requested
        if (include_barriers) {
          station_cols_to_select <- c("station_no", "cell_id", "x", "y", "raster_value",
                                     "cost_distance", "straight_distance", "tortuosity", "crosses_barrier")
        } else {
          station_cols_to_select <- c("station_no", "cell_id", "x", "y", "raster_value",
                                     "cost_distance", "straight_distance", "tortuosity")
        }

        current_prediction <- current_detections %>%
          dplyr::left_join(
            daily_de,
            by = c("station_id" = "station_no"),
            relationship = "many-to-many"
          ) %>%
          dplyr::left_join(
            station_distances_df %>%
              dplyr::select(dplyr::any_of(station_cols_to_select)),
            by = c("station_id" = "station_no", "cell_id"),
            relationship = "many-to-many"
          ) %>%
          dplyr::filter(!is.na(cell_id))

        # Apply barrier masking if enabled
        if (include_barriers && "crosses_barrier" %in% names(current_prediction)) {
          current_prediction <- current_prediction %>%
            dplyr::mutate(DE_pred = ifelse(crosses_barrier, 0, DE_pred))
        }

        prediction_data_list[[length(prediction_data_list) + 1]] <- current_prediction
      }
    }
    
    if (length(prediction_data_list) > 0) {
      prediction_data <- dplyr::bind_rows(prediction_data_list)
    } else {
      # Fallback to static DE if no temporal data available
      prediction_data <- detection_summary %>%
        dplyr::left_join(
          station_distances_df,
          by = c("station_id" = "station_no"),
          relationship = "many-to-many"
        ) %>%
        dplyr::filter(!is.na(cell_id))

      # Apply barrier masking if enabled
      if (include_barriers && "crosses_barrier" %in% names(prediction_data)) {
        prediction_data <- prediction_data %>%
          dplyr::mutate(DE_pred = ifelse(crosses_barrier, 0, DE_pred))
      }
    }
  } else {
    # Use static DE predictions from station_distances_df
    prediction_data <- detection_summary %>%
      dplyr::left_join(
        station_distances_df,
        by = c("station_id" = "station_no"),
        relationship = "many-to-many"
      ) %>%
      dplyr::filter(!is.na(cell_id))

    # Apply barrier masking if enabled
    if (include_barriers && "crosses_barrier" %in% names(prediction_data)) {
      prediction_data <- prediction_data %>%
        dplyr::mutate(DE_pred = ifelse(crosses_barrier, 0, DE_pred))
    }
  }
  
  # Select final columns and arrange
  prediction_data <- prediction_data %>%
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
        station_min = min(!!dplyr::sym(DE_col), na.rm = TRUE),
        station_max = max(!!dplyr::sym(DE_col), na.rm = TRUE),
        station_range = station_max - station_min,
        !!normalized_col := ifelse(station_range > 0 & is.finite(station_range),
                                   (!!dplyr::sym(DE_col) - station_min) / station_range,
                                   ifelse(is.finite(!!dplyr::sym(DE_col)), 0.5, 0))
      ) %>%
      dplyr::select(-station_min, -station_max, -station_range) %>%
      dplyr::ungroup()
  } else if (method == "z_score") {
    result <- data %>%
      dplyr::group_by(!!dplyr::sym(station_col)) %>%
      dplyr::mutate(
        station_mean = mean(!!dplyr::sym(DE_col), na.rm = TRUE),
        station_sd = stats::sd(!!dplyr::sym(DE_col), na.rm = TRUE),
        !!normalized_col := ifelse(is.finite(station_sd) & station_sd > 0,
                                   (!!dplyr::sym(DE_col) - station_mean) / station_sd,
                                   ifelse(is.finite(!!dplyr::sym(DE_col)), 0, 0))
      ) %>%
      dplyr::select(-station_mean, -station_sd) %>%
      dplyr::ungroup()
  } else if (method == "robust") {
    result <- data %>%
      dplyr::group_by(!!dplyr::sym(station_col)) %>%
      dplyr::mutate(
        station_median = stats::median(!!dplyr::sym(DE_col), na.rm = TRUE),
        station_mad = stats::mad(!!dplyr::sym(DE_col), na.rm = TRUE),
        !!normalized_col := ifelse(is.finite(station_mad) & station_mad > 0,
                                   (!!dplyr::sym(DE_col) - station_median) / station_mad,
                                   ifelse(is.finite(!!dplyr::sym(DE_col)), 0, 0))
      ) %>%
      dplyr::select(-station_median, -station_mad) %>%
      dplyr::ungroup()
  }

  return(result)
}

# Create dataset of non-detecting receivers for each fish-time combination
create_non_detections <- function(station_detections, 
                                  points_regular, 
                                  max_distance_from_detecting = 2000, 
                                  station_info = NULL, 
                                  use_temporal_filtering = FALSE, 
                                  time_aggregation = "seconds") {

  # Get all unique fish-time combinations from station_detections
  # Handle both simulation data (with step) and field data (without step)
  required_cols <- c("fish_id", "time_period")
  optional_cols <- c("step", "x", "y")
  available_cols <- intersect(optional_cols, names(station_detections))
  
  fish_time_data <- station_detections %>%
    dplyr::select(all_of(c(required_cols, available_cols))) %>%
    dplyr::distinct()

  # Get all available stations from points_regular
  # Extract coordinates before dropping geometry
  coords <- sf::st_coordinates(points_regular)
  all_stations_base <- points_regular %>%
    sf::st_drop_geometry() %>%
    dplyr::mutate(
      station_x = coords[, 1],
      station_y = coords[, 2]
    )
  
  # Conditionally include raster_value (depth) column if it exists
  if ("raster_value" %in% names(all_stations_base)) {
    all_stations <- all_stations_base %>%
      dplyr::select(station_id = point_id, station_x, station_y, depth = raster_value)
  } else {
    all_stations <- all_stations_base %>%
      dplyr::select(station_id = point_id, station_x, station_y) %>%
      dplyr::mutate(depth = NA_real_)  # Provide default depth value
  }
  
  # Filter stations by deployment dates if temporal filtering is enabled
  if (use_temporal_filtering && !is.null(station_info)) {
    # Add date information to fish_time_data
    if (time_aggregation == "seconds") {
      # Convert seconds to dates (assuming origin date as in main function)
      origin_date <- as.Date("2025-01-01")
      fish_time_data <- fish_time_data %>%
        dplyr::mutate(detection_date = origin_date + (time_period / 86400))
    } else {
      # Extract dates directly from POSIX time_period
      fish_time_data <- fish_time_data %>%
        dplyr::mutate(detection_date = as.Date(as.POSIXct(time_period, origin = "1970-01-01")))
    }
    
    # Create station-date combinations that are valid (station deployed on that date)
    valid_station_dates <- fish_time_data %>%
      dplyr::select(fish_id, time_period, detection_date) %>%
      dplyr::distinct() %>%
      dplyr::cross_join(station_info %>% dplyr::select(station_id, start_date, end_date)) %>%
      dplyr::filter(detection_date >= start_date & detection_date <= end_date) %>%
      dplyr::select(fish_id, time_period, station_id)
    
    # Filter all_stations to only include those that were deployed for the fish-time combinations
    all_stations_filtered <- all_stations %>%
      dplyr::inner_join(
        valid_station_dates %>% dplyr::select(station_id) %>% dplyr::distinct(),
        by = "station_id"
      )
    
    # Use filtered stations for further processing
    all_stations <- all_stations_filtered
  }

  # Get stations that DID detect for each fish-time
  detecting_stations <- station_detections %>%
    dplyr::select(fish_id, time_period, station_id) %>%
    dplyr::distinct()

  # Get coordinates of detecting stations for each fish-time combination
  detecting_station_coords <- detecting_stations %>%
    dplyr::left_join(all_stations, by = "station_id") %>%
    dplyr::select(fish_id, time_period, detecting_station_id = station_id,
                  detecting_x = station_x, detecting_y = station_y)

  # Filter non-detecting stations based on distance
  if (!is.null(max_distance_from_detecting) && max_distance_from_detecting > 0) {
    candidate_non_detections <- fish_time_data %>%
      dplyr::left_join(detecting_station_coords, by = c("fish_id", "time_period")) %>%
      dplyr::cross_join(all_stations %>%
                          dplyr::rename(candidate_station_id = station_id,
                                        candidate_x = station_x,
                                        candidate_y = station_y,
                                        candidate_depth = depth)) %>%
      dplyr::mutate(
        distance_between_stations = sqrt((candidate_x - detecting_x)^2 + (candidate_y - detecting_y)^2)
      ) %>%
      dplyr::filter(distance_between_stations <= max_distance_from_detecting) %>%
      dplyr::select(fish_id, time_period, station_id = candidate_station_id,
                    station_x = candidate_x, station_y = candidate_y, depth = candidate_depth,
                    any_of(c("step", "x", "y"))) %>%
      dplyr::distinct() %>%
      dplyr::mutate(
        distance_to_station = ifelse(
          "x" %in% names(.) & "y" %in% names(.),
          sqrt((x - station_x)^2 + (y - station_y)^2),
          NA_real_  # For field data without fish positions
        )
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
    dplyr::anti_join(detecting_stations, by = c("fish_id", "time_period", "station_id")) %>%
    dplyr::mutate(
      detected = 0,
      detection_prob = NA
    ) %>%
    dplyr::select(fish_id, any_of(c("step", "x", "y")), time_period, station_id,
                  distance_to_station, detection_prob, detected,
                  station_x, station_y, depth) %>%
    dplyr::arrange(fish_id, time_period, station_id)
  
  # Apply final temporal filtering if enabled
  if (use_temporal_filtering && !is.null(station_info)) {
    # Add detection dates to non_detections
    if (time_aggregation == "seconds") {
      origin_date <- as.Date("2025-01-01")
      non_detections <- non_detections %>%
        dplyr::mutate(detection_date = origin_date + (time_period / 86400))
    } else {
      non_detections <- non_detections %>%
        dplyr::mutate(detection_date = as.Date(as.POSIXct(time_period, origin = "1970-01-01")))
    }
    
    # Filter to only include stations that were deployed on the detection date
    non_detections <- non_detections %>%
      dplyr::left_join(
        station_info %>% dplyr::select(station_id, start_date, end_date),
        by = "station_id"
      ) %>%
      dplyr::filter(
        !is.na(start_date) & !is.na(end_date) &
        detection_date >= start_date & detection_date <= end_date
      ) %>%
      dplyr::select(-detection_date, -start_date, -end_date)
  }

  return(non_detections)
}

# Function to aggregate non-detections and generate cells for prediction
aggregate_non_detections <- function(non_detections,
                                     station_distances_df,
                                     daily_de_lookup = NULL,
                                     use_temporal_de = FALSE,
                                     time_aggregation = "seconds",
                                     include_barriers = FALSE) {

  # Aggregate non-detections by fish, time, and station
  non_detection_summary <- non_detections %>%
    dplyr::group_by(fish_id, time_period, station_id) %>%
    dplyr::summarise(
      n_detections = 0,
      first_detection_time = NA,
      last_detection_time = NA,
      mean_detection_prob = mean(detection_prob, na.rm = TRUE),
      total_distance = mean(distance_to_station, na.rm = TRUE),
      .groups = 'drop'
    )

  # Use temporal DE prediction if enabled, otherwise use static DE from station_distances_df
  if (use_temporal_de && !is.null(daily_de_lookup)) {
    # Get unique dates from non-detection data
    if (time_aggregation == "seconds") {
      # Convert seconds to dates (assuming origin date as in main function)
      origin_date <- as.Date("2025-01-01")
      non_detection_dates <- non_detection_summary %>%
        dplyr::mutate(detection_date = origin_date + (time_period / 86400)) %>%
        dplyr::select(fish_id, time_period, detection_date) %>%
        dplyr::distinct()
    } else {
      # Extract dates directly from POSIX time_period
      non_detection_dates <- non_detection_summary %>%
        dplyr::mutate(detection_date = as.Date(as.POSIXct(time_period, origin = "1970-01-01"))) %>%
        dplyr::select(fish_id, time_period, detection_date) %>%
        dplyr::distinct()
    }
    
    # Merge with temporal DE predictions
    prediction_data_list <- list()
    
    for (i in 1:nrow(non_detection_dates)) {
      current_date <- non_detection_dates$detection_date[i]
      current_fish <- non_detection_dates$fish_id[i]
      current_time <- non_detection_dates$time_period[i]
      
      # Get DE predictions for this date
      date_key <- as.character(current_date)
      
      if (date_key %in% names(daily_de_lookup)) {
        daily_de <- daily_de_lookup[[date_key]]
        
        # Get non-detection summary for this fish/time
        current_non_detections <- non_detection_summary %>%
          dplyr::filter(fish_id == current_fish, time_period == current_time)
        
        # Merge with temporal DE predictions
        # Select columns including barriers if requested
        if (include_barriers) {
          station_cols_to_select <- c("station_no", "cell_id", "x", "y", "raster_value",
                                     "cost_distance", "straight_distance", "tortuosity", "crosses_barrier")
        } else {
          station_cols_to_select <- c("station_no", "cell_id", "x", "y", "raster_value",
                                     "cost_distance", "straight_distance", "tortuosity")
        }

        current_prediction <- current_non_detections %>%
          dplyr::left_join(
            daily_de,
            by = c("station_id" = "station_no"),
            relationship = "many-to-many"
          ) %>%
          dplyr::left_join(
            station_distances_df %>%
              dplyr::select(dplyr::any_of(station_cols_to_select)),
            by = c("station_id" = "station_no", "cell_id"),
            relationship = "many-to-many"
          ) %>%
          dplyr::filter(!is.na(cell_id))

        # Apply barrier masking if enabled
        if (include_barriers && "crosses_barrier" %in% names(current_prediction)) {
          current_prediction <- current_prediction %>%
            dplyr::mutate(DE_pred = ifelse(crosses_barrier, 0, DE_pred))
        }

        prediction_data_list[[length(prediction_data_list) + 1]] <- current_prediction
      }
    }
    
    if (length(prediction_data_list) > 0) {
      prediction_data <- dplyr::bind_rows(prediction_data_list)
    } else {
      # Fallback to static DE if no temporal data available
      prediction_data <- non_detection_summary %>%
        dplyr::left_join(
          station_distances_df,
          by = c("station_id" = "station_no"),
          relationship = "many-to-many"
        ) %>%
        dplyr::filter(!is.na(cell_id))

      # Apply barrier masking if enabled
      if (include_barriers && "crosses_barrier" %in% names(prediction_data)) {
        prediction_data <- prediction_data %>%
          dplyr::mutate(DE_pred = ifelse(crosses_barrier, 0, DE_pred))
      }
    }
  } else {
    # Use static DE predictions from station_distances_df
    prediction_data <- non_detection_summary %>%
      dplyr::left_join(
        station_distances_df,
        by = c("station_id" = "station_no"),
        relationship = "many-to-many"
      ) %>%
      dplyr::filter(!is.na(cell_id))

    # Apply barrier masking if enabled
    if (include_barriers && "crosses_barrier" %in% names(prediction_data)) {
      prediction_data <- prediction_data %>%
        dplyr::mutate(DE_pred = ifelse(crosses_barrier, 0, DE_pred))
    }
  }
  
  # Select final columns and arrange
  prediction_data <- prediction_data %>%
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
aggregate_probability <- function(df, detection_weight = 0.5, non_detection_weight = 0.5,
                                  integration_method = "subtractive", normalize_method = "global") {

  # Validate weights for additive method
  if (integration_method == "additive") {
    if (abs(detection_weight + non_detection_weight - 1) > 1e-10) {
      stop("Weights must sum to 1 for additive integration")
    }
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

  # Apply normalization based on method
  if (normalize_method == "global") {
    combined_data <- combined_data %>%
      dplyr::mutate(
        det_min = min(weighted_mean_DE_normalized, na.rm = TRUE),
        det_max = max(weighted_mean_DE_normalized, na.rm = TRUE),
        det_range = det_max - det_min,
        weighted_mean_DE_normalized_scaled = ifelse(is.finite(det_range) & det_range > 0,
                                                    (weighted_mean_DE_normalized - det_min) / det_range,
                                                    ifelse(is.finite(weighted_mean_DE_normalized), 0.5, 0)),
        non_det_min = min(non_det_DE_normalized, na.rm = TRUE),
        non_det_max = max(non_det_DE_normalized, na.rm = TRUE),
        non_det_range = non_det_max - non_det_min,
        non_det_DE_normalized_scaled = ifelse(is.finite(non_det_range) & non_det_range > 0,
                                              (non_det_DE_normalized - non_det_min) / non_det_range,
                                              ifelse(is.finite(non_det_DE_normalized), 0.5, 0))
      ) %>%
      dplyr::select(-det_min, -det_max, -det_range, -non_det_min, -non_det_max, -non_det_range)
  } else if (normalize_method == "none") {
    # Skip second normalization - values already properly weighted
    combined_data <- combined_data %>%
      dplyr::mutate(
        weighted_mean_DE_normalized_scaled = weighted_mean_DE_normalized,
        non_det_DE_normalized_scaled = non_det_DE_normalized
      )
  } else {
    # Default case - no normalization
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
      non_det_DE_normalized_scaled = ifelse(is.na(non_det_DE_normalized_scaled), 0, non_det_DE_normalized_scaled)
    )

  if (integration_method == "subtractive") {
    # Detection field is the base; non-detection carves away from it
    # Cells with no detection evidence remain at 0
    result <- result %>%
      dplyr::mutate(
        integrated_prob = pmax(0, weighted_mean_DE_normalized_scaled -
          (non_det_DE_normalized_scaled * non_detection_weight))
      )
  } else if (integration_method == "multiplicative") {
    # Detection field scaled down proportionally by non-detection evidence
    # Smoother penalty than subtractive; naturally non-negative
    result <- result %>%
      dplyr::mutate(
        integrated_prob = weighted_mean_DE_normalized_scaled *
          (1 - non_det_DE_normalized_scaled * non_detection_weight)
      )
  } else {
    # Additive (original WADE formula)
    # Weighted sum of detection and inverted non-detection probabilities
    result <- result %>%
      dplyr::mutate(
        integrated_prob = (weighted_mean_DE_normalized_scaled * detection_weight) +
          ((1.0 - non_det_DE_normalized_scaled) * non_detection_weight)
      )
  }

  # Rescale integrated_prob to [0, 1] per fish/time period
  result <- result %>%
    dplyr::group_by(fish_id, time_period) %>%
    dplyr::mutate(
      ip_min = min(integrated_prob, na.rm = TRUE),
      ip_max = max(integrated_prob, na.rm = TRUE),
      ip_range = ip_max - ip_min,
      integrated_prob = ifelse(ip_range > 0,
        (integrated_prob - ip_min) / ip_range,
        integrated_prob)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-ip_min, -ip_max, -ip_range)

  result <- result %>%
    dplyr::select(
      fish_id, time_period, cell_id, x, y, detections,
      weighted_mean_DE_normalized, non_det_DE_normalized,
      weighted_mean_DE_normalized_scaled, non_det_DE_normalized_scaled, integrated_prob
    )

  return(result)
}
