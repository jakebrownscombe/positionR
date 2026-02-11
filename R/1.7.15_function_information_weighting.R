#' Calculate Station Detection Efficiency Characteristics
#'
#' Calculates robust station DE characteristics using the top percentile of DE values
#' for each station and temporal group. This provides a more stable estimate of station
#' capability than using mean or max values.
#'
#' @param detection_probs Data frame containing detection probabilities with columns
#'   for station_id, DE_pred, and time_period
#' @param percentile_cutoff Numeric between 0 and 1. The percentile cutoff for
#'   calculating robust station DE (default 0.95 for top 95% of values)
#' @param temporal_grouping Character. How to group temporally: "day", "hour", or "none"
#' @param time_aggregation Character. The time aggregation method being used
#'   ("seconds", "hour", "day", "month")
#'
#' @return Data frame with station characteristics including effective DE by temporal group
#' @export
calculate_station_de_characteristics <- function(detection_probs, 
                                                percentile_cutoff = 0.95,
                                                temporal_grouping = "day",
                                                time_aggregation = "seconds") {
  
  # Validate inputs
  if (percentile_cutoff <= 0 || percentile_cutoff > 1) {
    stop("percentile_cutoff must be between 0 and 1")
  }
  
  # Create temporal grouping column based on aggregation method
  if (temporal_grouping == "none") {
    # No temporal grouping - calculate single value per station
    detection_probs <- detection_probs %>%
      dplyr::mutate(temporal_group = 1)
  } else {
    # Create temporal grouping based on the time_period format
    # time_period comes from the main function's aggregation
    
    # Detect the format of time_period
    sample_time <- detection_probs$time_period[1]
    
    if (temporal_grouping == "day") {
      if (inherits(sample_time, c("Date", "POSIXt"))) {
        # Already date/datetime format
        detection_probs <- detection_probs %>%
          dplyr::mutate(temporal_group = as.Date(time_period))
      } else if (is.numeric(sample_time)) {
        if (time_aggregation == "seconds") {
          # Numeric seconds - convert to days
          detection_probs <- detection_probs %>%
            dplyr::mutate(temporal_group = floor(time_period / 86400))
        } else {
          # Numeric days or other - use as-is or convert with origin
          if (sample_time > 10000) {
            # Likely days since epoch
            detection_probs <- detection_probs %>%
              dplyr::mutate(temporal_group = as.Date(time_period, origin = "1970-01-01"))
          } else {
            # Small numbers - probably relative days
            detection_probs <- detection_probs %>%
              dplyr::mutate(temporal_group = time_period)
          }
        }
      } else {
        # Character or other - try to convert
        detection_probs <- detection_probs %>%
          dplyr::mutate(temporal_group = as.Date(time_period))
      }
    } else if (temporal_grouping == "hour") {
      if (inherits(sample_time, "POSIXt")) {
        detection_probs <- detection_probs %>%
          dplyr::mutate(temporal_group = lubridate::floor_date(time_period, "hour"))
      } else if (is.numeric(sample_time)) {
        if (time_aggregation == "seconds") {
          detection_probs <- detection_probs %>%
            dplyr::mutate(temporal_group = floor(time_period / 3600))
        } else {
          detection_probs <- detection_probs %>%
            dplyr::mutate(temporal_group = time_period)
        }
      }
    }
  }
  
  # Calculate station characteristics by temporal group
  station_chars <- detection_probs %>%
    dplyr::group_by(station_id, temporal_group) %>%
    dplyr::summarise(
      # Count of cells for this station
      n_cells = dplyr::n(),
      
      # Robust station DE using top percentile
      top_percentile_threshold = stats::quantile(DE_pred, percentile_cutoff, na.rm = TRUE),
      station_effective_DE = mean(DE_pred[DE_pred >= top_percentile_threshold], na.rm = TRUE),
      
      # Alternative metrics for diagnostics
      station_max_DE = max(DE_pred, na.rm = TRUE),
      station_mean_DE = mean(DE_pred, na.rm = TRUE),
      station_median_DE = stats::median(DE_pred, na.rm = TRUE),
      station_min_DE = min(DE_pred, na.rm = TRUE),
      station_range_DE = station_max_DE - station_min_DE,
      
      .groups = "drop"
    ) %>%
    # Handle any NA values in effective DE
    dplyr::mutate(
      station_effective_DE = ifelse(is.na(station_effective_DE) | !is.finite(station_effective_DE),
                                    station_mean_DE,  # Fallback to mean
                                    station_effective_DE)
    )
  
  return(station_chars)
}

#' Apply Information-Theoretic Weighting
#'
#' Applies information-theoretic weighting to detection or non-detection data.
#' For detections, weights are adjusted based on the station's effective DE to
#' account for station capability differences. For non-detections, raw DE values
#' are used as they already encode the strength of absence evidence.
#'
#' @param probs_data Data frame containing probability data with columns for
#'   station_id, DE_pred, and optionally n_detections
#' @param station_characteristics Data frame with station DE characteristics
#'   from calculate_station_de_characteristics()
#' @param weighting_type Character. Either "detection" or "non_detection" to
#'   specify the type of weighting to apply
#' @param dampening_factor Numeric. Optional dampening factor for extreme weights.
#'   Use sqrt for square root dampening (default 1 for no dampening)
#'
#' @return Data frame with added weighting columns
#' @export
apply_information_weighting <- function(probs_data, 
                                      station_characteristics,
                                      weighting_type = "detection",
                                      dampening_factor = 1) {
  
  # Validate weighting type
  if (!weighting_type %in% c("detection", "non_detection")) {
    stop("weighting_type must be 'detection' or 'non_detection'")
  }
  
  # Add temporal grouping to match station characteristics
  # Detect temporal grouping type from station_characteristics
  if ("temporal_group" %in% names(station_characteristics)) {
    if ("time_period" %in% names(probs_data)) {
      
      # Get the type of temporal_group from station_characteristics
      sample_station_temporal <- station_characteristics$temporal_group[1]
      sample_data_time <- probs_data$time_period[1]
      
      # Check for valid data
      if (is.na(sample_data_time) || is.null(sample_data_time) || nrow(probs_data) == 0) {
        # For empty or invalid data, match the station_characteristics temporal_group type
        warning("Invalid or empty time_period data - skipping temporal grouping")
        if (inherits(sample_station_temporal, "Date")) {
          # Match Date type
          probs_data <- probs_data %>% dplyr::mutate(temporal_group = as.Date("1970-01-01"))
        } else {
          # Match numeric type
          probs_data <- probs_data %>% dplyr::mutate(temporal_group = 1)
        }
      } else {
      
      if (inherits(sample_station_temporal, "Date")) {
        # Station characteristics use Date - convert probs_data to Date
        if (inherits(sample_data_time, c("Date", "POSIXt"))) {
          probs_data <- probs_data %>%
            dplyr::mutate(temporal_group = as.Date(time_period))
        } else if (is.numeric(sample_data_time)) {
          if (!is.na(sample_data_time) && sample_data_time > 10000) {
            # Large number - likely days since epoch
            probs_data <- probs_data %>%
              dplyr::mutate(temporal_group = as.Date(time_period, origin = "1970-01-01"))
          } else {
            # Small number - relative days, convert to Date
            probs_data <- probs_data %>%
              dplyr::mutate(temporal_group = as.Date(time_period, origin = "1970-01-01"))
          }
        }
      } else if (is.numeric(sample_station_temporal)) {
        # Station characteristics use numeric - convert probs_data to numeric
        if (inherits(sample_data_time, c("Date", "POSIXt"))) {
          probs_data <- probs_data %>%
            dplyr::mutate(temporal_group = as.numeric(as.Date(time_period)))
        } else {
          # Already numeric - use directly or convert
          if (!is.na(sample_data_time) && sample_data_time > 1000000) {
            # Likely seconds - convert to days
            probs_data <- probs_data %>%
              dplyr::mutate(temporal_group = floor(time_period / 86400))
          } else {
            # Already in right scale
            probs_data <- probs_data %>%
              dplyr::mutate(temporal_group = time_period)
          }
        }
      } else {
        # Other type - try to match
        probs_data <- probs_data %>%
          dplyr::mutate(temporal_group = time_period)
      }
      } # Close the else block for valid data
    }
  }
  
  # Join with station characteristics
  weighted_data <- probs_data %>%
    dplyr::left_join(station_characteristics, 
                     by = if("temporal_group" %in% names(probs_data)) 
                       c("station_id", "temporal_group") else "station_id")
  
  if (weighting_type == "detection") {
    # Information-theoretic weighting for detections
    weighted_data <- weighted_data %>%
      dplyr::mutate(
        # Detection multiplier based on station's effective DE
        # Low-DE stations get higher multipliers (more surprising detections)
        detection_multiplier = n_detections / (station_effective_DE + 0.01),
        
        # Apply optional dampening to prevent extreme weights
        detection_multiplier_dampened = if (dampening_factor == 1) {
          detection_multiplier
        } else if (dampening_factor == 0.5) {
          sqrt(detection_multiplier)  # Square root dampening
        } else {
          detection_multiplier^dampening_factor  # Custom power dampening
        },
        
        # Final weight: multiplier × spatial DE
        # This preserves within-station spatial gradients
        DE_pred_weighted = detection_multiplier_dampened * DE_pred,
        
        # Keep normalized column name for compatibility with existing code
        DE_pred_normalized = DE_pred_weighted
      )
    
  } else if (weighting_type == "non_detection") {
    # Raw DE weighting for non-detections (no information adjustment)
    # High-DE non-detection = strong absence evidence
    # Low-DE non-detection = weak absence evidence
    weighted_data <- weighted_data %>%
      dplyr::mutate(
        # Use raw DE - already encodes strength of absence evidence
        DE_pred_weighted = DE_pred,
        
        # For compatibility with existing code
        DE_pred_normalized = DE_pred
      )
  }
  
  # Add diagnostic information
  weighted_data <- weighted_data %>%
    dplyr::mutate(
      weighting_method_used = weighting_type,
      station_info_available = !is.na(station_effective_DE)
    )
  
  return(weighted_data)
}

#' Print Station Characteristics Summary
#'
#' Prints a summary of station DE characteristics for diagnostic purposes
#'
#' @param station_chars Data frame from calculate_station_de_characteristics()
#' @param n_stations Number of stations to display (default 10)
#'
#' @export
print_station_characteristics <- function(station_chars, n_stations = 10) {
  
  cat("\n=== Station DE Characteristics Summary ===\n")
  
  # Overall summary
  cat("\nNumber of stations:", length(unique(station_chars$station_id)))
  cat("\nNumber of temporal groups:", length(unique(station_chars$temporal_group)))
  
  # Summary statistics
  cat("\n\nEffective DE Statistics (top 95% mean):")
  cat("\n  Min:", round(min(station_chars$station_effective_DE, na.rm = TRUE), 3))
  cat("\n  Mean:", round(mean(station_chars$station_effective_DE, na.rm = TRUE), 3))
  cat("\n  Median:", round(median(station_chars$station_effective_DE, na.rm = TRUE), 3))
  cat("\n  Max:", round(max(station_chars$station_effective_DE, na.rm = TRUE), 3))
  
  # Station ranking
  cat("\n\nTop", n_stations, "stations by effective DE:\n")
  top_stations <- station_chars %>%
    dplyr::group_by(station_id) %>%
    dplyr::summarise(
      mean_effective_DE = mean(station_effective_DE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(mean_effective_DE)) %>%
    utils::head(n_stations)
  
  print(as.data.frame(top_stations))
  
  cat("\nBottom", n_stations, "stations by effective DE:\n")
  bottom_stations <- station_chars %>%
    dplyr::group_by(station_id) %>%
    dplyr::summarise(
      mean_effective_DE = mean(station_effective_DE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(mean_effective_DE) %>%
    utils::head(n_stations)
  
  print(as.data.frame(bottom_stations))
  
  # Information multiplier range
  cat("\n\nExpected detection multiplier range:")
  cat("\n  High-DE station (0.8):", round(1/0.8, 2), "x per detection")
  cat("\n  Low-DE station (0.2):", round(1/0.2, 2), "x per detection")
  cat("\n  Ratio:", round((1/0.2)/(1/0.8), 2), ":1\n")
}