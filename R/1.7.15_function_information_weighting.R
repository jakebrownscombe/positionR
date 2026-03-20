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

  dt <- data.table::as.data.table(detection_probs)

  # Create temporal grouping column based on aggregation method
  if (temporal_grouping == "none") {
    dt[, temporal_group := 1L]
  } else {
    sample_time <- dt$time_period[1]
    if (temporal_grouping == "day") {
      if (inherits(sample_time, c("Date", "POSIXt"))) {
        dt[, temporal_group := as.Date(time_period)]
      } else if (is.numeric(sample_time)) {
        if (time_aggregation == "seconds") {
          dt[, temporal_group := floor(time_period / 86400)]
        } else if (sample_time > 10000) {
          dt[, temporal_group := as.Date(time_period, origin = "1970-01-01")]
        } else {
          dt[, temporal_group := time_period]
        }
      } else {
        dt[, temporal_group := as.Date(time_period)]
      }
    } else if (temporal_grouping == "hour") {
      if (inherits(sample_time, "POSIXt")) {
        dt[, temporal_group := lubridate::floor_date(time_period, "hour")]
      } else if (is.numeric(sample_time)) {
        if (time_aggregation == "seconds") {
          dt[, temporal_group := floor(time_period / 3600)]
        } else {
          dt[, temporal_group := time_period]
        }
      }
    }
  }

  # Calculate station characteristics by temporal group
  station_chars <- dt[, {
    thresh <- stats::quantile(DE_pred, percentile_cutoff, na.rm = TRUE)
    eff_de <- mean(DE_pred[DE_pred >= thresh], na.rm = TRUE)
    max_de <- max(DE_pred, na.rm = TRUE)
    mean_de <- mean(DE_pred, na.rm = TRUE)
    med_de <- stats::median(DE_pred, na.rm = TRUE)
    min_de <- min(DE_pred, na.rm = TRUE)
    .(n_cells = .N,
      top_percentile_threshold = thresh,
      station_effective_DE = eff_de,
      station_max_DE = max_de,
      station_mean_DE = mean_de,
      station_median_DE = med_de,
      station_min_DE = min_de,
      station_range_DE = max_de - min_de)
  }, by = .(station_id, temporal_group)]

  # Handle NA values in effective DE
  station_chars[is.na(station_effective_DE) | !is.finite(station_effective_DE),
    station_effective_DE := station_mean_DE]

  return(as.data.frame(station_chars))
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

  dt <- data.table::as.data.table(probs_data)
  dt_sc <- data.table::as.data.table(station_characteristics)

  # Add temporal grouping to match station characteristics
  if ("temporal_group" %in% names(dt_sc) && "time_period" %in% names(dt)) {
    sample_station_temporal <- dt_sc$temporal_group[1]
    sample_data_time <- dt$time_period[1]

    if (is.na(sample_data_time) || is.null(sample_data_time) || nrow(dt) == 0) {
      warning("Invalid or empty time_period data - skipping temporal grouping")
      if (inherits(sample_station_temporal, "Date")) {
        dt[, temporal_group := as.Date("1970-01-01")]
      } else {
        dt[, temporal_group := 1L]
      }
    } else if (inherits(sample_station_temporal, "Date")) {
      if (inherits(sample_data_time, c("Date", "POSIXt"))) {
        dt[, temporal_group := as.Date(time_period)]
      } else if (is.numeric(sample_data_time)) {
        dt[, temporal_group := as.Date(time_period, origin = "1970-01-01")]
      }
    } else if (is.numeric(sample_station_temporal)) {
      if (inherits(sample_data_time, c("Date", "POSIXt"))) {
        dt[, temporal_group := as.numeric(as.Date(time_period))]
      } else if (!is.na(sample_data_time) && sample_data_time > 1000000) {
        dt[, temporal_group := floor(time_period / 86400)]
      } else {
        dt[, temporal_group := time_period]
      }
    } else {
      dt[, temporal_group := time_period]
    }
  }

  # Join with station characteristics
  join_cols <- if ("temporal_group" %in% names(dt)) c("station_id", "temporal_group") else "station_id"
  dt <- dt_sc[dt, on = join_cols]

  if (weighting_type == "detection") {
    dt[, detection_multiplier := n_detections / (station_effective_DE + 0.01)]
    if (dampening_factor == 1) {
      dt[, detection_multiplier_dampened := detection_multiplier]
    } else if (dampening_factor == 0.5) {
      dt[, detection_multiplier_dampened := sqrt(detection_multiplier)]
    } else {
      dt[, detection_multiplier_dampened := detection_multiplier^dampening_factor]
    }
    dt[, `:=`(DE_pred_weighted = detection_multiplier_dampened * DE_pred,
              DE_pred_normalized = detection_multiplier_dampened * DE_pred)]
  } else {
    dt[, `:=`(DE_pred_weighted = DE_pred,
              DE_pred_normalized = DE_pred)]
  }

  # Add diagnostic information
  dt[, `:=`(weighting_method_used = weighting_type,
            station_info_available = !is.na(station_effective_DE))]

  return(as.data.frame(dt))
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