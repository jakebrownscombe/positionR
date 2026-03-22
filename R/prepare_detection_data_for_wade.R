#' Prepare real fish detection data for WADE positioning
#'
#' Converts sparse real fish detection data into the time-aggregated format
#' required by calculate_fish_positions function. Creates detection summaries
#' by time period and validates station deployment consistency.
#'
#' @param fish_detections Data frame containing fish detection records with columns:
#'   fish_id, station_id, detection_timestamp_utc, and optionally species.
#' @param station_deployments Data frame or sf object containing station deployment
#'   information with columns: station_id, x, y, deploy_datetime_UTC, recover_datetime_UTC.
#' @param selected_fish_id Character. ID of the fish to process. Must exist in
#'   fish_detections$fish_id.
#' @param time_aggregation Character. Time period for aggregating detections.
#'   Options: "hour", "day", "week", "month". Default is "day".
#' @param start_time POSIXct. Start time for the analysis period. If NULL,
#'   uses the first detection time for the selected fish.
#' @param end_time POSIXct. End time for the analysis period. If NULL,
#'   uses the last detection time for the selected fish.
#'
#' @return A list containing:
#'   \item{station_detections}{Data frame with columns: path_id, datetime,
#'     station_id, detected, detection_prob, distance_to_station, station_x,
#'     station_y, n_detections. Compatible with calculate_fish_positions function.}
#'   \item{station_info}{Data frame with columns: station_id, x, y, start_date,
#'     end_date, depth_m (if available). Station information for active stations
#'     during analysis period, including deployment dates for temporal filtering.}
#'   \item{time_periods}{Data frame with time period information: time_period,
#'     time_period_label, n_detections.}
#'   \item{deployment_warnings}{Character vector of any deployment warnings.}
#'   \item{summary_stats}{List with analysis summary statistics.}
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Filters detections for the selected fish
#'   \item Determines analysis time period (start_time to end_time)
#'   \item Identifies stations active during the analysis period
#'   \item Warns if any stations are deployed/recovered during analysis
#'   \item Aggregates detections by time period (hour/day/week/month)
#'   \item Creates detection summary with presence/absence only for station-time
#'     combinations where the station was deployed (filters by deployment dates)
#'   \item Formats output compatible with calculate_fish_positions
#' }
#'
#' Note: The function now filters station-time combinations to only include periods
#' when each station was actually deployed. This prevents creation of invalid
#' combinations and reduces memory usage, particularly important for large datasets
#' or long analysis periods.
#'
#' The station_detections output structure:
#' - path_id: Always 1 (single fish track)
#' - datetime: Time period timestamp
#' - station_id: Station identifier (character, will be converted by WADE function)
#' - detected: 1 if fish detected at this station during time period, 0 otherwise
#' - detection_prob: 1.0 for detections, 0 for non-detections (field data certainty)
#' - distance_to_station: 0 for detections (fish at station), NA for non-detections
#' - station_x, station_y: Station coordinates
#' - n_detections: Number of detections in this time period at this station
#'
#' @examples
#' \dontrun{
#' # Prepare data for a specific fish with daily aggregation
#' wade_data <- prepare_detection_data_for_wade(
#'   fish_detections = stoney_fish_detections,
#'   station_deployments = stoney_rx_deploy,
#'   selected_fish_id = "Walleye-1512985",
#'   time_aggregation = "day"
#' )
#'
#' # Check for deployment warnings
#' if (length(wade_data$deployment_warnings) > 0) {
#'   cat("Deployment warnings:\n")
#'   cat(paste(wade_data$deployment_warnings, collapse = "\n"))
#' }
#'
#' # Use with WADE positioning
#' results <- calculate_fish_positions(
#'   station_detections = wade_data$station_detections,
#'   station_distances_df = station_distances,
#'   station_info = wade_data$station_info,
#'   # ... other parameters
#' )
#'
#' # Hourly aggregation for fine-scale analysis
#' wade_data_hourly <- prepare_detection_data_for_wade(
#'   fish_detections = stoney_fish_detections,
#'   station_deployments = stoney_rx_deploy,
#'   selected_fish_id = "Walleye-1512985",
#'   time_aggregation = "hour",
#'   start_time = as.POSIXct("2023-10-15 00:00:00", tz = "UTC"),
#'   end_time = as.POSIXct("2023-10-17 00:00:00", tz = "UTC")
#' )
#' }
#'
#' @seealso \code{\link{calculate_fish_positions}}
#'
#' @export
prepare_detection_data_for_wade <- function(fish_detections,
                                           station_deployments,
                                           selected_fish_id,
                                           time_aggregation = "day",
                                           start_time = NULL,
                                           end_time = NULL) {
  
  # Load required libraries
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for this function")
  }
  if (!requireNamespace("lubridate", quietly = TRUE)) {
    stop("Package 'lubridate' is required for this function")
  }
  
  # Validate inputs
  required_fish_cols <- c("fish_id", "station_id", "detection_timestamp_utc")
  missing_fish_cols <- setdiff(required_fish_cols, names(fish_detections))
  if (length(missing_fish_cols) > 0) {
    stop("Missing columns in fish_detections: ", paste(missing_fish_cols, collapse = ", "))
  }
  
  required_station_cols <- c("station_id", "x", "y", "deploy_datetime_UTC", "recover_datetime_UTC")
  missing_station_cols <- setdiff(required_station_cols, names(station_deployments))
  if (length(missing_station_cols) > 0) {
    stop("Missing columns in station_deployments: ", paste(missing_station_cols, collapse = ", "))
  }
  
  # Validate time_aggregation parameter
  valid_aggregations <- c("hour", "day", "week", "month")
  if (!time_aggregation %in% valid_aggregations) {
    stop("time_aggregation must be one of: ", paste(valid_aggregations, collapse = ", "))
  }
  
  # Check if selected fish exists
  if (!selected_fish_id %in% fish_detections$fish_id) {
    stop("Fish ID '", selected_fish_id, "' not found in fish_detections. ",
         "Available fish IDs: ", paste(head(unique(fish_detections$fish_id), 5), collapse = ", "))
  }
  
  cat("Preparing WADE detection data for fish:", selected_fish_id, "\n")
  
  # Filter detections for selected fish
  fish_data <- fish_detections %>%
    dplyr::filter(fish_id == selected_fish_id) %>%
    dplyr::arrange(detection_timestamp_utc)
  
  if (nrow(fish_data) == 0) {
    stop("No detections found for fish ID: ", selected_fish_id)
  }
  
  cat("Found", nrow(fish_data), "detections for", selected_fish_id, "\n")
  
  # Determine analysis time period
  if (is.null(start_time)) {
    start_time <- min(fish_data$detection_timestamp_utc, na.rm = TRUE)
    cat("Using first detection as start time:", format(start_time), "\n")
  }
  
  if (is.null(end_time)) {
    end_time <- max(fish_data$detection_timestamp_utc, na.rm = TRUE)
    cat("Using last detection as end time:", format(end_time), "\n")
  }
  
  # Ensure start_time < end_time
  if (start_time >= end_time) {
    stop("start_time must be before end_time")
  }
  
  analysis_duration <- as.numeric(difftime(end_time, start_time, units = "days"))
  cat("Analysis period:", format(start_time), "to", format(end_time), 
      "(", round(analysis_duration, 2), "days)\n")
  
  # Filter fish detections to analysis period
  fish_data <- fish_data %>%
    dplyr::filter(detection_timestamp_utc >= start_time & 
                  detection_timestamp_utc <= end_time)
  
  cat("Detections within analysis period:", nrow(fish_data), "\n")
  
  # Identify stations active during analysis period
  active_stations <- station_deployments %>%
    dplyr::filter(
      deploy_datetime_UTC <= end_time &
      recover_datetime_UTC >= start_time
    )
  
  cat("Stations active during analysis period:", nrow(active_stations), "\n")
  
  if (nrow(active_stations) == 0) {
    stop("No stations were active during the analysis period")
  }
  
  # Check for deployment changes during analysis period
  deployment_warnings <- character(0)
  
  # Check for deployments during analysis
  deployed_during <- station_deployments %>%
    dplyr::filter(
      deploy_datetime_UTC > start_time & 
      deploy_datetime_UTC < end_time
    )
  
  if (nrow(deployed_during) > 0) {
    warning_msg <- paste("WARNING:", nrow(deployed_during), 
                        "station(s) deployed during analysis period:",
                        paste(deployed_during$station_id, collapse = ", "))
    deployment_warnings <- c(deployment_warnings, warning_msg)
  }
  
  # Check for recoveries during analysis
  recovered_during <- station_deployments %>%
    dplyr::filter(
      recover_datetime_UTC > start_time & 
      recover_datetime_UTC < end_time
    )
  
  if (nrow(recovered_during) > 0) {
    warning_msg <- paste("WARNING:", nrow(recovered_during), 
                        "station(s) recovered during analysis period:",
                        paste(recovered_during$station_id, collapse = ", "))
    deployment_warnings <- c(deployment_warnings, warning_msg)
  }
  
  # Print deployment warnings
  if (length(deployment_warnings) > 0) {
    cat("\n=== DEPLOYMENT WARNINGS ===\n")
    for (warning in deployment_warnings) {
      cat(warning, "\n")
    }
    cat("Consider selecting a period with stable station deployment.\n\n")
  }
  
  # Create station info dataframe with deployment dates
  station_info <- active_stations %>%
    dplyr::mutate(
      start_date = as.Date(deploy_datetime_UTC),
      end_date = as.Date(recover_datetime_UTC)
    ) %>%
    dplyr::select(station_id, x, y, start_date, end_date, dplyr::any_of("depth_m")) %>%
    dplyr::arrange(station_id)

  # Convert to regular dataframe if sf object
  if ("sf" %in% class(station_info)) {
    station_info <- sf::st_drop_geometry(station_info)
  }
  
  # Create time-aggregated detection summary
  cat("Aggregating detections by", time_aggregation, "...\n")
  
  # Add time period column to fish data
  fish_data_aggregated <- fish_data %>%
    dplyr::mutate(
      time_period = lubridate::floor_date(detection_timestamp_utc, unit = time_aggregation)
    )
  
  # Summarize detections by time period and station
  detection_summary <- fish_data_aggregated %>%
    dplyr::group_by(time_period, station_id) %>%
    dplyr::summarise(
      n_detections = n(),
      first_detection = min(detection_timestamp_utc),
      last_detection = max(detection_timestamp_utc),
      .groups = "drop"
    )
  
  # Get all time periods in the analysis
  time_periods <- seq(
    from = lubridate::floor_date(start_time, unit = time_aggregation),
    to = lubridate::floor_date(end_time, unit = time_aggregation),
    by = time_aggregation
  )

  # Create only valid station-date combinations (where station was deployed)
  # This filters to only periods when each station was actually in the water
  all_combinations <- station_info %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      time_period = list(time_periods[
        as.Date(time_periods) >= start_date &
        as.Date(time_periods) <= end_date
      ])
    ) %>%
    tidyr::unnest(time_period) %>%
    dplyr::select(time_period, station_id) %>%
    dplyr::ungroup()
  
  # Join with detection summary and station coordinates
  station_detections <- all_combinations %>%
    dplyr::left_join(detection_summary, by = c("time_period", "station_id")) %>%
    dplyr::left_join(station_info, by = "station_id") %>%
    dplyr::mutate(
      # WADE-specific columns
      path_id = selected_fish_id,  # Use original fish ID
      datetime = time_period,  # Time period timestamp
      detected = ifelse(is.na(n_detections), 0, 1),  # 1 if detected, 0 if not
      detection_prob = ifelse(detected == 1, 1.0, 0.0),  # Field data certainty
      distance_to_station = ifelse(detected == 1, 0, NA_real_),  # At station if detected
      station_x = x,  # Station coordinates
      station_y = y,
      n_detections = ifelse(is.na(n_detections), 0, n_detections)  # Replace NA with 0
    ) %>%
    dplyr::select(path_id, datetime, station_id, detected, detection_prob, 
                  distance_to_station, station_x, station_y, n_detections) %>%
    dplyr::arrange(datetime, station_id)
  
  # Create time period summary
  time_period_summary <- station_detections %>%
    dplyr::group_by(datetime) %>%
    dplyr::summarise(
      time_period_label = format(datetime, "%Y-%m-%d"),
      n_detections = sum(n_detections),
      n_stations_detected = sum(detected),
      .groups = "drop"
    ) %>%
    dplyr::rename(time_period = datetime)
  
  cat("Created", nrow(time_period_summary), "time periods with", time_aggregation, "aggregation\n")
  
  # Calculate summary statistics
  total_detections <- sum(station_detections$n_detections)
  n_time_periods <- length(unique(station_detections$datetime))
  n_stations <- length(unique(station_detections$station_id))
  detection_rate <- sum(station_detections$detected) / nrow(station_detections)
  
  summary_stats <- list(
    fish_id = selected_fish_id,
    analysis_start = start_time,
    analysis_end = end_time,
    analysis_duration_days = analysis_duration,
    time_aggregation = time_aggregation,
    n_time_periods = n_time_periods,
    n_active_stations = n_stations,
    total_detections = total_detections,
    detection_rate = detection_rate,
    original_detections = nrow(fish_data)
  )
  
  # Print summary
  cat("\n=== WADE PREPARATION SUMMARY ===\n")
  cat("Fish ID:", summary_stats$fish_id, "\n")
  cat("Time period:", format(summary_stats$analysis_start), "to", format(summary_stats$analysis_end), "\n")
  cat("Duration:", round(summary_stats$analysis_duration_days, 2), "days\n")
  cat("Time aggregation:", summary_stats$time_aggregation, "\n")
  cat("Time periods:", summary_stats$n_time_periods, "\n")
  cat("Active stations:", summary_stats$n_active_stations, "\n")
  cat("Detection matrix rows:", nrow(station_detections), "\n")
  cat("Total detections:", summary_stats$total_detections, "from", summary_stats$original_detections, "original\n")
  cat("Detection rate:", round(summary_stats$detection_rate * 100, 2), "%\n")
  
  # Return results
  results <- list(
    station_detections = station_detections,
    station_info = station_info,
    time_periods = time_period_summary,
    deployment_warnings = deployment_warnings,
    summary_stats = summary_stats
  )
  
  return(results)
}