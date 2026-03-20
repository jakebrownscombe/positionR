#' Check if a straight line between two points crosses NA cells (barriers)
#'
#' @param x1 X coordinate of starting point
#' @param y1 Y coordinate of starting point
#' @param x2 X coordinate of ending point
#' @param y2 Y coordinate of ending point
#' @param raster RasterLayer to check for barriers
#' @param n_samples Number of points to sample along the line (default 50)
#'
#' @return Logical value: TRUE if line crosses any NA cells, FALSE otherwise
#'
#' @keywords internal
check_line_crosses_barrier <- function(x1, y1, x2, y2, raster, n_samples = 50) {
  # Create sequence of points along the line
  t_seq <- seq(0, 1, length.out = n_samples)

  # Linear interpolation between points
  x_points <- x1 + t_seq * (x2 - x1)
  y_points <- y1 + t_seq * (y2 - y1)

  # Create matrix of coordinates
  line_coords <- cbind(x_points, y_points)

  # Extract raster values at these points
  raster_values <- raster::extract(raster, line_coords)

  # Check if any values are NA (barrier encountered)
  any(is.na(raster_values))
}

#' Vectorized barrier crossing check for all cells from a single station
#'
#' Checks barrier crossings for multiple cell endpoints simultaneously using a
#' single raster::extract() call instead of one per cell.
#'
#' @param x1 Numeric. X coordinate of station
#' @param y1 Numeric. Y coordinate of station
#' @param x2 Numeric vector. X coordinates of all cells
#' @param y2 Numeric vector. Y coordinates of all cells
#' @param raster RasterLayer to check for barriers
#' @param n_samples Number of sample points per line (default 50)
#' @param na_mask Logical vector. TRUE for cells to skip (NA distance)
#'
#' @return Logical vector: TRUE if line crosses barrier, NA if masked
#' @keywords internal
check_line_crosses_barrier_vectorized <- function(x1, y1, x2, y2, raster,
                                                   n_samples = 50, na_mask = NULL) {
  n_cells <- length(x2)
  t_seq <- seq(0, 1, length.out = n_samples)

  # Determine which cells to actually check
  if (!is.null(na_mask)) {
    check_idx <- which(!na_mask)
  } else {
    check_idx <- seq_len(n_cells)
  }

  result <- rep(NA, n_cells)

  if (length(check_idx) == 0) return(result)

  # Generate all sample points for all lines at once
  # Each line gets n_samples points; total = length(check_idx) * n_samples
  x2_check <- x2[check_idx]
  y2_check <- y2[check_idx]

  # Vectorized interpolation: outer product approach
  # t_seq is length n_samples, x2_check is length n_check
  x_all <- rep(x1, length(check_idx) * n_samples) +
    rep(t_seq, length(check_idx)) * rep(x2_check - x1, each = n_samples)
  y_all <- rep(y1, length(check_idx) * n_samples) +
    rep(t_seq, length(check_idx)) * rep(y2_check - y1, each = n_samples)

  # Single raster::extract() call for all points
  all_coords <- cbind(x_all, y_all)
  all_values <- raster::extract(raster, all_coords)

  # Reshape: each line is n_samples consecutive values
  # Check if any value in each line's group is NA
  na_flags <- is.na(all_values)
  line_idx <- rep(seq_along(check_idx), each = n_samples)
  crosses <- tapply(na_flags, line_idx, any)

  result[check_idx] <- as.logical(crosses)
  return(result)
}

#' Calculate cost distances from receiver stations to all raster cells (Enhanced)
#'
#' Enhanced version that supports character station IDs and flexible column naming.
#' Computes both cost-weighted and straight-line distances from each receiver
#' station to all valid cells in a raster. Uses hybrid distance calculation that
#' employs straight-line distance in open water and least-cost paths around barriers,
#' eliminating grid artifacts while preserving meaningful tortuosity.
#'
#' @param raster A RasterLayer object representing the study area. Non-NA cells
#'   are treated as valid locations for distance calculations.
#' @param receiver_frame An sf object or data frame containing receiver station locations.
#'   Must have a column for station identification (default is 'point_id').
#' @param max_distance Numeric. Maximum distance (in map units) for calculations.
#'   Distances beyond this threshold are set to NA. Default is NULL (no limit).
#' @param station_col Character. Name of the column containing station identifiers.
#'   Default is "point_id". Can be any column name with unique station IDs.
#'
#' @return A data frame in long format with the following columns:
#'   \item{cell_id}{Unique identifier for each raster cell}
#'   \item{x}{X coordinate of the cell center}
#'   \item{y}{Y coordinate of the cell center}
#'   \item{raster_value}{Original value from the input raster}
#'   \item{station_no}{Station identifier from receiver_frame (preserves original type)}
#'   \item{cost_distance}{Hybrid distance: straight-line when no barriers present,
#'     least-cost when barriers encountered. Eliminates grid artifacts in open water.}
#'   \item{straight_distance}{Euclidean distance from station to cell}
#'   \item{tortuosity}{Ratio of cost distance to straight distance. Values ~1.0
#'     indicate open water paths, values >1.0 indicate barrier navigation.}
#'   \item{crosses_barrier}{Logical flag indicating whether the straight-line path
#'     between station and cell crosses NA cells (barriers)}
#'
#' @details
#' This enhanced version:
#' \itemize{
#'   \item Supports both numeric and character station IDs
#'   \item Allows flexible column naming via station_col parameter
#'   \item Preserves the original data type of station identifiers
#'   \item Maintains compatibility with existing workflows
#' }
#'
#' The function uses a hybrid distance calculation approach:
#' \enumerate{
#'   \item Creates a uniform cost surface from the raster (all valid cells = 1)
#'   \item Builds a transition matrix with 8-directional connectivity
#'   \item Calculates accumulated cost distances using accCost()
#'   \item Computes straight-line (Euclidean) distances
#'   \item Detects barrier crossings via ray-casting along straight-line paths
#'   \item Uses straight distance when no barriers present, cost distance when barriers encountered
#'   \item Converts results to long format for analysis
#' }
#'
#' The hybrid approach eliminates radial grid artifacts in open water while preserving
#' meaningful tortuosity around actual barriers. The output \code{cost_distance} column
#' contains straight-line distance for open water paths and least-cost distance for
#' barrier-crossing paths. Tortuosity values ~1.0 indicate open water, while values >1.0
#' indicate navigation around barriers.
#'
#' @examples
#' \dontrun{
#' # With numeric station IDs (default)
#' stations_numeric <- generate_random_points(depth_raster, n_points = 5, seed = 123)
#' distances_numeric <- calculate_station_distances(depth_raster, stations_numeric)
#'
#' # With character station IDs
#' stations_char <- stoney_rx_deploy  # Has character station_id column
#' distances_char <- calculate_station_distances(
#'   raster = depth_raster,
#'   receiver_frame = stations_char,
#'   max_distance = 3000,
#'   station_col = "station_id"
#' )
#'
#' # Verify station ID preservation
#' unique(distances_char$station_no)  # Should show character IDs
#' }
#'
#' @seealso \code{\link{generate_random_points}}, \code{\link{generate_spaced_points}}
#'
#' @importFrom stats na.omit
#' @export
calculate_station_distances <- function(raster,
                                        receiver_frame,
                                        max_distance = NULL,
                                        station_col = "point_id") {

  # Validate inputs
  if (!station_col %in% names(receiver_frame)) {
    stop("Column '", station_col, "' not found in receiver_frame. ",
         "Available columns: ", paste(names(receiver_frame), collapse = ", "))
  }
  
  # Check for unique station IDs
  station_ids <- receiver_frame[[station_col]]
  if (any(duplicated(station_ids))) {
    stop("Station IDs in column '", station_col, "' must be unique. ",
         "Found duplicates: ", paste(station_ids[duplicated(station_ids)], collapse = ", "))
  }
  
  # Create uniform cost surface
  cost_raster <- raster
  cost_raster[!is.na(cost_raster)] <- 1

  # Create transition matrix
  tr <- gdistance::transition(cost_raster, transitionFunction = mean, directions = 8)
  tr <- gdistance::geoCorrection(tr, type = "c")

  # Convert points to SpatialPoints if needed
  if ("sf" %in% class(receiver_frame)) {
    points_sp <- methods::as(receiver_frame, "Spatial")
  } else {
    # Assume it's a data frame with x,y coordinates
    coordinates_cols <- c("x", "y")
    if (!all(coordinates_cols %in% names(receiver_frame))) {
      # Try alternative column names
      if (all(c("station_x", "station_y") %in% names(receiver_frame))) {
        coordinates_cols <- c("station_x", "station_y")
      } else if (all(c("lon", "lat") %in% names(receiver_frame))) {
        coordinates_cols <- c("lon", "lat")
      } else {
        stop("Could not find coordinate columns. Expected 'x','y' or 'station_x','station_y' or 'lon','lat'")
      }
    }
    points_sp <- sp::SpatialPointsDataFrame(
      coords = receiver_frame[, coordinates_cols],
      data = receiver_frame,
      proj4string = sp::CRS(raster::projection(raster))
    )
  }

  # Get raster cell coordinates for the dataframe
  raster_coords <- raster::coordinates(cost_raster)
  cell_numbers <- 1:raster::ncell(cost_raster)
  valid_cells <- which(!is.na(raster::values(cost_raster)))

  # Initialize results dataframe
  results_df <- data.frame(
    cell_id = cell_numbers[valid_cells],
    x = raster_coords[valid_cells, 1],
    y = raster_coords[valid_cells, 2],
    raster_value = raster::values(raster)[valid_cells]  # Use original raster values
  )

  # Calculate cost distances from each station
  for (i in 1:nrow(receiver_frame)) {
    current_station_id <- station_ids[i]
    cat("\rCalculating distances for station", i, "of", nrow(receiver_frame),
        "(ID:", current_station_id, ")    ")
    flush.console()

    # Get single point
    single_point <- points_sp[i, ]
    station_coords <- raster::coordinates(single_point)

    # Calculate cost distance from this point using accCost
    cost_dist_raster <- gdistance::accCost(tr, single_point)

    # Extract distances for valid cells only
    cost_distances <- raster::values(cost_dist_raster)[valid_cells]

    # Calculate straight-line distances
    cell_coords <- raster_coords[valid_cells, ]
    straight_distances <- sqrt((cell_coords[,1] - station_coords[1])^2 +
                                 (cell_coords[,2] - station_coords[2])^2)

    # Replace infinite values with NA
    cost_distances[is.infinite(cost_distances)] <- NA

    # Apply maximum distance filter if specified
    if (!is.null(max_distance)) {
      cost_distances[cost_distances > max_distance] <- NA
      straight_distances[straight_distances > max_distance] <- NA
    }

    # Detect barrier crossings for all cells at once (vectorized)
    crosses_barrier <- check_line_crosses_barrier_vectorized(
      x1 = station_coords[1],
      y1 = station_coords[2],
      x2 = cell_coords[, 1],
      y2 = cell_coords[, 2],
      raster = raster,
      na_mask = is.na(straight_distances)
    )

    # Calculate hybrid distance: straight when possible, cost when barrier present
    hybrid_distances <- ifelse(
      crosses_barrier,
      cost_distances,      # Use cost distance when barrier encountered
      straight_distances   # Use straight distance in open water
    )

    # Add to results dataframe with sanitized column names
    # Replace any non-standard characters in station ID for column naming
    safe_station_id <- gsub("[^[:alnum:]_]", "_", as.character(current_station_id))
    cost_col_name <- paste0("cost_dist_station_", safe_station_id)
    straight_col_name <- paste0("straight_dist_station_", safe_station_id)
    barrier_col_name <- paste0("crosses_barrier_station_", safe_station_id)

    results_df[[cost_col_name]] <- hybrid_distances  # Now contains hybrid distance
    results_df[[straight_col_name]] <- straight_distances
    results_df[[barrier_col_name]] <- crosses_barrier
  }

  # Convert to long format using data.table::melt
  cat("\n")  # End the progress line

  dt_results <- data.table::as.data.table(results_df)

  # Identify column groups
  cost_cols <- grep("^cost_dist_station_", names(dt_results), value = TRUE)
  straight_cols <- grep("^straight_dist_station_", names(dt_results), value = TRUE)
  barrier_cols <- grep("^crosses_barrier_station_", names(dt_results), value = TRUE)

  # Create station ID mapping
  safe_names <- gsub("[^[:alnum:]_]", "_", as.character(station_ids))
  station_map <- data.table::data.table(
    safe_name = safe_names,
    station_no = station_ids
  )

  # Melt all three column groups separately
  id_cols <- c("cell_id", "x", "y", "raster_value")

  cost_long <- data.table::melt(dt_results, id.vars = id_cols,
    measure.vars = cost_cols, variable.name = "station_col", value.name = "cost_distance",
    variable.factor = FALSE)
  cost_long[, safe_name := gsub("cost_dist_station_", "", station_col)]
  cost_long <- station_map[cost_long, on = .(safe_name)]
  cost_long[, c("station_col", "safe_name") := NULL]

  straight_long <- data.table::melt(dt_results, id.vars = "cell_id",
    measure.vars = straight_cols, variable.name = "station_col", value.name = "straight_distance",
    variable.factor = FALSE)
  straight_long[, safe_name := gsub("straight_dist_station_", "", station_col)]
  straight_long <- station_map[straight_long, on = .(safe_name)]
  straight_long[, c("station_col", "safe_name") := NULL]

  barrier_long <- data.table::melt(dt_results, id.vars = "cell_id",
    measure.vars = barrier_cols, variable.name = "station_col", value.name = "crosses_barrier",
    variable.factor = FALSE)
  barrier_long[, safe_name := gsub("crosses_barrier_station_", "", station_col)]
  barrier_long <- station_map[barrier_long, on = .(safe_name)]
  barrier_long[, c("station_col", "safe_name") := NULL]

  # Combine all components
  final_dt <- cost_long[straight_long, on = .(cell_id, station_no)]
  final_dt <- final_dt[barrier_long, on = .(cell_id, station_no)]
  final_dt[, tortuosity := cost_distance / straight_distance]
  final_dt <- final_dt[!is.na(cost_distance)]
  final_dt <- final_dt[, .(cell_id, x, y, raster_value, station_no, cost_distance,
                            straight_distance, tortuosity, crosses_barrier)]
  data.table::setorderv(final_dt, c("station_no", "cell_id"))

  final_df <- as.data.frame(final_dt)

  cat("Distance calculations complete!\n")
  cat("Result contains", nrow(final_df), "station-cell combinations\n")
  cat("Station IDs (", class(final_df$station_no), "):",
      paste(head(unique(final_df$station_no), 5), collapse = ", "),
      if(length(unique(final_df$station_no)) > 5) "..." else "", "\n")

  return(final_df)
}