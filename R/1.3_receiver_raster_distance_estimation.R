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

#' Calculate cost distances from receiver stations to all raster cells
#'
#' Computes both cost-weighted and straight-line distances from each receiver
#' station to all valid cells in a raster. Uses hybrid distance calculation that
#' employs straight-line distance in open water and least-cost paths around barriers,
#' eliminating grid artifacts while preserving meaningful tortuosity.
#'
#' @param raster A RasterLayer object representing the study area. Non-NA cells
#'   are treated as valid locations for distance calculations.
#' @param receiver_frame An sf object containing receiver station locations.
#'   Must have a column for station identification (default is 'point_id').
#' @param max_distance Numeric. Maximum distance (in map units) for calculations.
#'   Distances beyond this threshold are set to NA. Default is NULL (no limit).
#' @param station_col Character. Name of the column containing station identifiers.
#'   Default is "point_id". Supports both numeric and character station IDs.
#'
#' @return A data frame in long format with the following columns:
#'   \item{cell_id}{Unique identifier for each raster cell}
#'   \item{x}{X coordinate of the cell center}
#'   \item{y}{Y coordinate of the cell center}
#'   \item{raster_value}{Original value from the input raster}
#'   \item{station_no}{Station identifier from receiver_frame$point_id}
#'   \item{cost_distance}{Hybrid distance: straight-line when no barriers present,
#'     least-cost when barriers encountered. Eliminates grid artifacts in open water.}
#'   \item{straight_distance}{Euclidean distance from station to cell}
#'   \item{tortuosity}{Ratio of cost distance to straight distance. Values ~1.0
#'     indicate open water paths, values >1.0 indicate barrier navigation.}
#'   \item{crosses_barrier}{Logical flag indicating whether the straight-line path
#'     between station and cell crosses NA cells (barriers)}
#'
#' @details
#' This function uses a hybrid distance calculation approach:
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
#' # Generate receiver stations with numeric IDs (default)
#' stations <- generate_random_points(depth_raster, n_points = 5, seed = 123)
#'
#' # Calculate distances with default point_id column
#' distances <- calculate_station_distances(depth_raster, stations)
#'
#' # Calculate distances with 1000m maximum
#' distances_limited <- calculate_station_distances(
#'   raster = depth_raster,
#'   receiver_frame = stations,
#'   max_distance = 1000
#' )
#'
#' # Use with character station IDs
#' # Assuming stoney_rx_deploy has a 'station_id' column with character IDs
#' distances_char <- calculate_station_distances(
#'   raster = depth_raster,
#'   receiver_frame = stoney_rx_deploy,
#'   max_distance = 3000,
#'   station_col = "station_id"  # Specify the column with station IDs
#' )
#'
#' # Analyze tortuosity patterns
#' library(dplyr)
#' tortuosity_summary <- distances %>%
#'   group_by(station_no) %>%
#'   summarise(
#'     mean_tortuosity = mean(tortuosity, na.rm = TRUE),
#'     max_cost_dist = max(cost_distance, na.rm = TRUE)
#'   )
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
  
  # Get station IDs from specified column
  station_ids <- receiver_frame[[station_col]]
  if (any(duplicated(station_ids))) {
    stop("Station IDs in column '", station_col, "' must be unique. ",
         "Found duplicates: ", paste(station_ids[duplicated(station_ids)], collapse = ", "))
  }

  # Create uniform cost surface
  cost_raster <- raster
  cost_raster[!is.na(cost_raster)] <- 1

  # Create transition matrix
  cat("Creating transition matrix...\n")
  tr <- gdistance::transition(cost_raster, transitionFunction = mean, directions = 8)
  tr <- gdistance::geoCorrection(tr, type = "c")

  # Convert points to SpatialPoints
  points_sp <- methods::as(receiver_frame, "Spatial")

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
    cat("Calculating distances for station", i, "of", nrow(receiver_frame),
        "(ID:", current_station_id, ")\n")

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

    # Detect barrier crossings for each cell
    cat("  Detecting barrier crossings...\n")
    crosses_barrier <- logical(length(valid_cells))
    for (j in seq_along(valid_cells)) {
      # Skip if distance is NA
      if (is.na(straight_distances[j])) {
        crosses_barrier[j] <- NA
        next
      }

      crosses_barrier[j] <- check_line_crosses_barrier(
        x1 = station_coords[1],
        y1 = station_coords[2],
        x2 = cell_coords[j, 1],
        y2 = cell_coords[j, 2],
        raster = raster
      )
    }

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

  # Convert to long format
  cat("Converting to long format...\n")

  # Separate cost, straight distance, and barrier columns
  cost_cols <- grep("^cost_dist_station_", names(results_df), value = TRUE)
  straight_cols <- grep("^straight_dist_station_", names(results_df), value = TRUE)
  barrier_cols <- grep("^crosses_barrier_station_", names(results_df), value = TRUE)

  # Create a mapping of safe column names back to original station IDs
  station_mapping <- data.frame(
    safe_name = gsub("[^[:alnum:]_]", "_", as.character(station_ids)),
    original_id = station_ids,
    stringsAsFactors = FALSE
  )

  # Pivot cost distances to long format (now contains hybrid distances)
  cost_long <- results_df %>%
    dplyr::select(cell_id, x, y, raster_value, dplyr::all_of(cost_cols)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(cost_cols),
                        names_to = "station_col",
                        values_to = "cost_distance") %>%
    dplyr::mutate(safe_name = gsub("cost_dist_station_", "", station_col)) %>%
    dplyr::left_join(station_mapping, by = "safe_name") %>%
    dplyr::rename(station_no = original_id) %>%
    dplyr::select(-station_col, -safe_name)

  # Pivot straight distances to long format
  straight_long <- results_df %>%
    dplyr::select(cell_id, dplyr::all_of(straight_cols)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(straight_cols),
                        names_to = "station_col",
                        values_to = "straight_distance") %>%
    dplyr::mutate(safe_name = gsub("straight_dist_station_", "", station_col)) %>%
    dplyr::left_join(station_mapping, by = "safe_name") %>%
    dplyr::rename(station_no = original_id) %>%
    dplyr::select(-station_col, -safe_name)

  # Pivot barrier flags to long format
  barrier_long <- results_df %>%
    dplyr::select(cell_id, dplyr::all_of(barrier_cols)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(barrier_cols),
                        names_to = "station_col",
                        values_to = "crosses_barrier") %>%
    dplyr::mutate(safe_name = gsub("crosses_barrier_station_", "", station_col)) %>%
    dplyr::left_join(station_mapping, by = "safe_name") %>%
    dplyr::rename(station_no = original_id) %>%
    dplyr::select(-station_col, -safe_name)

  # Combine all components
  final_df <- cost_long %>%
    dplyr::left_join(straight_long, by = c("cell_id", "station_no")) %>%
    dplyr::left_join(barrier_long, by = c("cell_id", "station_no")) %>%
    dplyr::mutate(tortuosity = cost_distance / straight_distance) %>%
    dplyr::filter(!is.na(cost_distance)) %>%
    dplyr::select(cell_id, x, y, raster_value, station_no, cost_distance,
                  straight_distance, tortuosity, crosses_barrier) %>%
    dplyr::arrange(station_no, cell_id)

  cat("Distance calculations complete!\n")
  cat("Result contains", nrow(final_df), "station-cell combinations\n")
  cat("Station IDs (", class(final_df$station_no), "):", 
      paste(head(unique(final_df$station_no), 5), collapse = ", "),
      if(length(unique(final_df$station_no)) > 5) "..." else "", "\n")

  return(final_df)
}
