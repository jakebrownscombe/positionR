#' Calculate cost distances from receiver stations to all raster cells
#'
#' Computes both cost-weighted and straight-line distances from each receiver
#' station to all valid cells in a raster. Uses least-cost path analysis to
#' account for landscape constraints on signal transmission.
#'
#' @param raster A RasterLayer object representing the study area. Non-NA cells
#'   are treated as valid locations for distance calculations.
#' @param receiver_frame An sf object containing receiver station locations.
#'   Must have a 'point_id' column for station identification.
#' @param max_distance Numeric. Maximum distance (in map units) for calculations.
#'   Distances beyond this threshold are set to NA. Default is NULL (no limit).
#'
#' @return A data frame in long format with the following columns:
#'   \item{cell_id}{Unique identifier for each raster cell}
#'   \item{x}{X coordinate of the cell center}
#'   \item{y}{Y coordinate of the cell center}
#'   \item{raster_value}{Original value from the input raster}
#'   \item{station_no}{Station identifier from receiver_frame$point_id}
#'   \item{cost_distance}{Least-cost distance from station to cell}
#'   \item{straight_distance}{Euclidean distance from station to cell}
#'   \item{tortuosity}{Ratio of cost distance to straight distance}
#'
#' @details
#' This function uses the gdistance package to perform least-cost path analysis.
#' The process involves:
#' \enumerate{
#'   \item Creating a uniform cost surface from the raster (all valid cells = 1)
#'   \item Building a transition matrix with 8-directional connectivity
#'   \item Calculating accumulated cost distances using accCost()
#'   \item Computing straight-line distances for comparison
#'   \item Converting results to long format for analysis
#' }
#'
#' The tortuosity metric (cost/straight distance) indicates how much the
#' least-cost path deviates from a straight line, with values > 1 indicating
#' increased path complexity.
#'
#' @examples
#' \dontrun{
#' # Generate receiver stations
#' stations <- generate_random_points(depth_raster, n_points = 5, seed = 123)
#'
#' # Calculate distances with no maximum limit
#' distances <- calculate_station_distances(depth_raster, stations)
#'
#' # Calculate distances with 1000m maximum
#' distances_limited <- calculate_station_distances(
#'   raster = depth_raster,
#'   receiver_frame = stations,
#'   max_distance = 1000
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
                                        max_distance = NULL) {

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
    cat("Calculating distances for station", i, "of", nrow(receiver_frame),
        "(ID:", receiver_frame$point_id[i], ")\n")

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

    # Add to results dataframe
    cost_col_name <- paste0("cost_dist_station_", receiver_frame$point_id[i])
    straight_col_name <- paste0("straight_dist_station_", receiver_frame$point_id[i])

    results_df[[cost_col_name]] <- cost_distances
    results_df[[straight_col_name]] <- straight_distances
  }

  # Convert to long format
  cat("Converting to long format...\n")

  # Separate cost and straight distance columns
  cost_cols <- grep("^cost_dist_station_", names(results_df), value = TRUE)
  straight_cols <- grep("^straight_dist_station_", names(results_df), value = TRUE)

  # Pivot cost distances to long format
  cost_long <- results_df %>%
    dplyr::select(cell_id, x, y, raster_value, dplyr::all_of(cost_cols)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(cost_cols),
                        names_to = "station_col",
                        values_to = "cost_distance") %>%
    dplyr::mutate(station_no = as.numeric(gsub("cost_dist_station_", "", station_col))) %>%
    dplyr::select(-station_col)

  # Pivot straight distances to long format
  straight_long <- results_df %>%
    dplyr::select(cell_id, dplyr::all_of(straight_cols)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(straight_cols),
                        names_to = "station_col",
                        values_to = "straight_distance") %>%
    dplyr::mutate(station_no = as.numeric(gsub("straight_dist_station_", "", station_col))) %>%
    dplyr::select(-station_col)

  # Combine cost and straight distances
  final_df <- cost_long %>%
    dplyr::left_join(straight_long, by = c("cell_id", "station_no")) %>%
    dplyr::mutate(tortuosity = cost_distance / straight_distance) %>%
    dplyr::select(cell_id, x, y, raster_value, station_no, cost_distance, straight_distance, tortuosity) %>%
    dplyr::arrange(station_no, cell_id)

  return(final_df)
}
