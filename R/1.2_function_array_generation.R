#' Generate exact number of regular points within spatial boundaries
#'
#' Creates a specified number of regularly spaced points within the boundaries
#' of a raster or polygon object. Points are arranged in a systematic grid pattern.
#'
#' @param input_obj A spatial object. Can be a RasterLayer, RasterBrick, RasterStack,
#'   sf/sfc object, or SpatialPolygons/SpatialPolygonsDataFrame.
#' @param n_points Integer. The exact number of points to generate.
#' @param seed Numeric. Optional random seed for reproducible results. Default is NULL.
#'
#' @return An sf object containing the generated points with columns:
#'   \item{x}{X coordinates}
#'   \item{y}{Y coordinates}
#'   \item{point_id}{Unique identifier for each point}
#'   \item{raster_value}{Raster values at point locations (only if input was raster)}
#'   The object also has a "method_info" attribute describing the spacing pattern.
#'
#' @details
#' This function uses sf::st_sample with type = "regular" to create evenly
#' distributed points. If the input is a raster, values are extracted at point
#' locations. The approximate spacing between points is calculated and stored
#' as an attribute.
#'
#' @examples
#' \dontrun{
#' # Generate 50 regular points in a raster
#' r <- raster::raster(matrix(1:100, 10, 10))
#' points <- generate_exact_regular_points(r, n_points = 50, seed = 123)
#'
#' # Generate points in a polygon
#' poly <- sf::st_sfc(sf::st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0)))))
#' points <- generate_exact_regular_points(poly, n_points = 25)
#' }
#'
#' @export
generate_exact_regular_points <- function(input_obj, n_points, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Prepare the polygon
  prep <- prepare_polygon(input_obj)
  polygon_sf <- prep$polygon

  # Generate exact number of regular points
  regular_points <- sf::st_sample(polygon_sf, size = n_points, type = "regular")
  regular_points_sf <- sf::st_sf(geometry = regular_points)

  # Extract coordinates
  coords <- sf::st_coordinates(regular_points_sf)
  regular_points_sf$x <- coords[,1]
  regular_points_sf$y <- coords[,2]

  # Add point ID
  regular_points_sf$point_id <- 1:nrow(regular_points_sf)

  # Extract values if input was raster
  if (prep$type == "raster") {
    points_sp <- methods::as(regular_points_sf, "Spatial")
    raster_values <- raster::extract(prep$original, points_sp)
    regular_points_sf$raster_value <- raster_values
  }

  # Calculate approximate spacing between points
  if (nrow(regular_points_sf) > 1) {
    # Calculate distances between adjacent points
    coords_matrix <- sf::st_coordinates(regular_points_sf)
    distances <- sqrt(diff(coords_matrix[,1])^2 + diff(coords_matrix[,2])^2)
    # Use median distance as representative spacing
    approx_spacing <- round(median(distances[distances > 0]), 0)

    # Add method info as attribute
    attr(regular_points_sf, "method_info") <- paste0(" (regular pattern, ~", approx_spacing, "m spacing)")
  } else {
    attr(regular_points_sf, "method_info") <- " (regular pattern)"
  }

  return(regular_points_sf)
}







#' Prepare polygon from various spatial input types
#'
#' This helper function converts different spatial object types (raster, sf, sp)
#' into a standardized sf polygon format for further processing. For rasters,
#' it creates a polygon only from areas with actual data (non-NA values).
#'
#' @param input_obj A spatial object. Can be a RasterLayer, RasterBrick, RasterStack,
#'   sf/sfc object, or SpatialPolygons/SpatialPolygonsDataFrame.
#'
#' @return A list with three elements:
#'   \item{polygon}{An sf polygon object representing the data boundary}
#'   \item{type}{Character string indicating input type ("raster" or "polygon")}
#'   \item{original}{The original input object}
#'
#' @details
#' For raster inputs, the function creates a mask of non-NA values and converts
#' only those areas to polygons using dissolve = TRUE. This ensures that points
#' are only generated within areas that contain actual data.
#'
#' @keywords internal
prepare_polygon <- function(input_obj) {
  if (inherits(input_obj, c("RasterLayer", "RasterBrick", "RasterStack"))) {
    # Convert raster to polygon - but only areas with actual data
    # First, create a mask of non-NA values
    raster_mask <- !is.na(input_obj)

    # Convert mask to polygon (this gives us only the data area)
    raster_poly <- raster::rasterToPolygons(raster_mask, fun = function(x) x == 1, dissolve = TRUE)
    polygon_sf <- sf::st_as_sf(raster_poly)

    # Keep only geometry to avoid attribute warnings
    polygon_sf <- sf::st_geometry(polygon_sf)
    polygon_sf <- sf::st_sf(geometry = polygon_sf)

    return(list(polygon = polygon_sf, type = "raster", original = input_obj))
  } else if (inherits(input_obj, c("sf", "sfc"))) {
    # Already a polygon
    return(list(polygon = input_obj, type = "polygon", original = input_obj))
  } else if (inherits(input_obj, c("SpatialPolygons", "SpatialPolygonsDataFrame"))) {
    # Convert sp to sf
    polygon_sf <- sf::st_as_sf(input_obj)
    return(list(polygon = polygon_sf, type = "polygon", original = input_obj))
  } else {
    stop("Input must be a raster (RasterLayer) or polygon (sf, sp)")
  }
}


#' Generate points with specified spacing within spatial boundaries
#'
#' Creates a grid of points with a user-defined spacing within the boundaries
#' of a raster or polygon object.
#'
#' @param input_obj A spatial object. Can be a RasterLayer, RasterBrick, RasterStack,
#'   sf/sfc object, or SpatialPolygons/SpatialPolygonsDataFrame.
#' @param spacing Numeric. Distance between points in map units (typically meters).
#' @param seed Numeric. Optional random seed for reproducible results. Default is NULL.
#'
#' @return An sf object containing the generated points with columns:
#'   \item{x}{X coordinates}
#'   \item{y}{Y coordinates}
#'   \item{point_id}{Unique identifier for each point}
#'   \item{raster_value}{Raster values at point locations (only if input was raster)}
#'   The object also has a "method_info" attribute describing the spacing.
#'
#' @details
#' This function creates a regular grid based on the extent of the input object,
#' then keeps only points that fall within the boundaries. The grid is offset
#' by spacing/2 from the extent edges to ensure even distribution.
#'
#' @examples
#' \dontrun{
#' # Generate points every 500 meters
#' r <- raster::raster(matrix(1:100, 10, 10))
#' points <- generate_spaced_points(r, spacing = 500, seed = 123)
#'
#' # Generate points in a polygon with 100m spacing
#' poly <- sf::st_sfc(sf::st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0)))))
#' points <- generate_spaced_points(poly, spacing = 100)
#' }
#'
#' @export
generate_spaced_points <- function(input_obj, spacing, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  # Prepare the polygon
  prep <- prepare_polygon(input_obj)
  polygon_sf <- prep$polygon
  # Get extent
  if (prep$type == "raster") {
    ext <- raster::extent(prep$original)
  } else {
    bbox <- sf::st_bbox(polygon_sf)
    ext <- raster::extent(bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"])
  }
  # Generate regular grid with specified spacing
  x_seq <- seq(ext@xmin + spacing/2, ext@xmax - spacing/2, by = spacing)
  y_seq <- seq(ext@ymin + spacing/2, ext@ymax - spacing/2, by = spacing)
  # Create all combinations
  grid_points <- expand.grid(x = x_seq, y = y_seq)
  # Convert to sf points
  grid_sf <- sf::st_as_sf(grid_points, coords = c("x", "y"), crs = sf::st_crs(polygon_sf))
  # Keep only points within the polygon (now this is the actual data area!)
  points_within <- sf::st_intersection(grid_sf, polygon_sf)
  # Extract coordinates
  coords <- sf::st_coordinates(points_within)
  points_within$x <- coords[,1]
  points_within$y <- coords[,2]
  # Add point ID
  points_within$point_id <- 1:nrow(points_within)
  # Extract values if input was raster
  if (prep$type == "raster") {
    points_sp <- methods::as(points_within, "Spatial")
    raster_values <- raster::extract(prep$original, points_sp)
    points_within$raster_value <- raster_values
  }
  # Add method info as attribute
  attr(points_within, "method_info") <- paste0(" (", spacing, "m spacing)")
  return(points_within)
}




#' Generate random points within spatial boundaries
#'
#' Creates a specified number of randomly distributed points within the boundaries
#' of a raster or polygon object.
#'
#' @param input_obj A spatial object. Can be a RasterLayer, RasterBrick, RasterStack,
#'   sf/sfc object, or SpatialPolygons/SpatialPolygonsDataFrame.
#' @param n_points Integer. The exact number of points to generate.
#' @param seed Numeric. Optional random seed for reproducible results. Default is NULL.
#'
#' @return An sf object containing the generated points with columns:
#'   \item{x}{X coordinates}
#'   \item{y}{Y coordinates}
#'   \item{point_id}{Unique identifier for each point}
#'   \item{raster_value}{Raster values at point locations (only if input was raster)}
#'   The object also has a "method_info" attribute indicating random distribution.
#'
#' @details
#' This function uses sf::st_sample with type = "random" to create randomly
#' distributed points. If the input is a raster, values are extracted at point
#' locations.
#'
#' @examples
#' \dontrun{
#' # Generate 100 random points in a raster
#' r <- raster::raster(matrix(1:100, 10, 10))
#' points <- generate_random_points(r, n_points = 100, seed = 456)
#'
#' # Generate random points in a polygon
#' poly <- sf::st_sfc(sf::st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0)))))
#' points <- generate_random_points(poly, n_points = 50)
#' }
#'
#' @export
generate_random_points <- function(input_obj, n_points, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Prepare the polygon
  prep <- prepare_polygon(input_obj)
  polygon_sf <- prep$polygon

  # Generate exact number of random points
  random_points <- sf::st_sample(polygon_sf, size = n_points, type = "random")
  random_points_sf <- sf::st_sf(geometry = random_points)

  # Extract coordinates
  coords <- sf::st_coordinates(random_points_sf)
  random_points_sf$x <- coords[,1]
  random_points_sf$y <- coords[,2]

  # Add point ID
  random_points_sf$point_id <- 1:nrow(random_points_sf)

  # Extract values if input was raster
  if (prep$type == "raster") {
    points_sp <- methods::as(random_points_sf, "Spatial")
    raster_values <- raster::extract(prep$original, points_sp)
    random_points_sf$raster_value <- raster_values
  }

  # Add method info as attribute
  attr(random_points_sf, "method_info") <- " (random)"

  return(random_points_sf)
}

#' Plot generated points on spatial input object
#'
#' Creates a visualization showing the generated points overlaid on the original
#' spatial object (raster or polygon). Handles different input types automatically.
#'
#' @param input_obj A spatial object. Can be a RasterLayer, RasterBrick, RasterStack,
#'   sf/sfc object, or SpatialPolygons/SpatialPolygonsDataFrame.
#' @param points_sf An sf object containing the points to plot, typically generated
#'   by one of the point generation functions.
#'
#' @return A ggplot2 object showing the spatial input with points overlaid.
#'   For rasters, displays as a heatmap with color scale. For polygons, displays
#'   as filled shapes with points.
#'
#' @details
#' The function automatically detects whether the input is a raster or polygon
#' and creates an appropriate visualization. Point count and generation method
#' information are included in the plot subtitle if available from the points_sf
#' attributes.
#'
#' @examples
#' \dontrun{
#' # Plot points on a raster
#' r <- raster::raster(matrix(1:100, 10, 10))
#' points <- generate_random_points(r, n_points = 50)
#' plot_points_on_input(r, points)
#'
#' # Plot points on a polygon
#' poly <- sf::st_sfc(sf::st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0)))))
#' points <- generate_spaced_points(poly, spacing = 0.2)
#' plot_points_on_input(poly, points)
#' }
#'
#' @export
plot_points_on_input <- function(input_obj, points_sf) {
  # Prepare the input
  prep <- prepare_polygon(input_obj)
  # Get number of points
  n_points <- nrow(points_sf)
  # Get method info from attributes
  method_info <- attr(points_sf, "method_info")
  if (is.null(method_info)) method_info <- ""
  # Create subtitle with point count and method info
  subtitle <- paste0(n_points, " points generated", method_info)

  if (prep$type == "raster") {
    # Plot raster with points
    raster_df <- raster::as.data.frame(prep$original, xy = TRUE)

    # Get raster extent for proper scaling
    raster_extent <- raster::extent(prep$original)

    p <- ggplot2::ggplot() +
      ggplot2::geom_raster(data = raster_df, ggplot2::aes(x = x, y = y, fill = layer)) +
      ggplot2::geom_sf(data = points_sf, color = "green", size = 2) +
      ggplot2::scale_fill_gradient(low = "blue4", high = "cornflowerblue", na.value = "transparent") +
      ggplot2::coord_sf(xlim = c(raster_extent@xmin, raster_extent@xmax),
                        ylim = c(raster_extent@ymin, raster_extent@ymax),
                        expand = FALSE) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Generated Points on Raster",
                    subtitle = subtitle,
                    fill = "Depth (m)")
  } else {
    # Plot polygon with points
    p <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = prep$polygon, fill = "lightblue", alpha = 0.5, color = "darkblue") +
      ggplot2::geom_sf(data = points_sf, color = "green", size = 2) +
      ggplot2::coord_sf(expand = FALSE) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Generated Points on Polygon",
                    subtitle = subtitle)
  }
  return(p)
}
