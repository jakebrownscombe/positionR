#' Prepare polygon from various spatial input types
#'
#' This helper function converts different spatial object types (raster, sf, sp)
#' into a standardized sf polygon format for further processing.
#'
#' @param input_obj A spatial object. Can be a RasterLayer, RasterBrick, RasterStack,
#'   sf/sfc object, or SpatialPolygons/SpatialPolygonsDataFrame.
#'
#' @return A list with three elements:
#'   \item{polygon}{An sf polygon object}
#'   \item{type}{Character string indicating input type ("raster" or "polygon")}
#'   \item{original}{The original input object}
#'
#' @details
#' For raster inputs, the function converts to polygons using dissolve = TRUE.
#' For polygon inputs (sf, sp), it ensures they are in sf format.
#'
#' @examples
#' \dontrun{
#' # With a raster
#' r <- raster::raster(matrix(1:100, 10, 10))
#' result <- prepare_polygon(r)
#'
#' # With an sf polygon
#' poly <- sf::st_sfc(sf::st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0)))))
#' result <- prepare_polygon(poly)
#' }
#'
#' @keywords internal
prepare_polygon <- function(input_obj) {
  if (inherits(input_obj, c("RasterLayer", "RasterBrick", "RasterStack"))) {
    # Convert raster to polygon
    raster_poly <- raster::rasterToPolygons(input_obj, dissolve = TRUE)
    polygon_sf <- sf::st_as_sf(raster_poly)
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
