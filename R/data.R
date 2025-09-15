#' Fish Detection Data from Stoney Lake
#'
#' A dataset containing acoustic telemetry detection data for walleye from 
#' Stoney Lake, Ontario, Canada. This dataset includes detections from tagged 
#' fish recorded by an acoustic receiver array.
#'
#' @format A data frame with 125,538 rows and 4 columns:
#' \describe{
#'   \item{species}{Character. Species name ("Walleye")}
#'   \item{fish_id}{Character. Unique identifier for individual fish (e.g., "Walleye-1512985")}
#'   \item{station_id}{Character. Acoustic receiver station identifier (e.g., "T4-1")}
#'   \item{detection_timestamp_utc}{POSIXct. Date and time of detection in UTC}
#' }
#'
#' @details
#' This dataset represents real acoustic telemetry data collected as part of 
#' fish movement studies in Stoney Lake. The data includes detections from 
#' multiple individual walleye across various receiver stations over the 
#' study period (2023). Each row represents a single detection event.
#'
#' @source Field data collected from Stoney Lake, Ontario, Canada
#' @name stoney_fish_detections
#' @docType data
#' @keywords datasets
#' @examples
#' \dontrun{
#' # Load the dataset
#' data("stoney_fish_detections")
#' 
#' # View dataset structure
#' str(stoney_fish_detections)
#' 
#' # Number of unique fish
#' length(unique(stoney_fish_detections$fish_id))
#' 
#' # Number of receiver stations
#' length(unique(stoney_fish_detections$station_id))
#' 
#' # Date range of detections
#' range(stoney_fish_detections$detection_timestamp_utc)
#' }
"stoney_fish_detections"

#' Acoustic Receiver Deployment Data from Stoney Lake
#'
#' A spatial dataset (sf object) containing information about acoustic receiver 
#' deployments in Stoney Lake, Ontario, Canada. This dataset includes receiver 
#' locations, deployment periods, and site characteristics.
#'
#' @format An sf data frame with 145 rows and 7 columns:
#' \describe{
#'   \item{station_id}{Character. Unique identifier for receiver station (e.g., "T1-1")}
#'   \item{x}{Numeric. Longitude coordinate (decimal degrees, WGS84)}
#'   \item{y}{Numeric. Latitude coordinate (decimal degrees, WGS84)}
#'   \item{depth_m}{Numeric. Water depth at receiver location (meters)}
#'   \item{deploy_datetime_UTC}{POSIXct. Date and time of receiver deployment (UTC)}
#'   \item{recover_datetime_UTC}{POSIXct. Date and time of receiver recovery (UTC)}
#'   \item{geometry}{sf geometry column with UTM Zone 17N coordinates}
#' }
#'
#' @details
#' This dataset contains spatial and temporal information for acoustic receivers 
#' deployed in Stoney Lake as part of fish telemetry studies. The receivers were 
#' deployed across multiple years (2022-2024) to monitor fish movements. 
#' The coordinate reference system is UTM Zone 17N (EPSG:32617).
#' 
#' Each receiver has coordinates in both geographic (decimal degrees) and 
#' projected (UTM) coordinate systems. Deployment periods vary by station, 
#' with some receivers being redeployed multiple times.
#'
#' @source Field data collected from Stoney Lake, Ontario, Canada
#' @name stoney_rx_deploy
#' @docType data
#' @keywords datasets spatial
#' @examples
#' \dontrun{
#' # Load the dataset
#' data("stoney_rx_deploy")
#' 
#' # View dataset structure
#' str(stoney_rx_deploy)
#' 
#' # Plot receiver locations
#' library(ggplot2)
#' library(sf)
#' ggplot(stoney_rx_deploy) +
#'   geom_sf(aes(color = depth_m)) +
#'   scale_color_viridis_c(name = "Depth (m)") +
#'   labs(title = "Acoustic Receiver Locations in Stoney Lake") +
#'   theme_minimal()
#' 
#' # Summary of deployment periods
#' stoney_rx_deploy$deployment_days <- as.numeric(
#'   stoney_rx_deploy$recover_datetime_UTC - stoney_rx_deploy$deploy_datetime_UTC
#' )
#' summary(stoney_rx_deploy$deployment_days)
#' 
#' # Depth distribution
#' hist(stoney_rx_deploy$depth_m, main = "Receiver Depth Distribution", 
#'      xlab = "Depth (m)")
#' }
"stoney_rx_deploy"