#' Launch Interactive Array Design Shiny Application
#'
#' Opens the Array DesignR Shiny application for interactive acoustic receiver
#' array design and optimization. The app provides tools for receiver placement,
#' detection range modeling, and array performance analysis.
#'
#' @details
#' The Array DesignR application includes four main modules:
#' \itemize{
#'   \item \strong{Array Design}: Interactive map for placing receivers manually
#'     (click-to-add) or generating regular arrays. Load custom raster files
#'     (.tif, .grd, .asc) or use the included sample depth data.
#'   \item \strong{Range Modeling}: Fit detection efficiency models using
#'     \code{\link{create_logistic_curve_depth}} with depth-dependent parameters
#'     (d50 and d95 values).
#'   \item \strong{Array Performance}: Quantitative assessment of detection
#'     coverage, sampling bias analysis, and performance metrics based on
#'     detection efficiency thresholds.
#'   \item \strong{Instructions}: Comprehensive user guide with methodology,
#'     best practices, and troubleshooting tips.
#' }
#'
#' Key features:
#' \itemize{
#'   \item Load bathymetry or habitat rasters for spatial context
#'   \item Place receivers manually by clicking or generate regular arrays
#'   \item Calculate cost-weighted distances to assess coverage
#'   \item Visualize detection efficiency surfaces using fitted models
#'   \item Export receiver coordinates as CSV for field deployment
#'   \item Analyze depth sampling bias and coverage thresholds
#' }
#'
#' @return Launches the Shiny application in the default web browser.
#'   No return value.
#'
#' @examples
#' \dontrun{
#' # Launch the array design application
#' run_array_design()
#' }
#'
#' @seealso \code{\link{generate_regular_points}},
#'   \code{\link{create_logistic_curve_depth}},
#'   \code{\link{calculate_station_distances}}
#'
#' @export
run_array_design <- function() {
  app_dir <- system.file("shiny-apps/array_design", package = "positionR")

  if (app_dir == "") {
    stop("Could not find array_design app. Please reinstall positionR package.")
  }

  # Check for required Shiny packages
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required to run this application. Install it with: install.packages('shiny')")
  }

  if (!requireNamespace("shinydashboard", quietly = TRUE)) {
    stop("Package 'shinydashboard' is required to run this application. Install it with: install.packages('shinydashboard')")
  }

  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required to run this application. Install it with: install.packages('plotly')")
  }

  if (!requireNamespace("DT", quietly = TRUE)) {
    stop("Package 'DT' is required to run this application. Install it with: install.packages('DT')")
  }

  shiny::runApp(app_dir, display.mode = "normal")
}
