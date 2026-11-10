#' positionR: tools for acoustic telemetry studies
#'
#' positionR is a suite of four R packages for acoustic telemetry studies. This
#' package installs them and lists them with [positionR_packages()]. It does
#' not attach them: load the one you need with `library()`.
#'
#' @keywords internal
"_PACKAGE"

# ---- Imports the code does not call ----
# positionR imports the four packages so that installing it installs them.
# R CMD check notes an Imports entry that no code refers to. This function
# names one export of each package, so the check finds every import used.
# Nothing calls it.
suite_imports <- function() {
  list(
    detqc::read_detections,
    detfield::de_field,
    arraydesign::place_receivers,
    arraydesignR::run_app
  )
}
