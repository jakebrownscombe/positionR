# ---- The suite ----
# One row per package, in the order of the suite list of every positionR
# README. The roles are that list's wording.
suite_table <- function() {
  data.frame(
    package = c("detqc", "detfield", "arraydesign", "arraydesignR"),
    role = c(
      "read, check and filter detections",
      "detection fields and positions",
      "design and test the performance of tracking systems",
      "arraydesign as a point-and-click application"
    ),
    stringsAsFactors = FALSE
  )
}

# Version of an installed package as a string, or NA when it is not installed.
# The argument is named lib.loc as in utils::packageVersion(), so it is
# exempt from the snake_case rule.
installed_version <- function(package, lib.loc = NULL) { # nolint: object_name_linter.
  path <- system.file(package = package, lib.loc = lib.loc)
  if (!nzchar(path)) return(NA_character_)
  as.character(utils::packageVersion(package, lib.loc = lib.loc))
}

#' The positionR packages
#'
#' Lists the four packages of positionR with what each does, the version
#' installed, and the address of its documentation site.
#'
#' @param lib.loc Library folders to look for the packages in, as in
#'   [utils::packageVersion()]. `NULL`, the default, searches the libraries
#'   of [.libPaths()].
#' @return A data frame of class `positionR_packages` with one row per
#'   package and four columns: `package`; `role`, what the package does;
#'   `version`, the installed version as a string, `NA` when the package is
#'   not installed; and `site`, the URL of its documentation.
#' @examples
#' positionR_packages()
#' @export
# The suite name in the function name is exempt from the snake_case rule.
positionR_packages <- function(lib.loc = NULL) { # nolint: object_name_linter.
  out <- suite_table()
  out$version <- vapply(out$package, installed_version, character(1),
                        lib.loc = lib.loc, USE.NAMES = FALSE)
  out$site <- paste0("https://jakebrownscombe.github.io/", out$package, "/")
  class(out) <- c("positionR_packages", "data.frame")
  out
}

#' @export
print.positionR_packages <- function(x, ...) {
  n_installed <- sum(!is.na(x$version))
  cat(sprintf("positionR: %d packages, %d installed\n", nrow(x),
              n_installed))
  for (i in seq_len(nrow(x))) {
    version <- if (is.na(x$version[i])) "not installed" else x$version[i]
    cat(sprintf("\n%s %s\n  %s\n  %s\n", x$package[i], version, x$role[i],
                x$site[i]))
  }
  invisible(x)
}
