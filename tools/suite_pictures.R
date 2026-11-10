#!/usr/bin/env Rscript
# ---- Make the four suite pictures of the README ----
# Usage: Rscript tools/suite_pictures.R   (from the package folder)
# The landing page shows one picture from each package of the suite. This
# script makes all four into man/figures from the installed packages and the
# sibling package folders at ../<package>/, so the page never links into
# another package's site:
#   suite_detqc.png        the review card picture of detqc's Get started
#                          article, from ../detqc/vignettes/figures/
#   suite_detfield.png     the result-positions figure of detfield's
#                          positioning article, run from the installed copy
#   suite_arraydesign.png  the coverage-map figure of arraydesign's README,
#                          run from ../arraydesign/README.Rmd
#   suite_arraydesignR.png the Map tab screenshot of the installed arraydesignR
# Each picture is at most 900 px wide and cut to 256 colours. Metadata such
# as the write time is stripped, so a rerun on the same inputs gives the same
# bytes. Rerun after a sibling package changes one of these sources, then
# rebuild the site so docs/reference/figures/ takes the new copies.

# ---- Libraries ----
suppressMessages({
  library(knitr)
  library(ggplot2)
})

# ---- Settings ----
suite_root <- normalizePath("..")
out_dir <- file.path("man", "figures")
width_px <- 900
stopifnot(dir.exists(out_dir))
cat("suite root:", suite_root, "\n")
for (pkg in c("detqc", "detfield", "arraydesign", "arraydesignR")) {
  cat(sprintf("%-12s %s from %s\n", pkg, as.character(utils::packageVersion(pkg)),
              dirname(find.package(pkg))))
}

# ---- Helpers ----
# Cut a picture to the README's width and 256 colours, strip its metadata and
# write it under man/figures.
write_picture <- function(img, target, source) {
  info_in <- magick::image_info(img)
  if (info_in$width > width_px) img <- magick::image_resize(img, paste0(width_px, "x"))
  img <- magick::image_quantize(img, max = 256, colorspace = "sRGB", dither = FALSE)
  img <- magick::image_strip(img)
  out <- file.path(out_dir, target)
  magick::image_write(img, out, format = "png")
  info_out <- magick::image_info(magick::image_read(out))
  cat(sprintf("%-24s from %s: %d x %d px in, %d x %d px out, %s bytes, md5 %s\n",
              target, source, info_in$width, info_in$height, info_out$width,
              info_out$height, format(file.size(out), big.mark = ","),
              unname(tools::md5sum(out))))
  invisible(out)
}

# Run an R Markdown file's code from chunk `from` (the first line when NULL)
# to the end of chunk `to` in a new environment. Returns the environment, with
# the last visible value of the code as `.value`.
run_chunks <- function(rmd, from, to) {
  code_file <- tempfile(fileext = ".R")
  invisible(purl(rmd, output = code_file, documentation = 1, quiet = TRUE))
  lines <- readLines(code_file)
  starts <- grep("^## ----", lines)
  first <- if (is.null(from)) 1L else grep(paste0("^## ----", from, "\\b"), lines)
  last_start <- grep(paste0("^## ----", to, "\\b"), lines)
  stopifnot(length(first) == 1, length(last_start) == 1)
  last <- c(starts[starts > last_start], length(lines) + 1)[1] - 1
  cat(sprintf("  %s: lines %d to %d of %d, %d chunks\n", basename(rmd), first, last,
              length(lines), sum(starts >= first & starts <= last)))
  env <- new.env()
  for (expr in parse(text = lines[first:last])) {
    val <- withVisible(eval(expr, envir = env))
    if (val$visible) env$.value <- val$value
  }
  env
}

# Save a ggplot at a size in inches and a resolution, and read it back as a
# picture. The size in inches sets the text size relative to the panel.
render_plot <- function(p, width_in, height_in, dpi) {
  stopifnot(inherits(p, "ggplot"))
  tmp <- tempfile(fileext = ".png")
  ggsave(tmp, p, width = width_in, height = height_in, dpi = dpi, bg = "white")
  magick::image_read(tmp)
}

# ---- detqc: The Review Card ----
# The card is an HTML widget; detqc's tools/vignette_screenshots.R takes the
# picture and the vignette keeps it under vignettes/figures/.
detqc_src <- file.path("detqc", "vignettes", "figures", "review_card_get_started.png")
write_picture(magick::image_read(file.path(suite_root, detqc_src)), "suite_detqc.png", detqc_src)

# ---- detfield: Positions Of One Animal Against Its Track ----
# The installed copy of the vignette matches the installed detfield it calls.
cat("detfield positioning.Rmd, to chunk result-positions:\n")
pos <- run_chunks(system.file("doc", "positioning.Rmd", package = "detfield"),
                  from = NULL, to = "result-positions")
with(pos, {
  cat("  animal:", animal, "\n")
  cat("  window:", format(window, "%Y-%m-%d %H:%M UTC"), "\n")
  cat("  hours in the window:", nrow(window_hours), "; single-receiver hours:", sum(single), "\n")
  cat("  receivers drawn:", nrow(window_rx), "; heard in the window:",
      sum(window_rx$detections > 0), "; deployments in the study:", nrow(receivers), "\n")
  cat("  median error, hours of 2 or more receivers (m):",
      round(median(window_hours$error_m[!single])), "\n")
  cat("  errors of the single-receiver hours (m):",
      round(sort(window_hours$error_m[single])), "\n")
  cat("  entries in the colour key:", length(.value$scales$get_scales("colour")$breaks), "\n")
})
write_picture(render_plot(pos$.value, 900 / 110, 560 / 110, dpi = 110), "suite_detfield.png",
              "detfield/doc/positioning.Rmd, chunk result-positions")

# ---- arraydesign: Coverage Of 20 Receivers On A Grid ----
# README.Rmd is not installed, so the chunk is read from the sibling folder
# and run against the installed arraydesign and detfield.
cat("arraydesign README.Rmd, chunk coverage-map:\n")
cov <- run_chunks(file.path(suite_root, "arraydesign", "README.Rmd"),
                  from = "coverage-map", to = "coverage-map")
with(cov, {
  p_any <- 1 - apply(1 - field$p, 1, prod)
  cat("  receivers placed:", nrow(receivers), "; water cells:", length(p_any), "\n")
  cat("  cells where at least one receiver hears with p >= 0.5:", sum(p_any >= 0.5), "\n")
  cat("  shape key:", paste(names(.value$scales$get_scales("shape")$palette(1)), collapse = ", "),
      "\n")
})
# The README chunk draws at 1600 / 150 by 7 inches; this keeps that size, at a
# resolution that gives 900 px across.
write_picture(render_plot(cov$.value, 1600 / 150, 7, dpi = 900 / (1600 / 150)),
              "suite_arraydesign.png", "arraydesign/README.Rmd, chunk coverage-map")

# ---- arraydesignR: The Map Tab ----
# arraydesignR's tools/app_shots.mjs takes this screenshot; the package
# installs it with its help figures.
app_map <- system.file("help", "figures", "app_map.png", package = "arraydesignR")
stopifnot(nzchar(app_map))
write_picture(magick::image_read(app_map), "suite_arraydesignR.png",
              "arraydesignR/help/figures/app_map.png (installed)")
