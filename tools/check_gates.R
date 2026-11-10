#!/usr/bin/env Rscript
# ---- Code gates for a positionR package ----
# Usage: Rscript tools/check_gates.R [package-dir] [--tarball=PATH]
#   Without --tarball the script builds one with R CMD build in a temporary
#   directory. Exit status 1 when any gate fails.
#
# Gates every package of the suite passes:
#   function_length  no function definition in R/ longer than 80 lines
#   export_formals   no exported function with more than 6 formals; '...'
#                    counts as one
#   tarball_size     the built tarball is at most 2,000,000 bytes, or the
#                    limit in the package's DESCRIPTION field
#                    Config/positionR/tarball_limit (bytes)
#   example_time     no Rd example takes more than 5 s, elapsed or CPU, or
#                    the seconds in POSITIONR_EXAMPLE_SECONDS
#   pages_cover      every export is on a site page; runs only in a checkout
#                    of the whole suite
# Two more enforce the suite's naming rules:
#   export_masking   no export masks a function in base R or in a package
#                    R attaches by default
#   uk_spelling      no US spelling in a function name or formal in R/
# One more for a Shiny app package:
#   app_size         app_ui() at most 40 lines and app_server() at most 30, so
#                    the modules hold the app; a package without either
#                    function passes with nothing to check
# One more for the figures, which are all drawn with ggplot2:
#   base_graphics    no call to a base-graphics function or a graphics
#                    device in R/, in the Rd examples or in the vignette code
#                    chunks; figures are saved with ggplot2::ggsave()
#
# A gate that cannot run fails and says why, because a gate that did not run
# shows nothing about the package. One exception: pages_cover reads the sites
# of every package of the suite, so on a checkout of one package
# (POSITIONR_LONE_CHECKOUT=1, which a package's CI workflow sets) it is
# reported as not run and does not fail the gates.

# ---- Settings ----
max_function_lines <- 80L
max_export_formals <- 6L
# The default; a package that needs more sets Config/positionR/tarball_limit
# in its DESCRIPTION, so this file stays the same in every package.
max_tarball_bytes <- 2e6
# POSITIONR_EXAMPLE_SECONDS sets another limit for a slower machine, such as a
# shared CI runner, where examples that take under 5 s here took up to 6 s.
max_example_seconds <- as.numeric(Sys.getenv("POSITIONR_EXAMPLE_SECONDS", "5"))
max_app_lines <- c(app_ui = 40L, app_server = 30L)
# POSITIONR_BASE_GRAPHICS=report turns the base_graphics gate into a report
# for a package that draws with base graphics.
base_graphics_mode <- Sys.getenv("POSITIONR_BASE_GRAPHICS", "fail")
device_calls <- c(
  "png", "pdf", "jpeg", "tiff", "bmp", "svg", "cairo_pdf", "cairo_ps", "postscript",
  "dev.off", "dev.new", "dev.set", "dev.cur", "dev.flush", "dev.hold"
)
# A bare plot() is a generic that may dispatch to a method returning a ggplot,
# so plot() is a hit only when called through one of these namespaces.
base_plot_namespaces <- c("graphics", "terra", "sf", "raster")
default_packages <- c("base", "stats", "utils", "graphics", "grDevices", "methods", "datasets")
# Telemetry packages a user is likely to attach beside this one. Their
# exports join the masking pool when they are installed; when one is not,
# the gate says so and does not fail on it.
neighbour_packages <- c("glatos", "actel")

# US forms, by stem (any ending) and by whole word. A name token is a part of
# the name between underscores or dots.
us_stems <- c(
  "normaliz", "optimiz", "summariz", "standardiz", "analyz", "initializ",
  "visualiz", "organiz", "recogniz", "minimiz", "maximiz", "customiz",
  "categoriz", "characteriz", "parameteriz", "randomiz", "regulariz",
  "discretiz", "rasteriz", "vectoriz", "prioritiz", "utiliz", "finaliz",
  "synchroniz", "serializ"
)
us_words <- c(
  "color", "colors", "colored", "center", "centers", "centered", "behavior",
  "behaviors", "behavioral", "gray", "modeling", "modeled", "labeled",
  "labeling", "favor", "neighbor", "neighbors", "meter", "meters", "liter",
  "liters", "catalog", "traveled", "canceled", "fiber", "defense"
)

# ---- Helpers ----
r_exe <- function(name) {
  file.path(R.home("bin"), if (.Platform$OS.type == "windows") paste0(name, ".exe") else name)
}

say <- function(...) cat(sprintf(...), "\n", sep = "")

# Gates that failed because they could not run say so in their detail. The
# summary names them apart, so "run" counts only the gates that tested the
# package.
did_not_run <- new.env(parent = emptyenv())

gate <- function(name, pass, detail, items = character()) {
  blocked <- isFALSE(pass) && grepl("gate did not run", detail, fixed = TRUE)
  label <- if (is.na(pass)) "SKIP" else if (isTRUE(pass)) "PASS" else if (blocked) "BLOCKED" else "FAIL"
  say("%s %-15s %s", label, name, detail)
  for (it in items) say("     %s", it)
  if (blocked) assign(name, TRUE, envir = did_not_run)
  if (is.na(pass)) NA else isTRUE(pass)
}

fmt_int <- function(x) formatC(x, format = "d", big.mark = ",")

# Top-level function definitions in R/: name, file, line span and formals.
function_defs <- function(pkg_dir) {
  files <- list.files(file.path(pkg_dir, "R"), pattern = "\\.[Rr]$", full.names = TRUE)
  rows <- list()
  for (f in files) {
    lines <- readLines(f, warn = FALSE, encoding = "UTF-8")
    exprs <- parse(text = lines, keep.source = TRUE, srcfile = srcfilecopy(f, lines))
    refs <- attr(exprs, "srcref")
    for (i in seq_along(exprs)) {
      e <- exprs[[i]]
      if (!is.call(e) || !is.name(e[[1]])) next
      if (!as.character(e[[1]]) %in% c("<-", "=", "<<-")) next
      rhs <- e[[3]]
      if (!is.call(rhs) || !identical(rhs[[1]], as.name("function"))) next
      rows[[length(rows) + 1]] <- data.frame(
        name = gsub("^[`'\"]|[`'\"]$", "", paste(deparse(e[[2]]), collapse = "")),
        file = file.path("R", basename(f)), first = refs[[i]][1], last = refs[[i]][3],
        n_formals = length(rhs[[2]]), formals = paste(names(rhs[[2]]), collapse = " "),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame(
    name = character(), file = character(), first = integer(), last = integer(),
    n_formals = integer(), formals = character(), stringsAsFactors = FALSE
  )
  out$lines <- out$last - out$first + 1L
  list(defs = out, n_files = length(files))
}

us_tokens <- function(x) {
  toks <- unlist(strsplit(tolower(x), "[._]+"))
  toks <- toks[nzchar(toks)]
  hit <- toks %in% us_words | vapply(toks, function(t) any(startsWith(t, us_stems)), logical(1))
  unique(toks[hit])
}

# ---- Gates on the source ----
gate_function_length <- function(defs, n_files) {
  long <- defs[defs$lines > max_function_lines, ]
  gate(
    "function_length", nrow(long) == 0,
    sprintf(
      "%d of %d functions over %d lines; longest %d (%d files in R/)",
      nrow(long), nrow(defs), max_function_lines, if (nrow(defs)) max(defs$lines) else 0L, n_files
    ),
    sprintf("%s, %s:%d-%d, %d lines", long$name, long$file, long$first, long$last, long$lines)
  )
}

gate_export_formals <- function(defs, exports, patterns) {
  missing <- setdiff(exports, defs$name)
  ex <- defs[defs$name %in% exports, ]
  wide <- ex[ex$n_formals > max_export_formals, ]
  items <- c(
    sprintf("%s has %d formals: %s", wide$name, wide$n_formals, wide$formals),
    sprintf("%s: no `name <- function(...)` definition in R/, formals not counted", missing),
    sprintf("exportPattern(%s) in NAMESPACE: exports cannot be listed", patterns)
  )
  gate(
    "export_formals", nrow(wide) == 0 && !length(missing) && !length(patterns),
    sprintf("%d of %d exports over %d formals", nrow(wide), length(exports), max_export_formals),
    items
  )
}

gate_export_masking <- function(exports) {
  for (p in default_packages) suppressPackageStartupMessages(library(p, character.only = TRUE))
  pool <- lapply(setNames(default_packages, default_packages), function(p) {
    if (p == "base") ls(baseenv(), all.names = TRUE) else ls(paste0("package:", p), all.names = TRUE)
  })
  installed <- neighbour_packages[vapply(neighbour_packages, requireNamespace, logical(1), quietly = TRUE)]
  for (p in installed) pool[[p]] <- getNamespaceExports(p)
  absent <- setdiff(neighbour_packages, installed)
  hits <- vapply(exports, function(x) {
    w <- names(pool)[vapply(pool, function(n) x %in% n, logical(1))]
    paste(w, collapse = ", ")
  }, "")
  masked <- exports[nzchar(hits)]
  gate(
    "export_masking", !length(masked),
    sprintf(
      "%d of %d exports mask one of the %s names in %s (R %s)%s",
      length(masked), length(exports), fmt_int(sum(lengths(pool))),
      paste(names(pool), collapse = ", "), getRversion(),
      if (length(absent)) sprintf("; absent from the library and so unchecked: %s", paste(absent, collapse = ", ")) else ""
    ),
    sprintf("%s masks %s::%s", masked, hits[nzchar(hits)], masked)
  )
}

gate_uk_spelling <- function(defs) {
  items <- character()
  n_formals <- 0L
  for (i in seq_len(nrow(defs))) {
    bad_name <- us_tokens(defs$name[i])
    fmls <- strsplit(defs$formals[i], " ", fixed = TRUE)[[1]]
    n_formals <- n_formals + length(fmls)
    bad_formals <- fmls[vapply(fmls, function(f) length(us_tokens(f)) > 0, logical(1))]
    if (length(bad_name)) {
      items <- c(items, sprintf("%s (%s:%d): name uses %s", defs$name[i], defs$file[i],
                                defs$first[i], paste(bad_name, collapse = ", ")))
    }
    if (length(bad_formals)) {
      items <- c(items, sprintf("%s (%s:%d): formal %s", defs$name[i], defs$file[i],
                                defs$first[i], paste(bad_formals, collapse = ", ")))
    }
  }
  gate(
    "uk_spelling", !length(items),
    sprintf("%d US spellings in %d function names and %d formals", length(items),
            nrow(defs), n_formals),
    items
  )
}

gate_app_size <- function(defs) {
  app <- defs[defs$name %in% names(max_app_lines), ]
  if (!nrow(app)) {
    return(gate("app_size", TRUE, "no app_ui() or app_server() in R/; nothing to check"))
  }
  limit <- unname(max_app_lines[app$name])
  over <- app$lines > limit
  gate(
    "app_size", !any(over),
    sprintf("%d of %d app functions over their limit: %s", sum(over), nrow(app),
            paste(sprintf("%s %d lines (limit %d)", app$name, app$lines, limit), collapse = ", ")),
    sprintf("%s, %s:%d-%d, %d lines; limit %d", app$name[over], app$file[over],
            app$first[over], app$last[over], app$lines[over], limit[over])
  )
}

# Calls of base-graphics functions in R code, read from the parser's tokens,
# so a string, a comment or a ggplot2 function such as geom_point() never
# matches. A bare call to a name the same code assigns or takes as a formal,
# or to a function in `own`, calls that object and not graphics. `offset` maps the lines of a code
# chunk back to its file.
base_graphics_hits <- function(code, where, offset = 0L, own = character()) {
  pd <- tryCatch(utils::getParseData(parse(text = code, keep.source = TRUE)),
                 error = function(e) NULL)
  if (is.null(pd)) return(data.frame(where = where, line = NA_integer_, call = "(could not parse)"))
  calls <- pd[pd$token == "SYMBOL_FUNCTION_CALL", ]
  pkgs <- pd[pd$token == "SYMBOL_PACKAGE", c("parent", "text")]
  ns <- pkgs$text[match(calls$parent, pkgs$parent)]
  term <- pd[pd$terminal, ]
  term <- term[order(term$line1, term$col1), ]
  nxt <- c(term$token[-1], "")
  assigned <- c(term$text[term$token == "SYMBOL" & nxt %in% c("LEFT_ASSIGN", "EQ_ASSIGN")],
                term$text[term$token == "SYMBOL_FORMALS"])
  local <- is.na(ns) & calls$text %in% c(assigned, own)
  drawing <- setdiff(union(getNamespaceExports("graphics"), device_calls), "plot")
  hit <- (calls$text %in% drawing & (is.na(ns) | ns != "ggplot2") & !local) |
    (!is.na(ns) & ns == "graphics") |
    (calls$text == "plot" & !is.na(ns) & ns %in% base_plot_namespaces)
  if (!any(hit)) return(data.frame(where = character(), line = integer(), call = character()))
  data.frame(where = where, line = calls$line1[hit] + offset,
             call = paste0(ifelse(is.na(ns[hit]), "", paste0(ns[hit], "::")), calls$text[hit], "()"))
}

# R code chunks of an R Markdown or Quarto file, each with the line it starts on.
rmd_chunks <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  open <- grep("^\\s*```+\\s*\\{[rR]([ ,}]|$)", lines)
  fence <- grep("^\\s*```+\\s*$", lines)
  lapply(open, function(o) {
    end <- fence[fence > o][1]
    if (is.na(end)) end <- length(lines) + 1L
    list(code = lines[seq_len(end - o - 1L) + o], offset = o)
  })
}

gate_base_graphics <- function(pkg_dir) {
  r_files <- list.files(file.path(pkg_dir, "R"), pattern = "\\.[Rr]$", full.names = TRUE)
  rd_files <- list.files(file.path(pkg_dir, "man"), pattern = "\\.Rd$", full.names = TRUE)
  vig_files <- list.files(file.path(pkg_dir, "vignettes"), pattern = "\\.(Rmd|qmd|R)$",
                          full.names = TRUE)
  rel <- function(f) sub(paste0("^", pkg_dir, "/"), "", f)
  own <- function_defs(pkg_dir)$defs$name
  hits <- list()
  for (f in r_files) {
    hits[[f]] <- base_graphics_hits(readLines(f, warn = FALSE, encoding = "UTF-8"), rel(f),
                                    own = own)
  }
  for (f in rd_files) {
    ex <- tempfile(fileext = ".R")
    tools::Rd2ex(f, ex, commentDontrun = FALSE, commentDonttest = FALSE)
    if (file.exists(ex)) hits[[f]] <- base_graphics_hits(readLines(ex, warn = FALSE), paste(rel(f), "example"))
  }
  for (f in vig_files) {
    if (grepl("\\.R$", f)) {
      hits[[f]] <- base_graphics_hits(readLines(f, warn = FALSE, encoding = "UTF-8"), rel(f))
      next
    }
    for (ch in rmd_chunks(f)) {
      hits[[length(hits) + 1]] <- base_graphics_hits(ch$code, rel(f), ch$offset)
    }
  }
  h <- do.call(rbind, hits)
  enforce <- identical(base_graphics_mode, "fail")
  detail <- sprintf(
    "%d base-graphics calls in %d of %d files (R/ %d, man/ %d, vignettes/ %d)%s",
    nrow(h), length(unique(h$where)), length(r_files) + length(rd_files) + length(vig_files),
    length(r_files), length(rd_files), length(vig_files),
    if (enforce) "" else "; reported only (POSITIONR_BASE_GRAPHICS=report)"
  )
  items <- sprintf("%s:%s %s", h$where, ifelse(is.na(h$line), "?", h$line), h$call)
  if (enforce || !nrow(h)) return(gate("base_graphics", !nrow(h), detail, items))
  say("%s %-15s %s", "WARN", "base_graphics", detail)
  for (it in items) say("     %s", it)
  structure(TRUE, reported = nrow(h))
}

# Every exported function appears on at least one page of the suite's sites,
# in a chunk that runs. The page checker runs from a checkout of the whole
# suite, where every package's site is at hand; it fails when an export is
# on no page and lists the functions on one page only.
gate_pages_cover <- function(pkg_dir) {
  script <- file.path(dirname(normalizePath(pkg_dir)), "tools", "check_pages.py")
  if (!file.exists(script)) {
    if (identical(Sys.getenv("POSITIONR_LONE_CHECKOUT"), "1")) {
      return(gate("pages_cover", NA, "not run: a checkout of one package has no suite page checker"))
    }
    return(gate("pages_cover", FALSE, "the suite page checker is not beside this package; gate did not run (POSITIONR_LONE_CHECKOUT=1 reports it as not run)"))
  }
  if (!nzchar(Sys.which("python3"))) {
    return(gate("pages_cover", FALSE, "python3 not found; gate did not run"))
  }
  pkg <- basename(normalizePath(pkg_dir))
  out <- suppressWarnings(system2("python3", c("-I", shQuote(script), pkg),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  summary <- grep(paste0("^", pkg, ":"), out, value = TRUE)
  gate("pages_cover", is.null(status) || status == 0L,
       if (length(summary)) summary[1] else "no summary line from the page checker",
       items = grep("FAIL|Error|Traceback", out, value = TRUE))
}

# ---- Gates on the built package ----
build_tarball <- function(pkg_dir, out_dir) {
  # R CMD build writes into the working directory, so it runs from out_dir.
  old <- setwd(out_dir)
  on.exit(setwd(old))
  log <- system2(r_exe("R"), c("CMD", "build", "--no-manual", shQuote(pkg_dir)),
                 stdout = TRUE, stderr = TRUE)
  tar <- list.files(out_dir, pattern = "\\.tar\\.gz$", full.names = TRUE)
  if (!length(tar)) {
    cat(tail(log, 20), sep = "\n")
    return(NA_character_)
  }
  tar[1]
}

tarball_limit <- function(pkg_dir) {
  field <- "Config/positionR/tarball_limit"
  value <- read.dcf(file.path(pkg_dir, "DESCRIPTION"), fields = field)[1, 1]
  if (is.na(value)) return(list(bytes = max_tarball_bytes, from = "the default"))
  bytes <- suppressWarnings(as.numeric(trimws(value)))
  if (!is.finite(bytes) || bytes <= 0) {
    return(list(bytes = NA_real_, from = sprintf("%s in DESCRIPTION reads '%s' and gives no byte count", field, value)))
  }
  list(bytes = bytes, from = sprintf("%s in DESCRIPTION", field))
}

gate_tarball_size <- function(tarball, pkg_dir) {
  if (is.na(tarball) || !file.exists(tarball)) {
    return(gate("tarball_size", FALSE, "no tarball: R CMD build failed; gate did not run"))
  }
  limit <- tarball_limit(pkg_dir)
  if (is.na(limit$bytes)) {
    return(gate("tarball_size", FALSE, sprintf("no limit: %s; gate did not run", limit$from)))
  }
  size <- file.size(tarball)
  gate(
    "tarball_size", size <= limit$bytes,
    sprintf("%s is %s bytes; limit %s, from %s", basename(tarball), fmt_int(size),
            fmt_int(limit$bytes), limit$from)
  )
}

time_one_example <- function(pkg, lib, ex_file, out_dir) {
  runner <- file.path(out_dir, paste0(basename(ex_file), ".run.R"))
  writeLines(c(
    sprintf("suppressPackageStartupMessages(library(%s, lib.loc = %s))", pkg, deparse(lib)),
    # A plotting example would otherwise write Rplots.pdf into the folder the
    # gates were started from, usually the package root.
    "grDevices::pdf(NULL)",
    "t0 <- proc.time()",
    sprintf("source(%s, echo = FALSE)", deparse(ex_file)),
    "t1 <- proc.time() - t0",
    "cpu <- sum(t1[c('user.self', 'sys.self', 'user.child', 'sys.child')], na.rm = TRUE)",
    "cat(sprintf('TIMING %.3f %.3f\\n', t1[['elapsed']], cpu))"
  ), runner)
  out <- suppressWarnings(system2(
    r_exe("Rscript"), c("--no-save", "--no-restore", "--no-init-file", shQuote(runner)),
    stdout = TRUE, stderr = TRUE
  ))
  tm <- grep("^TIMING ", out, value = TRUE)
  if (!length(tm)) return(c(elapsed = NA_real_, cpu = NA_real_))
  v <- as.numeric(strsplit(sub("^TIMING ", "", tm[1]), " ")[[1]])
  c(elapsed = v[1], cpu = v[2])
}

gate_example_time <- function(pkg, pkg_dir, tarball, out_dir) {
  rd_files <- list.files(file.path(pkg_dir, "man"), pattern = "\\.Rd$", full.names = TRUE)
  ex_files <- character()
  for (rd in rd_files) {
    ex <- file.path(out_dir, sub("\\.Rd$", ".R", basename(rd)))
    tools::Rd2ex(rd, ex, commentDontrun = TRUE, commentDonttest = FALSE)
    if (file.exists(ex) && file.size(ex) > 0) ex_files <- c(ex_files, ex)
  }
  denom <- sprintf("%d of %d Rd files carry examples", length(ex_files), length(rd_files))
  if (!length(ex_files)) {
    return(gate("example_time", TRUE, sprintf("0 examples to time; %s", denom)))
  }
  if (is.na(tarball)) return(gate("example_time", FALSE, "no tarball to install; gate did not run"))
  lib <- file.path(out_dir, "lib")
  dir.create(lib, showWarnings = FALSE)
  inst <- system2(r_exe("R"), c("CMD", "INSTALL", paste0("--library=", shQuote(lib)), shQuote(tarball)),
                  stdout = TRUE, stderr = TRUE)
  if (!dir.exists(file.path(lib, pkg))) {
    cat(tail(inst, 20), sep = "\n")
    return(gate("example_time", FALSE, "package did not install from the tarball; gate did not run"))
  }
  times <- t(vapply(ex_files, time_one_example, c(elapsed = 0, cpu = 0),
                    pkg = pkg, lib = lib, out_dir = out_dir))
  name <- sub("\\.R$", ".Rd", basename(ex_files))
  failed <- is.na(times[, "elapsed"])
  slow <- !failed & (times[, "elapsed"] > max_example_seconds | times[, "cpu"] > max_example_seconds)
  gate(
    "example_time", !any(slow) && !any(failed),
    sprintf("%d slow and %d failing of %d timed (%s); limit %g s; slowest %.2f s",
            sum(slow), sum(failed), length(ex_files), denom, max_example_seconds,
            max(c(0, times[!failed, ]))),
    c(sprintf("%s: %.2f s elapsed, %.2f s CPU", name[slow], times[slow, "elapsed"], times[slow, "cpu"]),
      sprintf("%s: example failed to run", name[failed]))
  )
}

# ---- Run ----
main <- function(args) {
  dirs <- args[!grepl("^--", args)]
  pkg_dir <- normalizePath(if (length(dirs)) dirs[1] else ".", mustWork = TRUE)
  tar_arg <- sub("^--tarball=", "", grep("^--tarball=", args, value = TRUE))
  desc <- read.dcf(file.path(pkg_dir, "DESCRIPTION"), fields = c("Package", "Version"))
  pkg <- desc[1, "Package"]
  say("gates: %s %s in %s", pkg, desc[1, "Version"], pkg_dir)
  out_dir <- tempfile("gates_")
  dir.create(out_dir)

  fd <- function_defs(pkg_dir)
  ns <- parseNamespaceFile(basename(pkg_dir), dirname(pkg_dir))
  exports <- ns$exports
  tarball <- if (length(tar_arg)) normalizePath(tar_arg[1]) else build_tarball(pkg_dir, out_dir)

  ok <- c(
    function_length = gate_function_length(fd$defs, fd$n_files),
    export_formals = gate_export_formals(fd$defs, exports, ns$exportPatterns),
    export_masking = gate_export_masking(exports),
    uk_spelling = gate_uk_spelling(fd$defs),
    app_size = gate_app_size(fd$defs),
    base_graphics = (bg <- gate_base_graphics(pkg_dir)),
    tarball_size = gate_tarball_size(tarball, pkg_dir),
    example_time = gate_example_time(pkg, pkg_dir, tarball, out_dir),
    pages_cover = gate_pages_cover(pkg_dir)
  )
  skipped <- names(ok)[is.na(ok)]
  ok <- ok[!is.na(ok)]
  not_ran <- intersect(names(ok), ls(did_not_run))
  ran <- ok[setdiff(names(ok), not_ran)]
  say("gates: %d run, %d passed, %d failed%s%s%s%s", length(ran), sum(ran), sum(!ran),
      if (any(!ran)) paste0(": ", paste(names(ran)[!ran], collapse = ", ")) else "",
      if (length(not_ran)) sprintf("; %d blocked, could not run: %s", length(not_ran),
                                   paste(not_ran, collapse = ", ")) else "",
      if (length(skipped)) sprintf("; %d not run: %s", length(skipped), paste(skipped, collapse = ", ")) else "",
      if (is.null(attr(bg, "reported"))) "" else
        sprintf("; base_graphics reported %d calls without failing", attr(bg, "reported")))
  invisible(all(ok))
}

if (sys.nframe() == 0L) {
  quit(status = if (main(commandArgs(trailingOnly = TRUE))) 0L else 1L)
}
