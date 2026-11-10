#!/usr/bin/env Rscript
# ---- Check the 1.x to 2.0 table in NEWS.md ----
# Usage: Rscript tools/check_news_map.R   (from the package folder)
# Reads the table under "Where the 1.x functions went" and tests two things
# against the NAMESPACE files in the suite folder: every export and dataset
# of positionR 1.7.0 is in the second column exactly once, and every function
# named in the third column is exported by the package it is listed under.
# Exit status 1 on any miss.

# ---- Settings ----
suite_root <- normalizePath("..")
v1_dir <- file.path(suite_root, "positionR")
v2_pkgs <- c("detqc", "detfield", "arraydesign", "arraydesignR")

# ---- Read the table ----
news <- readLines("NEWS.md", warn = FALSE)
rows <- grep("^\\| [A-Z]", news, value = TRUE)
rows <- rows[!grepl("^\\| Job in 1.x", rows)]
cells <- lapply(strsplit(rows, " \\| "), function(r) gsub("^\\| | \\|$", "", r))
ticked <- function(txt) {
  hits <- regmatches(txt, gregexpr("`[A-Za-z_.0-9]+(\\(\\))?`", txt))[[1]]
  sub("\\(\\)$", "", gsub("`", "", hits))
}
cat("table rows:", length(cells), "\n")

# ---- Column 2 against the 1.x exports and datasets ----
v1_ns <- readLines(file.path(v1_dir, "NAMESPACE"))
v1_exports <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", v1_ns, value = TRUE))
v1_data <- sub("^\"(.*)\"$", "\\1",
               grep("^\"", readLines(file.path(v1_dir, "R", "data.R"), warn = FALSE), value = TRUE))
v1_all <- c(v1_exports, v1_data)
in_table <- unlist(lapply(cells, function(r) ticked(r[2])))
cat("1.x exports:", length(v1_exports), "; 1.x datasets:", length(v1_data),
    "; names in column 2:", length(in_table), "\n")
missing_v1 <- setdiff(v1_all, in_table)
twice <- in_table[duplicated(in_table)]
unknown_v1 <- setdiff(in_table, v1_all)
cat("1.x names not in the table:", length(missing_v1), paste(missing_v1, collapse = " "), "\n")
cat("1.x names in the table twice:", length(twice), paste(twice, collapse = " "), "\n")
cat("column 2 names 1.x does not export:", length(unknown_v1),
    paste(unknown_v1, collapse = " "), "\n")

# ---- Column 3 against the 2.0 exports ----
v2_exports <- lapply(setNames(v2_pkgs, v2_pkgs), function(p) {
  ns <- readLines(file.path(suite_root, p, "NAMESPACE"))
  c(sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", ns, value = TRUE)),
    unique(sub("^S3method\\(([^,]+),.*$", "\\1", grep("^S3method\\(", ns, value = TRUE))))
})
n_checked <- 0
misses <- character()
for (r in cells) {
  # Each "(**pkg**)" closes the functions written before it
  parts <- strsplit(r[3], "\\(\\*\\*")[[1]]
  for (k in seq_len(length(parts) - 1)) {
    pkg <- sub("\\*\\*\\).*$", "", parts[k + 1])
    fns <- ticked(sub("^.*\\*\\*\\)", "", parts[k]))
    for (f in fns) {
      n_checked <- n_checked + 1
      if (!f %in% v2_exports[[pkg]]) misses <- c(misses, paste0(pkg, "::", f))
    }
  }
}
no_home <- sum(grepl("^No 2.0|^Partly", vapply(cells, `[`, "", 3)))
cat("2.0 names checked in column 3:", n_checked, "; not exported where listed:",
    length(misses), paste(misses, collapse = " "), "\n")
cat("rows with no 2.0 home or a partial one:", no_home, "of", length(cells), "\n")
quit(status = if (length(c(missing_v1, twice, unknown_v1, misses))) 1L else 0L)
