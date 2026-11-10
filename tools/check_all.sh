#!/usr/bin/env bash
# ---- Local check for the positionR umbrella package ----
# Usage: tools/check_all.sh   (from any folder)
# Steps: document, build, check, news, names. Every step runs even
# after a failure; the exit status is 1 when any step failed. Build and check
# output go to a temporary folder named on the last line. The check runs with
# --as-cran but with the remote CRAN checks off, so it makes no web calls for
# URLs or the CRAN database. The four suite packages must be installed.
# The names step checks how package names are written on the pages; it needs
# a checkout of the whole suite and fails outside one.

PKG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ROOT="$(dirname "$PKG_DIR")"
OUT="$(mktemp -d "${TMPDIR:-/tmp}/check_positionR.XXXXXX")"
FAILED=""

run_step() {
  echo "== $1"
  if ! "step_$1"; then FAILED="$FAILED $1"; echo "FAIL $1"; fi
}

step_document() { (cd "$PKG_DIR" && Rscript -e 'roxygen2::roxygenise()'); }
step_build() { (cd "$OUT" && R CMD build "$PKG_DIR"); }
step_check() {
  (cd "$OUT" && _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE \
    R CMD check --as-cran --no-manual positionR_*.tar.gz) | tail -20
  ! grep -q "ERROR\|WARNING" "$OUT/positionR.Rcheck/00check.log"
}
step_news() { (cd "$PKG_DIR" && Rscript tools/check_news_map.R); }
step_names() { (cd "$ROOT" && python3 -I tools/check_pkg_names.py "$(basename "$PKG_DIR")"); }

for s in document build check news names; do run_step "$s"; done
echo "output: $OUT"
if [ -n "$FAILED" ]; then echo "failed:$FAILED"; exit 1; fi
echo "all steps passed"
