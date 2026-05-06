#!/usr/bin/env Rscript
# ============================================================
# setup_r_packages.R
#
# Run this ONCE after creating the conda environment to ensure
# all R packages are installed. Conda-forge covers most of them,
# but this script catches any that didn't install cleanly and
# verifies everything is working.
#
# Usage (with the conda env active):
#   Rscript setup_r_packages.R
# ============================================================

cran_packages <- c(
  "optparse",       # CLI args for plot_genomic_region.R
  "ggplot2",        # plotting
  "gggenes",        # gene arrow diagrams
  "dplyr",          # data wrangling
  "shiny",          # Shiny app
  "DT",             # interactive tables
  "googlesheets4",  # Google Sheets access
  "rsconnect"       # deploy to shinyapps.io
)

cat("Checking and installing CRAN packages...\n\n")

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing: %s\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = FALSE)
  } else {
    cat(sprintf("  OK (already installed): %s\n", pkg))
  }
}

for (pkg in cran_packages) {
  install_if_missing(pkg)
}

cat("\n── Verification ─────────────────────────────────────────\n")
all_ok <- TRUE
for (pkg in cran_packages) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  cat(sprintf("  %-20s %s\n", pkg, if (ok) "✓" else "✗ FAILED"))
  if (!ok) all_ok <- FALSE
}

cat("\n")
if (all_ok) {
  cat("All packages installed successfully.\n")
} else {
  cat("Some packages failed — check errors above.\n")
  quit(status = 1)
}
