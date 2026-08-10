# Shared setup for the carbon-labeling replication package.

find_package_root <- function() {
  candidates <- unique(c(getwd(), dirname(getwd()), normalizePath(".", mustWork = FALSE)))
  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "data", "metadata_change.xlsx"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }
  stop("Could not find data/metadata_change.xlsx. Run scripts from the package root.")
}

root_dir <- find_package_root()
output_dir <- file.path(root_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(readxl)
  library(metafor)
})

metadata_change <- read_excel(file.path(root_dir, "data", "metadata_change.xlsx"))
data_regression_rob <- read_excel(file.path(root_dir, "data", "data_regression_rob.xlsx"))

metadata_change$sei <- sqrt(metadata_change$v)
metadata_change$z <- metadata_change$cohens_d / metadata_change$sei
if ("SampleSize" %in% names(metadata_change)) {
  metadata_change$sqrtninv <- sqrt(1 / metadata_change$SampleSize)
}
