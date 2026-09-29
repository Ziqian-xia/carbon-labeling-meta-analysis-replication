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

# Link risk-of-bias ratings to the canonical effect dataset by stable identifiers.
rob_key <- paste(data_regression_rob$Key, data_regression_rob$studyid, sep = "::")
meta_key <- paste(metadata_change$Key, metadata_change$studyid, sep = "::")
stopifnot(!anyDuplicated(rob_key), !anyDuplicated(meta_key),
          length(rob_key) == length(meta_key), setequal(rob_key, meta_key))
idx <- match(rob_key, meta_key)
data_regression_rob$cohens_d <- metadata_change$cohens_d[idx]
data_regression_rob$v <- metadata_change$v[idx]

metadata_change$sei <- sqrt(metadata_change$v)
metadata_change$z <- metadata_change$cohens_d / metadata_change$sei
if ("SampleSize" %in% names(metadata_change)) {
  metadata_change$sqrtninv <- sqrt(1 / metadata_change$SampleSize)
}

# The covariate workbook is a separate input for the Bayesian regressions.
# Fail early if its effect values drift from the canonical dataset.
regression_data <- read_excel(file.path(root_dir, "data", "data_regression_change.xlsx"))
reg_key <- paste(regression_data$Key, regression_data$studyid, sep = "::")
stopifnot(!anyDuplicated(reg_key), length(reg_key) == length(meta_key),
          setequal(reg_key, meta_key))
reg_idx <- match(reg_key, meta_key)
stopifnot(isTRUE(all.equal(regression_data$cohens_d,
                         metadata_change$cohens_d[reg_idx], tolerance = 1e-10)),
          isTRUE(all.equal(regression_data$v,
                         metadata_change$v[reg_idx], tolerance = 1e-10)))
