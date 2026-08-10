# One-command replication check for the carbon-labeling meta-analysis.
# Run from this folder with: Rscript --vanilla run_all.R

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
  setwd(dirname(script_path))
}

source("code/00_setup.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(weightr)
  library(zcurve)
})

check_close <- function(label, value, expected, tolerance) {
  ok <- isTRUE(abs(as.numeric(value) - expected) <= tolerance)
  data.frame(
    item = label,
    value = as.numeric(value),
    expected = expected,
    tolerance = tolerance,
    pass = ok
  )
}

d <- metadata_change[!is.na(metadata_change$SampleSize), ]

re <- rma(cohens_d, v, data = d, method = "REML")
re_pred <- predict(re)
mv <- rma.mv(cohens_d, v, random = list(~1 | studyid, ~1 | Key), data = d, method = "REML")
mv_pred <- predict(mv)

pet <- rma.mv(
  cohens_d, v,
  mods = ~ sei,
  random = list(~1 | studyid, ~1 | Key),
  data = d,
  method = "REML"
)

selection_cutpoints <- list(c(0.025), c(0.05), c(0.025, 0.975), c(0.05, 0.5, 0.95))
selection <- do.call(rbind, lapply(selection_cutpoints, function(st) {
  w <- weightfunct(d$cohens_d, d$v, steps = st)
  lrt <- 2 * (w[[1]]$value - w[[2]]$value)
  data.frame(
    cutpoints = paste(st, collapse = ", "),
    LRT = lrt,
    df = length(st),
    p = pchisq(lrt, length(st), lower.tail = FALSE)
  )
}))

set.seed(1)
zc <- zcurve(abs(d$z))
zcs <- summary(zc)
edr <- zcs$coefficients["EDR", "Estimate"]
soric_fdr <- 100 * (1 / edr - 1) * (.05 / .95)

pet_all <- rma.mv(
  cohens_d, v,
  mods = ~ sei + research_setting + scenario + TLL + product_category,
  random = list(~1 | studyid, ~1 | Key),
  data = d,
  method = "REML"
)
sqrt_n <- rma.mv(
  cohens_d, v,
  mods = ~ sqrtninv,
  random = list(~1 | studyid, ~1 | Key),
  data = d,
  method = "REML"
)

checks <- bind_rows(
  check_close("k effect sizes", nrow(d), 52, 0),
  check_close("articles", length(unique(d$Key)), 22, 0),
  check_close("multilevel d", mv$b[1], 0.133, 0.001),
  check_close("multilevel CI lower", mv$ci.lb, 0.029, 0.001),
  check_close("multilevel CI upper", mv$ci.ub, 0.236, 0.001),
  check_close("Q statistic", re$QE, 260.77, 0.02),
  check_close("PET intercept/lower-bound g", pet$b[1], 0.000, 0.001),
  check_close("PET CI lower", pet$ci.lb[1], -0.142, 0.001),
  check_close("PET CI upper", pet$ci.ub[1], 0.142, 0.001),
  check_close("RE prediction interval lower", re_pred$pi.lb, -0.347, 0.001),
  check_close("RE prediction interval upper", re_pred$pi.ub, 0.631, 0.001),
  check_close("selection-model smallest p", min(selection$p), 0.196, 0.001),
  check_close("z-curve ERR", zcs$coefficients["ERR", "Estimate"], 0.77, 0.005),
  check_close("z-curve EDR", edr, 0.73, 0.005),
  check_close("Soric implied FDR percent", soric_fdr, 1.9, 0.1),
  check_close("PET small-study slope", pet$b[2], 1.026, 0.001),
  check_close("PET all-moderator slope", pet_all$b[2], 1.058, 0.001),
  check_close("sqrt(1/n) slope", sqrt_n$b[2], 1.950, 0.001)
)

write.csv(checks, file.path(output_dir, "reproduction_check.csv"), row.names = FALSE)

source("code/15_si_tables_S7_S10.R")

summary_lines <- c(
  "Carbon-labeling replication check",
  sprintf("Package root: %s", root_dir),
  sprintf("k = %d effect sizes; articles = %d", nrow(d), length(unique(d$Key))),
  sprintf("Primary multilevel model: d = %.3f, 95%% CI [%.3f, %.3f]", mv$b[1], mv$ci.lb, mv$ci.ub),
  sprintf("Heterogeneity: Q(51) = %.2f, p < .001; RE I2 = %.1f%%", re$QE, re$I2),
  sprintf("PET lower-bound estimate: g = %.3f, 95%% CI [%.3f, %.3f]", pet$b[1], pet$ci.lb[1], pet$ci.ub[1]),
  sprintf("Prediction interval: RE [%.3f, %.3f]; multilevel [%.3f, %.3f]", re_pred$pi.lb, re_pred$pi.ub, mv_pred$pi.lb, mv_pred$pi.ub),
  sprintf("Three-parameter selection model: all LRTs non-significant; smallest p = %.3f", min(selection$p)),
  sprintf("z-curve: ERR = %.2f; EDR = %.2f; implied false-discovery rate = %.1f%%", zcs$coefficients["ERR", "Estimate"], edr, soric_fdr),
  sprintf("Small-study slope: PET b = %.3f, p = %.4f; all-moderator PET b = %.3f, p = %.4f; sqrt(1/n) b = %.3f, p = %.4f", pet$b[2], pet$pval[2], pet_all$b[2], pet_all$pval[2], sqrt_n$b[2], sqrt_n$pval[2]),
  sprintf("Automated numeric checks passed: %s", if (all(checks$pass)) "YES" else "NO")
)
writeLines(summary_lines, file.path(output_dir, "reproduction_summary.txt"))
print(checks)
cat(paste(summary_lines, collapse = "\n"), "\n")

if (!all(checks$pass)) {
  stop("At least one reproduction check failed. See output/reproduction_check.csv.")
}
