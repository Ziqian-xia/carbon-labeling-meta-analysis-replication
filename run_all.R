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

# Do not silently drop rows from a reproduction of the locked manuscript dataset.
stopifnot(nrow(metadata_change) == 52L,
          all(is.finite(metadata_change$cohens_d)),
          all(is.finite(metadata_change$v) & metadata_change$v > 0),
          all(is.finite(metadata_change$SampleSize) & metadata_change$SampleSize > 0))
d <- metadata_change

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
  check_close("multilevel d", mv$b[1], 0.187, 0.001),
  check_close("multilevel CI lower", mv$ci.lb, 0.070, 0.001),
  check_close("multilevel CI upper", mv$ci.ub, 0.305, 0.001),
  check_close("Q statistic", re$QE, 322.99, 0.02),
  check_close("PET-adjusted d", pet$b[1], 0.074, 0.001),
  check_close("PET CI lower", pet$ci.lb[1], -0.107, 0.001),
  check_close("PET CI upper", pet$ci.ub[1], 0.254, 0.001),
  check_close("RE prediction interval lower", re_pred$pi.lb, -0.345, 0.001),
  check_close("RE prediction interval upper", re_pred$pi.ub, 0.731, 0.001),
  check_close("selection-model smallest p", min(selection$p), 0.182, 0.001),
  check_close("z-curve ERR", zcs$coefficients["ERR", "Estimate"], 0.81, 0.005),
  check_close("z-curve EDR", edr, 0.78, 0.005),
  check_close("Soric implied FDR percent", soric_fdr, 1.5, 0.1),
  check_close("PET small-study slope", pet$b[2], 0.847, 0.001),
  check_close("PET all-moderator slope", pet_all$b[2], 0.938, 0.001),
  check_close("sqrt(1/n) slope", sqrt_n$b[2], 1.591, 0.001)
)

# Current manuscript subgroup tests and Table S5. These reproduce existing
# specifications; they do not validate extraction decisions or model assumptions.
moderators <- c("scenario", "research_setting", "product_category",
                "sample_characteristics", "TLL")
moderator_fits <- lapply(moderators, function(term) {
  rma.mv(cohens_d, v, mods = reformulate(term),
         random = list(~1 | studyid, ~1 | Key), data = d, method = "REML")
})
moderator_table <- data.frame(moderator = moderators,
  QM = vapply(moderator_fits, function(m) as.numeric(m$QM), numeric(1)),
  p = vapply(moderator_fits, function(m) as.numeric(m$QMp), numeric(1)))
moderator_table$q_BH <- p.adjust(moderator_table$p, "BH")
write.csv(moderator_table, file.path(output_dir, "moderator_tests.csv"), row.names = FALSE)
setting <- moderator_fits[[2]]
setting_pred <- predict(setting, newmods = matrix(c(0, 1), ncol = 1))
stopifnot(identical(levels(factor(d$research_setting)), c("Field", "Lab")))
checks <- bind_rows(checks,
  check_close("setting QM", setting$QM, 7.28, .01),
  check_close("setting p", setting$QMp, .007, .001),
  check_close("setting qBH", moderator_table$q_BH[2], .035, .001),
  check_close("Field subgroup predicted d", setting_pred$pred[1], .048, .001),
  check_close("Lab subgroup predicted d", setting_pred$pred[2], .323, .001))
rob_fits <- lapply(paste0("D", 1:5), function(term) {
  rma(cohens_d, v, mods = reformulate(term), data = data_regression_rob, method = "REML")
})
rob_table <- bind_rows(lapply(seq_along(rob_fits), function(i) {
  m <- rob_fits[[i]]
  data.frame(Model = paste0("D", i), Moderator = rownames(m$b),
    Estimate = as.numeric(m$b), SE = m$se, z = m$zval, p = m$pval,
    CI.lb = m$ci.lb, CI.ub = m$ci.ub, QM = as.numeric(m$QM), QM.p = as.numeric(m$QMp))
}))
write.csv(rob_table, file.path(output_dir, "Table_S5.csv"), row.names = FALSE)
rob_target <- c(.010, .067, .012, .074, .705)
for (i in seq_along(rob_fits)) {
  checks <- bind_rows(checks, check_close(paste0("S5 D", i, " omnibus p"),
    rob_fits[[i]]$QMp, rob_target[i], .001))
}
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
write.csv(data.frame(file = names(tools::md5sum(list.files("data", full.names = TRUE))),
  md5 = unname(tools::md5sum(list.files("data", full.names = TRUE)))),
  file.path(output_dir, "input_checksums.csv"), row.names = FALSE)

write.csv(checks, file.path(output_dir, "reproduction_check.csv"), row.names = FALSE)

source("code/15_si_tables_S7_S10.R")

summary_lines <- c(
  "Carbon-labeling replication check — manuscript numerical targets 2026-09-29",
  "Scope: primary/PET models, subgroup tests, S5 and S7-S10. Not a full manuscript validation.",
  "Not checked: Bayesian fits, extraction eligibility, TES assumptions, Figure S5 power inputs.",
  sprintf("Package root: %s", root_dir),
  sprintf("k = %d effect sizes; articles = %d", nrow(d), length(unique(d$Key))),
  sprintf("Primary multilevel model: d = %.3f, 95%% CI [%.3f, %.3f]", mv$b[1], mv$ci.lb, mv$ci.ub),
  sprintf("Heterogeneity: Q(51) = %.2f, p < .001; RE I2 = %.1f%%", re$QE, re$I2),
  sprintf("PET-adjusted estimate: d = %.3f, 95%% CI [%.3f, %.3f]", pet$b[1], pet$ci.lb[1], pet$ci.ub[1]),
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
