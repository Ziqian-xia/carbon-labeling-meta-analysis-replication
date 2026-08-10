# ------------------------------------------------------------------------------
# Script: 15_si_tables_S7_S10.R   (REPRODUCIBLE; seed=1 for z-curve)
# Purpose: Produce the contents of new Supplementary Tables S7-S10 added in the
#          manuscript, on the canonical k=52 dataset:
#            S7  three-parameter selection model (heterogeneity-robust bias test)
#            S8  z-curve analysis
#            S9  prediction intervals (overall + subgroups)
#            S10 precision-effect (small-study) tests, incl. SMD-artifact-robust
# Inputs:  data/metadata_change.xlsx
# Outputs: output/SI_S7.csv ... SI_S10.csv  + console summary
# ------------------------------------------------------------------------------
suppressMessages({library(readxl); library(metafor); library(weightr); library(zcurve)})
source("code/00_setup.R")
d <- metadata_change
d <- d[!is.na(d$SampleSize), ]
d$sei <- sqrt(d$v); d$z <- d$cohens_d/d$sei; d$sqrtninv <- sqrt(1/d$SampleSize)
dir.create(output_dir, showWarnings = FALSE)

## ---- S7: three-parameter selection model, cut-point sensitivity --------------
## NOTE: weightfunct() stores the NEGATIVE log-likelihood in $value, so the
## likelihood-ratio statistic is 2*(null$value - adjusted$value); getting this
## sign wrong yields a spurious negative LRT and p = 1.
S7 <- do.call(rbind, lapply(
  list(c(0.025), c(0.05), c(0.025, 0.975), c(0.05, 0.5, 0.95)),
  function(st) {
    w <- weightfunct(d$cohens_d, d$v, steps = st)
    un <- w[[1]]; ad <- w[[2]]
    lrt <- 2 * (un$value - ad$value)
    data.frame(cutpoints = paste(st, collapse = ", "),
               adj_mean  = round(ad$par[2], 3),
               omega     = paste(round(ad$par[-(1:2)], 2), collapse = ", "),
               LRT       = round(lrt, 3),
               df        = length(st),
               p         = round(pchisq(lrt, length(st), lower.tail = FALSE), 3))
  }))
write.csv(S7, file.path(output_dir, "SI_S7.csv"), row.names = FALSE)

## ---- S8: z-curve -------------------------------------------------------------
set.seed(1); zc <- zcurve(abs(d$z)); s <- summary(zc)
edr <- s$coefficients["EDR","Estimate"]
S8 <- data.frame(
  metric = c("Expected replication rate (ERR)","Expected discovery rate (EDR)",
             "Implied false-discovery rate (Soric)","Significant findings (|z|>1.96)"),
  estimate = c(sprintf("%.2f", s$coefficients["ERR","Estimate"]),
               sprintf("%.2f", edr),
               sprintf("%.1f%%", 100*(1/edr-1)*(.05/.95)),
               sprintf("%d of %d", sum(abs(d$z)>1.96), nrow(d))),
  ci95 = c(sprintf("%.2f-%.2f", s$coefficients["ERR","l.CI"], s$coefficients["ERR","u.CI"]),
           sprintf("%.2f-%.2f", s$coefficients["EDR","l.CI"], s$coefficients["EDR","u.CI"]), "", ""))
write.csv(S8, file.path(output_dir, "SI_S8.csv"), row.names = FALSE)

## ---- S9: prediction intervals (overall + subgroups) --------------------------
pirow <- function(lab, sub) {
  m <- rma(cohens_d, v, data = sub, method = "REML"); p <- predict(m)
  data.frame(group = lab, k = m$k, d = round(m$b,3),
             ci = sprintf("%.3f, %.3f", m$ci.lb, m$ci.ub),
             pi = sprintf("%.3f, %.3f", p$pi.lb, p$pi.ub),
             I2 = sprintf("%.1f%%", m$I2))
}
mv <- rma.mv(cohens_d, v, random=list(~1|studyid,~1|Key), data=d, method="REML"); pmv <- predict(mv)
S9 <- rbind(
  pirow("Overall (random-effects)", d),
  data.frame(group="Overall (three-level multilevel)", k=mv$k, d=round(mv$b,3),
             ci=sprintf("%.3f, %.3f", mv$ci.lb, mv$ci.ub),
             pi=sprintf("%.3f, %.3f", pmv$pi.lb, pmv$pi.ub), I2="-"),
  pirow("Online", subset(d, scenario=="Online")),
  pirow("Offline", subset(d, scenario=="Offline")),
  pirow("Lab", subset(d, research_setting=="Lab")),
  pirow("Field", subset(d, research_setting=="Field")))
write.csv(S9, file.path(output_dir, "SI_S9.csv"), row.names = FALSE)

## ---- S10: precision-effect (small-study) tests -------------------------------
p0 <- rma.mv(cohens_d, v, mods=~sei, random=list(~1|studyid,~1|Key), data=d, method="REML")
p1 <- rma.mv(cohens_d, v, mods=~sei+research_setting, random=list(~1|studyid,~1|Key), data=d, method="REML")
p2 <- rma.mv(cohens_d, v, mods=~sei+research_setting+scenario+TLL+product_category, random=list(~1|studyid,~1|Key), data=d, method="REML")
ps <- rma.mv(cohens_d, v, mods=~sqrtninv, random=list(~1|studyid,~1|Key), data=d, method="REML")
re <- rma(cohens_d, v, data=d, method="REML")
egg_se <- regtest(re, predictor="sei"); egg_vi <- regtest(re, predictor="vi")
S10 <- data.frame(
  specification = c("PET, predictor = SE (primary)",
                    "PET + research setting (lab/field)",
                    "PET + all moderators (setting, scenario, label, product)",
                    "Funnel test, predictor = sqrt(1/n)  [SMD-artifact-robust]",
                    "Egger, predictor = SE", "Egger, predictor = variance",
                    "Precision vs research setting (mean SE Lab vs Field)"),
  statistic = c(sprintf("b=%.3f", p0$b[2]), sprintf("b=%.3f", p1$b[2]), sprintf("b=%.3f", p2$b[2]),
                sprintf("b=%.3f", ps$b[2]), sprintf("z=%.2f", egg_se$zval), sprintf("z=%.2f", egg_vi$zval),
                sprintf("%.3f vs %.3f", mean(d$sei[d$research_setting=="Lab"]), mean(d$sei[d$research_setting=="Field"]))),
  p = c(round(p0$pval[2],3), round(p1$pval[2],3), round(p2$pval[2],3), round(ps$pval[2],3),
        round(egg_se$pval,3), round(egg_vi$pval,3), round(t.test(sei~research_setting, data=d)$p.value,3)),
  note = c("small-study effect present","survives adjustment","survives adjustment",
           "survives; not a metric artifact","", "","precision uncorrelated with setting"))
write.csv(S10, file.path(output_dir, "SI_S10.csv"), row.names = FALSE)

cat("Wrote SI_S7.csv, SI_S8.csv, SI_S9.csv, SI_S10.csv to output/\n")
print(S7); cat("\n"); print(S8); cat("\n"); print(S9); cat("\n"); print(S10)
cat("\n--- sessionInfo ---\n"); print(sessionInfo())
