# ------------------------------------------------------------------------------
# Script: 14_robma_bias_decomposition.R   (RoBMA 4.0.0 API; REPRODUCIBLE seed=1)
# Purpose: Resolve the tension between the submitted RoBMA bias BF and the null
#          3-parameter selection model. The submitted RoBMA (table.R) used
#          PET+PEESE bias priors ONLY (no selection models). Here we compare the
#          publication-bias inclusion BF under three bias ensembles:
#            A. PET+PEESE only (submitted spec)
#            B. DEFAULT ensemble (selection/weight-function models + PET/PEESE)
#            C. selection (weight-function) models only
# Inputs:  data/metadata_change.xlsx  (canonical k=52).
# Run:     PATH=/opt/homebrew/bin:$PATH DYLD_FALLBACK_LIBRARY_PATH=/opt/homebrew/lib Rscript ...
# Needs:   JAGS >= 4.3.1, rjags + RoBMA (4.0.0) built against it.
# ------------------------------------------------------------------------------
source("code/00_setup.R")
suppressWarnings(suppressMessages({library(readxl); library(RoBMA)}))
d <- metadata_change

bias_inclusion <- function(fit, label) {
  s <- summary(fit)
  cat("\n========== ", label, " ==========\n")
  print(s)                       # 'Components summary' includes Publication bias inclusion BF
  invisible(s)
}

## A. Submitted spec: PET + PEESE only ------------------------------------------
fitA <- RoBMA(yi = d$cohens_d, vi = d$v, measure = "SMD", cluster = d$Key,
              prior_bias = list(prior_PET  (distribution="Cauchy", parameters=list(0,1), prior_weights=1/2),
                                prior_PEESE(distribution="Cauchy", parameters=list(0,5), prior_weights=1/2)),
              parallel = FALSE, seed = 1)
bias_inclusion(fitA, "A. PET+PEESE-only (submitted spec)")

if (Sys.getenv("RUN_ROBMA_EXTENDED", unset = "0") != "1") {
  cat("\nSet RUN_ROBMA_EXTENDED=1 to run the default and selection-only RoBMA ensembles.\n")
  cat("\n========== sessionInfo ==========\n"); print(sessionInfo())
  quit(save = "no", status = 0)
}

## B. Default ensemble: selection models + PET/PEESE ----------------------------
fitB <- RoBMA(yi = d$cohens_d, vi = d$v, measure = "SMD", cluster = d$Key,
              parallel = FALSE, seed = 1)
bias_inclusion(fitB, "B. DEFAULT ensemble (selection + PET/PEESE)")

## C. Selection (weight-function) models only -----------------------------------
fitC <- RoBMA(yi = d$cohens_d, vi = d$v, measure = "SMD", cluster = d$Key,
              prior_bias = list(
                prior_weightfunction(distribution="one.sided", parameters=list(steps=c(0.05),      alpha=c(1,1)),   prior_weights=1/2),
                prior_weightfunction(distribution="one.sided", parameters=list(steps=c(0.05,0.10),  alpha=c(1,1,1)), prior_weights=1/2)),
              parallel = FALSE, seed = 1)
bias_inclusion(fitC, "C. selection-models-only")

cat("\n========== sessionInfo ==========\n"); print(sessionInfo())
