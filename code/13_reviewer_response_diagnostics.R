# ------------------------------------------------------------------------------
# Script: 13_reviewer_response_diagnostics.R
# Purpose: Reviewer-response diagnostics on the CANONICAL reviewed dataset
#          (PNAS Nexus "mar5" version, carbon-labeling, k=52, d=0.13, g=0.000).
# Inputs:  data/metadata_change.xlsx  (the OSF k=52 file — reproduces
#          d=0.133 [0.03,0.24], Q=260.77, PET g=0.000 [-0.142,0.142], slope 1.03).
#          NB: this is NOT the same as plot_v1108/data/metadata_change.xlsx (k=55).
# Addresses R1: prediction interval; small-study slope vs lab/field & online/
#          offline confounds; heterogeneity-robust bias tests (3PSM, caliper);
#          z-curve evidence that significant results are largely genuine.
# ------------------------------------------------------------------------------
suppressMessages({library(readxl); library(metafor); library(weightr); library(zcurve)})
source("code/00_setup.R")
d <- metadata_change
d$sei <- sqrt(d$v); d$z <- d$cohens_d/d$sei

## (1) Prediction interval -------------------------------------------------------
re <- rma(cohens_d, v, data=d, method="REML"); pr <- predict(re)
mv <- rma.mv(cohens_d, v, random=list(~1|studyid,~1|Key), data=d, method="REML"); prm <- predict(mv)
cat(sprintf("(1) PI  RE d=%.3f CI[%.3f,%.3f] PI[%.3f,%.3f] I2=%.0f%%\n",
            re$b,re$ci.lb,re$ci.ub,pr$pi.lb,pr$pi.ub,re$I2))
cat(sprintf("        MV d=%.3f CI[%.3f,%.3f] PI[%.3f,%.3f]\n\n", mv$b,mv$ci.lb,mv$ci.ub,prm$pi.lb,prm$pi.ub))

## (2) Is the small-study slope a measured-confound artifact? --------------------
cat(sprintf("(2) precision vs setting: sei Lab=%.3f Field=%.3f p=%.3f | Online=%.3f Offline=%.3f p=%.3f\n",
  mean(d$sei[d$research_setting=="Lab"]), mean(d$sei[d$research_setting=="Field"]), t.test(sei~research_setting,data=d)$p.value,
  mean(d$sei[d$scenario=="Online"]), mean(d$sei[d$scenario=="Offline"]), t.test(sei~scenario,data=d)$p.value))
p0 <- rma.mv(cohens_d,v,mods=~sei, random=list(~1|studyid,~1|Key), data=d, method="REML")
p1 <- rma.mv(cohens_d,v,mods=~sei+research_setting, random=list(~1|studyid,~1|Key), data=d, method="REML")
p2 <- rma.mv(cohens_d,v,mods=~sei+research_setting+scenario+TLL+product_category, random=list(~1|studyid,~1|Key), data=d, method="REML")
cat(sprintf("    PET sei-only : slope=%.3f p=%.4f  g=%.3f [%.3f,%.3f]\n", p0$b[2],p0$pval[2],p0$b[1],p0$ci.lb[1],p0$ci.ub[1]))
cat(sprintf("    + lab/field  : slope=%.3f p=%.4f\n", p1$b[2],p1$pval[2]))
cat(sprintf("    + all mods   : slope=%.3f p=%.4f\n\n", p2$b[2],p2$pval[2]))

## (3) Heterogeneity-robust bias tests ------------------------------------------
w <- weightfunct(d$cohens_d, d$v, steps=c(0.025,0.975)); lrt <- 2*(w[[1]]$value-w[[2]]$value)
jb <- sum(abs(d$z)>=1.76 & abs(d$z)<1.96); ja <- sum(abs(d$z)>=1.96 & abs(d$z)<2.16)
cat(sprintf("(3) 3PSM LRT X2(2)=%.2f p=%.3f | caliper below=%d above=%d p=%.3f\n\n",
            lrt, pchisq(lrt,2,lower.tail=FALSE), jb, ja, binom.test(ja, ja+jb, .5)$p.value))

## (4) z-curve -------------------------------------------------------------------
zc <- zcurve(abs(d$z)); s <- summary(zc)
edr <- s$coefficients["EDR","Estimate"]
cat(sprintf("(4) z-curve sig=%d/%d ERR=%.2f EDR=%.2f Soric-FDR=%.1f%%\n",
            sum(abs(d$z)>1.96), nrow(d), s$coefficients["ERR","Estimate"], edr, 100*(1/edr-1)*(.05/.95)))
