###Table 1：Summary of bias analyses
#1 Prepare the data
library(metafor)
library(clubSandwich)       # robust tests for mv models

# Copy the master metadata and add standard-error columns required downstream.
dat <- metadata_change
dat$sei <- sqrt(dat$v)      # for Egger
# Recalculate yi/vi explicitly to ensure metafor uses the intended columns.
dat <- escalc(measure    = "SMD",
              yi         = cohens_d,
              vi         = v,
              data       = dat,
              var.names  = c("yi", "vi"))

#2 Multilevel meta-analysis (base model)
res_mv <- rma.mv(
  yi     = cohens_d,
  V      = v,
  random = list(~1 | studyid,  # level-2
                ~1 | Key),     # level-3
  method = "REML",
  data   = dat
)

#Egger's test for publication bias
rm <- rma(yi=cohens_d, vi=v, data = metadata_change)
regtest(rm, model = "rma")

#the rank correlation test
ranktest(res_mv)

#Test of excess significance (TES)
dat_agg <- aggregate(dat, cluster = studyid, yi = "cohens_d", vi = "v",rho= 0.5)
tes(dat_agg$yi, dat_agg$vi, test = "chi2")


### BMA random
library(RoBMA)
fit_BMA_main <- RoBMA(d = metadata_change$cohens_d, v = metadata_change$v, study_ids = metadata_change$Key,
                      priors_bias= list(prior_PET(distribution = "Cauchy", parameters = list(0,1),  prior_weights = 1/2), 
                                        prior_PEESE(distribution = "Cauchy", parameters = list(0,5),  prior_weights = 1/2) ),
                      parallel = TRUE, seed = 1)
summary(fit_BMA_main)


