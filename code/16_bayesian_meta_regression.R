library(readxl)
library(RoBMA)

data_regression_change <- as.data.frame(read_excel("data/data_regression_change.xlsx"))

cols <- c("Key", "publication_year", "SampleSize", "sample_characteristics",
          "labeltype", "TLL", "product_category", "intervention_duration",
          "outcome_kind", "channel", "Country", "scenario", "direction",
          "iso3", "wtp_avg", "cohens_d", "v", "GDP_per_capita2023",
          "log_GDP_per_capita2023", "gccs_education_lev3")

dat <- na.omit(data_regression_change[, cols])
names(dat)[names(dat) == "cohens_d"] <- "d"
names(dat)[names(dat) == "v"] <- "var"
names(dat)[names(dat) == "SampleSize"] <- "n"
dat$SampleSize <- dat$n

predictors <- c("Publication Year" = "publication_year",
                "Willingness To Pay" = "wtp_avg",
                "Log GDP Per Capita (2023)" = "log_GDP_per_capita2023",
                "Sample Size" = "SampleSize",
                "Intervention Duration" = "intervention_duration")

results <- lapply(names(predictors), function(label) {
  term <- predictors[[label]]
  fit <- NoBMA.reg(as.formula(paste("~", term)), data = dat,
                   parallel = FALSE, seed = 1)
  s <- summary(fit, output_scale = "r")
  data.frame(Predictor = label,
             `P(incl)` = round(s$components_predictors[term, "post_prob"], 3),
             `BF(incl)` = round(s$components_predictors[term, "inclusion_BF"], 3),
             `β (Mean)` = round(s$estimates_predictors[term, "Mean"], 3),
             `95% Credible Interval` = sprintf("%.3f, %.3f",
                                              s$estimates_predictors[term, "0.025"],
                                              s$estimates_predictors[term, "0.975"]),
             N = nrow(dat), check.names = FALSE)
})

dir.create("output", showWarnings = FALSE)
write.csv(do.call(rbind, results),
          "output/Table_S4_bayesian_meta_regression_updated.csv", row.names = FALSE)
