# Reproduce existing TES and Figure S5 illustration without changing primary models.
# Run from package root: Rscript --vanilla code/17_illustrative_diagnostics.R
source("code/00_setup.R")
suppressPackageStartupMessages({library(metapower); library(ggplot2)})
d <- metadata_change
# Legacy TES specification: independent comparisons and a common true effect.
# studyid is unique for every row, so prior aggregation by studyid was a no-op.
tes_fit <- metafor::tes(d$cohens_d, d$v, tau2 = 0, test = "chi2")
write.csv(data.frame(observed = tes_fit$O, expected = tes_fit$E,
  common_effect = tes_fit$theta, tau2_assumed = 0, p = tes_fit$pval,
  effect_estimates = nrow(d), articles = length(unique(d$Key))),
  file.path(output_dir, "TES_common_effect.csv"), row.names = FALSE)
writeLines(c("TES under a common-effect assumption, treating comparisons as independent.",
  "This diagnostic does not account for between-effect heterogeneity or within-article dependence.",
  capture.output(print(tes_fit))), file.path(output_dir, "TES_common_effect.txt"))

# Preserve the existing numerical input construction for this illustration.
# The sample-size field mixes observation units. The resulting study_size is
# therefore only an assumed total N in a hypothetical equal-size design; it is
# not an estimate of the number of participants per included study or per arm.
k <- nrow(d)
assumed_n <- mean(d$SampleSize) / 2
re <- rma(cohens_d, v, data = d, method = "REML")
assumed_i2 <- re$tau2 / (re$tau2 + mean(d$v))
inputs <- data.frame(k = k, assumed_total_n_per_study = assumed_n,
  assumed_i2 = assumed_i2, alpha = .05)
write.csv(inputs, file.path(output_dir, "Figure_S5_assumptions.csv"), row.names = FALSE)
curve <- data.frame(d = seq(.01, 1, by = .01))
curve$power <- vapply(curve$d, function(es) {
  metapower::mpower(effect_size = es, study_size = assumed_n, k = k,
    i2 = assumed_i2, es_type = "d", p = .05)$power$random_power
}, numeric(1))
write.csv(curve, file.path(output_dir, "Figure_S5_values.csv"), row.names = FALSE)
i <- which(diff(curve$power >= .8) != 0)[1]
stopifnot(!is.na(i))
mde <- with(curve, d[i] + (.8 - power[i]) * (d[i+1] - d[i]) / (power[i+1] - power[i]))
p <- ggplot(curve, aes(d, power)) +
  geom_line(linewidth = 1.2, colour = "#1F4E79") +
  geom_hline(yintercept = .8, linetype = "dashed", colour = "darkred") +
  geom_vline(xintercept = mde, linetype = "dashed", colour = "#1F4E79") +
  annotate("point", x = mde, y = .8) +
  annotate("text", x = mde + .05, y = .75,
    label = sprintf("d = %.3f", mde), hjust = 0, size = 4) +
  scale_x_continuous(breaks = seq(0, 1, .1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Assumed effect size (Cohen's d)", y = "Statistical power",
    title = "Illustrative power under independence assumptions",
    subtitle = sprintf("k = %d; assumed total N per study = %.0f; assumed I² = %.1f%%",
      k, assumed_n, 100 * assumed_i2)) + theme_bw(base_size = 14)
ggsave(file.path(output_dir, "Figure_S5_power.png"), p,
  width = 10, height = 7.93, units = "in", dpi = 300, bg = "white")
