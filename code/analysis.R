### benjamini_fdr_adjustments
library(metafor)
library(dplyr)

# ---- Safety checks ----
required_cols <- c("cohens_d", "v", "studyid", "Key",
                   "scenario",      # online vs offline (please confirm your column name)
                   "research_setting",              # lab/field
                   "product_category",     # ready-to-eat, packaged/grocery, non-food
                   "sample_characteristics",# student vs general population
                   "TLL")                  # traffic-light vs other

missing_cols <- setdiff(required_cols, names(metadata_change))
if (length(missing_cols) > 0) {
  stop(paste0(
    "metadata_change is missing these required columns:\n  - ",
    paste(missing_cols, collapse = "\n  - "),
    "\n\nFix: rename your variables or edit required_cols to match your dataset."
  ))
}

# ---- Helper: run model + extract omnibus QM p-value ----
omnibus_p <- function(dat, moderator) {
  fml <- as.formula(paste0("~ ", moderator))
  m <- rma.mv(
    yi = cohens_d, V = v,
    mods = fml,
    random = list(~ 1 | studyid, ~ 1 | Key),
    data = dat,
    method = "REML"
  )
  
  # If moderator has only 1 level (or was dropped), no omnibus test is meaningful
  if (length(coef(m)) <= 1) return(NA_real_)
  
  p <- anova(m, btt = 2:length(coef(m)))$QMp
  as.numeric(p)
}

# ---- Run the 5 moderator analyses (each in its own model) ----
raw_p <- c(
  research_design = omnibus_p(metadata_change, "scenario"),
  setting         = omnibus_p(metadata_change, "research_setting"),
  product_type    = omnibus_p(metadata_change, "product_category"),
  sample          = omnibus_p(metadata_change, "sample_characteristics"),
  label_design    = omnibus_p(metadata_change, "TLL")
)

# Drop NA tests (in case any moderator had too few levels after filtering)
raw_p <- raw_p[!is.na(raw_p)]

# ---- BH adjustment + before/after summaries ----
q_BH <- p.adjust(raw_p, method = "BH")

alpha_raw <- 0.05
alpha_fdr <- 0.05

sig_raw <- raw_p < alpha_raw
sig_bh  <- q_BH  < alpha_fdr

summary_counts <- tibble(
  n_tests = length(raw_p),
  n_sig_raw = sum(sig_raw),
  n_sig_bh  = sum(sig_bh)
)

# ---- Output table (before + after) ----
out <- tibble(
  test = names(raw_p),
  p_raw = raw_p,
  q_BH  = q_BH,
  sig_raw_0.05 = sig_raw,
  sig_BH_0.05  = sig_bh
) %>%
  arrange(p_raw)

print(summary_counts)


knitr::kable(
  out %>%
    mutate(
      p_raw = signif(p_raw, 3),
      q_BH  = signif(q_BH,  3)
    ),
  caption = "Omnibus moderator tests(raw p-values and BH-adjusted q-values)."
)



### consider publication bias
rm_filter <- rma.mv(
  yi     = cohens_d,
  V      = v,
  random = list(~1 | studyid,  # level-2
                ~1 | Key),     # level-3
  method = "REML",
  data   = metadata_change)

# 1.Add the study-level standard error to data frame
metadata_change <- metadata_change |>
  dplyr::mutate(sei = sqrt(v))      # PET uses the SE; PEESE uses SE²

# 2.Fit the three-level PET model
pet_filter <- rma.mv(
  yi     = cohens_d,
  V      = v,
  mods   = ~ sei,                   # PET term
  random = list(~1 | studyid, ~1 | Key),
  method = "REML",
  data   = metadata_change
)
summary(pet_filter)

# 3.Fit the three-level PEESE model
peese_filter <- rma.mv(
  yi     = cohens_d,
  V      = v,
  mods   = ~ I(sei^2),              # PEESE term
  random = list(~1 | studyid, ~1 | Key),
  method = "REML",
  data   = metadata_change
)
summary(peese_filter)

# 4.Apply Stanley & Doucouliagos’ decision rule
if (pet_filter$pval[1] > .05) {           # PET nonsignificant
  adj_est <- pet_filter$b[1]
  adj_se  <- pet_filter$se[1]
  method  <- "PET"
} else {                           # PET significant → use PEESE
  adj_est <- peese_filter$b[1]
  adj_se  <- peese_filter$se[1]
  method  <- "PEESE"
}
adj_ci <- c(
  adj_est - qnorm(.975) * adj_se,
  adj_est + qnorm(.975) * adj_se
)
sprintf(
  "%s-corrected g = %.3f (95%% CI %.3f to %.3f)",
  method, adj_est, adj_ci[1], adj_ci[2]
)



###heterogeneity
library(metafor)


# Fit a random-effects model using REML
model <- rma(cohens_d,v,method="REML",data = metadata_change)
model <- rma.mv(yi=cohens_d,V=v,
                random = list(~ 1 | studyid,     # indicate level 2
                              ~ 1 | Key), 
                data = metadata_change)

# View the model results (includes τ², I², H², Q-test)
summary(model)

# Extract heterogeneity statistics explicitly
tau2 <- model$tau2        # between-study variance
I2   <- model$I2          # proportion of variance due to heterogeneity
H2   <- model$H2          # total/sampling variance ratio
Q    <- model$QE          # Q statistic
Q_p  <- model$QEp         # p-value for Q-test

heterogeneity_results <- list(
  tau2 = tau2,
  I2 = I2,
  H2 = H2,
  Q = Q,
  Q_p_value = Q_p
)

heterogeneity_results



