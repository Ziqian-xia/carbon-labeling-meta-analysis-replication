#figure s1
# Define model list
modellist <- c(rep('Random Effect', 4), rep('Multilevel Model', 2))
estlist <- c('DL', 'HS', 'ML', 'REML', 'ML', 'REML')
modolres <- data.frame(mlst = modellist, estlist)

# Sensitivity Analysis
Vmat0<-impute_covariance_matrix(metadata_change$v,
                                cluster = metadata_change$Key,
                                r=0.6)
modeldata<-data.frame()
rand1<-conf_int(obj=rma(cohens_d,v,method="DL",data = metadata_change),
                vcov = 'CR2',cluster=metadata_change$Key)
rand2<-conf_int(obj=rma(cohens_d,v,method="HS",data = metadata_change),
                vcov = 'CR2',cluster=metadata_change$Key)
rand3<-conf_int(obj=rma(cohens_d,v,method="ML",data = metadata_change),
                vcov = 'CR2',cluster=metadata_change$Key)
rand4<-conf_int(obj=rma(cohens_d,v,method="REML",data = metadata_change),
                vcov = 'CR2',cluster=metadata_change$Key)
mv1<-conf_int(obj=rma.mv(cohens_d ,
                         V = Vmat0 ,                 # sampling variances
                         random = list(~ 1 | studyid,     # indicate level 2
                                       ~ 1 | Key), # indicate level 3
                         data = metadata_change,
                         tdist = TRUE,
                         method = "ML"),vcov = 'CR2')
mv2<-conf_int(obj=rma.mv(cohens_d,
                         V=Vmat0,
                         method = 'REML',
                         random = list(~ 1 | studyid,     # indicate level 2
                                       ~ 1 | Key),
                         tdist = TRUE,
                         data=metadata_change),vcov = 'CR2')

modeldata <- bind_cols(
  modolres, bind_rows(rand1, rand2, rand3, rand4, mv1, mv2)
) %>%
  mutate(esttype = paste(mlst," with ", estlist," estimator"))  # Creating `esttype`

# Create a new column for faceting (Random Effects on top, Multilevel Models on bottom)
modeldata <- modeldata %>%
  mutate(model_type = factor(ifelse(mlst == "Random Effect", "Random Effects Model", "Multilevel Model"),
                             levels = c("Random Effects Model", "Multilevel Model")))  # Correct ordering
# Plot
s1 <- ggplot(modeldata, aes(y = beta, x = esttype, ymin = CI_L, ymax = CI_U)) +
  geom_pointrange(color = "#D74B4B", size = 1) +  # Red confidence interval lines
  geom_point(color = "#D74B4B") +  # Black dots for point estimates
  geom_vline(xintercept = 0, lty = 2) +  # Dotted line at x=0
  theme_bw() +
  coord_flip() +
  facet_wrap(~ model_type, ncol = 1, scales = "free_y") +  # Facet with no gap
  theme(
    plot.title = element_text(hjust = 0, size = 16, face = "bold"), 
    strip.text = element_text(size = 15, face = "bold"),  # Facet title styling
    axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.8,face = "bold"),
    legend.position = 'none',  # Remove legend
    text = element_text(size = 15,face = "bold"),
    axis.title.y = element_blank(),
    panel.spacing = unit(0, "lines")# Remove gap between facets
  ) +
  labs(y = "Cluster robust estimate (95% CI)", x = "Different model and estimation methods")

# Display the plot
print(s1)
ggsave(file.path(output_dir, "Figure_S1_sensitivity.png"), s1,
       width = 10, height = 7.93, units = "in", dpi = 300, bg = "white")

#figure s2
# Accumulate by publication year; studyid breaks ties within each year.
s2_data <- metadata_change[order(metadata_change$`Publication Year`, metadata_change$studyid), ]
m1 <- rma(yi = cohens_d, vi = v, data = s2_data, method = "REML")
s2_cumulative <- cumul(m1)
s2_points <- data.frame(year = s2_data$`Publication Year`, studyid = s2_data$studyid,
                        effect_size = s2_cumulative$estimate, I2 = s2_cumulative$I2)
s2_points$I2[1] <- NA_real_  # Heterogeneity is undefined for a single estimate.
s2 <- ggplot(s2_points, aes(x = effect_size, y = I2)) +
  geom_path(colour = "grey75", linewidth = 0.5, na.rm = TRUE) +
  geom_point(aes(colour = year), size = 2.5, na.rm = TRUE) +
  scale_colour_gradient(low = "grey75", high = "black", name = "Publication year",
                        breaks = c(2014, 2016, 2019, 2022, 2025)) +
  labs(x = "Cumulative effect size (Cohen's d)", y = "Heterogeneity (I², %)" ) +
  theme_bw(base_size = 14)
ggsave(file.path(output_dir, "Figure_S2_cumulative.png"), s2,
       width = 8, height = 6.9, units = "in", dpi = 300, bg = "white")
write.csv(s2_points, file.path(output_dir, "Figure_S2_values.csv"), row.names = FALSE)

#figure s3
rm <- rma(yi=cohens_d, vi=v, data = metadata_change)
png(file.path(output_dir, "Figure_S3_influence.png"),
    width = 2400, height = 1996, res = 300)
influence.rma.uni(rm) %>% plot()
dev.off()

#figure s4: subgroup analysis
metadata_change$scenario <- as.factor(metadata_change$scenario)
metadata_change$product_category <- as.factor(metadata_change$product_category)
metadata_change$sample_characteristics <- as.factor(metadata_change$sample_characteristics)
metadata_change$TLL <- as.factor(metadata_change$TLL)
metadata_change$research_setting <- as.factor(metadata_change$research_setting)

m1 <- rma.mv(yi = cohens_d, V = v, mods = ~ scenario, data = metadata_change,random = list(~ 1 | studyid,     # indicate level 2
                                                                                          ~ 1 | Key))
m2 <- rma.mv(yi = cohens_d, V = v, mods = ~ product_category, data = metadata_change,random = list(~ 1 | studyid,     # indicate level 2
                                                                                                   ~ 1 | Key))
m3 <- rma.mv(yi = cohens_d, V = v, mods = ~ sample_characteristics, data = metadata_change,random = list(~ 1 | studyid,     # indicate level 2
                                                                                                         ~ 1 | Key))
m4 <- rma.mv(yi = cohens_d, V = v, mods = ~ TLL, data = metadata_change,random = list(~ 1 | studyid,     # indicate level 2
                                                                                      ~ 1 | Key))
m5 <- rma.mv(yi = cohens_d, V = v, mods = ~ research_setting, data = metadata_change,random = list(~ 1 | studyid,     # indicate level 2
                                                                                      ~ 1 | Key))
# Predict each group's mean from the same moderator models.
# Non-intercept coefficients alone are contrasts, not subgroup means.
subgroup_predictions <- function(fit, variable) {
  lev <- levels(metadata_change[[variable]])
  nd <- setNames(data.frame(factor(lev, levels = lev)), variable)
  X <- model.matrix(reformulate(variable), nd)
  pr <- predict(fit, newmods = X[, -1, drop = FALSE])
  data.frame(factor = lev, effect_size = pr$pred, lower_ci = pr$ci.lb, upper_ci = pr$ci.ub)
}
subgroup_results1 <- subgroup_predictions(m1, "scenario")
subgroup_results2 <- subgroup_predictions(m2, "product_category")
subgroup_results3 <- subgroup_predictions(m3, "sample_characteristics")
subgroup_results4 <- subgroup_predictions(m4, "TLL")
subgroup_results5 <- subgroup_predictions(m5, "research_setting")

# Add a panel indicator
subgroup_results1$panel <- "Scenario"
subgroup_results2$panel <- "Type of\nProduct"
subgroup_results3$panel <- "Sample\nCharacteristics"
subgroup_results4$panel <- "Label\nDesign"
subgroup_results5$panel <- "Research\nSetting"



# Combine data
combined_results <- bind_rows( subgroup_results1,subgroup_results2,subgroup_results3,subgroup_results4,subgroup_results5)

# Sort data by effect size for better visualization
combined_results <- combined_results %>%
  arrange(desc(panel), effect_size)



###creat the plot
s4 <- ggplot(combined_results, aes(x = effect_size, y = reorder(factor, desc(panel) * effect_size))) +
  geom_point(aes(color = panel), size = 4, shape = 18) +  # Diamond shapes
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, color = panel), height = 0.15) +  # Slim error bars
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +  # Reference line at zero
  facet_grid(panel ~ ., scales = "free_y", space = "free") +  # Separate panels
  labs(title = "",x = "Model-estimated subgroup effect (Cohen's d; 95% CI)", y = "") +  # No title
  scale_color_manual(values = c("Scenario" = "#a6cee3", "Research\nSetting" = "#1f78b4", "Type of\nProduct" = "#b2df8a","Sample\nCharacteristics" = "#33a02c","Label\nDesign"="#fb9a99")) +  # Custom colors
  theme_minimal(base_size = 15) +  # Adjust base font size for readability
  theme(
    plot.title = element_text(face = "bold", hjust = 0),  # Left-align the title
    plot.title.position = "plot",  # Position title in the upper-left corner
    strip.background = element_rect(fill = "gray90", color = "black", linewidth = 1),  # Panel label background with border
    strip.text = element_text(face = "bold", hjust = 0.5),  # Centered panel labels
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  
    axis.text.y = element_text(size = 14),  # Adjust y-axis text size
    axis.title.x = element_text(size = 15, face = "bold"),  # Make x-axis label bolder
    plot.margin = margin(10, 10, 10, 10),  # Adjust plot margins
    legend.position = "none"  
  )
s4
ggsave(file.path(output_dir, "Figure_S4_subgroups.png"), s4,
       width = 10, height = 11.2, units = "in", dpi = 300, bg = "white")

write.csv(combined_results, file.path(output_dir, "Figure_S4_values.csv"), row.names = FALSE)

#figure s5: illustration only; the primary model is unchanged.
source("code/17_illustrative_diagnostics.R")

#table s3: Omnibus
library(dplyr)
library(emmeans)
library(openxlsx)

# Ensure categorical variables are factors
metadata_change <- metadata_change %>%
  mutate(
    channel = factor(channel),
    product_category = factor(product_category),
    sample_characteristics = factor(sample_characteristics),
    TLL = factor(TLL)
    
  )

# Omnibus ANOVA Test
omnibus_test <- aov(cohens_d ~ channel + product_category + sample_characteristics + TLL , data = metadata_change)
omnibus_summary <- summary(omnibus_test)

# Separate planned comparisons for each subgroup

# channel comparisons
channel_comparisons <- emmeans(omnibus_test, pairwise ~ channel)
channel_comparisons_summary <- summary(channel_comparisons$contrasts)

# product_category comparisons
product_category_comparisons <- emmeans(omnibus_test, pairwise ~ product_category)
product_category_comparisons_summary <- summary(product_category_comparisons$contrasts)

# Sample characteristics comparisons
sample_characteristics_comparisons <- emmeans(omnibus_test, pairwise ~ sample_characteristics)
sample_characteristics_comparisons_summary <- summary(sample_characteristics_comparisons$contrasts)

# TLL comparisons
TLL_comparisons <- emmeans(omnibus_test, pairwise ~ TLL)
TLL_comparisons_summary <- summary(TLL_comparisons$contrasts)


# Export all comparisons and omnibus results to Excel
wb <- createWorkbook()

# Add omnibus test summary
addWorksheet(wb, "Omnibus Test")
writeData(wb, "Omnibus Test", as.data.frame(summary(omnibus_test)[[1]]), rowNames = TRUE)

# Add channel comparisons
addWorksheet(wb, "channel Comparisons")
writeData(wb, "channel Comparisons", channel_comparisons_summary)

# Add type of product comparisons
addWorksheet(wb, "Type of Product Comparisons")
writeData(wb, "Type of Product Comparisons", product_category_comparisons_summary)

# Add sample characteristics comparisons
addWorksheet(wb, "Sample Char. Comparisons")
writeData(wb, "Sample Char. Comparisons", sample_characteristics_comparisons_summary)

# Add TLL comparisons
addWorksheet(wb, "TLL Comparisons")
writeData(wb, "TLL Comparisons", TLL_comparisons_summary)


# Save workbook
saveWorkbook(wb, "Planned_Comparisons_and_Omnibus.xlsx", overwrite = TRUE)



#table s5: risk-of-bias domains
library(metafor)
library(dplyr)
library(stringr)
library(officer)
library(flextable)

# 1) Helper: p-value formatting
fmt_p <- function(p) {
  dplyr::case_when(
    is.na(p)        ~ NA_character_,
    p < 0.0001      ~ "<0.0001",
    p < 0.001       ~ formatC(p, format = "f", digits = 4),
    TRUE            ~ formatC(p, format = "f", digits = 3)
  )
}

# 2) Helper: extract one rma() model to tidy rows
extract_rma_table <- function(fit, model_name) {
  
  # coef table
  ct <- as.data.frame(fit$beta)
  colnames(ct) <- "Estimate"
  ct$SE      <- fit$se
  ct$`z-value` <- fit$zval
  ct$`p-value` <- fit$pval
  ct$CI.lb   <- fit$ci.lb
  ct$CI.ub   <- fit$ci.ub
  ct$Moderator <- rownames(ct)
  
  # model-level QM
  qm_val <- unname(fit$QM)
  qm_p   <- unname(fit$QMp)
  
  ct %>%
    mutate(
      Model = model_name,
      QM    = qm_val,
      QM.p  = qm_p
    ) %>%
    select(Model, Moderator, Estimate, SE, `z-value`, `p-value`, CI.lb, CI.ub, QM, QM.p)
}

# 3) Fit models
fit_D1 <- rma(yi = cohens_d, vi = v, mods = ~ D1, data = data_regression_rob, method = "REML")
fit_D2 <- rma(yi = cohens_d, vi = v, mods = ~ D2, data = data_regression_rob, method = "REML")
fit_D3 <- rma(yi = cohens_d, vi = v, mods = ~ D3, data = data_regression_rob, method = "REML")
fit_D4 <- rma(yi = cohens_d, vi = v, mods = ~ D4, data = data_regression_rob, method = "REML")
fit_D5 <- rma(yi = cohens_d, vi = v, mods = ~ D5, data = data_regression_rob, method = "REML")

# 4) Build combined table
tab <- bind_rows(
  extract_rma_table(fit_D1, "D1"),
  extract_rma_table(fit_D2, "D2"),
  extract_rma_table(fit_D3, "D3"),
  extract_rma_table(fit_D4, "D4"),
  extract_rma_table(fit_D5, "D5")
) %>%
  # cleaner moderator names for publication
  mutate(
    Moderator = case_when(
      Moderator == "intrcpt" ~ "Intercept",
      TRUE ~ Moderator
    ),
    # numeric formatting
    Estimate = formatC(Estimate, format = "f", digits = 3),
    SE       = formatC(SE,       format = "f", digits = 3),
    `z-value`  = formatC(`z-value`,  format = "f", digits = 3),
    `p-value`  = fmt_p(as.numeric(`p-value`)),
    CI.lb    = formatC(CI.lb,    format = "f", digits = 3),
    CI.ub    = formatC(CI.ub,    format = "f", digits = 3),
    QM       = formatC(QM,       format = "f", digits = 3),
    QM.p     = fmt_p(as.numeric(QM.p))
  )

write.csv(tab, file.path(output_dir, "Table_S5_risk_of_bias.csv"), row.names = FALSE)

ft <- flextable(tab)

# Header labels
ft <- set_header_labels(
  ft,
  Model = "Model",
  Moderator = "Moderator",
  Estimate = "Estimate",
  SE = "SE",
  `z-value` = "z-value",
  `p-value` = "p-value",
  CI.lb = "CI.lb",
  CI.ub = "CI.ub",
  QM = "QM",
  QM.p = "QM.p"
)

ft <- theme_booktabs(ft)
ft <- font(ft, fontname = "Arial", part = "all")
ft <- fontsize(ft, size = 10, part = "all")
ft <- bold(ft, part = "header")
ft <- align(ft, align = "left",  j = c("Model", "Moderator"), part = "all")
ft <- align(ft, align = "right", j = c("Estimate","SE","z-value","p-value","CI.lb","CI.ub","QM","QM.p"), part = "all")
ft <- autofit(ft)

ft <- border_remove(ft)
ft <- hline_top(ft, border = fp_border(width = 1), part = "header")
ft <- hline(ft, border = fp_border(width = 0.75), part = "header")

# 5) Write to Word
doc <- read_docx()

doc <- body_add_par(doc, "Table X | Meta-regression results.", style = "Normal")
doc <- body_add_par(
  doc,
  "Mixed-effects meta-regression models (REML). QM and QM.p are omnibus tests of moderators within each model.",
  style = "Normal"
)
doc <- body_add_flextable(doc, ft)
doc <- body_add_par(doc, "", style = "Normal")

print(doc, target = "meta_regression_results_table.docx")

