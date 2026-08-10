#fig1.A
library(dplyr)
library(ggplot2)
library(scales)

# 1. Global settings
col_primary   <- "#4D7EA8"
col_primary_l <- "#A5BED6"
col_line      <- "black"

theme_p1a <- function(
    base_size = 10,                
    base_family = "sans",
    panel_border = TRUE
) {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Axis titles
      axis.title.x = element_text(size = base_size * 1.3, face = "bold",
                                  margin = margin(t = 3)),
      axis.title.y = element_text(size = base_size * 1.3, face = "bold",
                                  angle = 90,
                                  margin = margin(r = 4)),
      
      # Axis tick labels
      axis.text = element_text(size = base_size * 1.1, colour = "black"),
      
      # Panel label (A, B, C...)
      plot.tag = element_text(face = "bold",
                              size = base_size * 1.6),   
      plot.tag.position = c(0.01, 0.98),
      
      # Axes
      axis.ticks = element_line(linewidth = 0.5),        
      axis.ticks.length = unit(2, "mm"),                 
      axis.line = element_line(linewidth = 0.5),
      
      # Panel border
      panel.border = if (panel_border)
        element_rect(colour = col_line, fill = NA, linewidth = 0.5) else
          element_blank(),
      
      # Clean, publication-ready margins
      plot.margin = margin(3, 3, 3, 3, "mm"),
      
      plot.title = element_blank(),
      panel.grid = element_blank()
    )
}

# 2. Plotting function
make_figure1a <- function(metadata_change,
                          year_min = 2014,
                          year_max = 2025,
                          add_trend = TRUE) {
  
  metadata_change_dedup <- metadata_change %>%
    distinct(Key, .keep_all = TRUE)
  
  metadata_change_summary <- metadata_change_dedup %>%
    count(`Publication Year`, name = "n") %>%
    filter(`Publication Year` >= year_min,
           `Publication Year` <= year_max)
  
  all_years <- tibble(`Publication Year` = seq(year_min, year_max, 1))
  metadata_change_summary <- all_years %>%
    left_join(metadata_change_summary, by = "Publication Year") %>%
    mutate(n = if_else(is.na(n), 0L, n))
  
  p <- ggplot(metadata_change_summary,
              aes(x = `Publication Year`, y = n)) +
    geom_col(
      width     = 0.75,
      linewidth = 0.5,        # thicker bars
      colour    = col_line,
      fill      = col_primary_l
    ) +
    
    { if (add_trend) geom_line(linewidth = 0.7, colour = col_primary) } +
    { if (add_trend) geom_point(size = 2.2, colour = col_primary) } + 
    
    scale_x_continuous(
      breaks = seq(year_min, year_max, 1),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05)),
      limits = c(0, NA),
      breaks = pretty_breaks(5)
    ) +
    labs(
      tag = "A",
      x   = NULL,
      y   = "Articles included\nin meta-analysis"
    ) +
    theme_p1a(base_size = 10)   
  
  p
}

# 3. Create plot
p1a <- make_figure1a(
  metadata_change,
  year_min = 2014,
  year_max = 2025,
  add_trend = TRUE
)

p1a


#map
library(dplyr)
library(ggplot2)
library(scales)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# A specific theme
theme_map <- function(base_size = 12, base_family = "sans") {
  theme_p1a(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.title   = element_blank(),
      axis.text    = element_blank(),
      axis.ticks   = element_blank(),
      axis.line    = element_blank(),
      
      # Legend styling
      legend.title = element_text(size = base_size * 1.3, face = "bold"),
      legend.text  = element_text(size = base_size),
      legend.key.height = unit(3, "mm"),
      legend.key.width  = unit(4, "mm"),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.position   = c(0.05, 0.18),    
      legend.justification = c("left", "bottom"),
      
      panel.border = element_rect(colour = col_line, fill = NA, linewidth = 0.5),
      
      plot.margin = margin(3, 3, 3, 3, "mm"),
      
      plot.title = element_blank()
    )
}

# 1. Summarize data: number of studies per country
country_counts <- metadata_change %>%
  filter(!is.na(iso3)) %>%           # ensure iso3 not missing
  count(iso3, name = "Study_Count")

# 2. Load & merge world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

world_data <- world %>%
  left_join(country_counts, by = c("adm0_a3" = "iso3")) %>%
  mutate(
    Study_Count = if_else(is.na(Study_Count), 0L, Study_Count)
  )

max_count <- max(world_data$Study_Count, na.rm = TRUE)

# 3. Create the map
p1a_map <- ggplot(world_data) +
  geom_sf(aes(fill = Study_Count),
          colour = "grey30", size = 0.2) +  # subtle borders
  coord_sf(
    xlim   = c(-170, 190),   # trim off extreme edges
    ylim   = c(-60, 85),
    expand = FALSE
  ) +
  scale_fill_gradientn(
    colours = c("white", col_primary_l, col_primary),  # same palette family as p1a
    name    = "Number of studies",
    limits  = c(0, max_count),
    breaks  = pretty_breaks(n = 5),
    guide   = guide_colorbar(
      barheight    = unit(18, "mm"),
      barwidth     = unit(3, "mm"),
      frame.colour = "black",
      ticks.colour = "black"
    )
  ) +
  theme_map(base_size = 10)

p1a_map


#fig1.B
library(dplyr)
library(ggplot2)
library(scales)

# 1. Helpers for reordering within facets
reorder_within <- function(x, by, within, fun = mean, sep = "___") {
  paste(x, within, sep = sep) |>
    stats::reorder(by, FUN = fun)
}
scale_y_reordered <- function(sep = "___") {
  scale_y_discrete(labels = function(x) gsub(paste0(sep, ".*$"), "", x))
}

# 2. theme (lollipop version)
theme_nature_lollipop <- function(base_size = 12, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.title.x = element_text(size = base_size * 1.3, face = "bold",
                                  margin = margin(t = 4)),
      axis.title.y = element_blank(),
      axis.text.x  = element_text(size = base_size * 1.1),
      axis.text.y  = element_text(size = base_size * 1.0, hjust = 1),
      
      plot.tag = element_text(size = base_size * 1.6,
                              face = "bold",
                              hjust = 0, vjust = 1),
      plot.tag.position = c(0.01, 0.98),
      
      axis.ticks.length = unit(1.5, "mm"),
      axis.ticks = element_line(linewidth = 0.3),
      axis.line  = element_line(linewidth = 0.4),
      
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
      
      strip.background = element_blank(),
      strip.text = element_text(size = base_size * 1.2, face = "bold"),
      
      legend.position = "none",
      plot.margin = margin(3, 3, 3, 3, "mm")
    )
}

# 3. Variable labels & colours
variables <- c("sample_characteristics", "product_category", "TLL", "channel")

var_names <- c(
  "Sample Characteristics",
  "Product Type",
  "Label Design", 
  "Research Setting"
)
names(var_names) <- variables

# color for each panel (similar style to the example figure)
var_cols <- c(
  "Sample Characteristics" = "#4D7EA8",  
  "Product Type"           = "#4D7EA8",  
  "Label Design"           = "#4D7EA8",  
  "Research Setting"       = "#4D7EA8"
)

# 4. Data processing
list_dfs <- lapply(variables, function(v) {
  metadata_change %>%
    count(!!sym(v)) %>%
    filter(!is.na(!!sym(v))) %>%              # drop NA categories
    mutate(
      category   = as.character(!!sym(v)),    # the REAL category names
      proportion = n / sum(n),
      variable   = v
    )
})

df_long <- bind_rows(list_dfs) %>%
  mutate(
    variable_label = var_names[variable],
    variable_label = factor(variable_label, levels = var_names),
    prop           = proportion,
    prop_label     = scales::percent(prop, accuracy = 0.1)
  )

df_long <- df_long %>%
  mutate(category_ordered = reorder_within(category, prop, variable_label))

# 5. Final lollipop plot (p1b)
p1b <- ggplot(df_long,
              aes(x = prop,
                  y = category_ordered,
                  colour = variable_label)) +
  
  geom_segment(aes(x = 0, xend = prop,
                   y = category_ordered, yend = category_ordered),
               linewidth = 0.6,
               alpha = 0.7) +
  
  geom_point(aes(fill = variable_label),
             size = 3,
             shape = 21,
             stroke = 0.4) +
  
  geom_text(aes(label = prop_label),          # <-- ADD
            hjust = -0.2,                     # push text to the right of the point
            size = 3.3,
            colour = "black",
            inherit.aes = TRUE) +
  
  facet_wrap(~ variable_label,
             ncol = 2,
             scales = "free_y") +
  
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = scales::number_format(accuracy = 0.01),
    expand = expansion(mult = c(0, 0.12))
  ) +
  scale_y_reordered() +
  scale_colour_manual(values = var_cols) +
  scale_fill_manual(values = var_cols) +
  
  labs(
    tag = "B",
    x = "Proportion of studies",
    y = NULL
  ) +
  
  theme_nature_lollipop()

p1b




#fig2
library(dplyr)
library(stringr)
library(metafor)
library(ggplot2)
library(patchwork)

# 0) Define categories
lab_channels   <- c("Lab/Sim-store UI")
field_channels <- c("Offline foodservice", "Offline retail", "Online retail")

food_products     <- c("Packaged food & grocery", "Ready-to-eat meals")
non_food_products <- c("Non-food consumer goods")

# 1) Build labels + types + rank order
metadata_change_add <- metadata_change %>%
  mutate(
    FirstAuthor  = str_extract(Author, "^[^,]+"),
    ArticleLabel = paste0(FirstAuthor, " et al. (", `Publication Year`, ")"),
    
    Setting_Type = case_when(
      channel %in% lab_channels   ~ "Lab",
      channel %in% field_channels ~ "Field",
      TRUE ~ NA_character_
    ),
    
    Product_Type = case_when(
      product_category %in% food_products     ~ "Food",
      product_category %in% non_food_products ~ "Non-food",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(Key) %>%
  arrange(Key) %>%
  mutate(
    es_in_article = row_number(),
    k_in_article  = n(),
    StudyLabel = ifelse(
      k_in_article > 1,
      paste0(ArticleLabel, " — ES", es_in_article),
      ArticleLabel
    )
  ) %>%
  ungroup() %>%
  arrange(cohens_d)

# 2) Fit multilevel model (unchanged)
r <- rma.mv(
  yi = cohens_d, V = v,
  random = list(~ 1 | studyid, ~ 1 | Key),
  data = metadata_change_add,
  method = "REML"
)

# 3) Study CIs + formatted text
plot_studies <- metadata_change_add %>%
  mutate(
    sei   = sqrt(v),
    ci.lb = cohens_d - 1.96 * sei,
    ci.ub = cohens_d + 1.96 * sei,
    est_ci_text = sprintf("%.2f [%.2f, %.2f]", cohens_d, ci.lb, ci.ub)
  )

# Subtle metafor-like point-size scaling by inverse variance
rescale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (isTRUE(all.equal(rng[1], rng[2]))) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}
w <- 1 / plot_studies$v
plot_studies <- plot_studies %>%
  mutate(point_size = 2.0 + 2.5 * rescale01(w))

# Pooled
pred <- predict(r)
plot_pooled <- tibble(
  StudyLabel  = "Random-effects model",
  cohens_d    = as.numeric(pred$pred),
  ci.lb       = as.numeric(pred$ci.lb),
  ci.ub       = as.numeric(pred$ci.ub),
  est_ci_text = sprintf("%.2f [%.2f, %.2f]", pred$pred, pred$ci.lb, pred$ci.ub)
)

# 4) y positions
n_study <- nrow(plot_studies)

plot_studies <- plot_studies %>% mutate(y = rev(row_number()))
plot_pooled  <- plot_pooled  %>% mutate(y = 0)

y_min <- -2.2
y_max <- n_study + 1.2

# 5) X axis control
x_plot_min <- -1.5
x_plot_max <-  1.5
x_breaks   <- seq(x_plot_min, x_plot_max, by = 0.5)

# 6) Themes
base_fs <- 12

theme_table_col <- theme_void(base_size = base_fs) +
  theme(
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

theme_forest <- theme_classic(base_size = base_fs) +
  theme(
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid       = element_blank(),
    axis.line.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.text.y      = element_blank(),
    axis.title.y     = element_blank(),
    legend.position  = "bottom",
    legend.box       = "vertical",
    legend.title     = element_text(size = base_fs - 1, face = "bold"),
    legend.text      = element_text(size = base_fs - 1)
  )

# 7) LEFT: Study labels column (LEFT-aligned, tight to middle)
#    KEY: remove fixed x-limits so the label column doesn't waste width
p_left <- ggplot() +
  # header
  annotate("text",
           x = 0, y = n_study + 0.9,
           label = "Study",
           hjust = 0, fontface = "bold", size = 4) +
  
  # labels
  geom_text(
    data = plot_studies,
    aes(x = 0, y = y, label = StudyLabel),
    hjust = 0, size = 4
  ) +
  geom_text(
    data = plot_pooled,
    aes(x = 0, y = y, label = StudyLabel),
    hjust = 0, size = 4
  ) +
  
  # IMPORTANT: no x limits; tiny right padding only
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(
    limits = c(y_min, y_max),
    breaks = NULL,
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  
  theme_table_col +
  theme(
    # tight to middle: no right margin
    plot.margin = margin(t = 10, r = 0, b = 6, l = 10)
  )

# 8) MIDDLE: Forest panel (x-axis ONLY here)
p_mid <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6) +
  
  geom_errorbarh(
    data = plot_studies,
    aes(y = y, xmin = ci.lb, xmax = ci.ub,
        colour = Setting_Type, linetype = Setting_Type),
    height = 0.16, linewidth = 0.6
  ) +
  geom_point(
    data = plot_studies,
    aes(y = y, x = cohens_d,
        colour = Setting_Type, shape = Product_Type, size = point_size),
    stroke = 0.7
  ) +
  
  # pooled CI + diamond
  geom_errorbarh(
    data = plot_pooled,
    aes(y = y, xmin = ci.lb, xmax = ci.ub),
    height = 0.18, linewidth = 0.8, colour = "black"
  ) +
  geom_point(
    data = plot_pooled,
    aes(y = y, x = cohens_d),
    shape = 18, size = 4.6, colour = "black"
  ) +
  
  scale_x_continuous(
    limits = c(x_plot_min, x_plot_max),
    breaks = x_breaks
  ) +
  scale_y_continuous(
    limits = c(y_min, y_max),
    breaks = NULL,
    expand = c(0, 0)
  ) +
  
  scale_shape_manual(values = c("Food" = 16, "Non-food" = 15)) +
  scale_colour_manual(values = c("Field" = "#0072B2", "Lab" = "#E69F00")) +
  scale_linetype_manual(values = c("Field" = "dashed", "Lab" = "solid")) +
  scale_size_identity(guide = "none") +
  
  labs(
    x = "Effect size (Cohen's d)",
    colour = "Setting",
    linetype = "Setting",
    shape = "Product type"
  ) +
  
  theme_forest +
  theme(
    axis.title.x = element_text(size = base_fs),
    axis.text.x  = element_text(size = base_fs - 1),
    # tight to left + right panels
    plot.margin  = margin(t = 10, r = 0, b = 6, l = 0)
  )

# 9) RIGHT: numeric column (right-aligned, tight to middle)
p_right <- ggplot() +
  annotate("text",
           x = 1, y = n_study + 0.9,
           label = "Estimate [95% CI]",
           hjust = 1, fontface = "bold", size = 4) +
  
  geom_text(
    data = plot_studies,
    aes(x = 1, y = y, label = est_ci_text),
    hjust = 1, size = 4
  ) +
  geom_text(
    data = plot_pooled,
    aes(x = 1, y = y, label = est_ci_text),
    hjust = 1, size = 4
  ) +
  
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(
    limits = c(y_min, y_max),
    breaks = NULL,
    expand = c(0, 0)
  ) +
  theme_table_col +
  theme(
    plot.margin = margin(t = 10, r = 10, b = 6, l = 0)
  )

# 10) Combine (reduced gaps via margins + tuned widths)
final_plot <- (p_left + p_mid + p_right) +
  plot_layout(
    widths = c(0.36, 0.48, 0.16),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.justification = "center"
  )

final_plot

ggsave(
  "fig2.png",
  final_plot,
  width  = 240,
  height = 280,
  units  = "mm",
  dpi    = 600,
  bg     = "white",
  device = ragg::agg_png
)



#fig3
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


funnel(
  rm,
  xlab = "Effect size (Cohen's d)",
  ylab = "Standard error",
  level = c(90, 95, 99),
  shade = c("white", "gray90", "gray80"),
  refline = 0
)

legend(
  "topright",
  legend = c("p < 0.10", "p < 0.05", "p < 0.01"),
  fill = c("gray80", "gray90", "white"),
  cex = 0.8
)

plot(
  dat$sei,
  dat$cohens_d,
  xlab = "Standard error",
  ylab = "Effect size (Cohen's d)",
  pch = 19,
  cex = 0.7
)

abline(
  lm(cohens_d ~ sei, data = dat),
  lwd = 2
)

abline(h = 0, lty = 2)

