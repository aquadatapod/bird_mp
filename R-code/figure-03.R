# ==============================================================================
# FIGURE 3: WING-LENGTH THRESHOLD AND HABITAT PATTERNS IN MICROPLASTIC BURDEN
#
# Panel A: Wing-length threshold among terrestrial species
# Panel B: MP burden across wing length with tree-derived cutoff
# Panel C: Habitat-level MP burden normalized by body mass
# ==============================================================================

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(mgcv)
library(rstatix)
library(multcompView)

FIGURE_NUMBER <- "Figure3"

# ==============================================================================
# 1. HELPER FUNCTION
# ==============================================================================

get_mode <- function(x) {
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  names(sort(table(x), decreasing = TRUE))[1]
}

# ==============================================================================
# 2. LOAD FINAL ROW-LEVEL DATA
# ==============================================================================

tree_dat_raw <- read.csv(
  here("data", "processed", "row_level.csv"),
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
) %>%
  drop_na(
    Species,
    Mppersp,
    Beak.Length_Culmen,
    Wing.Length,
    Mass
  ) %>%
  mutate(
    across(
      c(Habitat, Migration, Trophic.Niche),
      as.factor
    )
  )

# ==============================================================================
# 3. SPECIES-LEVEL AGGREGATION
#
# This is the same species-level dataset used for the final decision tree.
# Continuous variables = median
# Categorical variables = mode
# ==============================================================================

tree_dat <- tree_dat_raw %>%
  group_by(Species) %>%
  summarise(
    Mppersp            = median(Mppersp, na.rm = TRUE),
    Beak.Length_Culmen = median(Beak.Length_Culmen, na.rm = TRUE),
    Wing.Length        = median(Wing.Length, na.rm = TRUE),
    Mass               = median(Mass, na.rm = TRUE),
    Habitat            = get_mode(Habitat),
    Migration          = get_mode(Migration),
    Trophic.Niche      = get_mode(Trophic.Niche),
    .groups            = "drop"
  ) %>%
  mutate(
    across(
      c(Habitat, Migration, Trophic.Niche),
      as.factor
    ),
    logMP   = log1p(Mppersp),
    logMass = log10(Mass),
    logBeak = log10(Beak.Length_Culmen),
    MP_per_g = Mppersp / Mass
  )

# ==============================================================================
# 4. DATA CONSISTENCY CHECKS
# ==============================================================================

cat("\n========================================\n")
cat("FIGURE 3 DATA CONSISTENCY CHECK\n")
cat("========================================\n")

cat("Raw complete-case records:", nrow(tree_dat_raw), "\n")
cat("Final species:", nrow(tree_dat), "\n")
cat("Unique species:", n_distinct(tree_dat$Species), "\n")
cat(
  "Duplicate species:",
  sum(duplicated(tree_dat$Species)),
  "\n"
)

cat("\nMissing values in final dataset:\n")
print(colSums(is.na(tree_dat)))

cat("\nSpecies by habitat:\n")
print(
  tree_dat %>%
    count(Habitat, sort = TRUE)
)

# The final decision-tree dataset should contain 101 species.
stopifnot(nrow(tree_dat) == 101)
stopifnot(n_distinct(tree_dat$Species) == 101)

# Save the exact species-level dataset used for Figure 3
write.csv(
  tree_dat,
  here(
    "data",
    "processed",
    "figure3_species_level.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 5. TREE-DERIVED WING-LENGTH THRESHOLD
# ==============================================================================

WING_THRESHOLD <- 163

# The 163-mm split occurs after the habitat split in the decision tree.
# Therefore, Panel A and Panel B use the terrestrial branch only.

aquatic_habitats <- c(
  "Aquatic",
  "Coastal",
  "Marine",
  "Wetland"
)

terrestrial_data <- tree_dat %>%
  filter(
    !Habitat %in% aquatic_habitats
  ) %>%
  mutate(
    Wing_group = factor(
      ifelse(
        Wing.Length < WING_THRESHOLD,
        "small",
        "large"
      ),
      levels = c("small", "large")
    )
  )

cat(
  "\nTerrestrial species:",
  nrow(terrestrial_data),
  "\n"
)

cat(
  "Wing < 163 mm:",
  sum(terrestrial_data$Wing_group == "small"),
  "\n"
)

cat(
  "Wing >= 163 mm:",
  sum(terrestrial_data$Wing_group == "large"),
  "\n"
)

# ==============================================================================
# 6. PANEL A: WING-LENGTH THRESHOLD COMPARISON
# ==============================================================================

n_small <- sum(
  terrestrial_data$Wing_group == "small"
)

n_large <- sum(
  terrestrial_data$Wing_group == "large"
)

wilcox_res <- terrestrial_data %>%
  wilcox_test(
    Mppersp ~ Wing_group
  )

effsize_res <- terrestrial_data %>%
  wilcox_effsize(
    Mppersp ~ Wing_group
  )

p_val_a <- wilcox_res$p
r_eff <- effsize_res$effsize

label_small <- sprintf(
  "<%d mm\n(n = %d)",
  WING_THRESHOLD,
  n_small
)

label_large <- sprintf(
  "\u2265%d mm\n(n = %d)",
  WING_THRESHOLD,
  n_large
)

pA <- ggplot(
  terrestrial_data,
  aes(
    x = Wing_group,
    y = Mppersp,
    fill = Wing_group
  )
) +
  
  geom_boxplot(
    alpha = 0.7,
    outlier.shape = NA,
    width = 0.5,
    linewidth = 0.3
  ) +
  
  geom_jitter(
    width = 0.10,
    alpha = 0.4,
    size = 0.9,
    color = "gray20"
  ) +
  
  scale_y_log10(
    breaks = c(
      1,
      10,
      100,
      1000
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "small" = "#A23B72",
      "large" = "#2E86AB"
    )
  ) +
  
  scale_x_discrete(
    labels = c(
      "small" = label_small,
      "large" = label_large
    )
  ) +
  
  labs(
    x = "",
    y = "MPs per individual",
    subtitle = sprintf(
      "Terrestrial species; Wilcoxon %s, r = %.2f",
      ifelse(
        p_val_a < 0.001,
        "p < 0.001",
        sprintf(
          "p = %.3f",
          p_val_a
        )
      ),
      r_eff
    )
  ) +
  
  theme_classic(
    base_size = 9
  ) +
  
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      size = 8,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      color = "black"
    ),
    axis.title.y = element_text(
      size = 9,
      color = "black"
    ),
    plot.subtitle = element_text(
      size = 7.5,
      face = "italic"
    )
  )

# ==============================================================================
# 7. PANEL B: MP BURDEN ACROSS WING LENGTH
# ==============================================================================

# Linear model
m_lin <- lm(
  log10(Mppersp + 1) ~ log10(Wing.Length),
  data = terrestrial_data
)

# GAM
m_gam <- mgcv::gam(
  log10(Mppersp + 1) ~
    s(
      log10(Wing.Length),
      k = 4
    ),
  data = terrestrial_data,
  method = "REML"
)

edf <- summary(m_gam)$s.table[
  1,
  "edf"
]

use_gam <- (
  AIC(m_lin) - AIC(m_gam) > 2 &&
    edf > 1.5
)

# Prediction grid
pred_grid <- data.frame(
  Wing.Length = seq(
    min(
      terrestrial_data$Wing.Length,
      na.rm = TRUE
    ),
    max(
      terrestrial_data$Wing.Length,
      na.rm = TRUE
    ),
    length.out = 200
  )
)

# Model prediction
if (use_gam) {
  
  pred <- predict(
    m_gam,
    newdata = pred_grid,
    se.fit = TRUE
  )
  
  model_label <- sprintf(
    "GAM: edf = %.1f",
    edf
  )
  
} else {
  
  pred <- predict(
    m_lin,
    newdata = pred_grid,
    se.fit = TRUE
  )
  
  model_label <- sprintf(
    "Linear: R\u00b2 = %.2f",
    summary(m_lin)$r.squared
  )
}

pred_grid$fit <- pred$fit

pred_grid$lwr <- (
  pred$fit -
    1.96 * pred$se.fit
)

pred_grid$upr <- (
  pred$fit +
    1.96 * pred$se.fit
)

# Panel B
pB <- ggplot(
  terrestrial_data,
  aes(
    x = Wing.Length,
    y = log10(Mppersp + 1)
  )
) +
  
  geom_ribbon(
    data = pred_grid,
    aes(
      x = Wing.Length,
      ymin = lwr,
      ymax = upr
    ),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.3
  ) +
  
  geom_line(
    data = pred_grid,
    aes(
      x = Wing.Length,
      y = fit
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.8
  ) +
  
  geom_vline(
    xintercept = WING_THRESHOLD,
    linetype = "dashed",
    color = "red",
    linewidth = 0.4
  ) +
  
  annotate(
    "text",
    x = WING_THRESHOLD,
    y = max(
      log10(
        terrestrial_data$Mppersp + 1
      ),
      na.rm = TRUE
    ),
    label = sprintf(
      "Tree split\n%d mm",
      WING_THRESHOLD
    ),
    size = 2.5,
    color = "red",
    hjust = -0.1,
    vjust = 1
  ) +
  
  geom_point(
    alpha = 0.7,
    size = 1.6
  ) +
  
  labs(
    x = "Wing length (mm)",
    y = expression(
      "log"[10] * "(MP + 1)"
    ),
    subtitle = sprintf(
      "%s; terrestrial species (n = %d)",
      model_label,
      nrow(terrestrial_data)
    )
  ) +
  
  theme_classic(
    base_size = 9
  ) +
  
  theme(
    axis.text = element_text(
      size = 8,
      color = "black"
    ),
    axis.title = element_text(
      size = 9,
      color = "black"
    ),
    plot.subtitle = element_text(
      size = 7.5,
      face = "italic"
    )
  )

# ==============================================================================
# 8. PANEL C: HABITAT-LEVEL MP BURDEN NORMALIZED BY BODY MASS
# ==============================================================================

MIN_SPECIES_PER_HABITAT <- 5

habitat_counts <- tree_dat %>%
  count(
    Habitat,
    sort = TRUE
  )

habitats_kept <- habitat_counts %>%
  filter(
    n >= MIN_SPECIES_PER_HABITAT
  ) %>%
  pull(Habitat)

habitat_data <- tree_dat %>%
  filter(
    Habitat %in% habitats_kept
  ) %>%
  mutate(
    Habitat = factor(
      Habitat,
      levels = habitats_kept
    )
  )

# Kruskal-Wallis test
kw_test <- kruskal.test(
  MP_per_g ~ Habitat,
  data = habitat_data
)

# Dunn post-hoc test
dunn_res <- habitat_data %>%
  dunn_test(
    MP_per_g ~ Habitat,
    p.adjust.method = "BH"
  )

# Pairwise adjusted p-values
pw <- dunn_res$p.adj

names(pw) <- paste(
  dunn_res$group1,
  dunn_res$group2,
  sep = "-"
)

# Compact letter display
letters_obj <- multcompView::multcompLetters(
  pw
)

letters_df <- data.frame(
  Habitat = names(
    letters_obj$Letters
  ),
  Letters = letters_obj$Letters,
  stringsAsFactors = FALSE
)

# Sample size
n_per_habitat <- habitat_data %>%
  count(Habitat)

# Letter positions
label_pos <- habitat_data %>%
  group_by(Habitat) %>%
  summarise(
    y = max(
      MP_per_g,
      na.rm = TRUE
    ) * 1.5,
    .groups = "drop"
  ) %>%
  left_join(
    letters_df,
    by = "Habitat"
  ) %>%
  left_join(
    n_per_habitat,
    by = "Habitat"
  )

# Panel C
pC <- ggplot(
  habitat_data,
  aes(
    x = Habitat,
    y = MP_per_g,
    fill = Habitat
  )
) +
  
  geom_boxplot(
    alpha = 0.7,
    outlier.shape = NA,
    width = 0.5,
    linewidth = 0.3
  ) +
  
  geom_jitter(
    width = 0.12,
    alpha = 0.3,
    size = 0.8,
    color = "gray20"
  ) +
  
  geom_text(
    data = label_pos,
    aes(
      x = Habitat,
      y = y,
      label = Letters
    ),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 3.2
  ) +
  
  scale_y_log10() +
  
  scale_x_discrete(
    labels = setNames(
      sprintf(
        "%s\n(n=%d)",
        label_pos$Habitat,
        label_pos$n
      ),
      label_pos$Habitat
    )
  ) +
  
  labs(
    x = "",
    y = "MPs per gram body mass",
    subtitle = sprintf(
      "Kruskal-Wallis %s; Dunn post-hoc (BH), n = %d species",
      ifelse(
        kw_test$p.value < 0.001,
        "p < 0.001",
        sprintf(
          "p = %.4f",
          kw_test$p.value
        )
      ),
      nrow(habitat_data)
    )
  ) +
  
  theme_classic(
    base_size = 9
  ) +
  
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      size = 7.5,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      color = "black"
    ),
    axis.title.y = element_text(
      size = 9,
      color = "black"
    ),
    plot.subtitle = element_text(
      size = 7.5,
      face = "italic"
    )
  )

# ==============================================================================
# 9. COMBINE FIGURE
# ==============================================================================

final_figure3 <- (
  pA | pB
) / pC +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(
        size = 11,
        face = "bold"
      )
    )
  )

# ==============================================================================
# 10. DISPLAY
# ==============================================================================

print(final_figure3)

# ==============================================================================
# 11. SAVE
# ==============================================================================

if (!dir.exists(here("output"))) {
  dir.create(
    here("output"),
    recursive = TRUE
  )
}

out_pdf <- here(
  "output",
  paste0(
    FIGURE_NUMBER,
    "_WingLength_Habitat.pdf"
  )
)

out_tiff <- here(
  "output",
  paste0(
    FIGURE_NUMBER,
    "_WingLength_Habitat.tiff"
  )
)

ggsave(
  filename = out_pdf,
  plot = final_figure3,
  width = 174,
  height = 160,
  units = "mm",
  dpi = 300
)

ggsave(
  filename = out_tiff,
  plot = final_figure3,
  width = 174,
  height = 160,
  units = "mm",
  dpi = 300,
  compression = "lzw"
)

cat("\n========================================\n")
cat("FIGURE 3 COMPLETE\n")
cat("========================================\n")
cat("Species-level dataset:", nrow(tree_dat), "\n")
cat("Terrestrial species:", nrow(terrestrial_data), "\n")
cat("Wing threshold:", WING_THRESHOLD, "mm\n")
cat("Panel C species:", nrow(habitat_data), "\n")
cat("PDF:", out_pdf, "\n")
cat("TIFF:", out_tiff, "\n")







################################################################################
# FIGURE 3 RESULTS VALUE EXTRACTION SCRIPT
################################################################################

cat("\n========================================================================\n")
cat("                  FIGURE 3 MANUSCRIPT DIAGNOSTICS & VALUES              \n")
cat("========================================================================\n\n")

# ------------------------------------------------------------------------------
# PANEL A & B: WING-LENGTH THRESHOLD & REGRESSIONS (TERRESTRIAL SPECIES)
# ------------------------------------------------------------------------------
cat("--- PANEL A & B: TERRESTRIAL SPECIES WING-LENGTH THRESHOLD ---\n")
cat(sprintf("Total Terrestrial Species (N): %d\n", nrow(terrestrial_data)))
cat(sprintf("Small Wing Group (< 163 mm): N = %d\n", n_small))
cat(sprintf("Large Wing Group (>= 163 mm): N = %d\n", n_large))

cat("\nWilcoxon Test (Small vs Large Wing Group):\n")
cat(sprintf("Wilcoxon W = %.2f, p = %.4e, Effect size r = %.2f\n", 
            wilcox_res$statistic, wilcox_res$p, effsize_res$effsize))

cat("\nWing Group Medians and IQRs (MPs per individual):\n")
terrestrial_data %>%
  group_by(Wing_group) %>%
  summarise(
    N = n(),
    Median_MP = median(Mppersp, na.rm = TRUE),
    IQR_MP = IQR(Mppersp, na.rm = TRUE)
  ) %>%
  print()

cat("\nWing Length Model Comparisons:\n")
cat(sprintf("Linear Fit: R2 = %.3f, p = %.4e (AIC: %.1f)\n", 
            summary(m_lin)$r.squared, summary(m_lin)$coefficients[2, 4], AIC(m_lin)))
cat(sprintf("GAM Fit: EDF = %.2f, Dev Expl = %.1f%% (AIC: %.1f)\n", 
            edf, summary(m_gam)$dev.expl * 100, AIC(m_gam)))

# ------------------------------------------------------------------------------
# PANEL C: HABITAT-LEVEL MP BURDEN NORMALIZED BY BODY MASS
# ------------------------------------------------------------------------------
cat("\n--- PANEL C: HABITAT-LEVEL MP BURDEN NORMALIZED BY BODY MASS ---\n")
cat(sprintf("Total Species Analyzed in Panel C (N): %d across %d Habitats\n", 
            nrow(habitat_data), n_distinct(habitat_data$Habitat)))
cat(sprintf("Kruskal-Wallis Chi2 = %.2f, df = %d, p = %.4e\n", 
            kw_test$statistic, kw_test$parameter, kw_test$p.value))

cat("\nHabitat Summary Statistics (MPs per gram body mass):\n")
habitat_summary_stats <- habitat_data %>%
  group_by(Habitat) %>%
  summarise(
    N = n(),
    Median_MP_per_g = median(MP_per_g, na.rm = TRUE),
    IQR_MP_per_g = IQR(MP_per_g, na.rm = TRUE)
  ) %>%
  left_join(letters_df, by = "Habitat") %>%
  arrange(desc(Median_MP_per_g))

print(habitat_summary_stats)

cat("\nDunn Post-Hoc Test Pairwise Comparisons (BH Adjusted):\n")
print(dunn_res %>% dplyr::select(group1, group2, p.adj, p.adj.signif))

cat("\n========================================================================\n")
