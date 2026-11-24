# ==============================================================================
# SCRIPT: FUNCTIONAL CLUSTERING, DECISION TREE MODELING & LOSO SENSITIVITY ANALYSIS
# Purpose: Reproducible workflow for avian microplastic trait predictors & validation
# ==============================================================================

# 1. Load Required Libraries
library(here)
library(tidyverse)
library(cluster)
library(rpart)
library(rpart.plot)

# Helper function for categorical mode
get_mode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

# Aggregation function from row-level records to species-level medians/modes
aggregate_species_data <- function(df_raw) {
  df_raw %>%
    group_by(Species) %>%
    summarise(
      Mppersp            = median(Mppersp, na.rm = TRUE),
      Beak.Length_Culmen = median(Beak.Length_Culmen, na.rm = TRUE),
      Wing.Length        = median(Wing.Length, na.rm = TRUE),
      Mass               = median(Mass, na.rm = TRUE),
      Habitat            = get_mode(Habitat),
      Migration          = get_mode(Migration),
      Trophic.Niche      = get_mode(Trophic.Niche),
      study_id           = get_mode(study_id),
      .groups            = "drop"
    ) %>%
    mutate(
      across(c(Habitat, Migration, Trophic.Niche), as.factor),
      logMP   = log1p(Mppersp),
      logMass = log10(Mass),
      logBeak = log10(Beak.Length_Culmen)
    )
}

# ==============================================================================
# 2. LOAD RAW DATA & BUILD BASELINE GLOBAL TREE
# ==============================================================================

tree_dat_raw <- read.csv(here("data/processed", "row_level.csv"), stringsAsFactors = FALSE) %>%
  drop_na(Species, Mppersp, Beak.Length_Culmen, Wing.Length, Mass, study_id, Habitat) %>%
  mutate(across(c(Habitat, Migration, Trophic.Niche, study_id), as.factor))

# Generate baseline species-level dataset
tree_dat <- aggregate_species_data(tree_dat_raw)

# Spearman Correlations & Gower Distance / PAM Clustering
cor_mass_sp <- cor.test(tree_dat$logMass, tree_dat$logMP, method = "spearman", exact = FALSE)
cor_beak_sp <- cor.test(tree_dat$logBeak, tree_dat$logMP, method = "spearman", exact = FALSE)

gower_dist <- daisy(
  tree_dat %>% dplyr::select(Beak.Length_Culmen, Habitat, Mass, Migration, Trophic.Niche, Wing.Length),
  metric = "gower"
)

pam_fit <- pam(gower_dist, k = 3)
tree_dat$cluster <- factor(pam_fit$clustering)
tree_dat$distinctiveness <- apply(as.matrix(gower_dist), 1, mean)

# Fit Baseline Global Decision Tree Model
tree_model <- rpart(
  logMP ~ distinctiveness + Beak.Length_Culmen + Wing.Length + Mass + Migration + Habitat + Trophic.Niche,
  data = tree_dat,
  method = "anova",
  control = rpart.control(minsplit = 10, minbucket = 5, cp = 0.005, maxdepth = 4, xval = 10)
)

# Robust Pruning
cp_table <- tree_model$cptable
if (nrow(cp_table) > 1) {
  min_idx <- which.min(cp_table[, "xerror"])
  best_cp <- cp_table[min_idx, "CP"]
  if (best_cp == cp_table[1, "CP"] && nrow(cp_table) > 1) {
    best_cp <- cp_table[2, "CP"]
  }
  final_tree <- prune(tree_model, cp = best_cp)
} else {
  final_tree <- tree_model
}


# ==============================================================================
# 3. LEAVE-ONE-STUDY-OUT (LOSO) SENSITIVITY ANALYSIS (TERRESTRIAL SUBSET)
# ==============================================================================

terrestrial_raw_base <- tree_dat_raw %>%
  filter(Habitat %in% c("Forest", "Grassland", "Human modified", "Shrubland", "Woodland"))

unique_studies <- unique(as.character(terrestrial_raw_base$study_id))
cat("Total independent studies in terrestrial dataset:", length(unique_studies), "\n")

tree_formula <- logMP ~ Beak.Length_Culmen + Wing.Length + Mass + Migration + Trophic.Niche

loso_terrestrial_results <- data.frame(
  excluded_study = unique_studies,
  root_variable  = NA_character_,
  threshold      = NA_real_,
  mean_diff_163  = NA_real_
)

for (i in seq_along(unique_studies)) {
  left_out <- unique_studies[i]
  
  # Exclude all rows from the left-out study
  raw_train <- terrestrial_raw_base %>% filter(study_id != left_out)
  
  # Aggregate to species level
  train_dat <- aggregate_species_data(raw_train)
  
  if(nrow(train_dat) < 15) next # Skip if too few species remain
  
  # Fit and prune tree on terrestrial subset
  fit <- tryCatch({
    t_raw <- rpart(
      formula = tree_formula, 
      data = train_dat, 
      method = "anova",
      control = rpart.control(minsplit = 10, minbucket = 5, cp = 0.005, maxdepth = 4, xval = 10)
    )
    cp_tab <- t_raw$cptable
    m_idx <- which.min(cp_tab[, "xerror"])
    b_cp <- cp_tab[m_idx, "CP"]
    if (length(b_cp) == 0 || is.na(b_cp)) b_cp <- cp_tab[1, "CP"]
    prune(t_raw, cp = b_cp)
  }, error = function(e) NULL)
  
  if (is.null(fit)) next
  
  root_var <- as.character(fit$frame$var[1])
  loso_terrestrial_results$root_variable[i] <- root_var
  
  # Record threshold if Wing.Length is selected at the root of the terrestrial tree
  if (root_var == "Wing.Length" && !is.null(fit$splits)) {
    loso_terrestrial_results$threshold[i] <- fit$splits[1, "index"]
  }
  
  # Calculate effect size at the fixed 163 mm threshold for this leave-one-out fold
  below <- train_dat$Wing.Length <= 163
  above <- train_dat$Wing.Length > 163
  
  if (sum(below) > 0 && sum(above) > 0) {
    loso_terrestrial_results$mean_diff_163[i] <- mean(train_dat$logMP[below]) - mean(train_dat$logMP[above])
  }
}

# ==============================================================================
# 4. SUMMARIZE & EXPORT RESULTS
# ==============================================================================

cat("\n===== TERRESTRIAL LEAVE-ONE-STUDY-OUT RESULTS =====\n")
print(table(loso_terrestrial_results$root_variable, useNA = "ifany"))

valid_thresh <- loso_terrestrial_results$threshold[!is.na(loso_terrestrial_results$threshold)]
if (length(valid_thresh) > 0) {
  cat("\nThreshold (mm) range across study omissions:", range(valid_thresh), "\n")
  cat("Median threshold:", median(valid_thresh), "\n")
}

valid_diffs <- loso_terrestrial_results$mean_diff_163[!is.na(loso_terrestrial_results$mean_diff_163)]
if (length(valid_diffs) > 0) {
  cat("\nEffect size difference range (<= 163 mm minus > 163 mm):", range(valid_diffs), "\n")
  cat("Mean effect size across omissions:", mean(valid_diffs), "\n")
}

# Export sensitivity output for supplementary archiving
write.csv(loso_terrestrial_results, here("output", "LOSO_Sensitivity_Analysis_Results.csv"), row.names = FALSE)
cat("\nWorkflow completed successfully. LOSO sensitivity table saved to output directory.\n")
# ==============================================================================
# GENERATE SUPPLEMENTARY TABLE: LOSO SENSITIVITY SUMMARY
# ==============================================================================

library(tidyverse)
library(here)

# Format the LOSO results dataframe cleanly for publication
loso_table_pub <- loso_terrestrial_results %>%
  mutate(
    Threshold_mm = ifelse(is.na(threshold), "Node Shift", round(threshold, 2)),
    Mean_LogMP_Diff = round(mean_diff_163, 3)
  ) %>%
  dplyr::select(
    Excluded_Study = excluded_study,
    Primary_Root_Split = root_variable,
    Threshold_Value_mm = Threshold_mm,
    Effect_Size_Diff = Mean_LogMP_Diff
  )

# View table in console
print(loso_table_pub)

# Save for Supplementary Material insertion
write.csv(loso_table_pub, here("output", "Table_S_LOSO_Sensitivity.csv"), row.names = FALSE)
cat("\nSaved Supplementary Table to output/Table_S_LOSO_Sensitivity.csv\n")




# ==============================================================================
# GENERATE SUPPLEMENTARY FIGURE: LOSO & BOOTSTRAP STABILITY PLOT
# ==============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)
library(here)

# Panel A: LOSO Threshold Distribution
fig_s_a <- ggplot(loso_terrestrial_results %>% filter(!is.na(threshold)), aes(x = threshold)) +
  geom_histogram(binwidth = 5, fill = "#588157", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 163.3, color = "#D62828", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Wing-length threshold (mm)",
    y = "Frequency (Study Omissions)",
    title = "A. Leave-One-Study-Out Stability",
    subtitle = "Median threshold = 163.3 mm (100% Wing.Length root split)"
  ) +
  theme_classic(base_size = 11) +
  theme(axis.title = element_text(face = "bold"), plot.title = element_text(face = "bold"))

# Panel B: LOSO Effect Size Distribution
fig_s_b <- ggplot(loso_terrestrial_results %>% filter(!is.na(mean_diff_163)), aes(x = mean_diff_163)) +
  geom_histogram(binwidth = 0.2, fill = "#4d7a9c", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 2.13, color = "#D62828", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Mean difference in logMP (<=163 mm minus >163 mm)",
    y = "Frequency (Study Omissions)",
    title = "B. Effect Size Robustness",
    subtitle = "Mean effect across omissions = +2.13"
  ) +
  theme_classic(base_size = 11) +
  theme(axis.title = element_text(face = "bold"), plot.title = element_text(face = "bold"))

# Combine panels
supp_fig <- (fig_s_a / fig_s_b) + plot_annotation(tag_levels = 'A')

# Save high-resolution PDF for submission
ggsave(here("output", "Figure_S_Sensitivity_Analysis.pdf"), supp_fig, width = 8, height = 9, dpi = 300, device = cairo_pdf)
cat("\nSaved Supplementary Figure to output/Figure_S_Sensitivity_Analysis.pdf\n")