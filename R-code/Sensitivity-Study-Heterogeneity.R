


# ==============================================================================
# SENSITIVITY ANALYSIS: EXCLUDING SINGLE-STUDY SPECIES
# ==============================================================================

# 1. Identify number of unique studies per species from the raw data
species_study_counts <- raw_data %>%
  drop_na(Species, Mppersp, study_id) %>%
  group_by(Species) %>%
  summarise(n_studies = n_distinct(study_id), .groups = "drop")

# 2. Merge study counts back into your trait_final or tree_dat dataset
# (Ensure your dataset has a 'Species' column)
sensitivity_dat <- trait_final %>%
  left_join(species_study_counts, by = "Species") %>%
  filter(n_studies > 1)

cat(sprintf("Total species before filtering: %d\n", nrow(trait_final)))
cat(sprintf("Species represented by > 1 study: %d\n", nrow(sensitivity_dat)))
cat(sprintf("Species represented by exactly 1 study: %d\n", nrow(trait_final) - nrow(sensitivity_dat)))

# 3. Re-run core regression model on multi-study species only
model_full <- lm(logMP ~ log10_Mass + log10_Wing, data = trait_final)
model_sens <- lm(logMP ~ log10_Mass + log10_Wing, data = sensitivity_dat)

cat("\n=== SENSITIVITY MODEL COMPARISON (Coefficients) ===\n")
cat("Full Dataset Model:\n")
print(summary(model_full)$coefficients)
cat("\nMulti-Study Only Model:\n")
print(summary(model_sens)$coefficients)
cat("===================================================\n")




# ==============================================================================
# SENSITIVITY ANALYSIS: STUDY REPRESENTATION & MIXED MODELS
# ==============================================================================

library(lme4)
library(lmerTest)

# 1. Count studies per species using raw row-level data
species_study_counts <- raw_data %>%
  drop_na(Species, Mppersp, study_id) %>%
  group_by(Species) %>%
  summarise(n_studies = n_distinct(study_id), .groups = "drop")

# 2. Merge study counts into row-level data to test via mixed models
row_level_sensitivity <- raw_data %>%
  left_join(species_study_counts, by = "Species") %>%
  mutate(
    logMP = log1p(Mppersp),
    logMass = log10(Mass),
    logWing = log10(Wing.Length),
    study_id = factor(study_id)
  ) %>%
  filter(!is.na(logMP), !is.na(logMass), !is.na(logWing))

# 3. Full mixed-effects model (all species, controlling for study random intercept)
model_mixed_full <- lmer(logMP ~ logMass + logWing + (1 | study_id), data = row_level_sensitivity)

# 4. Filtered mixed-effects model (excluding species represented by only 1 study)
row_level_multistudy <- row_level_sensitivity %>% filter(n_studies > 1)
model_mixed_multi <- lmer(logMP ~ logMass + logWing + (1 | study_id), data = row_level_multistudy)

cat("\n=== MIXED-EFFECTS SENSITIVITY COMPARISON (Accounting for ICC) ===\n")
cat("Full Dataset Mixed Model:\n")
print(summary(model_mixed_full)$coefficients)
cat("\nMulti-Study Only Mixed Model:\n")
print(summary(model_mixed_multi)$coefficients)
cat("================================================================\n")









# ==============================================================================
# SENSITIVITY ANALYSIS: STUDY HETEROGENEITY & MULTI-STUDY SPECIES
# ==============================================================================

library(tidyverse)
library(lme4)
library(lmerTest)

# 1. Count independent studies per species from raw data
species_study_counts <- raw_data %>%
  drop_na(Species, Mppersp, study_id) %>%
  group_by(Species) %>%
  summarise(n_studies = n_distinct(study_id), .groups = "drop")

cat("Total species:", nrow(species_study_counts), "\n")
cat("Species with > 1 study:", sum(species_study_counts$n_studies > 1), "\n")

# ------------------------------------------------------------------------------
# PART A: Species-Level Median Sensitivity (Paralleling Main Trait Analyses)
# ------------------------------------------------------------------------------
species_sens <- raw_data %>%
  drop_na(Species, Mppersp, Mass, Beak.Length_Culmen, Wing.Length, Habitat) %>%
  group_by(Species) %>%
  summarise(
    n_studies    = n_distinct(study_id),
    logMP_median = median(log1p(Mppersp), na.rm = TRUE),
    logMass      = median(log10(Mass), na.rm = TRUE),
    logBeak      = median(log10(Beak.Length_Culmen), na.rm = TRUE),
    logWing      = median(log10(Wing.Length), na.rm = TRUE),
    Habitat      = names(sort(table(Habitat), decreasing = TRUE))[1],
    .groups      = "drop"
  ) %>%
  filter(n_studies > 1)

cat("\n=== SPECIES-LEVEL MEDIAN SENSITIVITY (n_studies > 1) ===\n")
print(cor.test(species_sens$logMass, species_sens$logMP_median, method = "spearman", exact = FALSE))
print(cor.test(species_sens$logBeak, species_sens$logMP_median, method = "spearman", exact = FALSE))
print(summary(lm(logMP_median ~ Habitat, data = species_sens)))

# ------------------------------------------------------------------------------
# PART B: Row-Level Mixed Models (Separating Study vs. Species Variance)
# ------------------------------------------------------------------------------
row_level_sens <- raw_data %>%
  left_join(species_study_counts, by = "Species") %>%
  mutate(
    logMP    = log1p(Mppersp),
    logMass  = log10(Mass),
    logBeak  = log10(Beak.Length_Culmen),
    study_id = factor(study_id),
    Species  = factor(Species)
  ) %>%
  filter(!is.na(logMP), !is.na(logMass), !is.na(Habitat))

# Separate univariate models to avoid Mass-Wing collinearity, with crossed random effects
lmm_mass_full  <- lmer(logMP ~ logMass + (1 | study_id) + (1 | Species), data = row_level_sens)
lmm_mass_multi <- lmer(logMP ~ logMass + (1 | study_id) + (1 | Species), data = filter(row_level_sens, n_studies > 1))

lmm_hab_full   <- lmer(logMP ~ Habitat + (1 | study_id) + (1 | Species), data = row_level_sens)
lmm_hab_multi  <- lmer(logMP ~ Habitat + (1 | study_id) + (1 | Species), data = filter(row_level_sens, n_studies > 1))

cat("\n=== ROW-LEVEL LMM: BODY MASS (All vs. Multi-Study Species) ===\n")
print(summary(lmm_mass_full)$coefficients)
print(summary(lmm_mass_multi)$coefficients)

cat("\n=== ROW-LEVEL LMM: HABITAT (All vs. Multi-Study Species) ===\n")
print(anova(lmm_hab_full))
print(anova(lmm_hab_multi))