
# ==============================================================================
# TISSUE SENSITIVITY ANALYSIS
# ==============================================================================

library(tidyverse)

# 1. Inspect tissue types present in the raw dataset
cat("=== RAW TISSUE FREQUENCIES ===\n")
print(table(bird.dat$Tissue))
cat("==============================\n\n")

# 2. Filter dataset to include only comprehensive digestive tracts
# (Excluding isolated droppings, feathers/FP, or non-digestive samples)
comprehensive_tissues <- c("Crop", "Esophagus", "Gizzard", "GIT", "RP")

digestive_strict <- bird.dat %>%
  filter(Tissue %in% comprehensive_tissues) %>%
  drop_na(Species, Mppersp, Mass, Wing.Length)

# 3. Aggregate at the species level for the strict subset
species_strict <- digestive_strict %>%
  group_by(Species) %>%
  summarise(
    logMP   = log1p(median(Mppersp, na.rm = TRUE)),
    logMass = log10(median(Mass, na.rm = TRUE)),
    logWing = log10(median(Wing.Length, na.rm = TRUE)),
    .groups = "drop"
  )

# 4. Ensure tree_dat has the corresponding logWing variable mapped for fair comparison
tree_dat_fixed <- tree_dat %>%
  mutate(logWing = log10(Wing.Length))

# 5. Fit models on both datasets
model_all_tissues    <- lm(logMP ~ logMass + logWing, data = tree_dat_fixed)
model_strict_tissues <- lm(logMP ~ logMass + logWing, data = species_strict)

# 6. Print comparative model summaries
cat("=== TISSUE SENSITIVITY MODEL COMPARISON ===\n")
cat("Full Species Aggregate Model (All Tissues):\n")
print(summary(model_all_tissues)$coefficients)

cat("\nStrict Digestive Tract Only Model (Excluding Droppings/Other):\n")
print(summary(model_strict_tissues)$coefficients)
cat("===============================================\n")