# ==============================================================================
# VIF & COLLINEARITY CHECK FOR BODY MASS AND WING LENGTH
# ==============================================================================

library(car)

# 1. Fit standard linear model with both mass and wing length
vif_model <- lm(logMP ~ log10_Mass + log10_Wing, data = trait_final)

# 2. Calculate Variance Inflation Factors
vif_values <- vif(vif_model)
cat("\n=== VARIANCE INFLATION FACTORS (VIF) ===\n")
print(vif_values)

# 3. Calculate Spearman correlation between mass and wing length
mass_wing_cor <- cor.test(trait_final$log10_Mass, trait_final$log10_Wing, method = "spearman", exact = FALSE)
cat("\n=== SPEARMAN CORRELATION (Mass vs Wing Length) ===\n")
print(mass_wing_cor)








# ==============================================================================
# SENSITIVITY CHECK: FUNCTIONAL DISTINCTIVENESS VS UNDERLYING TRAITS 
# ==============================================================================

library(tidyverse)
library(car)

# 1. Model comparison
model_dist_only <- lm(logMP ~ distinctiveness, data = tree_dat)
model_controlled <- lm(logMP ~ distinctiveness + logMass + Wing.Length + Habitat, data = tree_dat)

cat("\n=== DISTINCTIVENESS MODEL COMPARISON ===\n")
cat("Distinctiveness alone:\n")
print(summary(model_dist_only)$coefficients)

cat("\nControlled for Mass, Wing Length, and Habitat:\n")
print(summary(model_controlled)$coefficients)
cat("===========================================\n")

# 2. Parus major & Influence check
parus_check <- tree_dat %>% 
  mutate(row_num = row_number()) %>%
  filter(str_detect(Species, "Parus major"))

cat("\n=== PARUS MAJOR INFLUENCE CHECK ===\n")
if(nrow(parus_check) > 0) {
  # Print cleanly using base selection to avoid pipe printing bugs
  print(as.data.frame(parus_check[, c("Species", "logMP", "distinctiveness", "logMass", "Wing.Length")]))
  
  # Calculate Cook's distance
  cooks_d <- cooks.distance(model_controlled)
  influential_points <- which(cooks_d > (4 / nrow(tree_dat)))
  cat("\nNumber of influential points (Cook's D > 4/n):", length(influential_points), "\n")
} else {
  cat("Parus major not found directly in tree_dat species names.\n")
}
cat("===================================\n")



