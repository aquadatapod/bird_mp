################################################################################
# PGLS ANALYSIS AND BRANCH-LENGTH SENSITIVITY
################################################################################

library(ape)
library(caper)

# 1. Verify that the pPCA scores contain the required variables

required_vars <- c(
  "species_std", "logMP", "PC1", "PC2",
  "Habitat", "Feeding_habit", "Migration"
)

missing_vars <- setdiff(required_vars, names(scores_df))

if (length(missing_vars) > 0) {
  stop(
    "The following variables are missing from scores_df: ",
    paste(missing_vars, collapse = ", "),
    ". Run the pPCA section first and check scores_df."
  )
}

# 2. Prepare one common dataset for all PGLS models

pgls_df <- as.data.frame(scores_df)

pgls_df$species_std_caper <- as.character(pgls_df$species_std)

pgls_df <- pgls_df[
  complete.cases(
    pgls_df[, c(
      "species_std_caper", "logMP", "PC1", "PC2",
      "Habitat", "Feeding_habit", "Migration"
    )]
  ),
  ,
  drop = FALSE
]

pgls_df$Habitat <- droplevels(as.factor(pgls_df$Habitat))
pgls_df$Feeding_habit <- droplevels(as.factor(pgls_df$Feeding_habit))
pgls_df$Migration <- droplevels(as.factor(pgls_df$Migration))

if (anyDuplicated(pgls_df$species_std_caper)) {
  stop("Duplicate species identifiers occur in pgls_df.")
}

if (nrow(pgls_df) < 5) {
  stop("Too few complete species remain for PGLS.")
}

if (nlevels(pgls_df$Habitat) < 2 ||
    nlevels(pgls_df$Feeding_habit) < 2 ||
    nlevels(pgls_df$Migration) < 2) {
  stop("At least one categorical predictor has fewer than two levels in the complete-case dataset.")
}

# 3. Match the phylogeny to the PGLS dataset

pgls_tree <- ape::keep.tip(
  tree,
  pgls_df$species_std_caper
)

pgls_df <- pgls_df[
  match(pgls_tree$tip.label, pgls_df$species_std_caper),
  ,
  drop = FALSE
]

rownames(pgls_df) <- pgls_df$species_std_caper

stopifnot(
  identical(pgls_tree$tip.label, pgls_df$species_std_caper)
)

cat("\nPGLS species included:", nrow(pgls_df), "\n")

# 4. Fit the three original PGLS models using the same species

comp_data <- caper::comparative.data(
  phy = pgls_tree,
  data = pgls_df,
  names.col = "species_std_caper",
  vcv = TRUE,
  warn.dropped = TRUE
)

pgls_morph <- caper::pgls(
  logMP ~ PC1 + PC2,
  data = comp_data,
  lambda = "ML"
)

pgls_eco <- caper::pgls(
  logMP ~ Habitat + Feeding_habit + Migration,
  data = comp_data,
  lambda = "ML"
)

pgls_glob <- caper::pgls(
  logMP ~ PC1 + PC2 + Habitat + Feeding_habit + Migration,
  data = comp_data,
  lambda = "ML"
)

get_lambda <- function(mod) {
  as.numeric(mod$param["lambda"])
}

pgls_aic_table <- data.frame(
  Model = c("Morphology", "Ecology", "Global"),
  AIC = c(
    AIC(pgls_morph),
    AIC(pgls_eco),
    AIC(pgls_glob)
  ),
  R_Squared = c(
    summary(pgls_morph)$r.squared,
    summary(pgls_eco)$r.squared,
    summary(pgls_glob)$r.squared
  ),
  Lambda = c(
    get_lambda(pgls_morph),
    get_lambda(pgls_eco),
    get_lambda(pgls_glob)
  )
)

cat("\n=== PGLS MODEL COMPARISON ===\n")
print(pgls_aic_table)

# 5. Define alternative branch-length assumptions

tree_original <- pgls_tree

tree_grafen_05 <- ape::compute.brlen(
  pgls_tree,
  method = "Grafen",
  power = 0.5
)

tree_grafen_2 <- ape::compute.brlen(
  pgls_tree,
  method = "Grafen",
  power = 2
)

tree_list <- list(
  Original = tree_original,
  Grafen_05 = tree_grafen_05,
  Grafen_2 = tree_grafen_2
)

# 6. Refit all three models under each branch-length assumption

model_list <- list(
  Morphology = logMP ~ PC1 + PC2,
  Ecology = logMP ~ Habitat + Feeding_habit + Migration,
  Global = logMP ~ PC1 + PC2 + Habitat + Feeding_habit + Migration
)

sensitivity_results <- list()
counter <- 1

for (tree_name in names(tree_list)) {
  
  comp_sens <- caper::comparative.data(
    phy = tree_list[[tree_name]],
    data = pgls_df,
    names.col = "species_std_caper",
    vcv = TRUE,
    warn.dropped = TRUE
  )
  
  for (model_name in names(model_list)) {
    
    fit <- caper::pgls(
      model_list[[model_name]],
      data = comp_sens,
      lambda = "ML"
    )
    
    coef_table <- as.data.frame(
      summary(fit)$coefficients
    )
    
    coef_table$Term <- rownames(coef_table)
    rownames(coef_table) <- NULL
    
    coef_table$Model <- model_name
    coef_table$Branch_lengths <- tree_name
    coef_table$AIC <- AIC(fit)
    coef_table$Lambda <- get_lambda(fit)
    
    sensitivity_results[[counter]] <- coef_table
    counter <- counter + 1
  }
}

# 7. Combine sensitivity results

sensitivity_table <- do.call(
  rbind,
  sensitivity_results
)

sensitivity_table <- sensitivity_table[
  ,
  c(
    "Model", "Branch_lengths", "Term",
    "Estimate", "Std. Error", "t value",
    "Pr(>|t|)", "AIC", "Lambda"
  )
]

cat("\n=== BRANCH-LENGTH SENSITIVITY RESULTS ===\n")
print(sensitivity_table)

# 8. Export results

if (!dir.exists(here::here("output"))) {
  dir.create(here::here("output"), recursive = TRUE)
}

write.csv(
  pgls_aic_table,
  here::here("output", "PGLS_Model_Comparison.csv"),
  row.names = FALSE
)

write.csv(
  sensitivity_table,
  here::here("output", "PGLS_Branch_Length_Sensitivity.csv"),
  row.names = FALSE
)

write.csv(
  summary(pgls_glob)$coefficients,
  here::here("output", "PGLS_Global_Coefficients.csv")
)

cat("\nPGLS and branch-length sensitivity analyses completed.\n")










################################################################################
# PGLS MODEL DIAGNOSTICS: LAMBDA PROFILES, CONFIDENCE INTERVALS AND RESIDUALS
################################################################################

# 1. Profile likelihood for Pagel's lambda
lambda_profile_morph <- caper::pgls.profile(
  pgls_morph, which = "lambda", N = 100
)

lambda_profile_eco <- caper::pgls.profile(
  pgls_eco, which = "lambda", N = 100
)

lambda_profile_global <- caper::pgls.profile(
  pgls_glob, which = "lambda", N = 100
)

# 2. Confidence intervals for lambda
lambda_ci_morph <- caper::pgls.confint(
  pgls_morph, which = "lambda", param.CI = 0.95
)

lambda_ci_eco <- caper::pgls.confint(
  pgls_eco, which = "lambda", param.CI = 0.95
)

lambda_ci_global <- caper::pgls.confint(
  pgls_glob, which = "lambda", param.CI = 0.95
)

# 3. Extract profile-likelihood results
extract_lambda_profile <- function(profile, model_name) {
  
  data.frame(
    Model = model_name,
    Lambda_ML = profile$opt,
    Lower_Bound = profile$bounds.val[1],
    Lower_Bound_p = profile$bounds.p[1],
    Upper_Bound = profile$bounds.val[2],
    Upper_Bound_p = profile$bounds.p[2],
    CI_Lower = profile$ci.val[1],
    CI_Upper = profile$ci.val[2]
  )
}

lambda_diagnostics <- bind_rows(
  extract_lambda_profile(lambda_ci_morph, "Morphology"),
  extract_lambda_profile(lambda_ci_eco, "Ecology"),
  extract_lambda_profile(lambda_ci_global, "Global")
)

cat("\n=== PAGEL'S LAMBDA PROFILE-LIKELIHOOD RESULTS ===\n")
print(lambda_diagnostics, digits = 5)


# Plot likelihood profiles
pdf(
  here::here("output", "PGLS_Lambda_Profiles.pdf"),
  width = 12,
  height = 5
)

par(mfrow = c(1, 3))

plot(
  lambda_profile_morph,
  main = "Morphology model"
)

plot(
  lambda_profile_eco,
  main = "Ecology model"
)

plot(
  lambda_profile_global,
  main = "Global model"
)

dev.off()


# 5. Residual diagnostic plots
pdf(
  here::here("output", "PGLS_Residual_Diagnostics.pdf"),
  width = 10,
  height = 8
)

par(mfrow = c(3, 2))

models <- list(
  Morphology = pgls_morph,
  Ecology = pgls_eco,
  Global = pgls_glob
)

for (model_name in names(models)) {
  
  model <- models[[model_name]]
  
  plot(
    model$fitted,
    model$res,
    xlab = "Fitted values",
    ylab = "Residuals",
    main = paste(model_name, "- Residuals vs fitted")
  )
  
  abline(h = 0, lty = 2)
  
  qqnorm(
    model$res,
    main = paste(model_name, "- Normal Q-Q")
  )
  
  qqline(model$res)
}

dev.off()

# 6. Export results
write.csv(
  lambda_diagnostics,
  here::here("output", "PGLS_Lambda_Profile_Results.csv"),
  row.names = FALSE
)

cat("\nPGLS diagnostic files saved in the output directory.\n")



str(lambda_profile_morph)
str(lambda_profile_eco)
str(lambda_profile_global)