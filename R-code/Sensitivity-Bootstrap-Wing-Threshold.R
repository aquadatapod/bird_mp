# ==============================================================================
# SENSITIVITY & BOOTSTRAP STABILITY FOR WING-LENGTH THRESHOLD 
# ==============================================================================

set.seed(42)
n_boot <- 1000
threshold_splits <- numeric(n_boot)

for(i in 1:n_boot) {
  boot_dat <- tree_dat %>% sample_frac(replace = TRUE)
  boot_tree <- rpart(
    logMP ~ distinctiveness + Beak.Length_Culmen + Wing.Length + Mass + Migration + Habitat + Trophic.Niche,
    data = boot_dat,
    method = "anova",
    control = rpart.control(minsplit = 10, minbucket = 5, cp = 0.005, maxdepth = 4)
  )
  
  # Extract primary split variable and split point if Wing.Length is root split
  splits <- boot_tree$splits
  if(!is.null(splits) && "Wing.Length" %in% rownames(splits)) {
    threshold_splits[i] <- splits["Wing.Length", "index"]
  } else {
    threshold_splits[i] <- NA
  }
}

valid_thresholds <- na.omit(threshold_splits)
cat("\n=== BOOTSTRAP THRESHOLD STABILITY (1000 iterations) ===\n")
cat("Selection frequency of Wing Length as primary split:", length(valid_thresholds)/n_boot * 100, "%\n")
cat("Mean threshold split value:", mean(valid_thresholds), "mm\n")
cat("95% Confidence Interval:", quantile(valid_thresholds, c(0.025, 0.975)), "mm\n")


################################################################################
# Revised bootstrap

library(dplyr)
library(rpart)

terrestrial_dat <- tree_dat %>%
  filter(
    Habitat %in% c(
      "Forest",
      "Grassland",
      "Human modified",
      "Shrubland",
      "Woodland"
    )
  ) %>%
  droplevels() %>%
  filter(
    !is.na(logMP),
    !is.na(Wing.Length)
  )

cat("Terrestrial sample size:", nrow(terrestrial_dat), "\n")
print(table(terrestrial_dat$Habitat))

cat("\nWing length groups around 163 mm:\n")
print(table(
  ifelse(
    terrestrial_dat$Wing.Length <= 163,
    "At or below 163 mm",
    "Above 163 mm"
  )
))


########
########
set.seed(42)

n_boot <- 1000

predictors <- c(
  "distinctiveness",
  "Beak.Length_Culmen",
  "Wing.Length",
  "Mass",
  "Migration",
  "Trophic.Niche"
)

tree_formula <- logMP ~
  distinctiveness +
  Beak.Length_Culmen +
  Wing.Length +
  Mass +
  Migration +
  Trophic.Niche

boot_results <- data.frame(
  iteration = seq_len(n_boot),
  root_variable = NA_character_,
  threshold = NA_real_,
  n_below = NA_integer_,
  n_above = NA_integer_,
  mean_below = NA_real_,
  mean_above = NA_real_,
  difference_logMP = NA_real_
)

for (i in seq_len(n_boot)) {
  
  boot_dat <- terrestrial_dat[
    sample(
      seq_len(nrow(terrestrial_dat)),
      size = nrow(terrestrial_dat),
      replace = TRUE
    ),
  ]
  
  boot_tree <- tryCatch(
    rpart(
      formula = tree_formula,
      data = boot_dat,
      method = "anova",
      control = rpart.control(
        minsplit = 10,
        minbucket = 5,
        cp = 0.005,
        maxdepth = 4,
        xval = 10
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(boot_tree)) next
  
  cp_table <- boot_tree$cptable
  
  if (nrow(cp_table) > 1) {
    
    min_idx <- which.min(cp_table[, "xerror"])
    best_cp <- cp_table[min_idx, "CP"]
    
    if (min_idx == 1) {
      best_cp <- cp_table[2, "CP"]
    }
    
    final_boot_tree <- prune(
      boot_tree,
      cp = best_cp
    )
    
  } else {
    
    final_boot_tree <- boot_tree
  }
  
  root_var <- as.character(final_boot_tree$frame$var[1])
  
  boot_results$root_variable[i] <- root_var
  
  if (root_var == "Wing.Length") {
    
    threshold <- final_boot_tree$splits[1, "index"]
    
    boot_results$threshold[i] <- threshold
    
    below <- boot_dat$Wing.Length <= threshold
    above <- boot_dat$Wing.Length > threshold
    
    if (sum(below) > 0 && sum(above) > 0) {
      
      boot_results$n_below[i] <- sum(below)
      boot_results$n_above[i] <- sum(above)
      
      boot_results$mean_below[i] <-
        mean(boot_dat$logMP[below])
      
      boot_results$mean_above[i] <-
        mean(boot_dat$logMP[above])
      
      boot_results$difference_logMP[i] <-
        boot_results$mean_below[i] -
        boot_results$mean_above[i]
    }
  }
  
  # Evaluate the original 163 mm threshold independently
  below_original <- boot_dat$Wing.Length <= 163
  above_original <- boot_dat$Wing.Length > 163
  
  if (sum(below_original) > 0 &&
      sum(above_original) > 0) {
    
    boot_results$mean_below[i] <-
      mean(boot_dat$logMP[below_original])
    
    boot_results$mean_above[i] <-
      mean(boot_dat$logMP[above_original])
    
    # Store the fixed-threshold effect separately below
  }
}

# Calculate the effect at  the original 163 mm threshold

set.seed(42)

fixed_threshold_effects <- rep(NA_real_, n_boot)

for (i in seq_len(n_boot)) {
  
  boot_dat <- terrestrial_dat[
    sample(
      seq_len(nrow(terrestrial_dat)),
      size = nrow(terrestrial_dat),
      replace = TRUE
    ),
  ]
  
  below <- boot_dat$Wing.Length <= 163
  above <- boot_dat$Wing.Length > 163
  
  if (sum(below) > 0 && sum(above) > 0) {
    
    fixed_threshold_effects[i] <-
      mean(boot_dat$logMP[below]) -
      mean(boot_dat$logMP[above])
  }
}

valid_effects <- fixed_threshold_effects[
  is.finite(fixed_threshold_effects)
]

# Summaries of the boostrap

valid_thresholds <- boot_results$threshold[
  is.finite(boot_results$threshold)
]

root_frequency <- mean(
  boot_results$root_variable == "Wing.Length",
  na.rm = TRUE
) * 100

cat("\n===== TERRESTRIAL BOOTSTRAP RESULTS =====\n")

cat("\nWing length selected as root split:",
    round(root_frequency, 2), "%\n")

cat("\nNumber of trees selecting wing length:",
    length(valid_thresholds), "of", n_boot, "\n")

if (length(valid_thresholds) > 0) {
  
  cat("\nBootstrap threshold distribution (mm):\n")
  
  print(
    quantile(
      valid_thresholds,
      probs = c(0.025, 0.50, 0.975)
    )
  )
}

cat("\nValid fixed-threshold effect estimates:",
    length(valid_effects), "of", n_boot, "\n")

if (length(valid_effects) > 0) {
  
  cat("\nMean difference in logMP (at or below 163 mm minus above):\n")
  print(mean(valid_effects))
  
  cat("\n95% bootstrap interval:\n")
  print(
    quantile(
      valid_effects,
      probs = c(0.025, 0.50, 0.975)
    )
  )
  
  cat("\nIterations with higher mean logMP at or below 163 mm:\n")
  print(mean(valid_effects > 0) * 100)
}
