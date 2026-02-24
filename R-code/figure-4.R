## Functional traits and clustering
library(cluster)
library(ape)
library(phytools)

## Tree-based modelling
library(rpart)
library(rpart.plot)
library(tidyverse)
library(dplyr)
library(ggplot2)

# -----------------------------
# Load and prepare data
# -----------------------------

tree_dat_raw <- read.csv("data_tree.csv") %>%
  drop_na() %>%
  mutate(
    across(
      c(Habitat, Migration, Trophic.Niche, Species),
      as.factor
    )
  )

# -----------------------------
# Collapse to species-level means
# -----------------------------

tree_dat <- tree_dat_raw %>%
  group_by(Species) %>%
  summarise(
    Mppersp = mean(Mppersp, na.rm = TRUE),
    Beak.Length_Culmen = mean(Beak.Length_Culmen, na.rm = TRUE),
    Wing.Length = mean(Wing.Length, na.rm = TRUE),
    Mass = mean(Mass, na.rm = TRUE),
    Habitat = first(Habitat),
    Migration = first(Migration),
    Trophic.Niche = first(Trophic.Niche)
  ) %>%
  ungroup()

# Log-transform MP response
tree_dat <- tree_dat %>%
  mutate(logMP = log1p(Mppersp))

# -----------------------------
# Functional clustering
# -----------------------------

gower_dist <- daisy(
  tree_dat %>%
    select(
      Beak.Length_Culmen,
      Habitat,
      Mass,
      Migration,
      Trophic.Niche,
      Wing.Length
    ),
  metric = "gower"
)

pam_fit <- pam(gower_dist, k = 3)
tree_dat$cluster <- factor(pam_fit$clustering)

# Functional distinctiveness (species-level)
tree_dat$distinctiveness <- apply(
  as.matrix(gower_dist),
  1,
  mean
)

# -----------------------------
# Decision tree model
# -----------------------------

tree_model_01 <- rpart(
  logMP ~ distinctiveness +
    Beak.Length_Culmen +
    Wing.Length +
    Mass +
    Migration,
  data = tree_dat,
  method = "anova",
  control = rpart.control(minsplit = 5, cp = 0.01)
)

# -----------------------------
# Visualization
# -----------------------------

rpart.plot(
  tree_model_01,
  box.palette = "-RdBu",
  branch.col = "grey40",
  main = "Species-level trait predictors of microplastic accumulation"
)

printcp(tree_model_01)




################################################################################

# ---------------------------------------
# 1. Basic data check
# ---------------------------------------

cat("Number of species:", nrow(tree_dat), "\n")

# Quick extreme value check (informational only)
cat("Beak >100 mm:", sum(tree_dat$Beak.Length_Culmen > 100), "\n")
cat("Mass >2000 g:", sum(tree_dat$Mass > 2000), "\n")

# Ensure Migration is factor
tree_dat$Migration <- as.factor(tree_dat$Migration)

# ---------------------------------------
# 2. Build cross-validated decision tree
# ---------------------------------------

tree_model <- rpart(
  logMP ~ distinctiveness +
    Beak.Length_Culmen +
    Wing.Length +
    Mass +
    Migration,
  data = tree_dat,
  method = "anova",
  control = rpart.control(
    minsplit = 30,
    minbucket = 15,
    cp = 0.01,
    maxdepth = 4,
    xval = 10
  )
)



rpart.plot(
  tree_model,
  type = 4,
  extra = 101,
  box.palette = "-RdBu",
  shadow.col = "gray",
  nn = FALSE,
  fallen.leaves = TRUE,
  main = "Trait Predictors of Microplastic Accumulation",
  sub = paste("Based on", nrow(tree_dat), "species")
)

# ---------------------------------------
# 3. Select optimal tree (cross-validation)
# ---------------------------------------

cp_table <- tree_model$cptable
best_cp <- cp_table[which.min(cp_table[, "xerror"]), "CP"]

final_tree <- prune(tree_model, cp = best_cp)

cat("Optimal CP selected:", best_cp, "\n")

# ---------------------------------------
# 4. Plot final tree
# ---------------------------------------

pdf("Decision_Tree_Final.pdf", width = 10, height = 8)

rpart.plot(
  final_tree,
  type = 4,
  extra = 101,
  box.palette = "-RdBu",
  shadow.col = "gray",
  nn = FALSE,
  fallen.leaves = TRUE,
  main = "Trait Predictors of Microplastic Accumulation",
  sub = paste("Based on", nrow(tree_dat), "species")
)

dev.off()

cat("✓ Final tree saved: Decision_Tree_Final.pdf\n")

# ---------------------------------------
# 5. Extract terminal node summaries
# ---------------------------------------

frame <- final_tree$frame
terminal_nodes <- frame[frame$var == "<leaf>", ]

terminal_nodes <- terminal_nodes[order(terminal_nodes$yval, decreasing = TRUE), ]

cat("\nTerminal Profiles:\n")

for(i in 1:nrow(terminal_nodes)) {
  
  node_id <- as.numeric(rownames(terminal_nodes)[i])
  node_info <- terminal_nodes[i, ]
  
  # Get decision path
  rule <- path.rpart(final_tree, node = node_id, print.it = FALSE)[[1]]
  conditions <- paste(rule[-1], collapse = " AND ")
  
  # Species in node
  node_rows <- final_tree$where == node_id
  node_data <- tree_dat[node_rows, ]
  
  cat("\nProfile", i, "\n")
  cat(" Conditions:", conditions, "\n")
  cat(" Species (n):", node_info$n, "\n")
  cat(" Mean log(MP+1):", round(node_info$yval, 2), "\n")
  cat(" Mean MP:", round(mean(node_data$Mppersp), 1), "\n")
}



# -----------------------------
# 6. Representative species per terminal node (median-based)
# -----------------------------

library(dplyr)

# Prepare empty dataframe
representatives <- data.frame()

frame <- final_tree$frame
terminal_nodes <- frame[frame$var == "<leaf>", ]
terminal_nodes <- terminal_nodes[order(terminal_nodes$yval, decreasing = TRUE), ]

for(i in 1:nrow(terminal_nodes)) {
  
  node_id <- as.numeric(rownames(terminal_nodes)[i])
  
  # Rows belonging to this node
  node_rows <- final_tree$where == node_id
  node_data <- tree_dat[node_rows, ]
  
  # Median logMP for the node
  node_median_logMP <- median(node_data$logMP, na.rm = TRUE)
  
  # Back-transform to original MP
  node_median_MP <- expm1(node_median_logMP)
  
  # Find species closest to median
  node_data$dist_to_median <- abs(node_data$logMP - node_median_logMP)
  rep_species <- node_data %>%
    arrange(dist_to_median) %>%
    slice_head(n = 3)  # top 3 closest to median
  
  # Save representative species
  representatives <- rbind(
    representatives,
    data.frame(
      Profile = i,
      Species = rep_species$Species,
      Median_logMP = round(node_median_logMP, 2),
      Median_MP = round(node_median_MP, 1),
      Sample_Size = nrow(node_data)
    )
  )
}

# View the table
print(representatives)

# Save to CSV for publication/annotation
write.csv(representatives, "Decision_Tree_Node_Representatives.csv", row.names = FALSE)

cat("✓ Representative species table saved: Decision_Tree_Node_Representatives.csv\n")




# ---------------------------------------
# 6. Save outputs
# ---------------------------------------

write.csv(terminal_nodes,
          "Decision_Tree_Terminal_Nodes.csv",
          row.names = TRUE)

saveRDS(final_tree, "Decision_Tree_Model.rds")

cat("\n✓ Analysis complete\n")



# Get terminal nodes
frame <- final_tree$frame
terminal_nodes <- frame[frame$var == "<leaf>", ]

# Create empty list to store species per node
node_species_list <- list()

for(i in 1:nrow(terminal_nodes)) {
  
  node_id <- as.numeric(rownames(terminal_nodes)[i])
  
  # Rows belonging to this node
  node_rows <- final_tree$where == node_id
  node_data <- tree_dat[node_rows, ]
  
  node_species_list[[paste0("Node_", i)]] <- unique(node_data$Species)
}

node_species_list


# Get terminal nodes
frame <- final_tree$frame
terminal_nodes <- frame[frame$var == "<leaf>", ]

representatives <- data.frame()

for(j in 1:3) {
  
  node_id <- as.numeric(rownames(terminal_nodes)[i])
  node_rows <- final_tree$where == node_id
  node_data <- tree_dat[node_rows, ]
  
  # Median instead of mean
  node_median <- median(node_data$logMP)
  
  # Distance to median
  node_data$distance_to_median <- abs(node_data$logMP - node_median)
  
  # Species closest to median
  rep_species <- node_data[order(node_data$distance_to_median), ][1:3, ]
  
  representatives <- rbind(representatives,
                           data.frame(
                             Profile = i,
                             Species = as.character(rep_species$Species[j]),
                             Median_logMP = round(node_median, 2),
                             Median_MP = round(median(node_data$Mppersp), 1),
                             Sample_Size = nrow(node_data)
                           ))
}

representatives


# Get terminal nodes
frame <- final_tree$frame
terminal_nodes <- frame[frame$var == "<leaf>", ]

representatives <- data.frame()

for(i in 1:nrow(terminal_nodes)) {
  
  node_id <- as.numeric(rownames(terminal_nodes)[i])
  node_rows <- final_tree$where == node_id
  node_data <- tree_dat[node_rows, ]
  
  # Median logMP for the node
  node_median_log <- median(node_data$logMP)
  
  # Distance to median
  node_data$dist <- abs(node_data$logMP - node_median_log)
  
  # Get top 3 species closest to median
  top3 <- node_data[order(node_data$dist), ][1:3, ]
  
  # Add each species with its OWN MP value
  for(j in 1:3) {
    representatives <- rbind(representatives,
                             data.frame(
                               Profile = i,
                               Rank = j,
                               Node_ID = node_id,
                               Species = as.character(top3$Species[j]),
                               Median_logMP = round(node_median_log, 2),
                               Species_MP = round(top3$Mppersp[j], 1),  # ✅ Individual species MP
                               Node_Median_MP = round(median(node_data$Mppersp), 1),
                               Sample_Size = nrow(node_data)
                             ))
  }
}

representatives

