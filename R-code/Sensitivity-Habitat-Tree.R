# ==============================================================================
# FINE-SCALE HABITAT SENSITIVITY TREE
# ==============================================================================

library(rpart)
library(rpart.plot)
library(here)

# Run decision tree using the full multi-category Habitat factor
# (without collapsing into binary aquatic vs terrestrial)
fine_habitat_tree <- rpart(
  logMP ~ distinctiveness + Beak.Length_Culmen + Wing.Length + Mass + Migration + Habitat + Trophic.Niche,
  data = tree_dat,
  method = "anova",
  control = rpart.control(minsplit = 10, minbucket = 5, cp = 0.005, maxdepth = 4, xval = 10)
)

# Print complexity parameter table and variable importance to verify stability
cat("\n=== FINE-SCALE HABITAT SENSITIVITY CP TABLE ===\n")
print(fine_habitat_tree$cptable)

cat("\n=== VARIABLE IMPORTANCE (Fine-Scale Tree) ===\n")
print(fine_habitat_tree$variable.importance)

# Export standalone Fine-Scale Decision Tree Plot PDF for Supplementary Material
pdf(here("output", "Figure_S_Fine_Habitat_Tree.pdf"), width = 10, height = 7)
rpart.plot(
  fine_habitat_tree,
  type = 4,
  extra = 101,
  box.palette = "-RdBu",
  shadow.col = "gray",
  nn = FALSE,
  fallen.leaves = TRUE,
  main = "Fine-Scale Habitat Decision Tree Sensitivity Analysis",
  sub = paste("Multi-category habitat classification (N =", nrow(tree_dat), "species)")
)
dev.off()

cat("\nFine-scale habitat sensitivity tree generated and saved to /output/Figure_S_Fine_Habitat_Tree.pdf\n")
