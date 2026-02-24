################################################################################
# 1. DATA AGGREGATION (Handling 174 rows -> 94 species)
################################################################################
library(rotl)
library(phytools)
library(dplyr)

################################################################################
# 1. DATA PREPARATION AND pPCA ANALYSIS
################################################################################

# Read trait data
trait_data <- read.csv("ppca.csv", stringsAsFactors = FALSE)

# We use dplyr to aggregate so we don't lose data
# This averages duplicates and keeps one row per species
trait_clean <- trait_data %>%
  group_by(Species) %>%
  summarize(
    Wing.Length = mean(Wing.Length, na.rm = TRUE),
    Mass = mean(Mass, na.rm = TRUE),
    Mppersp = mean(Mppersp, na.rm = TRUE)
  ) %>%
  filter(!is.na(Mass) & !is.na(Wing.Length)) %>%
  as.data.frame()

rownames(trait_clean) <- trait_clean$Species

################################################################################
# 2. REAL PHYLOGENY BUILDING
################################################################################

# Match names to Open Tree of Life
taxa_matches <- tnrs_match_names(names = trait_clean$Species)
taxa_found <- taxa_matches[!is.na(taxa_matches$ott_id), ]

# Retrieve the tree
tree <- tol_induced_subtree(ott_ids = taxa_found$ott_id)

# Map labels back to original names
tree_ott_ids <- gsub(".*_ott", "", tree$tip.label)
mapping <- taxa_found$search_string[match(tree_ott_ids, as.character(taxa_found$ott_id))]
tree$tip.label <- mapping

# Standardize names for perfect matching
tree$tip.label <- tolower(gsub(" ", "", tree$tip.label))
rownames(trait_clean) <- tolower(gsub(" ", "", rownames(trait_clean)))

# Intersect to find species present in both data and tree
common_species <- intersect(tree$tip.label, rownames(trait_clean))
tree <- keep.tip(tree, common_species)

# Final formatting
tree <- ape::compute.brlen(tree)
tree <- phytools::force.ultrametric(tree, method="extend")

################################################################################
# 3. pPCA ANALYSIS (On Aggregated Data)
################################################################################

# Align morphological data to the tree
morpho_final <- trait_clean[tree$tip.label, c("Mass", "Wing.Length")]
trait_scaled <- scale(log(morpho_final))

# Run pPCA
ppca_result <- phyl.pca(tree, trait_scaled, method = "lambda", mode = "corr")

# Extract scores for plotting
scores_df <- as.data.frame(ppca_result$S)
scores_df$Species <- tree$tip.label
# Add back the averaged MP data for these specific species
scores_df$MP <- trait_clean[tree$tip.label, "Mppersp"]
scores_df$logMP <- log1p(scores_df$MP)

cat("Analysis complete. Final N =", length(tree$tip.label), "unique species.\n")
cat("Phylogenetic Signal (Lambda):", ppca_result$lambda, "\n")


# 1. Calculate Loadings and variance for the plot
eigenvalues <- diag(ppca_result$Eval)
percent_var <- eigenvalues / sum(eigenvalues) * 100
morpho_loadings <- as.data.frame(ppca_result$L[, 1:2])
colnames(morpho_loadings) <- c("PC1", "PC2")
morpho_loadings$Trait <- rownames(morpho_loadings)

# 2. Calculate MP Correlation Vector
mp_cor_pc1 <- cor(scores_df$logMP, scores_df$PC1, use = "complete.obs")
mp_cor_pc2 <- cor(scores_df$logMP, scores_df$PC2, use = "complete.obs")
mp_vector <- data.frame(PC1 = mp_cor_pc1, PC2 = mp_cor_pc2, Trait = "MP Burden")


library(ggrepel)
# 3. Enhanced Biplot
p_pca <- ggplot() +
  # Reference Circle
  ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = 1), color = "gray90", linetype = "dashed") +
  # Species Points
  geom_point(data = scores_df, aes(x = PC1, y = PC2, color = logMP, size = MP), alpha = 0.7) +
  # Morphological Vectors (Blue)
  geom_segment(data = morpho_loadings, aes(x = 0, y = 0, xend = PC1*3, yend = PC2*3), 
               arrow = arrow(length = unit(0.2, "cm")), color = "#1A3A5A", linewidth = 1) +
  # MP Vector (Red)
  geom_segment(data = mp_vector, aes(x = 0, y = 0, xend = PC1*3, yend = PC2*3), 
               arrow = arrow(length = unit(0.2, "cm")), color = "red", linewidth = 1.2) +
  # Labels
  geom_text_repel(data = morpho_loadings, aes(x = PC1*3.2, y = PC2*3.2, label = Trait), color = "#1A3A5A", fontface = "bold") +
  geom_text_repel(data = mp_vector, aes(x = PC1*3.2, y = PC2*3.2, label = Trait), color = "red", fontface = "bold") +
  # Styling
  scale_color_viridis_c(option = "plasma", name = "log(MP + 1)") +
  labs(x = paste0("pPC1 (", round(percent_var[1], 1), "%)"),
       y = paste0("pPC2 (", round(percent_var[2], 1), "%)"),
       title = "Phylogenetic PCA: Morphology and Microplastic Burden",
       subtitle = paste0("N = 91 species, Pagel's Lambda = ", round(ppca_result$lambda, 3))) +
  theme_minimal() +
  coord_fixed()


# Save as PDF (vector, for publication/editing)
ggsave("pPCA.pdf",
       plot = p_pca,
       width = 15,
       height = 12)







################################################################################
# FINAL ANALYSIS: pPCA WITH SPECIES LABELING & MP GRADIENTS
################################################################################
library(ggplot2)
library(ggrepel)
library(ggforce)
library(viridis)

# 1. Calculate Variance and Loadings
eigenvalues <- diag(ppca_result$Eval)
percent_var <- eigenvalues / sum(eigenvalues) * 100

# Extract Morphological Loadings
morpho_loadings <- as.data.frame(ppca_result$L[, 1:2])
colnames(morpho_loadings) <- c("PC1", "PC2")
morpho_loadings$Trait <- c("Mass", "Wing Length")

# 2. Identify Top 15 Species (Above Median & High Burden)
mp_median <- median(scores_df$MP, na.rm = TRUE)
scores_df$Category <- ifelse(scores_df$MP > mp_median, "Above Median", "Below Median")

# Sort and pick top 15 species for labeling
top_15_species <- scores_df %>%
  filter(Category == "Above Median") %>%
  arrange(desc(MP)) %>%
  slice(1:15)

# 3. Calculate Microplastic Correlation Vector
# This shows the direction of MP increase in the morphospace
mp_cor_pc1 <- cor(scores_df$logMP, scores_df$PC1, use = "complete.obs")
mp_cor_pc2 <- cor(scores_df$logMP, scores_df$PC2, use = "complete.obs")
mp_vector <- data.frame(PC1 = mp_cor_pc1, PC2 = mp_cor_pc2, Trait = "Microplastic Burden")

# Scaling factor for vectors to make them visible on the plot
v_scale <- 2.8 
library(ggforce)

# 4. FINAL PUBLICATION PLOT
final_plot <- ggplot() +
  # Reference Unit Circle (Shows strength of correlation)
  geom_circle(aes(x0 = 0, y0 = 0, r = 1 * v_scale), 
              color = "gray90", linetype = "dotted", linewidth = 0.5) +
  
  # All Species Points
  geom_point(data = scores_df, 
             aes(x = PC1, y = PC2, color = logMP, size = MP), 
             alpha = 0.6) +
  
  # Highlight Top 15 Species with a ring
  geom_point(data = top_15_species, 
             aes(x = PC1, y = PC2), 
             shape = 1, size = 4, color = "black", stroke = 0.8) +
  
  # Morphological vectors (Deep Blue)
  geom_segment(data = morpho_loadings, 
               aes(x = 0, y = 0, xend = PC1 * v_scale, yend = PC2 * v_scale),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"), 
               color = "#003366", linewidth = 1.2) +
  
  # Microplastic Vector (Bright Red)
  geom_segment(data = mp_vector, 
               aes(x = 0, y = 0, xend = PC1 * v_scale, yend = PC2 * v_scale),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"), 
               color = "#CC0000", linewidth = 1.5) +
  
  # Labels for Vectors
  geom_label_repel(data = morpho_loadings, 
                   aes(x = PC1 * v_scale, y = PC2 * v_scale, label = Trait),
                   fill = "white", color = "#003366", fontface = "bold", size = 5) +
  
  geom_label_repel(data = mp_vector, 
                   aes(x = PC1 * v_scale, y = PC2 * v_scale, label = Trait),
                   fill = "#FFF0F0", color = "#CC0000", fontface = "bold", size = 5) +
  
  # Labels for Top 15 Species
  geom_text_repel(data = top_15_species, 
                  aes(x = PC1, y = PC2, label = Species),
                  size = 3.5, fontface = "italic", color = "black",
                  box.padding = 0.5, point.padding = 0.5, max.overlaps = Inf) +
  
  # Aesthetics & Scales
  scale_color_viridis_c(option = "magma", name = expression(log[10](MP + 1))) +
  scale_size_continuous(name = "MP Burden", range = c(2, 10)) +
  
  # Final Formatting for Science
  labs(x = paste0("Phylogenetic PC1 (", round(percent_var[1], 1), "% variance)"),
       y = paste0("Phylogenetic PC2 (", round(percent_var[2], 1), "% variance)"),
       caption = paste0("N = 94 species | Pagel's Lambda = ", round(ppca_result$lambda, 4))) +
  theme_minimal(base_family = "sans") +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "right",
    plot.caption = element_text(hjust = 0, face = "italic", size = 10),
    aspect.ratio = 1
  ) +
  coord_fixed()


final_plot 

# Save for publication
ggsave("pPCA_Science.pdf", final_plot, width = 10, height = 10, device = cairo_pdf)


########################################################################################
# Supplementary
# Verifying species and its relation with other traits


# Check where your species actually fall
species_positions <- scores_df %>%
  select(Species, PC1, PC2, logMP) %>%
  mutate(
    PC1_category = case_when(
      PC1 > quantile(PC1, 0.75) ~ "High PC1 (+)",
      PC1 < quantile(PC1, 0.25) ~ "Low PC1 (-)",
      TRUE ~ "Middle PC1 (~0)"
    )
  ) %>%
  arrange(desc(logMP))

# Print top MP accumulators with their PC1 positions
top_mp_species <- species_positions %>%
  filter(logMP > median(logMP)) %>%
  head(20)

print(top_mp_species)








# Check beak lengths of high PC1 birds------------------------------------------


# First, let's properly merge pPCA scores with original data
# Assuming scores_df has PC1, PC2, Species columns from earlier


# 1. Standardize Species names in trait_clean (morphology)
trait_clean$Species <- tolower(gsub(" ", "", trait_clean$Species))

# 2. Standardize Species names in scores_df (pPCA results)
scores_df$Species <- tolower(gsub(" ", "", scores_df$Species))

# 3. Create the list of species to check in the SAME standardized format
high_pc1_species_clean <- tolower(gsub(" ", "", c(
  "Parus major", "Catharus guttatus", "Pycnonotus sinensis", 
  "Motacilla alba", "Copsychus saularis", "Turdus philomelos"
)))

# 4. Perform the Merge
tree_dat_with_pc <- trait_clean %>%
  select(Species, Mass, Wing.Length) %>% 
  left_join(scores_df %>% select(Species, PC1, PC2, MP, logMP), by = "Species")



cat("\n=== CRITICAL DIAGNOSTIC ===\n")

# 1. PC1 Loadings
cat("\n1. PC1 LOADINGS (Morphological direction):\n")
# Use the loadings we calculated from the ppca_result
loadings_df <- as.data.frame(ppca_result$L[, 1:2])
print(loadings_df)

# 2. Beak Length Distribution
cat("\n2. BEAK LENGTH DISTRIBUTION (Context for 31mm):\n")
mass_stats <- tree_dat_with_pc %>%
  summarise(
    N = n(),
    Mean_mass = mean(Mass, na.rm = TRUE),
    Median_mass = median(Mass, na.rm = TRUE),
    Q25 = quantile(Mass, 0.25, na.rm = TRUE),
    Q75 = quantile(Mass, 0.75, na.rm = TRUE),
    Min = min(Mass, na.rm = TRUE),
    Max = max(Mass, na.rm = TRUE),
    `Prop <164gm` = mean(Mass > 164, na.rm = TRUE)
  )
print(mass_stats) #------------------------------------------------------------- 1) beak_stats

# 3. Check High PC1 Species (Standardized for the join)
high_pc1_species_clean <- tolower(gsub(" ", "", c(
  "Parus major", "Catharus guttatus", "Pycnonotus sinensis", 
  "Motacilla alba", "Copsychus saularis", "Turdus philomelos"
)))

cat("\n3. BEAK LENGTHS OF HIGH PC1 SPECIES:\n")
mass_check <- tree_dat_with_pc %>%
  filter(Species %in% high_pc1_species_clean) %>%
  select(Species, Mass, logMP, PC1) %>%
  arrange(desc(PC1))
print(mass_check) #------------------------------------------------------------- 2) beak check

# 4. Decision Tree Birds
cat("\n4. DECISION TREE'S '<164' BIRDS IN PC1 SPACE:\n")
tree_threshold_birds <- tree_dat_with_pc %>%
  filter(Mass > 164) %>%
  select(Species, Mass, logMP, PC1) %>%
  arrange(desc(logMP)) %>%
  head(10)
print(tree_threshold_birds) # -------------------------------------------------- 3) tree_threshold_birds

# 5. Correlation
cat("\n5. CORRELATION: mass vs PC1:\n")
cor_test <- cor.test(tree_dat_with_pc$Mass, tree_dat_with_pc$PC1, use = "complete.obs")
print(cor_test) #--------------------------------------------------------------- 4) cor_test


# 7. Statistical test
cat("\n7. STATISTICAL COMPARISON:\n")
model1 <- lm(logMP ~ Mass, data = tree_dat_with_pc)
model2 <- lm(logMP ~ PC1, data = tree_dat_with_pc)
model3 <- lm(logMP ~ Mass + PC1, data = tree_dat_with_pc)

cat("Model 1 (Beak only): R² =", round(summary(model1)$r.squared, 3), "\n")
cat("Model 2 (PC1 only): R² =", round(summary(model2)$r.squared, 3), "\n")
cat("Model 3 (Beak + PC1): R² =", round(summary(model3)$r.squared, 3), "\n")

# 8. Check decision tree structure again
cat("\n8. DECISION TREE STRUCTURE:\n")
print(final_tree)
rpart.plot(final_tree, main = "Current Decision Tree")

# Save results for discussion
diagnostic_results <- list(
  mass_stats = mass_stats,
  high_pc1_beaks = mass_check,
  tree_threshold_birds = tree_threshold_birds,
  correlation = cor_test,
  models = list(mass_only = model1, pc1_only = model2, combined = model3)
)

saveRDS(diagnostic_results, "contradiction_diagnostic.rds")





# Visualization-----------------------------------------------------------------
# 2 Statistical comparision
library(ggplot2)
library(patchwork)

# 6. Visualize Patterns
# Panel 1: Beak vs PC1
p1 <- ggplot(tree_dat_with_pc, aes(x = Mass, y = PC1)) +
  geom_point(aes(color = logMP), size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", color = "red", se = FALSE, linetype = "dashed") +
  geom_vline(xintercept = 164, linetype = "dotted", color = "blue", linewidth = 1) +
  scale_color_viridis_c(option = "plasma") +
  labs(title = "Mass vs PC1", x = "Mass (gm)", y = "PC1 Score") +
  theme_minimal()

# Panel 2: The U-Shape Pattern (Quadratic Fit)
p2 <- ggplot(tree_dat_with_pc, aes(x = Mass, y = logMP)) +
  geom_point(aes(color = PC1), size = 3, alpha = 0.7) +
  geom_vline(xintercept = 164, linetype = "dotted", color = "blue", linewidth = 1) +
  # Use a quadratic formula to show the U-shape: y ~ poly(x, 2)
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "darkred", size = 1.2) +
  scale_color_viridis_c() +
  labs(title = "U-Shaped Accumulation Pattern", x = "Mass (gm)", y = "log(MP+1)") +
  theme_minimal()

# Panel 3: PC1 vs logMP
p3 <- ggplot(tree_dat_with_pc, aes(x = PC1, y = logMP)) +
  geom_point(aes(color = Mass >= 164), size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  scale_color_manual(values = c("gray", "darkblue"), name = "Mass<164gm") +
  labs(title = "PC1 vs MP Accumulation", x = "PC1 Score", y = "log(MP+1)") +
  theme_minimal()

combined <- (p1 | p2 | p3) + plot_annotation(
  #title = "Mechanistic Resolution: Morphological Drivers of MP Accumulation",
  #subtitle = "Note how both extremes of the beak length spectrum (Short & Long) exhibit high MP values."
)

print(combined)
ggsave("Figure_Diagnostic_Final.pdf", combined, width = 16, height = 5)


##################################################################################
# Use comparision model
# Linear vs quadratic model
m_lin <- lm(logMP ~ Mass, data = tree_dat_with_pc)
m_quad <- lm(logMP ~ Mass + I(Mass^2), data = tree_dat_with_pc)

anova(m_lin, m_quad)
summary(m_quad)


library(mgcv)

m_gam <- gam(
  logMP ~ s(Mass, k = 4),
  data = tree_dat_with_pc,
  method = "REML"
)

summary(m_gam)
plot(m_gam, shade = TRUE, residuals = TRUE)


AIC(m_lin, m_gam)


# Results
# Microplastic (MP) accumulation exhibited a significant non-linear association
# with beak length, outperforming the simple linear framework. A quadratic model 
# provided a significantly superior fit to the data 
# ($\Delta RSS = 47.3, F_{1,88} = 19.1, p < 0.001$), characterized by a 
# negative linear term ($\beta_1 = -0.05, p < 0.001$) and a positive quadratic
# term ($\beta_2 = 1.57 \times 10^{-4}, p < 0.001$). This confirms a convex 
# (U-shaped) relationship where MP burdens are elevated at both morphological 
# extremes (abbreviated and elongated beaks), explaining approximately 19.6% 
# of the total variance ($R^2_{adj} = 0.196$).

# Interpretation
# The transition from a univariate "U-shape" to a monotonic increase in
# phylogenetic principal component space (pPC1) provides a critical mechanistic
# insight. It suggests that the apparent non-linearity of a single trait—beak
# length—is a projection of a more complex, multidimensional feeding morphology. 
# When evolutionarily conserved traits (body size, wing length, and beak shape) 
# are integrated into a single pPC1 axis, MP burden increases linearly. 
# This indicates that while "short beaks" and "long beaks" both face high risk,
# they represent the same high-exposure evolutionary trajectory when 
# multidimensional morphospace is considered.


# 3. Statistical Comparison (Table Generation)

# Figure caption. Mechanistic drivers of microplastic (MP) accumulation across avian morphospace. (Left) Correlation between beak length and Phylogenetic Principal Component 1 (pPC1). Species with abbreviated beaks occupy high pPC1 space, while elongated beaks correspond to lower pPC1 scores. (Center) Quadratic regression reveals a U-shaped accumulation pattern ($y \sim \text{poly}(x, 2)$), where MP burden is maximized at both morphological extremes. The dashed blue line indicates the 31mm threshold identified by non-phylogenetic decision tree models. (Right) Linear regression of pPC1 against $\log_{10}$-transformed MP burden ($log(MP+1)$), highlighting that while pPC1 captures a significant portion of the variance, it primarily identifies the risk associated with small-bodied, short-beaked taxa (High pPC1/Blue points).

cat("\n7. STATISTICAL COMPARISON (Linear vs Quadratic):\n")
model_linear <- lm(logMP ~ Mass, data = tree_dat_with_pc)
model_quad   <- lm(logMP ~ poly(Mass, 2), data = tree_dat_with_pc) # The U-shape
model_pc1    <- lm(logMP ~ PC1, data = tree_dat_with_pc)

summary(model_linear)
summary(model_quad)
summary(model_pc1 )

cat("Model 1 (Linear Beak): R² =", round(summary(model_linear)$r.squared, 3), "\n")
cat("Model 2 (Quadratic Beak/U-Shape): R² =", round(summary(model_quad)$r.squared, 3), "\n")
cat("Model 3 (PC1 Score): R² =", round(summary(model_pc1)$r.squared, 3), "\n")

# Interpretation Tip: If Model 2 R-squared is higher than Model 1, 
# you have statistically proven the U-shaped relationship.







################################################################################
# PGLS analysis
################################################################################
# PGLS ANALYSIS: COMBINING MORPHOLOGY & ECOLOGY
################################################################################

library(caper)
library(dplyr)
library(ape)

# ==============================================================================
# 1. PREPARE ECOLOGICAL DATA
# ==============================================================================
cat("\n=== 1. PREPARING ECOLOGICAL DATA ===\n")

# Aggregating trait_data to ensure ONE row per species
# We take the 'first' non-NA value for categorical traits
eco_traits_clean <- trait_data %>%
  group_by(Species) %>%
  summarise(
    Feeding_habit = first(na.omit(Feeding_habit)),
    Habitat = first(na.omit(Habitat)),
    Migration = first(na.omit(Migration)),
    Species_type = first(na.omit(Species_type))
  ) %>%
  # Standardize Species names to match tree (lowercase, no spaces/punct)
  mutate(Species_std = tolower(gsub("[[:punct:][:space:]]", "", Species))) %>%
  filter(!is.na(Species_std) & Species_std != "")

# ==============================================================================
# 2. MERGE WITH pPCA SCORES
# ==============================================================================
cat("=== 2. MERGING WITH pPCA SCORES ===\n")

# Prepare scores_df (from your previous pPCA step)
scores_clean <- scores_df %>%
  mutate(Species_std = tolower(gsub("[[:punct:][:space:]]", "", Species))) %>%
  select(Species_std, PC1, PC2, logMP, MP)

# Merge everything into one master dataframe
pgls_data <- inner_join(scores_clean, eco_traits_clean, by = "Species_std") %>%
  # Convert character columns to factors for PGLS
  mutate(across(c(Feeding_habit, Habitat, Migration, Species_type), as.factor)) %>%
  # Remove rows with NA in key columns (PGLS cannot handle NAs)
  filter(!is.na(logMP), !is.na(PC1), !is.na(Habitat)) %>%
  # Ensure distinct rows
  distinct(Species_std, .keep_all = TRUE) %>%
  as.data.frame()

# Set row names to species names (required for caper)
rownames(pgls_data) <- pgls_data$Species_std

cat("Data ready. N =", nrow(pgls_data), "species.\n")

# ==============================================================================
# 3. SYNCHRONIZE TREE AND DATA
# ==============================================================================
cat("=== 3. SYNCHRONIZING TREE ===\n")

# Standardize Tree Tip Labels to match data exactly
tree$tip.label <- tolower(gsub("[[:punct:][:space:]]", "", tree$tip.label))

# Find the intersection
common_species <- intersect(tree$tip.label, pgls_data$Species_std)

if(length(common_species) == 0) {
  stop("CRITICAL ERROR: No matching species found between tree and data!")
} 

# Prune tree and data to match exactly
pgls_tree <- keep.tip(tree, common_species)
pgls_data_final <- pgls_data[common_species, ]

cat("Synchronization successful. Analyzing", length(pgls_tree$tip.label), "species.\n")

# ==============================================================================
# 4. CREATE COMPARATIVE DATA OBJECT
# ==============================================================================

comp_data <- comparative.data(
  phy = pgls_tree,
  data = pgls_data_final,
  names.col = "Species_std",
  vcv = TRUE,
  na.omit = FALSE, # We already cleaned NAs
  warn.dropped = TRUE
)

# ==============================================================================
# 5. RUN PGLS MODELS
# ==============================================================================
cat("\n=== 5. RUNNING PGLS MODELS ===\n")

# Model A: Morphology Only (Testing your PC1/PC2 hypothesis)
model_morph <- pgls(logMP ~ PC1 + PC2, 
                    data = comp_data, lambda = "ML")

# Model B: Ecology Only (Testing Habitat/Feeding)
model_eco <- pgls(logMP ~ Habitat + Feeding_habit + Migration, 
                  data = comp_data, lambda = "ML")

# Model C: The "Global" Model (Morphology + Ecology combined)
# This tests if Habitat matters *after* accounting for beak shape
model_global <- pgls(logMP ~ PC1 + PC2 + Habitat + Feeding_habit + Migration, 
                     data = comp_data, lambda = "ML")

# ==============================================================================
# 6. RESULTS & COMPARISON
# ==============================================================================

# 1. Compare Models (AIC)
aic_table <- data.frame(
  Model = c("Morphology Only", "Ecology Only", "Global (Morph+Eco)"),
  AIC = c(AIC(model_morph), AIC(model_eco), AIC(model_global)),
  R_Squared = c(summary(model_morph)$r.squared, 
                summary(model_eco)$r.squared, 
                summary(model_global)$r.squared)
)
print(aic_table)

# 2. Detailed Summary of the Best Model (Usually Global)
cat("\n=== SUMMARY OF GLOBAL MODEL ===\n")
summary(model_global)

# 3. ANOVA Table (To check significance of categorical variables)
cat("\n=== ANOVA TABLE (Significance of Factors) ===\n")
anova(model_global)

# ==============================================================================
# 7. EXPORT RESULTS
# ==============================================================================

# Save coefficients to CSV
coef_df <- as.data.frame(summary(model_global)$coefficients)
write.csv(coef_df, "PGLS_Global_Coefficients.csv")

cat("\nAnalysis Complete. Results saved to 'PGLS_Global_Coefficients.csv'.\n")
