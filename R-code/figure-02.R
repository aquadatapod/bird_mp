# ==============================================================================
# FIGURE 2: TISSUE-LEVEL & TRAIT-BASED MICROPLASTIC ACCUMULATION
# Outputs Panels B (Tissue comparisons), C (Body mass), D (Wing length), 
# and Panel E (Decision Tree Analysis)
# ==============================================================================

# 1. Load Required Libraries
# ------------------------------------------------------------------------------
library(here)
library(tidyverse)
library(ggrepel)
library(rstatix)
library(multcompView)
library(readxl)
library(patchwork)
library(mgcv)
library(car)
library(cluster)
library(rpart)
library(rpart.plot)

# 2. Shared Formatting & Theme
# ------------------------------------------------------------------------------
mp_low <- c(
  "#5e9cc2", "#8fc6a4", "#588157", "#d9b650",
  "#d96c50", "#4d7a9c", "#6a8d92", "#b68b6d",
  "#7aa5b5", "#8b9a6e", "#A23B72", "#2E86AB"
)

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(size = 8, color = "black"),
      axis.title = element_text(size = 9, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.3),
      plot.margin = margin(10, 10, 10, 10),
      plot.title = element_text(size = 9, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 8, face = "italic", color = "gray30")
    )
}

# ==============================================================================
# SECTION 1: PANELS B, C, AND D (Tissue Comparison & Trait Regressions)
# ==============================================================================

bird.dat <- read_excel(here("data/processed", "row_level_02.xlsx"))

# Clean tissue dataset
df_tissue <- bird.dat %>%
  filter(
    !is.na(Mppersp), !is.na(Tissue),
    !is.na(Longitude), !is.na(Latitude)
  ) %>%
  mutate(
    Tissue = factor(
      Tissue,
      levels = c("Crop", "Esophagus", "FP", "Gizzard", "GIT", "RP", "Dropping", "Other")
    )
  ) %>%
  filter(!is.na(Tissue)) %>%
  droplevels()

# Panel B Statistics & Letters
kw_test <- kruskal_test(df_tissue, Mppersp ~ Tissue)
dunn_res <- df_tissue %>% dunn_test(Mppersp ~ Tissue, p.adjust.method = "BH")

pw <- dunn_res$p.adj
names(pw) <- paste(dunn_res$group1, dunn_res$group2, sep = "-")
letters_raw <- multcompView::multcompLetters(pw)$Letters

letters_df <- data.frame(
  Tissue = names(letters_raw),
  Letters = letters_raw,
  stringsAsFactors = FALSE
)

y_pos <- df_tissue %>%
  group_by(Tissue) %>%
  summarise(y = max(Mppersp, na.rm = TRUE) * 1.5, .groups = "drop")

letters_df <- left_join(letters_df, y_pos, by = "Tissue")
n_df <- df_tissue %>% group_by(Tissue) %>% summarise(n = n(), .groups = "drop")

# Build Panel B
fig2b <- ggplot(df_tissue, aes(x = Tissue, y = Mppersp, fill = Tissue)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5, size = 0.3) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.0, color = "gray30") +
  geom_text(
    data = letters_df, aes(x = Tissue, y = y, label = Letters),
    inherit.aes = FALSE, size = 3.5, fontface = "bold"
  ) +
  geom_text(
    data = n_df, aes(x = Tissue, y = -0.5, label = paste0("N=", n)),
    inherit.aes = FALSE, size = 2.5, color = "gray20"
  ) +
  scale_fill_manual(values = rep(mp_low, length.out = length(levels(df_tissue$Tissue)))) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 1, 5, 10, 20, 50, 100, 200, 400, 800),
    labels = c("0", "1", "5", "10", "20", "50", "100", "200", "400", "800")
  ) +
  labs(
    x = "", y = expression("Microplastics per individual ("*MP~ind^{-1}*")"),
    title = paste0("Kruskal-Wallis: p = ", format.pval(kw_test$p, digits = 2))
  ) +
  theme_pub() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))


################################################################################

library(ggplot2)
library(dplyr)

# 1. Prepare clean tissue data and explicitly set the factor levels in your preferred order
df_main_tissue <- df_tissue %>%
  mutate(Tissue_Clean = case_when(
    Tissue %in% c("Crop", "Esophagus", "Gizzard", "GIT") ~ as.character(Tissue),
    TRUE ~ "External/Peripheral"
  )) %>%
  mutate(Tissue_Clean = factor(Tissue_Clean, levels = c("Esophagus", "Crop", "Gizzard", "GIT", "External/Peripheral")))

# 2. Calculate sample sizes per category for annotation
tissue_counts <- df_main_tissue %>%
  group_by(Tissue_Clean) %>%
  summarise(n = n(), .groups = "drop")

# 3. Generate Elsevier/Environmental Pollution-styled plot (using x = Tissue_Clean directly)
fig2b_main_styled <- ggplot(df_main_tissue, aes(x = Tissue_Clean, y = Mppersp, fill = Tissue_Clean)) +
  # Clean boxplots with subtle outlines
  geom_boxplot(alpha = 0.65, outlier.shape = 21, outlier.size = 1.2, outlier.alpha = 0.4, width = 0.45, linewidth = 0.4) +
  # Jittered data points for transparency
  geom_jitter(width = 0.12, alpha = 0.3, size = 0.8, color = "gray30") +
  # Sample size annotations positioned neatly near the baseline
  geom_text(
    data = tissue_counts,
    aes(x = Tissue_Clean, y = -0.3, label = paste0("n = ", n)),
    inherit.aes = FALSE, size = 3, color = "black", fontface = "italic"
  ) +
  # Professional palette matching journal standards
  scale_fill_viridis_d(option = "mako", begin = 0.3, end = 0.85) +
  # Pseudo-log scale matching your statistical backbone
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 1, 5, 10, 50, 100, 400, 800),
    labels = c("0", "1", "5", "10", "50", "100", "400", "800")
  ) +
  labs(
    x = "Gastrointestinal Compartment / Tissue", 
    y = expression(bold("Microplastics per individual") ~ bold((MP~ind^{-1}))),
    subtitle = paste0("Kruskal-Wallis test: p = ", format.pval(kw_test$p, digits = 2))
  ) +
  # Classic journal theme with light y-axis grid lines for data tracking
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 9.5, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9.5, face = "italic", color = "gray35"),
    plot.margin = margin(10, 10, 10, 10)
  )

print(fig2b_main_styled)

################################################################################





# Species-level Aggregation for Panels C & D
digestive_tissues <- c("Crop", "Esophagus", "Gizzard", "GIT", "RP")
digestive <- bird.dat %>%
  filter(
    Tissue %in% digestive_tissues,
    !is.na(Mppersp), !is.na(Mass), !is.na(Beak.Length_Culmen)
  ) %>%
  mutate(
    logMP   = log1p(Mppersp),
    logMass = log(Mass),
    logBeak = log(Beak.Length_Culmen)
  )

species_dig <- digestive %>%
  group_by(Species) %>%
  summarise(
    n_records  = n(),
    logMass    = first(logMass),
    logBeak    = first(logBeak),
    logMP_median = median(logMP, na.rm = TRUE),
    .groups    = "drop"
  )

cor_mass_sp <- cor.test(species_dig$logMass, species_dig$logMP_median, method = "spearman")
cor_beak_sp <- cor.test(species_dig$logBeak, species_dig$logMP_median, method = "spearman")

lm_mass <- lm(logMP_median ~ logMass, data = species_dig)
r2_mass <- summary(lm_mass)$r.squared

gam_beak <- mgcv::gam(logMP_median ~ s(logBeak, k = 4), data = species_dig, method = "REML")
edf_beak <- summary(gam_beak)$s.table[1, "edf"]
dev_beak <- summary(gam_beak)$dev.expl * 100

labeled_species <- species_dig %>%
  filter(Species %in% c("Gallus gallus", "Tachycineta bicolor"))

# Build Panel C
fig2c <- ggplot(species_dig, aes(x = logMass, y = logMP_median)) +
  annotate("rect", xmin = 7.5, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray90", alpha = 0.5) +
  geom_rug(sides = "b", alpha = 0.5, color = "gray30") +
  geom_point(color = "#4d7a9c", alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", color = "black", fill = "gray70", size = 0.8) +
  geom_text_repel(
    data = labeled_species, aes(label = paste0(Species, " (n=", n_records, ")")),
    size = 2.8, fontface = "italic", box.padding = 0.5
  ) +
  labs(
    x = expression("log"[e]*"(Body mass, g)"),
    y = expression("median log(MP + 1) in digestive tract"),
    title = sprintf("Linear fit: R² = %.2f\nrho = %.2f, p < 0.001 (n = %d species)", 
                    r2_mass, cor_mass_sp$estimate, nrow(species_dig))
  ) +
  theme_pub()

# Build Panel D
fig2d <- ggplot(species_dig, aes(x = logBeak, y = logMP_median)) +
  annotate("rect", xmin = 4.2, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray90", alpha = 0.5) +
  geom_rug(sides = "b", alpha = 0.5, color = "gray30") +
  geom_point(color = "#4d7a9c", alpha = 0.7, size = 2) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4), color = "black", fill = "gray70", size = 0.8) +
  geom_text_repel(
    data = labeled_species, aes(label = paste0(Species, " (n=", n_records, ")")),
    size = 2.8, fontface = "italic", box.padding = 0.5
  ) +
  labs(
    x = expression("log"[e]*"(Beak length, mm)"),
    y = expression("median log(MP + 1) in digestive tract"),
    title = sprintf("GAM: edf = %.1f, dev. expl. = %.0f%%\nrho = %.2f, p < 0.001 (n = %d species)", 
                    edf_beak, dev_beak, cor_beak_sp$estimate, nrow(species_dig))
  ) +
  theme_pub()

# ==============================================================================
# SECTION 2: PANEL E (Functional Clustering & Decision Tree Modeling)
# ==============================================================================

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

# 1. Load Data
tree_dat_raw <- read.csv(here("data/processed", "row_level.csv"), stringsAsFactors = FALSE) %>%
  drop_na(Species, Mppersp, Beak.Length_Culmen, Wing.Length, Mass) %>%
  mutate(across(c(Habitat, Migration, Trophic.Niche), as.factor))

# 2. Species-level median aggregation
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
    across(c(Habitat, Migration, Trophic.Niche), as.factor),
    logMP   = log1p(Mppersp),
    logMass = log10(Mass),
    logBeak = log10(Beak.Length_Culmen)
  )

# ------------------------------------------------------------------------------
# Spearman Correlation Warning (Ties in Data)
# Add `exact = FALSE` to perform asymptotic p-value estimation
# ------------------------------------------------------------------------------
cor_mass_sp <- cor.test(tree_dat$logMass, tree_dat$logMP, method = "spearman", exact = FALSE)
cor_beak_sp <- cor.test(tree_dat$logBeak, tree_dat$logMP, method = "spearman", exact = FALSE)

# ------------------------------------------------------------------------------
# Decision Tree Over-Pruning
# Adjust control parameters for species-level dataset (N ~ 90-100)
# ------------------------------------------------------------------------------
gower_dist <- daisy(
  tree_dat %>% dplyr::select(Beak.Length_Culmen, Habitat, Mass, Migration, Trophic.Niche, Wing.Length),
  metric = "gower"
)

pam_fit <- pam(gower_dist, k = 3)
tree_dat$cluster <- factor(pam_fit$clustering)
tree_dat$distinctiveness <- apply(as.matrix(gower_dist), 1, mean)

# Relaxed split controls suited for species-level sample sizes:
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
  
  # Prevent pruning back to root node
  if (best_cp == cp_table[1, "CP"] && nrow(cp_table) > 1) {
    best_cp <- cp_table[2, "CP"]
  }
  final_tree <- prune(tree_model, cp = best_cp)
} else {
  final_tree <- tree_model
}

# ------------------------------------------------------------------------------
# Node Representative Extraction & File Exports
# ------------------------------------------------------------------------------
if (!dir.exists(here("output"))) {
  dir.create(here("output"), recursive = TRUE)
}

frame <- final_tree$frame
terminal_nodes <- frame[frame$var == "<leaf>", ]

representatives <- data.frame()
for (i in 1:nrow(terminal_nodes)) {
  node_id   <- as.numeric(rownames(terminal_nodes)[i])
  node_rows <- final_tree$where == node_id
  node_data <- tree_dat[node_rows, ]
  
  if (nrow(node_data) > 0) {
    node_median_log <- median(node_data$logMP, na.rm = TRUE)
    node_data$dist  <- abs(node_data$logMP - node_median_log)
    top3 <- node_data[order(node_data$dist), ][1:min(3, nrow(node_data)), ]
    
    for (j in 1:nrow(top3)) {
      representatives <- rbind(
        representatives,
        data.frame(
          Profile        = i,
          Rank           = j,
          Node_ID        = node_id,
          Species        = as.character(top3$Species[j]),
          Median_logMP   = round(node_median_log, 2),
          Species_MP     = round(top3$Mppersp[j], 1),
          Node_Median_MP = round(median(node_data$Mppersp, na.rm = TRUE), 1),
          Sample_Size    = nrow(node_data)
        )
      )
    }
  }
}

write.csv(representatives, here("output", "Decision_Tree_Node_Representatives.csv"), row.names = FALSE)
saveRDS(final_tree, here("output", "Decision_Tree_Model.rds"))

# Export standalone Decision Tree Plot PDF
pdf(here("output", "Figure_2e_decision_tree.pdf"), width = 10, height = 7)
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

# ==============================================================================
# SECTION 3: ASSEMBLE COMBINED FIGURE
# ==============================================================================

# Combined panel B, C, D
fig2_bcd <- (fig2b_main_styled | fig2c | fig2d) + 
  plot_annotation(tag_levels = list(c("b", "c", "d"))) &
  theme(plot.tag = element_text(size = 12, face = "bold"))

ggsave(here("output", "Figure_2_BCD.pdf"), fig2_bcd, width = 34, height = 13, units = "cm", device = "pdf")
ggsave(here("output", "Figure_2_BCD.tiff"), fig2_bcd, width = 28, height = 9, units = "cm", dpi = 300, compression = "lzw")

print(fig2_bcd)
cat("\n Processing complete. Outputs saved to /output folder.\n")





# Extract node assignments and terminal-node summary
node_summary <- tree_dat %>%
  mutate(node = as.integer(final_tree$where)) %>%
  group_by(node) %>%
  summarise(
    n_species = n(),
    percent_species = round(100 * n_species / nrow(tree_dat), 1),
    
    # Response summary
    median_logMP = median(logMP, na.rm = TRUE),
    median_MP = expm1(median_logMP),
    
    # Key morphological variable used by the tree
    median_Wing_Length = median(Wing.Length, na.rm = TRUE),
    
    # Supporting body-size information
    median_Mass = median(Mass, na.rm = TRUE),
    
    # Habitat information
    top_habitats = paste(
      head(sort(table(Habitat), decreasing = TRUE), 3),
      collapse = ", "
    ),
    
    # Representative species
    example_species = paste(
      head(Species, 3),
      collapse = "; "
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    median_logMP = round(median_logMP, 2),
    median_MP = round(median_MP, 1),
    median_Wing_Length = round(median_Wing_Length, 1),
    median_Mass = round(median_Mass, 1)
  )

print(node_summary)
write.csv(node_summary, here("output", "Decision_Tree_Node_Descriptions.csv"), row.names = FALSE)



# log MP to original MP values per individuals
library(dplyr)
library(rpart)
library(rpart.plot)

# ------------------------------------------------------------------------------
# 1. Back-Transformation & Summary Function
# ------------------------------------------------------------------------------
# NOTE: Check if you used log1p(x) [natural log, ln(x + 1)] or log10(x + 1).
# We calculate both below so you can verify which transformation matches your pipeline.

node_summary_counts <- node_summary %>%
  mutate(
    # Natural log back-transformation: exp(x) - 1
    median_MP_natural = exp(median_logMP) - 1,
    
    # Base-10 log back-transformation: 10^(x) - 1
    median_MP_log10 = (10^median_logMP) - 1
  )

# Display the calculated raw counts
print(node_summary_counts)









################################################################################
# PHYLOGENETIC PCA (pPCA) & PGLS ANALYSIS WORKFLOW--------------------------------Figure 2f
################################################################################

# --- 0. LIBRARIES ---
library(here)
library(tidyverse)
library(rotl)
library(phytools)
library(ape)
library(caper)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(viridis)
library(patchwork)
library(mgcv)

# Helper function for categorical mode
get_mode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

################################################################################
# 1. DATA PREPARATION & SPECIES MEDIAN AGGREGATION
################################################################################

raw_data <- read.csv(
  here("data", "processed", "row_level_02.csv"),
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

# Species-level median aggregation
trait_clean <- raw_data %>%
  drop_na(Species, Mass, Wing.Length, Mppersp) %>%
  group_by(Species) %>%
  summarise(
    Mass               = median(Mass, na.rm = TRUE),
    Wing.Length        = median(Wing.Length, na.rm = TRUE),
    Mppersp            = median(Mppersp, na.rm = TRUE),
    Feeding_habit      = get_mode(Trophic.Niche),
    Habitat            = get_mode(Habitat),
    Migration          = get_mode(Migration),
    Species_type       = get_mode(Species_type),
    .groups            = "drop"
  ) %>%
  mutate(
    # Standardize species key: lower case, single underscore separator
    species_std = tolower(gsub("[[:punct:][:space:]]+", "_", trimws(Species))),
    log10_Mass  = log10(Mass),
    log10_Wing  = log10(Wing.Length),
    logMP       = log1p(Mppersp)
  ) %>%
  filter(!is.na(Mass) & !is.na(Wing.Length) & Mass > 0 & Wing.Length > 0)

################################################################################
# 2. PHYLOGENY RECONSTRUCTION (Corrected)
################################################################################

# Step 1: Match taxonomic names (specify context_name if needed)
taxa_matches <- tnrs_match_names(
  names = trait_clean$Species, 
  context_name = "Animals"  # Restricts matching to animal taxonomy
)

# Filter for successfully matched taxa
taxa_found <- taxa_matches[!is.na(taxa_matches$ott_id), ]

# Step 2: Retrieve the induced subtree (WITHOUT default_hierarchy)
tree_raw <- tol_induced_subtree(ott_ids = taxa_found$ott_id)

# Map labels back to standardized species names
tree_ott_ids <- gsub(".*_ott", "", tree_raw$tip.label)
taxa_found$ott_id_str <- as.character(taxa_found$ott_id)

mapping <- taxa_found %>%
  dplyr::select(ott_id_str, search_string) %>%
  distinct(ott_id_str, .keep_all = TRUE)

matched_indices <- match(tree_ott_ids, mapping$ott_id_str)
tree_raw$tip.label <- mapping$search_string[matched_indices]

# Standardize tip labels for matching with data
tree_raw$tip.label <- tolower(gsub("[[:punct:][:space:]]+", "_", trimws(tree_raw$tip.label)))

# Prune common species and format branch lengths
common_species <- intersect(tree_raw$tip.label, trait_clean$species_std)
tree <- keep.tip(tree_raw, common_species)
tree <- ape::compute.brlen(tree)
tree <- phytools::force.ultrametric(tree, method = "extend")

# Align dataframe strictly with tree tip labels
trait_final <- trait_clean %>%
  filter(species_std %in% tree$tip.label) %>%
  arrange(match(species_std, tree$tip.label)) %>%
  as.data.frame()

rownames(trait_final) <- trait_final$species_std

################################################################################
# 3. PHYLOGENETIC PCA (pPCA)
################################################################################

# Morphological matrix (log10-transformed)
morpho_mat <- trait_final[, c("log10_Mass", "log10_Wing")]

# Run pPCA (Pagel's Lambda model)
ppca_result <- phyl.pca(tree, morpho_mat, method = "lambda", mode = "corr")

# Extract variance explained
eigenvalues <- diag(ppca_result$Eval)
percent_var <- (eigenvalues / sum(eigenvalues)) * 100

# Assemble species pPCA scores
scores_df <- as.data.frame(ppca_result$S)
colnames(scores_df) <- c("PC1", "PC2")
scores_df$species_std <- rownames(scores_df)
scores_df <- left_join(scores_df, trait_final, by = "species_std")

# Calculate Morphological Loadings
morpho_loadings <- as.data.frame(ppca_result$L[, 1:2])
colnames(morpho_loadings) <- c("PC1", "PC2")
morpho_loadings$Trait <- c("Body Mass", "Wing Length")

# Calculate Microplastic Burden Correlation Vector
mp_cor_pc1 <- cor(scores_df$logMP, scores_df$PC1, use = "complete.obs")
mp_cor_pc2 <- cor(scores_df$logMP, scores_df$PC2, use = "complete.obs")
mp_vector  <- data.frame(PC1 = mp_cor_pc1, PC2 = mp_cor_pc2, Trait = "MP Burden")

################################################################################
# 4. PUBLICATION-READY pPCA BIPLOT
################################################################################

# Ensure output directory exists before saving anything
if (!dir.exists(here("output"))) {
  dir.create(here("output"), recursive = TRUE)
}

# Identify Top 10 High-Burden Accumulators for Labeling (Keeps plot clean)
top_species <- scores_df %>%
  arrange(desc(Mppersp)) %>%
  slice(1:10) %>%
  mutate(Display_Name = gsub("_", " ", stringr::str_to_sentence(species_std)))

# Scale vectors to extend past the species cloud for clear reading
v_scale <- max(abs(c(scores_df$PC1, scores_df$PC2))) * 0.85

# Adjust loading label endpoints slightly outward to prevent label overlap
morpho_loadings_plot <- morpho_loadings %>%
  mutate(
    x_end = PC1 * v_scale,
    y_end = PC2 * v_scale,
    # Manual slight offset for Wing Length vs Body Mass text
    x_label = x_end * 1.15,
    y_label = y_end * 1.15
  )

mp_vector_plot <- mp_vector %>%
  mutate(
    x_end = PC1 * v_scale,
    y_end = PC2 * v_scale,
    x_label = x_end * 1.15,
    y_label = y_end * 1.15
  )

pPCA_plot <- ggplot() +
  # Reference Circle
  geom_circle(aes(x0 = 0, y0 = 0, r = v_scale), color = "grey80", linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = 0, color = "grey90", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey90", linewidth = 0.3) +
  
  # All Species Points
  geom_point(
    data = scores_df,
    aes(x = PC1, y = PC2, color = logMP, size = Mppersp),
    alpha = 0.65
  ) +
  
  # High-Burden Species Highlights
  geom_point(
    data = top_species,
    aes(x = PC1, y = PC2),
    shape = 21, stroke = 0.8, color = "black", fill = NA, size = 4.5
  ) +
  
  # Morphological Loadings Vectors (Deep Blue)
  geom_segment(
    data = morpho_loadings_plot,
    aes(x = 0, y = 0, xend = x_end, yend = y_end),
    arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
    color = "#0D3B66", linewidth = 1.1
  ) +
  
  # Microplastic Correlation Vector (Bright Red)
  geom_segment(
    data = mp_vector_plot,
    aes(x = 0, y = 0, xend = x_end, yend = y_end),
    arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
    color = "#D62828", linewidth = 1.2
  ) +
  
  # Morphological Vector Labels (Clean text without solid box blocking data)
  geom_text_repel(
    data = morpho_loadings_plot,
    aes(x = x_label, y = y_label, label = Trait),
    color = "#0D3B66", fontface = "bold", size = 4.5,
    box.padding = 0.3, segment.color = NA
  ) +
  
  # MP Vector Label
  geom_text_repel(
    data = mp_vector_plot,
    aes(x = x_label, y = y_label, label = Trait),
    color = "#D62828", fontface = "bold", size = 4.5,
    box.padding = 0.3, segment.color = NA
  ) +
  
  # Species Name Labels (Italic, repelled cleanly)
  geom_text_repel(
    data = top_species,
    aes(x = PC1, y = PC2, label = Display_Name),
    size = 3.2, fontface = "italic", color = "grey15",
    box.padding = 0.5, point.padding = 0.4,
    min.segment.length = 0, segment.color = "grey50", segment.size = 0.3,
    max.overlaps = 15
  ) +
  
  # Scales & Colors
  scale_color_viridis_c(option = "magma", name = expression(log[10](MP + 1))) +
  scale_size_continuous(name = expression("Median MPs g"^-1), range = c(2, 7)) +
  
  # Axis Titles
  labs(
    x = paste0("Phylogenetic PC1 (", round(percent_var[1], 1), "% variance)"),
    y = paste0("Phylogenetic PC2 (", round(percent_var[2], 1), "% variance)"),
    caption = paste0("N = ", length(tree$tip.label), " species | Pagel's \u03BB = ", round(ppca_result$lambda, 3))
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold", size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    plot.caption = element_text(hjust = 0, face = "italic", color = "grey30"),
    aspect.ratio = 1
  ) +
  coord_fixed()

print(pPCA_plot)

# Save high-resolution PDF
ggsave(here("output", "Figure_pPCA_Clean.pdf"), pPCA_plot, width = 9, height = 8, device = cairo_pdf)


###################################################################################
###############################################################################

################################################################################
# 6. PHYLOGENETIC GENERALIZED LEAST SQUARES (PGLS)
################################################################################

# Prepare dataset for caper
pgls_df <- scores_df %>%
  mutate(
    across(c(Feeding_habit, Habitat, Migration, Species_type), as.factor),
    species_std_caper = species_std
  ) %>%
  filter(!is.na(logMP) & !is.na(PC1) & !is.na(Habitat)) %>%
  as.data.frame()

rownames(pgls_df) <- pgls_df$species_std_caper

# Synchronize Tree
pgls_tree <- keep.tip(tree, pgls_df$species_std_caper)

# Build Comparative Data Object
comp_data <- comparative.data(
  phy = pgls_tree,
  data = pgls_df,
  names.col = "species_std_caper",
  vcv = TRUE,
  warn.dropped = TRUE
)

# Fit PGLS Models
pgls_morph <- pgls(logMP ~ PC1 + PC2, data = comp_data, lambda = "ML")
pgls_eco   <- pgls(logMP ~ Habitat + Feeding_habit + Migration, data = comp_data, lambda = "ML")
pgls_glob  <- pgls(logMP ~ PC1 + PC2 + Habitat + Feeding_habit + Migration, data = comp_data, lambda = "ML")

# Extract Lambda safely
get_lambda <- function(mod) {
  as.numeric(mod$param["lambda"])
}

# Fixed Model Selection Table (No parenthesis/nesting errors)
pgls_aic_table <- data.frame(
  Model = c("Morphology (pPC1 + pPC2)", 
            "Ecology (Habitat + Niche + Migration)", 
            "Global (Morphology + Ecology)"),
  AIC = c(AIC(pgls_morph), 
          AIC(pgls_eco), 
          AIC(pgls_glob)),
  R_Squared = c(summary(pgls_morph)$r.squared, 
                summary(pgls_eco)$r.squared, 
                summary(pgls_glob)$r.squared),
  Lambda = c(get_lambda(pgls_morph), 
             get_lambda(pgls_eco), 
             get_lambda(pgls_glob))
)

cat("\n=== PGLS MODEL SELECTION TABLE ===\n")
print(pgls_aic_table)

# Write output safely (Directory verified above)
write.csv(summary(pgls_glob)$coefficients, here("output", "PGLS_Global_Coefficients.csv"))
write.csv(pgls_aic_table, here("output", "PGLS_Model_Comparison.csv"))

cat("\nSuccessfully saved PGLS model output and clean pPCA figure to 'output/' directory.\n")


summary(pgls_glob)



################################################################################
# SUPPORTING INFORMATION (SI) 
# Creates: (Diagnostic Panels)
################################################################################

library(here)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(mgcv)

# 1. Ensure SI Directory Exists
si_dir <- here("output", "SI")
if (!dir.exists(si_dir)) {
  dir.create(si_dir, recursive = TRUE)
}

# ==============================================================================
# FIGURE S1: 3-PANEL MECHANISTIC & DIAGNOSTIC PLOT
# ==============================================================================

# Panel A: Body Mass vs. pPC1 Score
p1 <- ggplot(scores_df, aes(x = Mass, y = PC1)) +
  geom_point(aes(color = logMP), size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", color = "#D62828", se = FALSE, linetype = "dashed") +
  scale_x_log10() +
  scale_color_viridis_c(option = "magma", name = expression(log[10](MP + 1))) +
  labs(
    title = "A. Body Mass vs. pPC1 Space",
    x = "Body Mass (g, log10 scale)",
    y = "Phylogenetic PC1 Score"
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11))

# Panel B: Body Mass vs. Microplastic Burden (Quadratic Fit)
p2 <- ggplot(scores_df, aes(x = log10_Mass, y = logMP)) +
  geom_point(aes(color = PC1), size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "#0D3B66", linewidth = 1) +
  scale_color_viridis_c(option = "viridis", name = "pPC1 Score") +
  labs(
    title = "B. Quadratic Accumulation Pattern",
    x = expression(log[10]*"(Body Mass, g)"),
    y = expression(log[10]*"(MP Burden + 1)")
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11))

# Panel C: pPC1 vs. Microplastic Burden (Linear Trajectory)
p3 <- ggplot(scores_df, aes(x = PC1, y = logMP)) +
  geom_point(color = "grey30", size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", color = "#D62828", se = TRUE, fill = "grey85") +
  labs(
    title = "C. Linear pPC1 Risk Trajectory",
    x = "Phylogenetic PC1 Score",
    y = expression(log[10]*"(MP Burden + 1)")
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11))

# Combine Panels via Patchwork
fig_s1_combined <- (p1 | p2 | p3) + 
  plot_annotation(
    title = "Figure S1. Diagnostic and mechanistic drivers of avian microplastic accumulation.",
    subtitle = "Comparing univariate non-linearity across body mass vs. linear trajectories in pPC space.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(face = "italic", size = 10, color = "grey30")
    )
  )

# Export Figure S1 PDF
ggsave(
  filename = file.path(si_dir, "Figure_S1_Diagnostics.pdf"),
  plot = fig_s1_combined,
  width = 15,
  height = 4.8,
  device = cairo_pdf
)

cat("Saved: Figure S1 (Figure_S1_Diagnostics.pdf)\n")

# ==============================================================================
# SI TABLES EXPORT
# ==============================================================================

# Table 1: Full pPCA Loadings Matrix & Eigenvalues
table_s1_loadings <- as.data.frame(ppca_result$L) %>%
  rownames_to_column(var = "Trait") %>%
  mutate(
    PC1_Var_Explained = percent_var[1],
    PC2_Var_Explained = percent_var[2],
    Pagels_Lambda = ppca_result$lambda
  )

write.csv(table_s1_loadings, file.path(si_dir, "Table_S1_pPCA_Loadings.csv"), row.names = FALSE)
cat("Saved: Table S1 (Table_S1_pPCA_Loadings.csv)\n")

# Table 2: Linear vs Quadratic vs GAM Model Comparison Diagnostics
model_lin  <- lm(logMP ~ log10_Mass, data = scores_df)
model_quad <- lm(logMP ~ poly(log10_Mass, 2), data = scores_df)
model_pc1  <- lm(logMP ~ PC1, data = scores_df)
gam_mod    <- gam(logMP ~ s(log10_Mass, k = 4), data = scores_df, method = "REML")

table_s2_diagnostics <- data.frame(
  Model_Specification = c("Linear (log10 Mass)", "Quadratic (log10 Mass^2)", "Linear (pPC1)", "GAM s(log10 Mass)"),
  AIC = c(AIC(model_lin), AIC(model_quad), AIC(model_pc1), AIC(gam_mod)),
  R_Squared_Adj = c(
    summary(model_lin)$adj.r.squared,
    summary(model_quad)$adj.r.squared,
    summary(model_pc1)$adj.r.squared,
    summary(gam_mod)$r.sq
  )
)

write.csv(table_s2_diagnostics, file.path(si_dir, "Table_S2_Univariate_Model_Comparison.csv"), row.names = FALSE)
cat("Saved: Table S2 (Table_S2_Univariate_Model_Comparison.csv)\n")

# Table 3: PGLS Model Selection Table
write.csv(pgls_aic_table, file.path(si_dir, "Table_S3_PGLS_Model_Selection.csv"), row.names = FALSE)
cat("Saved: Table S3 (Table_S3_PGLS_Model_Selection.csv)\n")

cat("\nAll Supporting Information outputs exported successfully to:", si_dir, "\n")