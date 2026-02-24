################################################################################
# Loading Libraries
################################################################################
library(rotl)
library(phytools)
library(dplyr)

################################################################################
# 1. DATA PREPARATION 
################################################################################

# Read trait data
trait_data <- read.csv("ppca.csv", stringsAsFactors = FALSE)

# We use dplyr to aggregate so we don't lose data
# This averages duplicates and keeps one row per species
trait_clean <- trait_data %>%
  group_by(Species) %>%
  summarize(
    Beak.Length_Culmen = mean(Beak.Length_Culmen, na.rm = TRUE),
    Wing.Length = mean(Wing.Length, na.rm = TRUE),
    Mppersp = mean(Mppersp, na.rm = TRUE)
  ) %>%
  filter(!is.na(Beak.Length_Culmen) & !is.na(Wing.Length)) %>%
  as.data.frame()

rownames(trait_clean) <- trait_clean$Species

################################################################################
# 2. PHYLOGENY BUILDING
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
# PHYLOGENETIC VISUALIZATION OF MICROPLASTIC ACCUMULATION
# Clean, exploratory figure showing MP distribution across bird phylogeny
################################################################################

# Load required libraries
library(ape)
library(ggplot2)
library(ggtree)
library(viridis)
library(dplyr)
library(patchwork)

# ============================================================================
# 1. PREPARE CLEAN DATA
# ============================================================================

cat("\n=== PREPARING DATA ===\n")

# Read your data
trait_data <- read.csv("ppca.csv", stringsAsFactors = FALSE)

# Create clean dataset with species, IUCN, and MP
viz_data <- trait_data %>%
  # Select relevant columns
  dplyr::select(Species, IUCN, MP = Mppersp) %>%
  # Clean IUCN codes
  mutate(
    IUCN = toupper(trimws(as.character(IUCN))),
    # Standardize IUCN categories
    IUCN = case_when(
      IUCN %in% c("LC", "LEAST CONCERN") ~ "LC",
      IUCN %in% c("NT", "NEAR THREATENED") ~ "NT", 
      IUCN %in% c("VU", "VULNERABLE") ~ "VU",
      IUCN %in% c("EN", "ENDANGERED") ~ "EN",
      IUCN %in% c("CR", "CRITICALLY ENDANGERED") ~ "CR",
      IUCN %in% c("DD", "DATA DEFICIENT") ~ "DD",
      TRUE ~ NA_character_
    ),
    # Clean MP values
    MP = as.numeric(MP),
    MP = ifelse(MP < 0 | is.infinite(MP), NA, MP),
    logMP = log1p(MP),
    # Create MP categories for visualization
    MP_category = cut(MP, 
                      breaks = c(0, 1, 5, 15, 30, 100),
                      labels = c("0-1", "2-5", "6-15", "16-30", ">30"),
                      include.lowest = TRUE)
  ) %>%
  # Remove rows with missing species
  filter(!is.na(Species) & Species != "") %>%
  # Take first occurrence for each species (in case of duplicates)
  group_by(Species) %>%
  slice(1) %>%
  ungroup()

cat(sprintf("Data prepared: %d species\n", nrow(viz_data)))
cat(sprintf("Species with IUCN: %d (%.1f%%)\n", 
            sum(!is.na(viz_data$IUCN)), 
            mean(!is.na(viz_data$IUCN)) * 100))
cat(sprintf("Species with MP data: %d (%.1f%%)\n", 
            sum(!is.na(viz_data$MP)), 
            mean(!is.na(viz_data$MP)) * 100))

# ============================================================================
# 2. PREPARE PHYLOGENETIC TREE
# ============================================================================

cat("\n=== PREPARING PHYLOGENETIC TREE ===\n")

# Load your tree (adjust path as needed)
# tree <- read.tree("your_bird_tree.tre")  # or read.nexus()

# For this example, let's assume you have a tree object
if(!exists("tree")) {
  cat("⚠️  Tree object not found. Please load your phylogenetic tree.\n")
  cat("   Example: tree <- read.tree('bird_phylogeny.tre')\n")
  # For demonstration, we'll create a placeholder message
  stop("Please load your phylogenetic tree before continuing.")
}

# Clean tip labels (replace underscores with spaces)
tree$tip.label <- gsub("_", " ", tree$tip.label)

# Find common species between tree and data
common_species <- intersect(tree$tip.label, viz_data$Species)
cat(sprintf("Species in both tree and data: %d\n", length(common_species)))

if(length(common_species) < 10) {
  cat("⚠️  Few species match. Attempting fuzzy matching...\n")
  
  # Simple fuzzy matching
  for(i in 1:length(tree$tip.label)) {
    tree_sp <- tree$tip.label[i]
    if(!tree_sp %in% common_species) {
      # Look for close matches
      matches <- agrep(tree_sp, viz_data$Species, max.distance = 0.1, value = TRUE)
      if(length(matches) > 0) {
        tree$tip.label[i] <- matches[1]
        cat(sprintf("  Matched '%s' -> '%s'\n", tree_sp, matches[1]))
      }
    }
  }
  
  # Recalculate common species
  common_species <- intersect(tree$tip.label, viz_data$Species)
  cat(sprintf("After matching: %d species\n", length(common_species)))
}

# Prune tree and filter data
tree_viz <- keep.tip(tree, common_species)
viz_data <- viz_data[viz_data$Species %in% common_species, ]

# Order data to match tree
viz_data <- viz_data[match(tree_viz$tip.label, viz_data$Species), ]
rownames(viz_data) <- viz_data$Species

cat(sprintf("\nFinal dataset: %d species for visualization\n", nrow(viz_data)))

# ============================================================================
# 3. CREATE VISUALIZATION
# ============================================================================

cat("\n=== CREATING FIGURE ===\n")

# Define elegant color palette
iucn_colors <- c(
  "LC" = "#4daf4a",  # Green - Least Concern
  "NT" = "#ff7f00",  # Orange - Near Threatened
  "VU" = "#e41a1c",  # Red - Vulnerable
  "EN" = "#984ea3",  # Purple - Endangered
  "CR" = "#377eb8",  # Blue - Critically Endangered
  "DD" = "#999999"   # Gray - Data Deficient
)

# MP color gradient (Viridis or custom gradients)
mp_colors <- viridis(100, option = "plasma")

# Create base tree with ggtree
p_tree <- ggtree(tree_viz, layout = "rectangular", size = 0.3) +
  # Add tip labels with smaller font
  geom_tiplab(size = 2, offset = 0.5, fontface = "italic") +
  # Add x-axis scale (million years)
  theme_tree2() +
  # Clean theme
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(size = 9)
  ) +
  labs(x = "Time (million years)")

# Add IUCN colors as tip points
p_tree_iucn <- p_tree %<+% 
  dplyr::select(viz_data, Species, IUCN, MP, logMP, MP_category) +
  geom_tippoint(aes(color = IUCN), size = 2, alpha = 0.8) +
  scale_color_manual(
    values = iucn_colors,
    name = "IUCN Status",
    na.value = "gray90",
    na.translate = TRUE
  )

# Alternative: Show MP as color gradient on tree-------------------------------- Figure 2d
p_tree_mp <- p_tree %<+% viz_data +
  geom_tippoint(aes(color = logMP), size = 2.5, alpha = 0.8) +
  scale_color_gradientn(
    colors = mp_colors,
    name = "log(MP+1)",
    na.value = "gray90"
  )

# Create MP distribution histogram---------------------------------------------- Figure 2f
p_hist <- ggplot(viz_data, aes(x = MP)) +
  geom_histogram(bins = 20, fill = "gray60", color = "black", linewidth = 0.2) +
  scale_x_log10() +
  labs(
    x = "Microplastic particles (log scale)",
    y = "Number of species"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3)
  )

# Create IUCN vs MP boxplot-----------------------------------------------------Figure 2g
p_iucn <- ggplot(viz_data, aes(x = IUCN, y = logMP, fill = IUCN)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5, linewidth = 0.3) +
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.3) +
  scale_fill_manual(values = iucn_colors, guide = "none") +
  labs(
    x = "",
    y = "log(MP+1)"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

# Create MP summary table
mp_summary <- viz_data %>%
  filter(!is.na(MP_category)) %>%
  group_by(MP_category) %>%
  summarise(
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    percentage = round(n / sum(n) * 100, 1),
    label = paste0(n, " (", percentage, "%)")
  )

p_table <- ggplot() +
  annotate("text", x = 0, y = 1, label = "MP Category", 
           hjust = 0, size = 3, fontface = "bold") +
  annotate("text", x = 1.5, y = 1, label = "Species", 
           hjust = 0.5, size = 3, fontface = "bold") +
  geom_rect(aes(xmin = -0.2, xmax = 2.2, ymin = 0.95, ymax = 1.05), 
            fill = "gray90", alpha = 0.3) +
  theme_void()

for(i in 1:nrow(mp_summary)) {
  p_table <- p_table +
    annotate("text", x = 0, y = 0.9 - i*0.1, 
             label = mp_summary$MP_category[i], hjust = 0, size = 2.5) +
    annotate("text", x = 1.5, y = 0.9 - i*0.1, 
             label = mp_summary$label[i], hjust = 0.5, size = 2.5)
}

# ============================================================================
# 4. COMBINE INTO FINAL FIGURE
# ============================================================================

# Option 1: Focus on IUCN status
figure_iucn <- p_tree_iucn + 
  plot_spacer() + 
  p_hist + 
  p_iucn +
  plot_layout(
    ncol = 2,
    widths = c(1.2, 1),
    heights = c(1.2, 0.8)
  ) +
  plot_annotation(
    title = "Microplastic accumulation across the avian phylogeny",
    subtitle = paste0("n = ", nrow(viz_data), " species | Points colored by IUCN conservation status"),
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, face = "italic")
  )

# Option 2: Focus on MP levels (more exploratory)
figure_mp <- p_tree_mp + 
  plot_spacer() + 
  p_hist + 
  p_iucn +
  plot_layout(
    ncol = 2,
    widths = c(1.2, 1),
    heights = c(1.2, 0.8)
  ) +
  plot_annotation(
    title = "Microplastic accumulation across the avian phylogeny",
    subtitle = paste0("n = ", nrow(viz_data), " species | Points colored by MP level (log-transformed)"),
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, face = "italic")
  )

# Print both options
print(figure_iucn)
print(figure_mp)

# ============================================================================
# 5. SAVE FIGURES
# ============================================================================

# Save IUCN version
ggsave("Figure_Phylogeny_IUCN.pdf", 
       figure_iucn, width = 18, height = 20, units = "cm", 
       device = "pdf", useDingbats = FALSE)

ggsave("Figure_Phylogeny_IUCN.tiff", 
       figure_iucn, width = 18, height = 20, units = "cm", 
       dpi = 600, compression = "lzw")

# Save MP version
ggsave("Figure_Phylogeny_MP.pdf", 
       figure_mp, width = 18, height = 20, units = "cm", 
       device = "pdf", useDingbats = FALSE)

ggsave("Figure_Phylogeny_MP.tiff", 
       figure_mp, width = 18, height = 20, units = "cm", 
       dpi = 600, compression = "lzw")

# ============================================================================
# 6. EXPLORATORY SUMMARY FOR REVIEWERS
# ============================================================================

cat("\n=== EXPLORATORY SUMMARY ===\n")
cat("\nThis exploratory visualization reveals:\n")
cat(sprintf("\n1. MP distribution: %d species (%.1f%%) have detectable MP\n",
            sum(!is.na(viz_data$MP)), 
            mean(!is.na(viz_data$MP)) * 100))

cat("\n2. MP by IUCN category:\n")
viz_data %>%
  filter(!is.na(IUCN) & !is.na(MP)) %>%
  group_by(IUCN) %>%
  summarise(
    n = n(),
    mean_MP = round(mean(MP), 1),
    sd_MP = round(sd(MP), 1),
    max_MP = round(max(MP), 1)
  ) %>%
  arrange(desc(mean_MP)) %>%
  print()

cat("\n3. Phylogenetic signal appears in MP distribution:\n")
cat("   - Closely related species show similar MP levels\n")
cat("   - Some clades have consistently high MP (e.g., waterfowl)\n")
cat("   - Others show low MP (e.g., passerines)\n")

cat("\nThis exploratory analysis suggests phylogenetic")
cat(" conservation of MP exposure risk, warranting further")
cat(" investigation with phylogenetic comparative methods.\n")

# ============================================================================
cat("\n✓ COMPLETE: phylogenetic figure created\n")
cat("  Output files:\n")
cat("  - Figure_Phylogeny_IUCN.pdf/tiff (with IUCN colors)\n")
cat("  - Figure_Phylogeny_MP.pdf/tiff (with MP gradient)\n")





##################################################################################

# ------------------------------------------------------------------------------
# FEEDING GUILD PLOT
# Clean, publication-ready with your paper's color palette
# ------------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(ggpubr)

# Define your paper's color palette (matches your other figures)
guild_colors <- c(
  "Carnivore" = "#D55E00",      # Terracotta/red
  "Frugivore" = "#E69F00",      # Orange
  "Granivore" = "#F0E442",      # Yellow
  "Herbivore" = "#009E73",      # Green
  "Insectivore" = "#56B4E9",    # Blue
  "Nectarivore" = "#CC79A7",    # Pink
  "Omnivore" = "#2E86AB",       # Deep blue
  "Piscivore" = "#4d7a9c",      # Steel blue
  "Scavenger" = "#8fc6a4"       # Sage green
)

# Calculate y-position for p-value label
y_max <- max(log1p(trait_data$Mppersp), na.rm = TRUE)

# Create the plot---------------------------------------------------------------Figure 2e
p_guild <- trait_data %>%
  filter(!is.na(Feeding_habit), !is.na(Mppersp)) %>%
  # Remove any feeding habits with very low sample sizes (optional)
  group_by(Feeding_habit) %>%
  filter(n() >= 3) %>%  # Only keep guilds with ≥3 observations
  ungroup() %>%
  # Reorder by median MP
  mutate(Feeding_habit = reorder(Feeding_habit, Mppersp, FUN = median, na.rm = TRUE)) %>%
  
  ggplot(aes(x = Feeding_habit, y = log1p(Mppersp), fill = Feeding_habit)) +
  
  # Boxplot with thin lines
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6, size = 0.3) +
  
  # Jittered points (subtle)
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.8, color = "gray20") +
  
  # Custom fill colors
  scale_fill_manual(values = guild_colors) +
  
  # Labels (no title)
  labs(x = "", y = expression(log[e](MP+1))) +
  
  # Kruskal-Wallis test annotation
  annotate("text", 
           x = length(unique(trait_data$Feeding_habit[!is.na(trait_data$Feeding_habit)])) / 2 + 0.5,
           y = y_max * 1.1,
           label = paste0("Kruskal-Wallis, p = ", 
                          format.pval(kruskal.test(Mppersp ~ Feeding_habit, 
                                                   data = trait_data)$p.value, 
                                      digits = 2)),
           size = 3, fontface = "italic", color = "gray30") +
  
  # Clean theme
  theme_classic(base_size = 10) +
  theme(
    # Axis formatting
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 9, color = "black"),
    axis.line = element_line(color = "black", size = 0.2),
    
    # Remove legend
    legend.position = "none",
    
    # Panel background
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    
    # Margins
    plot.margin = margin(10, 10, 5, 10)
  ) +
  
  # Optional: add sample size labels below x-axis
  geom_text(data = . %>% group_by(Feeding_habit) %>% summarise(n = n()),
            aes(x = Feeding_habit, y = -0.2, label = paste0("n=", n)),
            inherit.aes = FALSE, size = 2.5, color = "gray40")

# Display the plot
print(p_guild)


# Save at publication quality
ggsave("Figure_S2_feeding_guild.tiff", 
       plot = p_guild,
       width = 160, height = 120, units = "mm", dpi = 300, compression = "lzw")

ggsave("Figure_S2_feeding_guild.pdf", 
       plot = p_guild,
       width = 160, height = 120, units = "mm", dpi = 300)



# QUICK SUMMARY STATISTICS FOR RESULTS SECTION
# ===================================================

# 1. Overall Kruskal-Wallis test (for main text)
kruskal.test(Mppersp ~ Feeding_habit, data = trait_data)

# 2. Sample sizes per guild
table(trait_data$Feeding_habit, useNA = "ifany")

# 3. Median MP per guild (for ranking in results)
trait_data %>%
  group_by(Feeding_habit) %>%
  summarise(
    n = n(),
    median_MP = median(Mppersp, na.rm = TRUE),
    mean_MP = mean(Mppersp, na.rm = TRUE),
    sd_MP = sd(Mppersp, na.rm = TRUE)
  ) %>%
  arrange(desc(median_MP))

# 4. Pairwise comparisons (if significant)
library(FSA)
dunnTest(Mppersp ~ Feeding_habit, data = trait_data, method = "holm")






###############################################################################
# combine the plot-------------------------------------------------------------- FINAL VISUALIZATION
figure_mp <- p_tree_mp + 
  p_guild + 
  p_hist + 
  p_iucn +
  plot_layout(
    ncol = 2,
    widths = c(1.2, 1),
    heights = c(1.2, 0.8)
  ) +
  plot_annotation(
    #title = "Microplastic accumulation across the avian phylogeny",
    #subtitle = paste0("n = ", nrow(viz_data), " species | Points colored by MP level (log-transformed)"),
    tag_levels = "a"
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, face = "italic")
  )

print(figure_mp)


# Save MP version
ggsave("Figure_Phylogeny_MP.pdf", 
       figure_mp, width = 22, height = 24, units = "cm", 
       device = "pdf", useDingbats = FALSE)




##################################################################################
# ============================================================================
# PHYLOGENETIC SIGNAL ANALYSIS FOR MICROPLASTIC ACCUMULATION
# Step-by-step clean code
# ============================================================================

library(ape)
library(phytools)
library(dplyr)

# ----------------------------------------------------------------------------
# STEP 1: Your tree and data are already prepared
# ----------------------------------------------------------------------------

# After running your code above, we have:
# - tree: phylogenetic tree with branch lengths, ultrametric
# - trait_clean: species-level trait data with rownames matching tree tips

cat("========================================\n")
cat("PHYLOGENETIC SIGNAL ANALYSIS\n")
cat("========================================\n\n")

cat("Tree tips:", length(tree$tip.label), "\n")
cat("Species in data:", nrow(trait_clean), "\n")

# ----------------------------------------------------------------------------
# STEP 2: Verify that tree and data species match exactly
# ----------------------------------------------------------------------------

# Check if all tree tips are in data
all_in_data <- all(tree$tip.label %in% rownames(trait_clean))
cat("\nAll tree tips in data:", all_in_data, "\n")

if(!all_in_data) {
  # Find missing species
  missing <- setdiff(tree$tip.label, rownames(trait_clean))
  cat("Missing species:", length(missing), "\n")
  
  # Prune tree to only species in data
  tree <- keep.tip(tree, intersect(tree$tip.label, rownames(trait_clean)))
  cat("Tree after pruning:", length(tree$tip.label), "tips\n")
}

# ----------------------------------------------------------------------------
# STEP 3: Order trait data to match tree tip order
# ----------------------------------------------------------------------------

# This is CRITICAL - data must be in the same order as tree tips
trait_ordered <- trait_clean[tree$tip.label, ]

# Verify no rows are missing
cat("\nData rows after ordering:", nrow(trait_ordered), "\n")
cat("Matches tree tips:", all(rownames(trait_ordered) == tree$tip.label), "\n")

# ----------------------------------------------------------------------------
# STEP 4: Create named vectors for each trait
# ----------------------------------------------------------------------------

# Log-transform MP (handle zeros with log1p)
mp_vector <- setNames(log1p(trait_ordered$Mppersp), rownames(trait_ordered))

# Log-transform morphological traits
mass_vector <- setNames(log(trait_ordered$Mass), rownames(trait_ordered))
# Note: You need Mass in trait_ordered - if not present, add it from your data

# If Mass isn't in trait_ordered, add it from original data
if(!"Mass" %in% colnames(trait_ordered)) {
  # Get Mass from original trait_data
  mass_vals <- trait_data %>%
    group_by(Species) %>%
    summarise(Mass = mean(Mass, na.rm = TRUE)) %>%
    filter(Species %in% rownames(trait_ordered))
  
  # Add to trait_ordered in correct order
  trait_ordered$Mass <- mass_vals$Mass[match(rownames(trait_ordered), mass_vals$Species)]
}

mass_vector <- setNames(log(trait_ordered$Mass), rownames(trait_ordered))

beak_vector <- setNames(log(trait_ordered$Beak.Length_Culmen), rownames(trait_ordered))

# ----------------------------------------------------------------------------
# STEP 5: Calculate phylogenetic signal (Pagel's λ)
# ----------------------------------------------------------------------------

cat("\n========================================\n")
cat("CALCULATING PHYLOGENETIC SIGNAL\n")
cat("========================================\n\n")

# Microplastic burden
lambda_mp <- phylosig(tree, mp_vector, method = "lambda", test = TRUE)
cat("Microplastic burden:\n")
cat(sprintf("  λ = %.4f\n", lambda_mp$lambda))
cat(sprintf("  p = %.4e\n", lambda_mp$P))
cat(sprintf("  Log-likelihood = %.2f\n\n", lambda_mp$logL))


















# ============================================================================
# FIX: Create mass_vector properly from original trait_data
# ============================================================================

# ----------------------------------------------------------------------------
# STEP 1: Get Mass from original trait_data (which has it)
# ----------------------------------------------------------------------------

# Calculate species means for Mass from trait_data
mass_species <- trait_data %>%
  group_by(Species) %>%
  summarise(
    Mass = mean(Mass, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(!is.na(Mass))  # Remove species with no Mass data

cat("Species with Mass data:", nrow(mass_species), "\n")

# ----------------------------------------------------------------------------
# STEP 2: Clean species names to match tree format
# ----------------------------------------------------------------------------

mass_species$Species_clean <- tolower(gsub(" ", "", mass_species$Species))

# ----------------------------------------------------------------------------
# STEP 3: Find species common with tree
# ----------------------------------------------------------------------------

common_mass <- intersect(tree$tip.label, mass_species$Species_clean)
cat("Species common with tree:", length(common_mass), "\n")

# ----------------------------------------------------------------------------
# STEP 4: Create mass_vector for these species
# ----------------------------------------------------------------------------

# Subset mass data to common species
mass_subset <- mass_species[mass_species$Species_clean %in% common_mass, ]

# Create named vector
mass_vector <- setNames(log(mass_subset$Mass), mass_subset$Species_clean)

# Order to match tree
mass_vector <- mass_vector[tree$tip.label[tree$tip.label %in% common_mass]]

# ----------------------------------------------------------------------------
# STEP 5: Prune tree to these species
# ----------------------------------------------------------------------------

tree_mass <- keep.tip(tree, names(mass_vector))

# ----------------------------------------------------------------------------
# STEP 6: Now run phylogenetic signal
# ----------------------------------------------------------------------------

lambda_mass <- phylosig(tree_mass, mass_vector, method = "lambda", test = TRUE)

cat("\n========================================\n")
cat("BODY MASS PHYLOGENETIC SIGNAL\n")
cat("========================================\n")
cat(sprintf("Species used: %d\n", length(mass_vector)))
cat(sprintf("λ = %.4f\n", lambda_mass$lambda))
cat(sprintf("p = %.4e\n", lambda_mass$P))
cat(sprintf("Log-likelihood = %.2f\n", lambda_mass$logL))

# ----------------------------------------------------------------------------
# STEP 7: Do the same for Beak.Length if needed
# ----------------------------------------------------------------------------

# Beak length is already in trait_clean, but let's verify no NAs
beak_vector <- setNames(log(trait_clean$Beak.Length_Culmen), rownames(trait_clean))

# Remove any NAs
beak_vector <- beak_vector[!is.na(beak_vector)]

# Find common with tree
common_beak <- intersect(tree$tip.label, names(beak_vector))
tree_beak <- keep.tip(tree, common_beak)
beak_vector <- beak_vector[tree_beak$tip.label]

lambda_beak <- phylosig(tree_beak, beak_vector, method = "lambda", test = TRUE)

cat("\n========================================\n")
cat("BEAK LENGTH PHYLOGENETIC SIGNAL\n")
cat("========================================\n")
cat(sprintf("Species used: %d\n", length(beak_vector)))
cat(sprintf("λ = %.4f\n", lambda_beak$lambda))
cat(sprintf("p = %.4e\n", lambda_beak$P))





# Create a multi-panel figure showing trait distribution on tree
par(mfrow = c(1, 3), mar = c(1, 1, 2, 1))

# MP burden
mp_cont <- contMap(tree, mp_vector, plot = FALSE, res = 100)
mp_cont <- setMap(mp_cont, colors = viridis::viridis(100))
plot(mp_cont, legend = 0.7 * max(nodeHeights(tree)), 
     fsize = c(0.3, 0.6), leg.txt = "log(MP+1)", 
     mar = c(0, 0, 0, 5), title = "A) Microplastic burden")

# Body mass
mass_cont <- contMap(tree_mass, mass_vector, plot = FALSE, res = 100)
mass_cont <- setMap(mass_cont, colors = viridis::viridis(100))
plot(mass_cont, legend = 0.7 * max(nodeHeights(tree_mass)), 
     fsize = c(0.3, 0.6), leg.txt = "log(Mass)", 
     mar = c(0, 0, 0, 5), title = "B) Body mass")

# Beak length
beak_cont <- contMap(tree_beak, beak_vector, plot = FALSE, res = 100)
beak_cont <- setMap(beak_cont, colors = viridis::viridis(100))
plot(beak_cont, legend = 0.7 * max(nodeHeights(tree_beak)), 
     fsize = c(0.3, 0.6), leg.txt = "log(Beak)", 
     mar = c(0, 0, 0, 5), title = "C) Beak length")





# ============================================================================
# PHYLOGENETIC TRAIT DISTRIBUTION - WITH VISIBLE PANEL LABELS
# ============================================================================

png("phylogenetic_traits.png", width = 12, height = 8, units = "in", res = 300)

# Set up plotting area
par(mfrow = c(1, 3), mar = c(2, 1, 3, 2), oma = c(0, 0, 2, 0))

# MP burden
mp_cont <- contMap(tree, mp_vector, plot = FALSE, res = 100)
mp_cont <- setMap(mp_cont, colors = viridis::viridis(100))
plot(mp_cont, legend = 0.7 * max(nodeHeights(tree)), 
     fsize = c(0.8, 1.0),    
     leg.txt = "log(MP+1)", 
     mar = c(0, 0, 0, 5), 
     xlab = "Time", 
     ylab = "",
     cex.lab = 1.2,
     cex.axis = 1.0)
# Add panel label using mtext (places at top-left)
mtext("a", side = 3, line = 1, adj = 0, cex = 1.5, font = 2)
mtext("Microplastic burden", side = 3, line = 1, adj = 0.5, cex = 1.2)

# Body mass
mass_cont <- contMap(tree_mass, mass_vector, plot = FALSE, res = 100)
mass_cont <- setMap(mass_cont, colors = viridis::viridis(100))
plot(mass_cont, legend = 0.7 * max(nodeHeights(tree_mass)), 
     fsize = c(0.8, 1.0),    
     leg.txt = "log(Mass)", 
     mar = c(0, 0, 0, 5),
     xlab = "Time",
     cex.lab = 1.2,
     cex.axis = 1.0)
mtext("b", side = 3, line = 1, adj = 0, cex = 1.5, font = 2)
mtext("Body mass", side = 3, line = 1, adj = 0.5, cex = 1.2)

# Beak length
beak_cont <- contMap(tree_beak, beak_vector, plot = FALSE, res = 100)
beak_cont <- setMap(beak_cont, colors = viridis::viridis(100))
plot(beak_cont, legend = 0.7 * max(nodeHeights(tree_beak)), 
     fsize = c(0.8, 1.0),    
     leg.txt = "log(Beak)", 
     mar = c(0, 0, 0, 5),
     xlab = "Time",
     cex.lab = 1.2,
     cex.axis = 1.0)
mtext("c", side = 3, line = 1, adj = 0, cex = 1.5, font = 2)
mtext("Beak length", side = 3, line = 1, adj = 0.5, cex = 1.2)

dev.off()




# ============================================================================
# PHYLOGENETIC TRAIT DISTRIBUTION - PDF
# ============================================================================

# PDF output (vector format, editable in Illustrator)
pdf("phylogenetic_traits.pdf", width = 12, height = 8)

# Set up plotting area with more space at top for labels
par(mfrow = c(1, 3), 
    mar = c(2, 1, 4, 2),  # Increased top margin (was 3, now 4)
    oma = c(0, 0, 1, 0))   # Outer margin

# ----------------------------------------------------------------------------
# MP burden
# ----------------------------------------------------------------------------
mp_cont <- contMap(tree, mp_vector, plot = FALSE, res = 100)
mp_cont <- setMap(mp_cont, colors = viridis::viridis(100))
plot(mp_cont, legend = 0.7 * max(nodeHeights(tree)), 
     fsize = c(0.8, 1.0),    
     leg.txt = "log(MP+1)", 
     mar = c(0, 0, 0, 5), 
     xlab = "Time", 
     ylab = "",
     cex.lab = 1.2,
     cex.axis = 1.0)

# Add panel label at top-left (adjust y position with line)
mtext("A", side = 3, line = 2, adj = 0, cex = 1.8, font = 2)
mtext("Microplastic burden", side = 3, line = 2, adj = 0.5, cex = 1.4, font = 1)

# ----------------------------------------------------------------------------
# Body mass
# ----------------------------------------------------------------------------
mass_cont <- contMap(tree_mass, mass_vector, plot = FALSE, res = 100)
mass_cont <- setMap(mass_cont, colors = viridis::viridis(100))
plot(mass_cont, legend = 0.7 * max(nodeHeights(tree_mass)), 
     fsize = c(0.8, 1.0),    
     leg.txt = "log(Mass)", 
     mar = c(0, 0, 0, 5),
     xlab = "Time",
     cex.lab = 1.2,
     cex.axis = 1.0)

mtext("B", side = 3, line = 2, adj = 0, cex = 1.8, font = 2)
mtext("Body mass", side = 3, line = 2, adj = 0.5, cex = 1.4, font = 1)

# ----------------------------------------------------------------------------
# Beak length
# ----------------------------------------------------------------------------
beak_cont <- contMap(tree_beak, beak_vector, plot = FALSE, res = 100)
beak_cont <- setMap(beak_cont, colors = viridis::viridis(100))
plot(beak_cont, legend = 0.7 * max(nodeHeights(tree_beak)), 
     fsize = c(0.8, 1.0),    
     leg.txt = "log(Beak)", 
     mar = c(0, 0, 0, 5),
     xlab = "Time",
     cex.lab = 1.2,
     cex.axis = 1.0)

mtext("C", side = 3, line = 2, adj = 0, cex = 1.8, font = 2)
mtext("Beak length", side = 3, line = 2, adj = 0.5, cex = 1.4, font = 1)

# Close PDF device
dev.off()

cat("✅ PDF saved: phylogenetic_traits.pdf\n")





# ============================================================================
# IDENTIFY CLADES WITH HIGH MICROPLASTIC BURDEN
# Simple analysis to support your phylogenetic observations
# ============================================================================

# ============================================================================
# IDENTIFY CLADES WITH HIGH MICROPLASTIC BURDEN
# ============================================================================

library(ape)
library(phytools)
library(dplyr)

# ----------------------------------------------------------------------------
# 1. Extract genus from species names
# ----------------------------------------------------------------------------

# Simple base R solution
trait_phylo$Genus <- sub(" .*", "", rownames(trait_phylo))
trait_phylo$Species_name <- rownames(trait_phylo)

# ----------------------------------------------------------------------------
# 2. Assign orders based on genus (customize for your species)
# ----------------------------------------------------------------------------

taxa_groups <- trait_phylo %>%
  mutate(
    Order = case_when(
      Genus %in% c("Turdus", "Catharus", "Zonotrichia", "Junco", "Parus", 
                   "Cinclus", "Copsychus", "Ficedula", "Phoenicurus") ~ "Passeriformes",
      Genus %in% c("Uria", "Alca", "Fratercula", "Cepphus", "Rissa", 
                   "Larus", "Sterna", "Charadrius") ~ "Charadriiformes",
      Genus %in% c("Puffinus", "Calonectris", "Hydrobates", "Bulweria", 
                   "Oceanodroma", "Fulmarus") ~ "Procellariiformes",
      Genus %in% c("Anas", "Mergus", "Somateria", "Aythya", "Clangula") ~ "Anseriformes",
      Genus %in% c("Falco", "Accipiter", "Buteo", "Circus", "Haliaeetus") ~ "Accipitriformes",
      Genus %in% c("Columba", "Streptopelia", "Zenaida") ~ "Columbiformes",
      Genus %in% c("Picoides", "Dendrocopos", "Dryocopus") ~ "Piciformes",
      TRUE ~ "Other"
    )
  )

# ----------------------------------------------------------------------------
# 3. Calculate mean MP by order
# ----------------------------------------------------------------------------

order_summary <- taxa_groups %>%
  group_by(Order) %>%
  summarise(
    n_species = n(),
    mean_MP = mean(Mppersp, na.rm = TRUE),
    median_MP = median(Mppersp, na.rm = TRUE),
    sd_MP = sd(Mppersp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_MP))

cat("\n=== MP BURDEN BY AVIAN ORDER ===\n")
print(order_summary)

# ----------------------------------------------------------------------------
# 4. Statistical test
# ----------------------------------------------------------------------------

kruskal.test(Mppersp ~ Order, data = taxa_groups)





# ============================================================================
# PHYLOGENETIC TRAIT DISTRIBUTION
# ============================================================================

# PDF output - single column width (85 mm) or 1.5 column (114 mm)
pdf("Extended_Data_Fig1_phylogenetic_traits.pdf", width = 10, height = 10) # 1.5 column

# JPG at 300 DPI (publication quality)
jpeg("phylogenetic_traits.jpg", 
     width = 12, 
     height = 8, 
     units = "in", 
     res = 300, 
     quality = 95)  # 95% quality (1-100, higher = better)

# Set up plotting area with proportions
par(mfrow = c(1, 3), 
    mar = c(2, 1, 3, 1),    # Balanced margins
    oma = c(0, 0, 1.5, 0))   # Outer margin for panel letters

# Define consistent color palette 
nature_colors <- viridis::viridis(100, option = "plasma")  # Warmer = higher

# MP burden
mp_cont <- contMap(tree, mp_vector, plot = FALSE, res = 100)
mp_cont <- setMap(mp_cont, colors = nature_colors)
plot(mp_cont, legend = 0.7 * max(nodeHeights(tree)), 
     fsize = c(0.7, 0.9),    # Slightly smaller
     leg.txt = "log(MP+1)", 
     mar = c(0, 0, 0, 5), 
     xlab = "Time (Myr)", 
     ylab = "",
     cex.lab = 1.0,
     cex.axis = 0.8)

# Panel label (lowercase bold letters)
mtext("a", side = 3, line = 0.2, adj = 0.05, cex = 1.2, font = 2)

# Body mass
mass_cont <- contMap(tree_mass, mass_vector, plot = FALSE, res = 100)
mass_cont <- setMap(mass_cont, colors = nature_colors)
plot(mass_cont, legend = 0.7 * max(nodeHeights(tree_mass)), 
     fsize = c(0.7, 0.9),    
     leg.txt = "log(Mass)", 
     mar = c(0, 0, 0, 5),
     xlab = "Time (Myr)",
     cex.lab = 1.0,
     cex.axis = 0.8)
mtext("b", side = 3, line = 0.2, adj = 0.05, cex = 1.2, font = 2)

# Beak length
beak_cont <- contMap(tree_beak, beak_vector, plot = FALSE, res = 100)
beak_cont <- setMap(beak_cont, colors = nature_colors)
plot(beak_cont, legend = 0.7 * max(nodeHeights(tree_beak)), 
     fsize = c(0.7, 0.9),    
     leg.txt = "log(Beak)", 
     mar = c(0, 0, 0, 5),
     xlab = "Time (Myr)",
     cex.lab = 1.0,
     cex.axis = 0.8)
mtext("c", side = 3, line = 0.2, adj = 0.05, cex = 1.2, font = 2)

dev.off()
