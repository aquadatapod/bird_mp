################################################################################
# FIGURE 1: GEOGRAPHIC, PHYLOGENETIC & FUNCTIONAL DIVERSITY
#
# Generates Panels a-i :
# 1. Latitudinal pseudo-replication (site-level means + mixed-model callout)
# 2. Species-level functional uniqueness (FUni) and originality (FOrig)
# 3. Canonical species-level IUCN & Feeding Guild Kruskal-Wallis tests
################################################################################

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(sf)
library(rnaturalearth)
library(viridis)
library(patchwork)
library(FD)
library(cluster)
library(ape)
library(rotl)
library(phytools)
library(ggtree)

# ==============================================================================
# 0. THEME DEFINITION
# ==============================================================================
theme_nature <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.3),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.text = element_text(color = "black", size = base_size - 1),
      axis.title = element_text(color = "black", size = base_size),
      plot.tag = element_text(size = 11, face = "bold"),
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text = element_text(size = base_size - 2),
      plot.margin = margin(2, 2, 2, 2)
    )
}

# ==============================================================================
# 1. DATA PREPARATION (RAW & CANONICAL SPECIES LEVEL)
# ==============================================================================
cat("\n=== 1. LOADING AND PREPARING DATA ===\n")

row_level <- read.csv(here("data", "processed", "row_level_02.csv"),
                      stringsAsFactors = FALSE, na.strings = c("", "NA"))

# Clean known Order typos
row_level <- row_level %>%
  mutate(Order = case_when(
    Order == "Charadriformes" ~ "Charadriiformes",
    Order == "Phalacrocorcidae" ~ "Suliformes",
    TRUE ~ Order
  ))

# Raw records dataset for geography
df_clean <- row_level %>%
  filter(!is.na(Mppersp), !is.na(Latitude), !is.na(Longitude)) %>%
  mutate(logMP = log1p(Mppersp))

# Canonical Species-Level Table (Avoids Pseudo-Replication)
species_phylo <- row_level %>%
  group_by(Species) %>%
  summarise(
    n_records           = n(),
    Mass                = median(Mass, na.rm = TRUE),
    Beak.Length_Culmen  = median(Beak.Length_Culmen, na.rm = TRUE),
    Beak.Width          = median(Beak.Width, na.rm = TRUE),
    Wing.Length         = median(Wing.Length, na.rm = TRUE),
    MP                  = median(Mppersp, na.rm = TRUE),
    Latitude            = median(Latitude, na.rm = TRUE),
    Longitude           = median(Longitude, na.rm = TRUE),
    Order               = names(sort(table(Order), decreasing = TRUE))[1],
    IUCN                = names(sort(table(IUCN), decreasing = TRUE))[1],
    Feeding_habit       = names(sort(table(Feeding_habit), decreasing = TRUE))[1],
    Habitat             = names(sort(table(Habitat), decreasing = TRUE))[1],
    Migration           = names(sort(table(Migration), decreasing = TRUE))[1],
    .groups             = "drop"
  ) %>%
  filter(!is.na(MP)) %>%
  mutate(
    logMP = log1p(MP),
    IUCN = factor(IUCN, levels = c("LC", "NT", "VU", "EN", "CR", "DD")),
    Species_std = tolower(gsub("[[:punct:][:space:]]", "", Species))
  )

cat(sprintf("Record-level observations: %d\n", nrow(df_clean)))
cat(sprintf("Canonical species-level records: %d\n", nrow(species_phylo)))




# ==============================================================================
# 2. PANELS A, B, C, D — GEOGRAPHIC & LATITUDINAL PATTERNS
# ==============================================================================
cat("\n=== 2. GENERATING PANELS A, B, C, D ===\n")

# Panel A: Conceptual diagram placeholder (assembled via image or layout)
# Panel B: World Map
world <- ne_countries(scale = "medium", returnclass = "sf")

p_b_map <- ggplot() +
  # Base map - very subtle
  geom_sf(data = world, fill = "#f0f0f0", color = "#cccccc", size = 0.2) +
  # Points
  geom_point(data = df_clean, 
             aes(x = Longitude, y = Latitude, 
                 color = logMP, size = Mppersp),
             alpha = 0.6) +
  # Scales
  scale_color_viridis_c(option = "plasma", 
                        name = "log(MP+1)",
                        guide = guide_colorbar(barwidth = 0.8, barheight = 4)) +
  scale_size_continuous(name = "MP particles",
                        range = c(0.8, 5),
                        breaks = c(1, 10, 50, 100, 500),
                        guide = guide_legend(override.aes = list(alpha = 1))) +
  # Labels - minimal
  labs(tag = "a") +
  # Theme - clean
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.tag = element_text(size = 14, face = "bold"),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 10),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 90), expand = FALSE)

# Panel C: Continent Boxplot
continent_order <- df_clean %>%
  group_by(Continent) %>%
  summarise(med = median(Mppersp, na.rm = TRUE)) %>%
  arrange(desc(med)) %>% pull(Continent)

df_clean$Continent <- factor(df_clean$Continent, levels = continent_order)
continent_sizes <- df_clean %>% group_by(Continent) %>% summarise(n = n(), .groups = 'drop')

# pseudo-log scale for raw Mppersp
library(scales)

p_c_continent <- ggplot(df_clean, aes(x = Continent, y = Mppersp)) +
  geom_boxplot(fill = "#d9d9d9", outlier.shape = NA, width = 0.5, linewidth = 0.3, color = "#4d4d4d") +
  geom_jitter(width = 0.15, alpha = 0.2, size = 0.6, color = "#4d4d4d") +
  scale_y_continuous(
    trans = "log1p",
    breaks = c(0, 1, 10, 100, 1000),
    labels = c("0", "1", "10", "100", "1000")
  ) +
  geom_text(data = continent_sizes, aes(x = Continent, y = 0, label = paste0("N=", n)),
            size = 2, color = "gray30", vjust = 1.5, inherit.aes = FALSE) +
  labs(x = "", y = "MP particles / individual", tag = "c") +
  theme_nature()

ggsave(here("output", "continet.pdf"), p_c_continent, width = 14, height = 5, units = "cm", device = "pdf")

# Panel D: Latitudinal Pattern Corrected (Site-Level Means + Mixed Model Annotation)
site_level <- df_clean %>%
  mutate(site_id = paste(round(Latitude, 4), round(Longitude, 4), sep = "_")) %>%
  group_by(site_id) %>%
  summarise(Latitude = first(Latitude), Mppersp_median = median(Mppersp), n_records = n(), .groups = "drop") %>%
  mutate(oversampled = n_records >= 5)

p_d_lat <- ggplot() +
  geom_point(data = df_clean, aes(x = Latitude, y = Mppersp), color = "grey70", alpha = 0.35, size = 0.8) +
  geom_point(data = site_level, aes(x = Latitude, y = Mppersp_median, size = n_records, color = oversampled), alpha = 0.85) +
  scale_color_manual(values = c("FALSE" = "#2E6F95", "TRUE" = "#D7263D"), guide = "none") +
  scale_size_continuous(range = c(1.2, 4), guide = "none") +
  geom_vline(xintercept = c(23.5, -23.5), linetype = "dotted", color = "#4d4d4d", linewidth = 0.3) +
  annotate("text", x = 0, y = 300, label = "Equator", size = 2.3, color = "#4d4d4d") +
  annotate("label", x = -55, y = 600, label = "Mixed model:\nNo significant latitudinal effect\n(\u03b2 = -0.004, p = 0.606)",
           size = 2.1, hjust = 0, fill = "white", label.size = 0.15) +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000), labels = c("0.1", "1", "10", "100", "1000")) +
  labs(x = "Latitude", y = "MP particles", tag = "d") +
  theme_nature()

ggsave(here("output", "latitude.pdf"), p_d_lat, width = 14, height = 5, units = "cm", device = "pdf")

# ==============================================================================
# 3. PANELS E, F, G — PHYLOGENY, GUILD & IUCN
# ==============================================================================
cat("\n=== 3. GENERATING PANELS E, F, G ===\n")

# Taxonomic Match & Tree Generation
taxa_matches <- tnrs_match_names(names = species_phylo$Species)
taxa_found <- taxa_matches[!is.na(taxa_matches$ott_id), ]
tree <- tol_induced_subtree(ott_ids = taxa_found$ott_id)
tree_ott_ids <- gsub(".*_ott", "", tree$tip.label)
tree$tip.label <- taxa_found$search_string[match(tree_ott_ids, as.character(taxa_found$ott_id))]
tree$tip.label <- tolower(gsub("[[:punct:][:space:]]", "", tree$tip.label))

species_phylo_tree <- species_phylo %>% filter(Species_std %in% tree$tip.label)
tree_viz <- keep.tip(tree, species_phylo_tree$Species_std)
tree_viz <- compute.brlen(tree_viz, method = "Grafen")

# Panel E: Tree Plot
plot_data <- species_phylo_tree %>% dplyr::select(Species_std, logMP)
p_e_tree <- ggtree(tree_viz, layout = "rectangular", size = 0.25) %<+% plot_data +
  geom_tippoint(aes(color = logMP), size = 1.5, alpha = 0.85) +
  scale_color_viridis_c(option = "plasma", name = "log(MP+1)") +
  theme_tree2() +
  labs(tag = "e")

###############################################################################
# ==============================================================================
# 3. PANELS E, F, G — PHYLOGENY, GUILD & IUCN (Cleaned Tree Code)
# ==============================================================================
cat("\n=== 3. GENERATING PANELS E, F, G ===\n")

# Taxonomic Match & Tree Generation
taxa_matches <- tnrs_match_names(names = species_phylo$Species)
taxa_found <- taxa_matches[!is.na(taxa_matches$ott_id), ]
tree <- tol_induced_subtree(ott_ids = taxa_found$ott_id)
tree_ott_ids <- gsub(".*_ott", "", tree$tip.label)
tree$tip.label <- taxa_found$search_string[match(tree_ott_ids, as.character(taxa_found$ott_id))]
tree$tip.label <- tolower(gsub("[[:punct:][:space:]]", "", tree$tip.label))
species_phylo_tree <- species_phylo %>% filter(Species_std %in% tree$tip.label)
tree_viz <- keep.tip(tree, species_phylo_tree$Species_std)
tree_viz <- compute.brlen(tree_viz, method = "Grafen")

# Panel E: Simplified Tree Plot (Fixed ggtree parameters)
plot_data <- species_phylo_tree %>% dplyr::select(Species_std, logMP)

p_e_tree <- ggtree(tree_viz, layout = "rectangular") %<+% plot_data +
  geom_tree(linewidth = 0.2) +
  geom_tippoint(aes(color = logMP), size = 1.2, alpha = 0.9) +
  scale_color_viridis_c(option = "plasma", name = "log(MP+1)") +
  theme_tree2() +
  labs(tag = "e", subtitle = "Phylogenetic distribution of MP load") +
  theme(
    legend.position = "right",
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.2, "cm"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8, face = "bold"),
    plot.subtitle = element_text(size = 8, face = "italic")
  )

print(p_e_tree)
ggsave(here("output", "p_tree.pdf"), p_e_tree, width = 5, height = 7.5, units = "cm", device = "pdf")

##############################################################################################

# Panel F: Feeding Guild
guild_data <- species_phylo %>%
  filter(!is.na(Feeding_habit)) %>%
  group_by(Feeding_habit) %>% filter(n() >= 3) %>% ungroup() %>%
  mutate(Feeding_habit = reorder(Feeding_habit, logMP, FUN = median))

kw_guild <- kruskal.test(logMP ~ Feeding_habit, data = guild_data)

p_f_guild <- ggplot(guild_data, aes(x = Feeding_habit, y = logMP)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, linewidth = 0.3, fill = "#d9d9d9") +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1, color = "gray20") +
  scale_fill_viridis_d(option = "turbo", guide = "none") +
  labs(x = "", y = expression(log[e](MP+1)), tag = "f",
       subtitle = sprintf("Kruskal-Wallis p = %.2e", kw_guild$p.value)) +
  theme_nature() +
  theme(plot.subtitle = element_text(size = 6.5, face = "italic"))

ggsave(here("output", "guild.pdf"), p_f_guild, width = 12, height = 7.5, units = "cm", device = "pdf")
ggsave(here("output", "funi.pdf"), p_h_funi, width = 12, height = 7.5, units = "cm", device = "pdf")


# Panel G: IUCN Status
iucn_data <- species_phylo %>%
  filter(!is.na(IUCN)) %>%
  group_by(IUCN) %>% filter(n() >= 2) %>% ungroup()

kw_iucn <- kruskal.test(logMP ~ IUCN, data = iucn_data)
iucn_colors <- c("LC"="#4daf4a","NT"="#ff7f00","VU"="#e41a1c","EN"="#984ea3","CR"="#377eb8")

p_g_iucn <- ggplot(iucn_data, aes(x = IUCN, y = logMP)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, linewidth = 0.3, fill = "#d9d9d9") +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1, color = "gray20") +
  scale_fill_manual(values = iucn_colors, guide = "none") +
  labs(x = "", y = expression(log[e](MP+1)), tag = "g",
       subtitle = sprintf("Kruskal-Wallis p = %.3f", kw_iucn$p.value)) +
  theme_nature() +
  theme(plot.subtitle = element_text(size = 6.5, face = "italic"))

ggsave(here("output", "iucn.pdf"), p_g_iucn, width = 9, height = 8, units = "cm", device = "pdf")

# ==============================================================================
# 4. PANELS H & I — FUNCTIONAL UNIQUENESS AND ORIGINALITY
# ==============================================================================
cat("\n=== 4. GENERATING PANELS H & I ===\n")

# 1. Prepare continuous traits (log-transformed) and factor traits
fd_input_data <- species_phylo %>%
  filter(
    !is.na(Mass), !is.na(Beak.Length_Culmen), !is.na(Beak.Width), !is.na(Wing.Length),
    !is.na(Feeding_habit), !is.na(Habitat), !is.na(Migration)
  ) %>%
  mutate(
    log_Mass = log(Mass),
    log_BeakLength = log(Beak.Length_Culmen),
    log_BeakWidth = log(Beak.Width),
    log_WingLength = log(Wing.Length),
    Feeding_habit = as.factor(Feeding_habit),
    Habitat = as.factor(Habitat),
    Migration = as.factor(Migration)
  )

# 2. Extract matrix and set row names safely
trait_matrix <- fd_input_data %>%
  dplyr::select(log_Mass, log_BeakLength, log_BeakWidth, log_WingLength, 
                Feeding_habit, Habitat, Migration) %>%
  as.data.frame()

rownames(trait_matrix) <- fd_input_data$Species

# 3. Compute Gower Distance & Metrics once
gower_dist <- daisy(trait_matrix, metric = "gower")
pcoa <- cmdscale(gower_dist, k = 2, eig = TRUE)

dist_mat <- as.matrix(gower_dist)
diag(dist_mat) <- NA
funi <- apply(dist_mat, 1, mean, na.rm = TRUE)

centroid <- colMeans(pcoa$points)
forig <- apply(pcoa$points, 1, function(x) sqrt(sum((x - centroid)^2)))

fd_df <- data.frame(
  Species = rownames(trait_matrix),
  FUni = funi,
  FOrig = forig,
  logMP = fd_input_data$logMP
)

cor_funi  <- cor.test(fd_df$FUni,  fd_df$logMP, method = "spearman", exact = FALSE)
cor_forig <- cor.test(fd_df$FOrig, fd_df$logMP, method = "spearman", exact = FALSE)

# Panel H: Functional Uniqueness
p_h_funi <- ggplot(fd_df, aes(x = logMP, y = FUni)) +
  geom_point(color = "#2E6F95", alpha = 0.6, size = 1.5) +
  geom_smooth(method = "loess", color = "black", linewidth = 0.6, se = TRUE, alpha = 0.15) +
  geom_text_repel(data = fd_df %>% slice_max(FUni, n = 2), aes(label = Species),
                  size = 2.2, fontface = "italic", box.padding = 0.3) +
  labs(x = "log(Microplastic particles + 1)", y = "Functional uniqueness", tag = "h",
       subtitle = sprintf("rho = %.2f, p = %.3f (N = %d species)", 
                          cor_funi$estimate, cor_funi$p.value, nrow(fd_df))) +
  theme_nature() +
  theme(plot.subtitle = element_text(size = 6.5, face = "italic"))

ggsave(here("output", "funi.pdf"), p_h_funi, width = 12, height = 7.5, units = "cm", device = "pdf")

# Panel I: Functional Originality
p_i_forig <- ggplot(fd_df, aes(x = logMP, y = FOrig)) +
  geom_point(color = "#2E6F95", alpha = 0.6, size = 1.5) +
  geom_smooth(method = "loess", color = "black", linewidth = 0.6, se = TRUE, alpha = 0.15) +
  geom_text_repel(data = fd_df %>% slice_max(FOrig, n = 3), aes(label = Species),
                  size = 2.2, fontface = "italic", box.padding = 0.3) +
  labs(x = "log(Microplastic particles + 1)", y = "Functional originality", tag = "i",
       subtitle = sprintf("rho = %.2f, p = %.3f (N = %d species)", 
                          cor_forig$estimate, cor_forig$p.value, nrow(fd_df))) +
  theme_nature() +
  theme(plot.subtitle = element_text(size = 6.5, face = "italic"))

ggsave(here("output", "forig.pdf"), p_i_forig, width = 12, height = 7.5, units = "cm", device = "pdf")

# ==============================================================================
# 5. ASSEMBLE COMPLETE MULTI-PANEL FIGURE (PANELS A–I)
# ==============================================================================
cat("\n=== 5. ASSEMBLING FINAL FIGURE ===\n")

# Top section: Map + Boxplots
top_row <- p_b_map | (p_c_continent / p_d_lat)

# Middle section: Tree + Guilds
mid_row <- p_e_tree | (p_f_guild / p_g_iucn)

# Bottom section: Functional Diversity
bot_row <- p_h_funi | p_i_forig

# Combine using patchwork
final_figure <- top_row / mid_row / bot_row +
  plot_layout(heights = c(1.2, 1.2, 0.8))

print(final_figure)

# Save Outputs
ggsave(here("output", "Figure_Composite_Panels_a_i.pdf"), final_figure,
       width = 24, height = 28, units = "cm", device = "pdf", useDingbats = FALSE)

ggsave(here("output", "Figure_Composite_Panels_a_i.tiff"), final_figure,
       width = 24, height = 28, units = "cm", dpi = 600, compression = "lzw")

cat("\n Composite figure successfully saved to outputs/\n")
































################################################################################
# RESULTS VALUE EXTRACTION & OUTPUT
# Extracts aligned raw and log-transformed metrics.
################################################################################

cat("\n========================================================================\n")
cat("                  EXTRACTED RESULTS VALUES                                \n")
cat("========================================================================\n\n")

# ------------------------------------------------------------------------------
# 1. GEOGRAPHIC PATTERNS (PANELS B, C, D)
# ------------------------------------------------------------------------------
cat("--- 1. GEOGRAPHIC PATTERNS ---\n")
cat(sprintf("Total Record Observations (n): %d\n", nrow(df_clean)))
cat(sprintf("Total Species Analyzed (N): %d\n", length(unique(df_clean$Species))))
cat(sprintf("Percentage of records < 20 MP: %.1f%%\n", mean(df_clean$Mppersp < 20, na.rm = TRUE) * 100))

# Continent Summaries: Reports Raw Median [IQR] and Transformed log1p Median
cat("\nContinent Microplastic Burdens:\n")
df_clean %>%
  group_by(Continent) %>%
  summarise(
    N = n(),
    Median_Raw = median(Mppersp, na.rm = TRUE),
    IQR_Raw = IQR(Mppersp, na.rm = TRUE),
    Median_log1p = log1p(median(Mppersp, na.rm = TRUE))
  ) %>%
  arrange(desc(Median_Raw)) %>%
  mutate(
    Formatted_Text = sprintf("%s (N=%d): raw median = %.1f MP/individual [IQR: %.1f], log1p median = %.2f",
                             Continent, N, Median_Raw, IQR_Raw, Median_log1p)
  ) %>%
  pull(Formatted_Text) %>%
  walk(~ cat(.x, "\n"))

# Spearman Correlation with Latitude
spearman_lat <- cor.test(df_clean$Latitude, df_clean$Mppersp, method = "spearman", exact = FALSE)
cat(sprintf("\nLatitude Spearman Rank Test: Rho = %.3f, p = %.4f\n", 
            spearman_lat$estimate, spearman_lat$p.value))


# ------------------------------------------------------------------------------
# 2. PHYLOGENETIC & ECOLOGICAL DRIVERS (PANELS E, F, G)
# ------------------------------------------------------------------------------
cat("\n--- 2. PHYLOGENETIC & ECOLOGICAL DRIVERS ---\n")
cat(sprintf("Phylogenetic Species Pool (N): %d species\n", nrow(species_phylo)))

# Taxonomic Order Comparison
kw_order <- kruskal.test(logMP ~ Order, data = species_phylo)
cat(sprintf("\nOrder Kruskal-Wallis Test: Chi2 = %.2f, df = %d, p = %.4e\n", 
            kw_order$statistic, kw_order$parameter, kw_order$p.value))

cat("\nTop Taxonomic Order Medians:\n")
species_phylo %>%
  group_by(Order) %>%
  summarise(
    N = n(),
    Median_Raw = median(MP, na.rm = TRUE),
    Median_log1p = median(logMP, na.rm = TRUE)
  ) %>%
  filter(N >= 3) %>%
  arrange(desc(Median_Raw)) %>%
  print()

# Feeding Guild Comparison
guild_data <- species_phylo %>% 
  filter(!is.na(Feeding_habit)) %>% 
  group_by(Feeding_habit) %>% 
  filter(n() >= 3) %>% 
  ungroup()

kw_guild <- kruskal.test(logMP ~ Feeding_habit, data = guild_data)
cat(sprintf("\nFeeding Guild Kruskal-Wallis Test (Panel F): Chi2 = %.2f, df = %d, p = %.4e\n", 
            kw_guild$statistic, kw_guild$parameter, kw_guild$p.value))

cat("\nFeeding Guild Summary:\n")
guild_data %>%
  group_by(Feeding_habit) %>%
  summarise(
    N = n(),
    Median_Raw = median(MP, na.rm = TRUE),
    IQR_Raw = IQR(MP, na.rm = TRUE),
    Median_log1p = median(logMP, na.rm = TRUE)
  ) %>%
  arrange(desc(Median_Raw)) %>%
  print()

# IUCN Category Comparison (Switching to Median/IQR for Non-Parametric Consistency)
iucn_data <- species_phylo %>% 
  filter(!is.na(IUCN)) %>% 
  group_by(IUCN) %>% 
  filter(n() >= 2) %>% 
  ungroup()

kw_iucn <- kruskal.test(logMP ~ IUCN, data = iucn_data)
cat(sprintf("\nIUCN Status Kruskal-Wallis Test (Panel G): Chi2 = %.2f, df = %d, p = %.4f\n", 
            kw_iucn$statistic, kw_iucn$parameter, kw_iucn$p.value))

cat("\nIUCN Category Medians:\n")
iucn_data %>%
  group_by(IUCN) %>%
  summarise(
    N = n(),
    Median_Raw = median(MP, na.rm = TRUE),
    IQR_Raw = IQR(MP, na.rm = TRUE),
    Median_log1p = median(logMP, na.rm = TRUE)
  ) %>%
  print()


# ------------------------------------------------------------------------------
# 3. FUNCTIONAL DIVERSITY RELATIONSHIPS (PANELS H & I)
# ------------------------------------------------------------------------------
cat("\n--- 3. FUNCTIONAL DIVERSITY ---\n")
# Leverages fd_df computed in Section 4 to prevent redundant calculations
cat(sprintf("Functional Diversity Species Subset: N = %d species\n", nrow(fd_df)))

cor_funi  <- cor.test(fd_df$FUni,  fd_df$logMP, method = "spearman", exact = FALSE)
cor_forig <- cor.test(fd_df$FOrig, fd_df$logMP, method = "spearman", exact = FALSE)

cat(sprintf("Functional Uniqueness (Panel H)  : Spearman Rho = %.3f, p = %.4f\n", 
            cor_funi$estimate, cor_funi$p.value))
cat(sprintf("Functional Originality (Panel I) : Spearman Rho = %.3f, p = %.4f\n", 
            cor_forig$estimate, cor_forig$p.value))

# Key Species Callouts for Manuscript Text
top_funi <- fd_df %>% slice_max(FUni, n = 3) %>% pull(Species) %>% paste(collapse = ", ")
top_forig <- fd_df %>% slice_max(FOrig, n = 3) %>% pull(Species) %>% paste(collapse = ", ")

cat(sprintf("\nHighest FUni Species: %s\n", top_funi))
cat(sprintf("Highest FOrig Species: %s\n", top_forig))

cat("\n========================================================================\n")



























