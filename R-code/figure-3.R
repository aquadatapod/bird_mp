# Load required libraries
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(rstatix)
library(multcompView)
library(readxl)
library(patchwork)  # For consistent theming

# ------------------------------------------------------------------------------
# Define consistent color palettes to match Figures 1 & 2
# ------------------------------------------------------------------------------

# Default color palette for tissues (matches your previous style)
mp_low <- c(
  "#5e9cc2", "#8fc6a4", "#588157", "#d9b650",
  "#d96c50", "#4d7a9c", "#6a8d92", "#b68b6d",
  "#7aa5b5", "#8b9a6e", "#A23B72", "#2E86AB"
)

# Alternative grayscale for print (uncomment if needed)
# mp_low <- rep(c("#f0f0f0", "#d9d9d9", "#bdbdbd", "#969696", "#737373"), 2)

# ------------------------------------------------------------------------------
# Load and prepare data
# ------------------------------------------------------------------------------

# Load bird microplastic data
bird.dat <- read_excel(here("data", "map-bird-mp.xlsx"),
  sheet = "bird"
)

# Initial inspection
summary(bird.dat$Mppersp)
table(bird.dat$Mppersp == 0)

# Data filtering and factor ordering
df_tissue <- bird.dat %>%
  filter(
    !is.na(Mppersp),
    !is.na(Tissue),
    !is.na(Longitude),
    !is.na(Latitude)
  ) %>%
  mutate(
    Tissue = factor(
      Tissue,
      levels = c(
        "Crop", "Esophagus", "FP", "Gizzard",
        "GIT", "RP", "Dropping", "Other",
        "Preen oil", "Serum"
      )
    )
  ) %>%
  droplevels()

# Sample sizes per tissue
n_df <- df_tissue %>%
  group_by(Tissue) %>%
  summarise(n = n(), .groups = "drop")

# ------------------------------------------------------------------------------
# Statistical testing (matching your approach)
# ------------------------------------------------------------------------------

# Kruskal-Wallis test
kw_test <- kruskal_test(df_tissue, Mppersp ~ Tissue)

# Dunn post-hoc test with BH correction (more conservative than Bonferroni)
dunn_res <- df_tissue %>%
  dunn_test(Mppersp ~ Tissue, p.adjust.method = "BH")

print(dunn_res, n = 28)

# ------------------------------------------------------------------------------
# Create compact letter display for significance
# ------------------------------------------------------------------------------

# Extract adjusted p-values and create comparison names
pw <- dunn_res$p.adj
names(pw) <- paste(dunn_res$group1, dunn_res$group2, sep = "-")

# Get compact letters
letters_raw <- multcompView::multcompLetters(pw)$Letters

# Create data frame with letters
letters_df <- data.frame(
  Tissue = names(letters_raw),
  Letters = letters_raw,
  stringsAsFactors = FALSE
)

# Calculate y-position for letters (above the highest point)
y_pos <- df_tissue %>%
  group_by(Tissue) %>%
  summarise(y = max(Mppersp, na.rm = TRUE) * 1.15, .groups = "drop")

# Merge letters with y-positions
letters_df <- left_join(letters_df, y_pos, by = "Tissue")

# ------------------------------------------------------------------------------
# FIGURE 3: Tissue-level MP accumulation
# Styled to match Figures 1 & 2
# ------------------------------------------------------------------------------

fig3_tissue <- ggplot(df_tissue, aes(x = Tissue, y = Mppersp, fill = Tissue)) +
  
  # Boxplot (matching style of Figures 1 & 2)
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5, size = 0.3) +
  
  # Jittered points
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.2, color = "gray20") +
  
  # Significance letters (bold, matching Figure 1)
  geom_text(
    data = letters_df,
    aes(x = Tissue, y = y, label = Letters),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold",
    color = "black"
  ) +
  
  # Sample sizes below x-axis (matching Figure 1)
  geom_text(
    data = n_df,
    aes(x = Tissue, y = -1.5, label = paste0("n=", n)),
    inherit.aes = FALSE,
    size = 3,
    color = "gray30"
  ) +
  
  # Fill colors (using your mp_low palette, recycled if needed)
  scale_fill_manual(
    values = rep(mp_low, length.out = length(levels(df_tissue$Tissue)))
  ) +
  
  # Y-axis with pseudo-log transformation (handles zeros well)
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 1, 5, 10, 20, 50, 100, 200, 400, 800),
    labels = c("0", "1", "5", "10", "20", "50", "100", "200", "400", "800")
  ) +
  
  # Labels
  labs(
    x = "",
    y = expression("Microplastics per individual" ~ (MP~ind^{-1})),
    title = NULL,
    subtitle = NULL
  ) +
  
  # Annotation for Kruskal-Wallis result
  annotate(
    "text",
    x = 0.6,
    y = max(df_tissue$Mppersp, na.rm = TRUE) * 1.35,
    label = paste0(
      "Kruskal-Wallis: p = ", 
      format.pval(kw_test$p, digits = 2, eps = 0.001)
    ),
    hjust = 0,
    size = 3.5,
    fontface = "italic",
    color = "gray20"
  ) +
  
  # Theme - matching Figures 1 & 2 exactly
  theme_classic(base_size = 10) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      size = 8, 
      color = "black",
      face = "plain"
    ),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 9, color = "black"),
    axis.line = element_line(color = "black", size = 0.2),
    plot.margin = margin(15, 10, 20, 10),
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9, face = "italic", hjust = 0.5)
  )

# View the plot
print(fig3_tissue)

# ------------------------------------------------------------------------------
# Save the figure (matching dimensions from Figures 1 & 2)
# ------------------------------------------------------------------------------

# For main text (single column width - 85mm)
ggsave(
  "Figure_3_tissue_main.tiff",
  plot = fig3_tissue,
  width = 85,
  height = 85,
  units = "mm",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  "Figure_3_tissue_main.pdf",
  plot = fig3_tissue,
  width = 120,
  height = 120,
  units = "mm",
  dpi = 300
)

# For supplementary (wider version if needed)
ggsave(
  "Figure_S2_tissue_detail.tiff",
  plot = fig3_tissue,
  width = 174,  # 1.5 column width
  height = 120,
  units = "mm",
  dpi = 300,
  compression = "lzw"
)








################################################################################ 
# ============================================================================
# LINK TISSUE ACCUMULATION TO FUNCTIONAL TRAITS
# Using your df_tissue with correct tissue categories
# ============================================================================

library(dplyr)
library(ggplot2)
library(FSA)
library(ggpubr)

# ----------------------------------------------------------------------------
# 1. Define tissue categories based on your table
# ----------------------------------------------------------------------------

tissue_categories <- df_tissue %>%
  mutate(
    Tissue_Group = case_when(
      Tissue %in% c("Crop", "Esophagus", "Gizzard", "GIT", "RP") ~ "Digestive",
      Tissue %in% c("FP", "Dropping") ~ "Excreta/Other",
      Tissue %in% c("Other") ~ "Peripheral",
      TRUE ~ "Other"
    ),
    logMP = log1p(Mppersp)
  )

# Verify categorization
table(tissue_categories$Tissue, tissue_categories$Tissue_Group)

# ----------------------------------------------------------------------------
# 2. Digestive tissues only - by Feeding Guild (Feeding_habit)
# ----------------------------------------------------------------------------

digestive <- tissue_categories %>%
  filter(Tissue_Group == "Digestive", !is.na(Feeding_habit))

cat("\n=== DIGESTIVE MP BY FEEDING GUILD ===\n")
cat("Sample sizes:\n")
print(table(digestive$Feeding_habit))

# Kruskal-Wallis test
kw_guild <- kruskal.test(logMP ~ Feeding_habit, data = digestive)
cat(sprintf("\nKruskal-Wallis: χ² = %.2f, df = %d, p = %.4f\n", 
            kw_guild$statistic, kw_guild$parameter, kw_guild$p.value))

# Summary statistics
guild_summary <- digestive %>%
  group_by(Feeding_habit) %>%
  summarise(
    n = n(),
    mean_MP = round(mean(Mppersp, na.rm = TRUE), 1),
    median_MP = round(median(Mppersp, na.rm = TRUE), 1),
    sd_MP = round(sd(Mppersp, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(median_MP))

print(guild_summary)

# Pairwise comparisons if significant
if(kw_guild$p.value < 0.05) {
  dunn_guild <- dunnTest(logMP ~ Feeding_habit, data = digestive, method = "holm")
  cat("\nPairwise comparisons (Holm-adjusted):\n")
  print(dunn_guild)
}

# ----------------------------------------------------------------------------
# 3. Body mass correlation with digestive MP
# ----------------------------------------------------------------------------

digestive_mass <- digestive %>%
  filter(!is.na(Mass), Mass > 0) %>%
  mutate(logMass = log(Mass))

cor_mass <- cor.test(digestive_mass$logMP, digestive_mass$logMass, 
                     method = "spearman")

cat("\n=== BODY MASS vs DIGESTIVE MP ===\n")
cat(sprintf("Spearman's ρ = %.3f, p = %.4f\n", cor_mass$estimate, cor_mass$p.value))
cat(sprintf("n = %d\n", nrow(digestive_mass)))

# ----------------------------------------------------------------------------
# 4. Beak length correlation with digestive MP
# ----------------------------------------------------------------------------

digestive_beak <- digestive %>%
  filter(!is.na(Beak.Length_Culmen), Beak.Length_Culmen > 0) %>%
  mutate(logBeak = log(Beak.Length_Culmen))

cor_beak <- cor.test(digestive_beak$logMP, digestive_beak$logBeak, 
                     method = "spearman")

cat("\n=== BEAK LENGTH vs DIGESTIVE MP ===\n")
cat(sprintf("Spearman's ρ = %.3f, p = %.4f\n", cor_beak$estimate, cor_beak$p.value))
cat(sprintf("n = %d\n", nrow(digestive_beak)))

# ----------------------------------------------------------------------------
# 5. Visualizations
# ----------------------------------------------------------------------------

# Boxplot by guild
p1 <- ggplot(digestive, aes(x = reorder(Feeding_habit, logMP, FUN = median), 
                            y = logMP, fill = Feeding_habit)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  labs(x = "Feeding guild", y = "log(MP+1) in digestive tract",
       title = "Digestive MP by feeding guild") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  stat_compare_means(method = "kruskal.test", label.y = max(digestive$logMP) * 1.1)

# Body mass correlation
p2 <- ggplot(digestive_mass, aes(x = logMass, y = logMP)) +
  geom_point(aes(color = Feeding_habit), alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  labs(x = "log(Body mass)", y = "log(MP+1) in digestive tract",
       title = "Body mass predicts digestive MP") +
  theme_classic() +
  theme(legend.position = "bottom")

# Beak length correlation
p3 <- ggplot(digestive_beak, aes(x = logBeak, y = logMP)) +
  geom_point(aes(color = Feeding_habit), alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  labs(x = "log(Beak length)", y = "log(MP+1) in digestive tract",
       title = "Beak length predicts digestive MP") +
  theme_classic() +
  theme(legend.position = "bottom")

# Save plots
ggsave("tissue_by_guild.pdf", p1, width = 8, height = 6)
ggsave("tissue_mass_correlation.pdf", p2, width = 7, height = 5)
ggsave("tissue_beak_correlation.pdf", p3, width = 7, height = 5)

# ============================================================================
# TISSUE PLOTS
# Clean, simple, publication-ready
# ============================================================================

library(ggplot2)
library(dplyr)
library(ggpubr)
library(patchwork)

# Define color palette (consistent with your paper)
guild_colors <- c(
  "Carnivore" = "#D55E00",      # Terracotta
  "Herbivore" = "#009E73",       # Green
  "Insectivore" = "#56B4E9",     # Blue
  "Omnivore" = "#2E86AB",        # Deep blue
  "Piscivore" = "#4d7a9c"        # Steel blue
  # Add other guilds as needed
)

# ----------------------------------------------------------------------------
# PLOT 1: Digestive MP by Feeding Guild (Boxplot)
# ----------------------------------------------------------------------------

# Calculate y-position for Kruskal-Wallis label
y_max_guild <- max(digestive$logMP, na.rm = TRUE)

p1 <- ggplot(digestive, aes(x = reorder(Feeding_habit, logMP, FUN = median), 
                            y = logMP)) +
  
  # Boxplot with thin lines, no fill
  geom_boxplot(aes(fill = Feeding_habit), 
               alpha = 0.7, 
               outlier.shape = NA, 
               width = 0.5,
               size = 0.3) +
  
  # Jittered points (subtle)
  geom_jitter(width = 0.15, 
              alpha = 0.3, 
              size = 0.8, 
              color = "gray20") +
  
  # Custom fill colors
  scale_fill_manual(values = guild_colors) +
  
  # Labels (no title)
  labs(x = "", 
       y = expression("log"[e] * "(MP+1) in digestive tract")) +
  
  # Clean theme
  theme_classic(base_size = 10) +
  theme(
    # Axis formatting
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 9, color = "black"),
    axis.line = element_line(color = "black", size = 0.2),
    
    # No legend
    legend.position = "none",
    
    # Panel background
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    
    # Margins
    plot.margin = margin(10, 10, 5, 10)
  ) +
  
  # Add Kruskal-Wallis result (subtle, in corner)
  annotate("text", 
           x = length(unique(digestive$Feeding_habit)) - 0.5,
           y = y_max_guild * 0.9,
           label = paste0("Kruskal-Wallis, p = ", 
                          format.pval(kw_guild$p.value, digits = 2)),
           size = 3, 
           fontface = "italic", 
           color = "gray40",
           hjust = 1) +
  
  # Panel label
  annotate("text", 
           x = 0.5, 
           y = y_max_guild * 1.05,
           label = "a", 
           size = 5, 
           fontface = "bold", 
           color = "black")

# ----------------------------------------------------------------------------
# PLOT 2: Body Mass Correlation--------------------------------------------------------- Figure 3c
# ----------------------------------------------------------------------------

y_max_mass <- max(digestive_mass$logMP, na.rm = TRUE)

p2 <- ggplot(digestive_mass, aes(x = logMass, y = logMP)) +
  
  # Points with guild colors
  geom_point(aes(color = Feeding_habit), 
             alpha = 0.6, 
             size = 1.5) +
  
  # Simple regression line (no confidence band for cleaner look)
  #geom_smooth(method = "lm", color = "black", se = FALSE, size = 0.5) +
  geom_smooth(method = "loess", color = "black", se = TRUE)+
  
  # Custom colors
  scale_color_manual(values = guild_colors, name = "Feeding guild") +
  
  # Labels
  labs(x = expression("log"[e] * "(Body mass, g)"), 
       y = expression("log"[e] * "(MP+1) in digestive tract")) +
  
  # Clean theme
  theme_classic(base_size = 10) +
  theme(
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 9, color = "black"),
    axis.line = element_line(color = "black", size = 0.2),
    
    # Legend at bottom, smaller
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8, face = "bold"),
    legend.key.size = unit(0.4, "cm"),
    
    # Margins
    plot.margin = margin(10, 10, 5, 10)
  ) +
  
  # Add correlation result
  annotate("text", 
           x = min(digestive_mass$logMass) + 0.5,
           y = y_max_mass * 0.95,
           label = paste0("ρ = ", round(cor_mass$estimate, 2), 
                          ", p = ", format.pval(cor_mass$p.value, digits = 2)),
           size = 3, 
           fontface = "italic", 
           color = "gray40",
           hjust = 0) +
  
  # Panel label
  annotate("text", 
           x = min(digestive_mass$logMass) - 0.1,
           y = y_max_mass * 1.05,
           label = "b", 
           size = 5, 
           fontface = "bold", 
           color = "black")

# ----------------------------------------------------------------------------
# PLOT 3: Beak Length Correlation--------------------------------------------------------Figure 3d
# ----------------------------------------------------------------------------

y_max_beak <- max(digestive_beak$logMP, na.rm = TRUE)

p3 <- ggplot(digestive_beak, aes(x = logBeak, y = logMP)) +
  
  geom_point(aes(color = Feeding_habit), 
             alpha = 0.6, 
             size = 1.5) +
  
  #geom_smooth(method = "lm", color = "black", se = FALSE, size = 0.5) +
  geom_smooth(method = "loess", color = "black", se = TRUE)+
  
  scale_color_manual(values = guild_colors, name = "Feeding guild") +
  
  labs(x = expression("log"[e] * "(Beak length, mm)"), 
       y = expression("log"[e] * "(MP+1) in digestive tract")) +
  
  theme_classic(base_size = 10) +
  theme(
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 9, color = "black"),
    axis.line = element_line(color = "black", size = 0.2),
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8, face = "bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(10, 10, 5, 10)
  ) +
  
  annotate("text", 
           x = min(digestive_beak$logBeak) + 0.2,
           y = y_max_beak * 0.95,
           label = paste0("ρ = ", round(cor_beak$estimate, 2), 
                          ", p = ", format.pval(cor_beak$p.value, digits = 2)),
           size = 3, 
           fontface = "italic", 
           color = "gray40",
           hjust = 0) +
  
  annotate("text", 
           x = min(digestive_beak$logBeak) - 0.1,
           y = y_max_beak * 1.05,
           label = "c", 
           size = 5, 
           fontface = "bold", 
           color = "black")

# ----------------------------------------------------------------------------
# COMBINE PLOTS
# ----------------------------------------------------------------------------

# Arrange in a single row
tissue_figure <- p2 | p3

# Save
ggsave("Figure_tissue_linkage.pdf", 
       plot = tissue_figure,
       width = 10, height = 5, units = "in", dpi = 300)

ggsave("Figure_tissue_linkage.tiff", 
       plot = tissue_figure,
       width = 12, height = 4.5, units = "in", dpi = 300, compression = "lzw")

# Display
print(tissue_figure)













# ============================================================================
# SLOPE CHART - Connect body mass to MP for each guild
# Shows trend lines per guild
# ============================================================================

# Create bins for visualization
digestive_binned <- digestive_mass %>%
  mutate(Mass_bin = cut(logMass, breaks = 5))

p_slope <- ggplot(digestive_mass, aes(x = logMass, y = logMP, color = Feeding_habit)) +
  
  # Add trend lines per guild
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2, size = 0.8) +
  
  # Points
  geom_point(alpha = 0.5, size = 1.5) +
  
  scale_color_manual(values = guild_colors, name = "Feeding guild") +
  
  labs(x = "log(Body mass)", y = "log(MP+1) in digestive tract") +
  
  theme_classic() +
  theme(legend.position = "bottom")


# ============================================================================
# BIVARIATE HEATMAP - Shows where most data cluster
# Great for highlighting high-density regions
# ============================================================================

library(MASS)

p_heatmap <- ggplot(digestive_mass, aes(x = logMass, y = logMP)) +
  
  # 2D density contour
  stat_density_2d(aes(fill = after_stat(density)), 
                  geom = "raster", contour = FALSE, alpha = 0.8) +
  
  # Add points on top
  geom_point(alpha = 0.3, size = 1, color = "white") +
  
  # Add regression line
  geom_smooth(method = "lm", color = "red", se = FALSE, size = 0.8) +
  
  scale_fill_viridis_c(option = "plasma", name = "Density") +
  
  labs(x = "log(Body mass)", y = "log(MP+1) in digestive tract") +
  
  theme_classic() +
  
  # Highlight the high-density region where small birds have high MP
  annotate("text", x = 4, y = 3.5, 
           label = "High-risk zone", 
           color = "white", fontface = "bold", size = 4)




# ----------------------------------------------------------------------------
# 6. Digestive vs Peripheral/Excreta comparison
# ----------------------------------------------------------------------------

# Add Excreta to comparison
tissue_all <- tissue_categories %>%
  filter(!is.na(Feeding_habit), Tissue_Group %in% c("Digestive", "Excreta/Other"))

# Test if digestive > excreta
wt_excreta <- wilcox.test(logMP ~ Tissue_Group, data = tissue_all)
cat("\n=== DIGESTIVE vs EXCRETA MP ===\n")
cat(sprintf("Wilcoxon p = %.4f\n", wt_excreta$p.value))

# Summary by guild and tissue
tissue_guild_summary <- tissue_all %>%
  group_by(Feeding_habit, Tissue_Group) %>%
  summarise(
    n = n(),
    median_MP = round(median(Mppersp, na.rm = TRUE), 1),
    mean_MP = round(mean(Mppersp, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(Feeding_habit, desc(Tissue_Group))

print(tissue_guild_summary)

# ----------------------------------------------------------------------------
# 7. Create final summary table for paper
# ----------------------------------------------------------------------------

final_table <- digestive %>%
  group_by(Feeding_habit) %>%
  summarise(
    n_digestive = n(),
    Median_digestive_MP = round(median(Mppersp, na.rm = TRUE), 1),
    Mean_digestive_MP = round(mean(Mppersp, na.rm = TRUE), 1),
    SD_digestive_MP = round(sd(Mppersp, na.rm = TRUE), 1),
    Mass_correlation = round(cor(logMP, log(Mass), method = "spearman"), 2),
    Mass_p = cor.test(logMP, log(Mass), method = "spearman")$p.value,
    .groups = "drop"
  )

print(final_table)

# Save
write.csv(final_table, "tissue_trait_linkage.csv", row.names = FALSE)

cat("\n✅ Analysis complete! Results saved.\n")






















# ------------------------------------------------------------------------------
# OPTIONAL: Create a compact summary table of tissue statistics
# ------------------------------------------------------------------------------

tissue_summary <- df_tissue %>%
  group_by(Tissue) %>%
  summarise(
    N = n(),
    Mean = round(mean(Mppersp, na.rm = TRUE), 1),
    SD = round(sd(Mppersp, na.rm = TRUE), 1),
    Median = round(median(Mppersp, na.rm = TRUE), 1),
    IQR = round(IQR(Mppersp, na.rm = TRUE), 1),
    Min = min(Mppersp, na.rm = TRUE),
    Max = max(Mppersp, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(Median))

# Print table
knitr::kable(tissue_summary, 
             caption = "Table S2. Microplastic accumulation by tissue type",
             digits = 1)

# Save as CSV
write.csv(tissue_summary, "Table_S2_tissue_summary.csv", row.names = FALSE)



################################################################################
# Calculate mean, median, SD, and other summary statistics for each tissue
tissue_stats <- df_tissue %>%
  group_by(Tissue) %>%
  summarise(
    n = n(),
    mean_MP = mean(Mppersp, na.rm = TRUE),
    sd_MP = sd(Mppersp, na.rm = TRUE),
    median_MP = median(Mppersp, na.rm = TRUE),
    IQR_MP = IQR(Mppersp, na.rm = TRUE),
    min_MP = min(Mppersp, na.rm = TRUE),
    max_MP = max(Mppersp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_MP))  # Sort by mean MP (highest to lowest)

# Print the complete statistics table
cat("\n=== TISSUE-SPECIFIC MP ACCUMULATION STATISTICS ===\n")
print(tissue_stats)

# Format for manuscript: "mean ± SD" format
tissue_stats <- tissue_stats %>%
  mutate(
    mean_sd = sprintf("%.1f ± %.1f", mean_MP, sd_MP),
    median_iqr = sprintf("%.1f (%.1f-%.1f)", 
                         median_MP, 
                         quantile(df_tissue$Mppersp[df_tissue$Tissue == Tissue], 0.25, na.rm = TRUE),
                         quantile(df_tissue$Mppersp[df_tissue$Tissue == Tissue], 0.75, na.rm = TRUE))
  )

# Create a clean table for manuscript
manuscript_table <- tissue_stats %>%
  select(Tissue, n, mean_sd, median_iqr)

cat("\n=== FOR MANUSCRIPT (mean ± SD format) ===\n")
print(manuscript_table)

# Alternative: Create a publication-ready table
library(kableExtra)
manuscript_table_formatted <- manuscript_table %>%
  kbl(caption = "Microplastic accumulation across avian tissues",
      col.names = c("Tissue", "n", "Mean ± SD (particles)", "Median (IQR)"),
      align = c("l", "r", "r", "r")) %>%
  kable_classic(full_width = FALSE) %>%
  row_spec(1, bold = TRUE, color = "white", background = "#2E5A87")

# Save as HTML/PDF for supplementary materials
save_kable(manuscript_table_formatted, "Table_Tissue_MP_Statistics.html")

# For your manuscript text, extract specific tissues:
cat("\n=== SPECIFIC VALUES FOR MANUSCRIPT TEXT ===\n")

# Get GIT statistics
git_stats <- tissue_stats %>% filter(Tissue == "GIT")
cat(sprintf("GIT: %.1f ± %.1f particles (n = %d)\n", 
            git_stats$mean_MP, git_stats$sd_MP, git_stats$n))

# Get dropping statistics
drop_stats <- tissue_stats %>% filter(Tissue == "Dropping")
cat(sprintf("Droppings: %.1f ± %.1f particles (n = %d)\n", 
            drop_stats$mean_MP, drop_stats$sd_MP, drop_stats$n))

# Get crop statistics
crop_stats <- tissue_stats %>% filter(Tissue == "Crop")
cat(sprintf("Crop: %.1f ± %.1f particles (n = %d)\n", 
            crop_stats$mean_MP, crop_stats$sd_MP, crop_stats$n))

# Get gizzard statistics
giz_stats <- tissue_stats %>% filter(Tissue == "Gizzard")
cat(sprintf("Gizzard: %.1f ± %.1f particles (n = %d)\n", 
            giz_stats$mean_MP, giz_stats$sd_MP, giz_stats$n))

# Calculate the fold-difference between highest and lowest
highest <- max(tissue_stats$mean_MP, na.rm = TRUE)
lowest <- min(tissue_stats$mean_MP[tissue_stats$mean_MP > 0], na.rm = TRUE)
fold_difference <- highest / lowest

cat(sprintf("\nFold difference (highest/lowest): %.0f-fold\n", fold_difference))

# Save the complete statistics
write.csv(tissue_stats, "tissue_mp_statistics_complete.csv", row.names = FALSE)
cat("\n✓ Complete statistics saved as 'tissue_mp_statistics_complete.csv'\n")













##################################################################################
# Supplementary material

# Create Table S1: Data Characteristics
table_s1 <- bird.dat %>%
  summarise(
    `Total Observations` = n(),
    `Species with MP Data` = n_distinct(Species[!is.na(Mppersp)]),
    `Range of MP per Specimen` = paste0(min(Mppersp, na.rm = TRUE), " - ", 
                                        max(Mppersp, na.rm = TRUE)),
    `Median MP per Specimen` = median(Mppersp, na.rm = TRUE),
    `Zero MP Observations` = sum(Mppersp == 0, na.rm = TRUE),
    `Proportion Zero` = round(mean(Mppersp == 0, na.rm = TRUE) * 100, 1),
    `Skewness` = moments::skewness(Mppersp, na.rm = TRUE),
    `Kurtosis` = moments::kurtosis(Mppersp, na.rm = TRUE)
  ) %>%
  t() %>%  # Transpose to vertical format
  as.data.frame() %>%
  tibble::rownames_to_column("Statistic")

# Save as CSV
write.csv(table_s1, "Table_S1_Data_Characteristics.csv", row.names = FALSE)

# Print for manuscript
print(table_s1)



# First, inspect your data for blank values
print("=== CHECKING FOR BLANK VALUES ===")
print(paste("Rows with blank Species:", sum(bird.dat$Species == "" | is.na(bird.dat$Species))))
print(paste("Rows with blank Mppersp:", sum(bird.dat$Mppersp == "")))
print(paste("Rows with NA Mppersp:", sum(is.na(bird.dat$Mppersp))))

# Convert blanks to NA for consistent handling
bird.dat_clean <- bird.dat %>%
  mutate(
    # Convert empty strings to NA
    Species = ifelse(Species == "", NA, Species),
    Mppersp = ifelse(Mppersp == "", NA, Mppersp),
    # Convert Mppersp to numeric if it's character
    Mppersp = as.numeric(Mppersp)
  )

# Now create Table S1 with proper handling
table_s1 <- bird.dat_clean %>%
  summarise(
    # Total rows in dataset
    `Total Rows in Dataset` = n(),
    
    # Species with any data (non-blank, non-NA)
    `Species with Any Data` = n_distinct(Species[!is.na(Species)]),
    
    # Species with MP data (non-NA Mppersp)
    `Species with MP Data` = n_distinct(Species[!is.na(Mppersp)]),
    
    # Observations with MP data
    `Observations with MP Data` = sum(!is.na(Mppersp)),
    
    # Range of MP (only non-NA values)
    `Range of MP per Specimen` = ifelse(
      sum(!is.na(Mppersp)) > 0,
      paste0(min(Mppersp, na.rm = TRUE), " - ", max(Mppersp, na.rm = TRUE)),
      "No MP data"
    ),
    
    # Central tendency measures
    `Median MP per Specimen` = median(Mppersp, na.rm = TRUE),
    `Mean MP per Specimen` = mean(Mppersp, na.rm = TRUE),
    `SD of MP per Specimen` = sd(Mppersp, na.rm = TRUE),
    
    # Zero values (actual zeros, not missing)
    `Zero MP Observations` = sum(Mppersp == 0, na.rm = TRUE),
    
    # Proportion of zeros among observations WITH MP data
    `Proportion Zero (of MP data)` = ifelse(
      sum(!is.na(Mppersp)) > 0,
      round(sum(Mppersp == 0, na.rm = TRUE) / sum(!is.na(Mppersp)) * 100, 1),
      NA
    ),
    
    # Proportion of missing MP data
    `Proportion Missing MP Data` = round(sum(is.na(Mppersp)) / n() * 100, 1),
    
    # Distribution properties (only if sufficient data)
    `Skewness` = ifelse(
      sum(!is.na(Mppersp)) > 2,
      round(moments::skewness(Mppersp, na.rm = TRUE), 2),
      NA
    ),
    `Kurtosis` = ifelse(
      sum(!is.na(Mppersp)) > 2,
      round(moments::kurtosis(Mppersp, na.rm = TRUE), 2),
      NA
    )
  ) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Statistic")

# Print results
print(table_s1)






# Calculate Goodman's simultaneous multinomial confidence intervals
library(DescTools)

# Sample sizes per tissue
n_df <- df_tissue %>%
  group_by(Tissue) %>%
  summarise(n = n(), .groups = "drop")

# Total sample size
N <- sum(n_df$n)

# Calculate proportions
n_df <- n_df %>%
  mutate(
    proportion = n / N,
    percentage = proportion * 100
  )

# Calculate Goodman's 95% simultaneous CIs
# For multinomial proportions with k categories
goodman_ci <- MultinomCI(n_df$n, conf.level = 0.95, method = "goodman")

# Add to dataframe
n_df <- n_df %>%
  mutate(
    CI_lower = goodman_ci[, 2] * 100,  # Convert to percentage
    CI_upper = goodman_ci[, 3] * 100,
    CI_formatted = sprintf("%.1f (%.1f–%.1f)", percentage, CI_lower, CI_upper)
  )

# Print results
print(n_df[, c("Tissue", "n", "percentage", "CI_lower", "CI_upper", "CI_formatted")])

# Create a plot to show sampling representation
library(ggplot2)

fig_sampling <- ggplot(n_df, aes(x = reorder(Tissue, -percentage), y = percentage)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                width = 0.2, color = "darkred", linewidth = 0.8) +
  geom_text(aes(label = paste0("n=", n)), 
            vjust = -0.5, size = 3.5) +
  labs(x = "Tissue Type", y = "Sampling Proportion (%)",
       title = "Tissue Sampling Representation with 95% Goodman CIs",
       subtitle = "Error bars show simultaneous confidence intervals") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Figure_Sampling_Representation.pdf", fig_sampling, width = 10, height = 6)

# Save as table
write.csv(n_df, "Table_Tissue_Sampling.csv", row.names = FALSE)