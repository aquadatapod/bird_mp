################################################################################
# GEOGRAPHICAL DISTRIBUTION OF MICROPLASTIC BURDEN IN BIRDS
################################################################################

library(readxl)
library(dplyr)
library(ggplot2)
library(viridis)
library(sf)
library(rnaturalearth)
library(patchwork)
library(here)

# ============================================================================
# 1. LOAD AND PREPARE DATA
# ============================================================================

cat("\n=== LOADING GEOGRAPHICAL DATA ===\n")

# Load your Excel file
df <- read_excel(here("data", "map-bird-mp.xlsx"), sheet = "bird")

# Clean and prepare
df_clean <- df %>%
  filter(!is.na(Mppersp), !is.na(Latitude), !is.na(Longitude)) %>%
  mutate(
    logMP = log1p(Mppersp),
    MP_category = case_when(
      Mppersp == 0 ~ "0",
      Mppersp < 5 ~ "1-5",
      Mppersp < 20 ~ "6-20",
      Mppersp < 50 ~ "21-50",
      TRUE ~ ">50"
    ),
    MP_category = factor(MP_category, 
                         levels = c("0", "1-5", "6-20", "21-50", ">50"))
  )

cat(sprintf("Total observations: %d\n", nrow(df_clean)))
cat(sprintf("Unique species: %d\n", n_distinct(df_clean$Species)))

# ============================================================================
# 2. PERCENTAGE SUMMARIES (for your co-author's suggested text)
# ============================================================================

cat("\n=== MICROPLASTIC DISTRIBUTION SUMMARY ===\n")

# Overall distribution
total_n <- nrow(df_clean)

# Below 20 MP threshold
below_20 <- df_clean %>% filter(Mppersp < 20) %>% nrow()
pct_below_20 <- (below_20 / total_n) * 100

cat(sprintf("\nOverall MP distribution:\n"))
cat(sprintf("  Observations below 20 MP: %d (%.1f%%)\n", below_20, pct_below_20))
cat(sprintf("  Observations above 20 MP: %d (%.1f%%)\n", 
            total_n - below_20, 100 - pct_below_20))

# By category
cat("\nMP category breakdown:\n")
df_clean %>%
  group_by(MP_category) %>%
  summarise(
    n = n(),
    percentage = n() / total_n * 100
  ) %>%
  arrange(MP_category) %>%
  mutate(
    display = sprintf("  %s: %d (%.1f%%)", MP_category, n, percentage)
  ) %>%
  pull(display) %>%
  cat(sep = "\n")

# Continent-wise summary
cat("\n=== CONTINENT-WISE SUMMARY ===\n")

continent_summary <- df_clean %>%
  group_by(Continent) %>%
  summarise(
    n = n(),
    mean_MP = mean(Mppersp, na.rm = TRUE),
    sd_MP = sd(Mppersp, na.rm = TRUE),
    median_MP = median(Mppersp, na.rm = TRUE),
    pct_below_20 = sum(Mppersp < 20) / n() * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(mean_MP))

print(continent_summary)

# Save summary
write.csv(continent_summary, "continent_summary.csv", row.names = FALSE)



################################################################################
# GEOGRAPHICAL VISUALIZATION
# Map + Continent boxplot + Latitudinal distribution
################################################################################

library(ggplot2)
library(dplyr)
library(sf)
library(rnaturalearth)
library(patchwork)
library(viridis)

# ============================================================================
# 1. PREPARE DATA
# ============================================================================

# Get world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Calculate latitudinal distribution summary
lat_summary <- df_clean %>%
  mutate(
    lat_band = cut(Latitude, 
                   breaks = seq(-60, 90, by = 15),
                   labels = c("60°S-45°S", "45°S-30°S", "30°S-15°S", 
                              "15°S-0°", "0°-15°N", "15°N-30°N", 
                              "30°N-45°N", "45°N-60°N", "60°N-75°N", "75°N-90°N"))
  ) %>%
  group_by(lat_band) %>%
  summarise(
    n = n(),
    median_MP = median(Mppersp, na.rm = TRUE),
    mean_MP = mean(Mppersp, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(!is.na(lat_band))

# ============================================================================
# 2. PANEL A: WORLD MAP 
# ============================================================================

p_map <- ggplot() +
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

# ============================================================================
# 3. PANEL B: CONTINENT BOXPLOT (Simplified)
# ============================================================================

# Reorder continents by median MP
continent_order <- df_clean %>%
  group_by(Continent) %>%
  summarise(median_MP = median(Mppersp, na.rm = TRUE)) %>%
  arrange(desc(median_MP)) %>%
  pull(Continent)

df_clean$Continent <- factor(df_clean$Continent, levels = continent_order)

# Continent summary for sample sizes
continent_sizes <- df_clean %>%
  group_by(Continent) %>%
  summarise(n = n(), .groups = 'drop')

p_continent <- ggplot(df_clean, aes(x = Continent, y = Mppersp)) +
  geom_boxplot(fill = "#d9d9d9", 
               outlier.shape = NA, 
               width = 0.5, 
               size = 0.3,
               color = "#4d4d4d") +
  geom_jitter(width = 0.15, 
              alpha = 0.15, 
              size = 0.5, 
              color = "#4d4d4d") +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000),
                labels = c("0.1", "1", "10", "100", "1000"),
                expand = expansion(mult = c(0.05, 0.1))) +
  labs(x = "", y = "MP particles", tag = "b") +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 7, color = "black"),
    axis.line = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    plot.tag = element_text(size = 14, face = "bold"),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  geom_text(data = continent_sizes,
            aes(x = Continent, y = 0.03, label = paste0("n=", n)),
            size = 2, color = "gray30", vjust = 1.5, inherit.aes = FALSE)

# ============================================================================
# 4. PANEL C: LATITUDINAL DISTRIBUTION
# ============================================================================

p_lat <- ggplot(df_clean, aes(x = Latitude, y = Mppersp)) +
  # Scatter points
  geom_point(aes(color = logMP), alpha = 0.5, size = 1) +
  # Smoothed trend
  geom_smooth(method = "loess", color = "#d95f02", se = TRUE, alpha = 0.2, size = 0.8) +
  # Vertical lines for tropics
  geom_vline(xintercept = c(23.5, -23.5), linetype = "dotted", color = "#4d4d4d", size = 0.2) +
  annotate("text", x = 0, y = 500, label = "Equator", size = 2, color = "#4d4d4d") +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000),
                labels = c("0.1", "1", "10", "100", "1000")) +
  scale_color_viridis_c(option = "plasma", guide = "none") +
  labs(x = "Latitude", y = "MP particles", tag = "c") +
  theme_classic(base_size = 8) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 7, color = "black"),
    axis.line = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    plot.tag = element_text(size = 14, face = "bold"),
    plot.margin = margin(2, 2, 2, 2)
  )

# ============================================================================
# 5. COMBINE ALL PANELS
# ============================================================================

# Layout: Map on top, continent and latitude below
final_figure <- p_map / (p_continent | p_lat) +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    #title = "Geographical patterns of microplastic burden in birds",
    #subtitle = paste0("n = ", nrow(df_clean), " observations"),
    theme = theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5, face = "italic", color = "gray30")
    )
  )

# Save 
ggsave("Figure_Geography_Nature.pdf", 
       final_figure, 
       width = 30,  # cm (2-column width)
       height = 15, # cm
       units = "cm",
       dpi = 300,
       device = "pdf")

ggsave("Figure_Geography_Nature.tiff", 
       final_figure, 
       width = 18, 
       height = 16, 
       units = "cm", 
       dpi = 600, 
       compression = "lzw")

# ============================================================================
# 6. STATISTICAL SUMMARY FOR TEXT
# ============================================================================

cat("\n=== LATITUDINAL PATTERNS ===\n")
cat("\nBy latitude band:\n")
print(lat_summary)

# Test for latitudinal trend
cor_test <- cor.test(df_clean$Latitude, log1p(df_clean$Mppersp), method = "spearman")
cat("\nLatitudinal correlation:\n")
cat(sprintf("Spearman's ρ = %.2f, p = %.3f\n", cor_test$estimate, cor_test$p.value))




#################################################################################
################################################################################

# ============================================================================
# 5. PRINT STATISTICS
# ============================================================================

cat("\n=== MAP SUMMARY ===\n")
cat(sprintf("Dimensions: 30cm × 15cm\n"))
cat(sprintf("Total observations: %d\n", nrow(df_clean)))
cat(sprintf("Latitude range: %.1f° to %.1f°\n", 
            min(df_clean$Latitude, na.rm = TRUE),
            max(df_clean$Latitude, na.rm = TRUE)))
cat(sprintf("Longitude range: %.1f° to %.1f°\n", 
            min(df_clean$Longitude, na.rm = TRUE),
            max(df_clean$Longitude, na.rm = TRUE)))

cat("\n✓ Map saved as:\n")
cat("  - Figure_Map_Nature.pdf (vector)\n")
cat("  - Figure_Map_Nature.tiff (300 dpi, compressed)\n")
