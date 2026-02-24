








# Load required libraries
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggpubr)



# Read data
fd_raw <- read.csv("fd-ind.csv", stringsAsFactors = FALSE, na.strings = c("", "NA"))

# Calculate MP load normalized by body mass (this is your "ingestion proxy")
fd_raw$MP_per_g <- fd_raw$Mppersp / fd_raw$Mass

# For log-transformation (log-transform allometric data)
fd_raw$log_MP <- log10(fd_raw$Mppersp + 1)  # +1 to handle zeros
fd_raw$log_Mass <- log10(fd_raw$Mass)
fd_raw$log_Beak <- log10(fd_raw$Beak_Length_Culmen)

# Check your new data structure
summary(fd_raw$MP_per_g)


# ------------------------------------------------------------------------------
# PANEL A: Body mass threshold effect with sample sizes-------------------------Figure 6a
# ------------------------------------------------------------------------------

# Calculate sample sizes
n_small <- sum(fd_raw$Mass < 164, na.rm = TRUE)
n_large <- sum(fd_raw$Mass >= 164, na.rm = TRUE)

# Calculate medians for annotation
median_small <- median(fd_raw$Mppersp[fd_raw$Mass < 164], na.rm = TRUE)
median_large <- median(fd_raw$Mppersp[fd_raw$Mass >= 164], na.rm = TRUE)
fold_diff <- median_large / median_small

pA <- ggplot(fd_raw, aes(x = Mass < 164, y = Mppersp, fill = Mass < 164)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5, size = 0.3) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 0.8, color = "gray20") +
  scale_y_log10(breaks = c(1, 10, 100, 1000),
                labels = c("1", "10", "100", "1000")) +
  scale_fill_manual(values = c("FALSE" = "#2E86AB", "TRUE" = "#A23B72")) +
  labs(x = "", y = "MPs per individual") +
  scale_x_discrete(labels = c("FALSE" = paste0("≥164g\n(n = ", n_large, ")"),
                              "TRUE" = paste0("<164g\n(n = ", n_small, ")"))) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.title.y = element_text(size = 8, color = "black"),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  annotate("text", x = 1.5, y = max(fd_raw$Mppersp, na.rm = TRUE) * 0.7, 
           label = paste0("p = 0.006\nd = -0.76 [-1.11, -0.40]"), 
           size = 2.5, hjust = 0.5, fontface = "italic")

# ------------------------------------------------------------------------------
# PANEL B: Allometric relationship with threshold line--------------------------Figure 6b
# ------------------------------------------------------------------------------

pB <- ggplot(fd_raw, aes(x = Mass, y = Mppersp, color = Habitat)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_vline(xintercept = 164, linetype = "dashed", color = "red", size = 0.3) +
  scale_x_log10(breaks = c(10, 100, 1000, 10000),
                labels = c("10", "100", "1000", "10000")) +
  scale_y_log10(breaks = c(1, 10, 100, 1000),
                labels = c("1", "10", "100", "1000")) +
  scale_color_manual(values = c(
    "Aquatic" = "#5e9cc2", "Coastal" = "#8fc6a4", "Forest" = "#588157",
    "Grassland" = "#d9b650", "Human modified" = "#d96c50", "Marine" = "#4d7a9c",
    "Reverine" = "#6a8d92", "Shrubland" = "#b68b6d", "Wetland" = "#7aa5b5",
    "Woodland" = "#8b9a6e"
  )) +
  labs(x = "Body mass (g)", y = "MPs per individual", color = "Habitat") +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 8, color = "black"),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

# ------------------------------------------------------------------------------
# PANEL C: Habitat ranking (clean version with letters for significance)--------Figure 6c
# ------------------------------------------------------------------------------

# Prepare habitat data for panel C
habitat_c <- fd_raw %>%
  filter(Habitat %in% c("Forest", "Shrubland", "Coastal", "Marine")) %>%
  mutate(Habitat = factor(Habitat, 
                          levels = c("Marine", "Coastal", "Forest", "Shrubland")))

# Sample sizes for panel C
sample_c <- habitat_c %>%
  group_by(Habitat) %>%
  summarise(n = n(), .groups = 'drop') %>%
  mutate(label = paste0(Habitat, "\n(n=", n, ")"))

pC <- ggplot(habitat_c, aes(x = Habitat, y = MP_per_g, fill = Habitat)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.4, size = 0.3) +
  geom_jitter(width = 0.1, alpha = 0.25, size = 0.7, color = "gray20") +
  scale_y_log10(limits = c(0.001, 8),
                breaks = c(0.001, 0.01, 0.1, 1, 10),
                labels = c("0.001", "0.01", "0.1", "1", "10")) +
  scale_fill_manual(values = c(
    "Marine" = "#4d7a9c",    # blue
    "Coastal" = "#8fc6a4",   # green
    "Forest" = "#588157",    # forest green
    "Shrubland" = "#b68b6d"  # brown
  )) +
  labs(x = "", y = "MPs per gram") +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.title.y = element_text(size = 8, color = "black"),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  scale_x_discrete(labels = setNames(sample_c$label, sample_c$Habitat)) +
  
  # Add significance letters
  annotate("text", x = 1, y = 0.3, label = "a", size = 3.5, fontface = "bold") +
  annotate("text", x = 2, y = 0.8, label = "b", size = 3.5, fontface = "bold") +
  annotate("text", x = 3, y = 2.5, label = "c", size = 3.5, fontface = "bold") +
  annotate("text", x = 4, y = 3.5, label = "c", size = 3.5, fontface = "bold")

# ------------------------------------------------------------------------------
# COMBINE ALL PANELS
# ------------------------------------------------------------------------------

final_figure <- (pA | pB) / pC +
  plot_annotation(
    tag_levels = 'a',
    theme = theme(plot.tag = element_text(size = 10, face = "bold"))
  ) &
  theme(plot.tag.position = c(0, 0.95))



# Save main figure
ggsave("Figure_1_main.tiff", 
       plot = final_figure,
       width = 174, 
       height = 160, 
       units = "mm", 
       dpi = 300, 
       compression = "lzw")

ggsave("Figure_1_main.pdf", 
       plot = final_figure,
       width = 174, 
       height = 160, 
       units = "mm", 
       dpi = 300)

# Display
print(final_figure)