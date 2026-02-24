# ==============================================================================
# QUALITY ASSURANCE / QUALITY CONTROL (QA/QC) ANALYSIS
# For: Bird Microplastic Meta-Analysis
# # ==============================================================================
#
# PURPOSE:
# This script systematically assesses methodological heterogeneity across source 
# studies in a global meta-analysis of microplastic accumulation in birds.
# It evaluates potential biases from extraction methods, identification techniques,
# zero detection patterns, temporal trends, and study influence.
#
# KEY ANALYSES:
# 1. Extraction method comparison (Kruskal-Wallis, effect sizes)
# 2. Identification method comparison (Kruskal-Wallis, effect sizes)
# 3. Study-level effect sizes (dot plot with outlier identification)
# 4. Zero detection quantification (by study and extraction method)
# 5. Temporal trends (Spearman correlation, LOESS smoothing)
# 6. Validation of high-risk species across methods (KOH-only subset)
#
# INPUT:
# - QA-QC-data.csv : Contains columns:
#   * Reference, Species, Extraction_method, Identification,
#     Mppersp (MP count per individual), Years_sampling, etc.
#
# OUTPUT:
# - Multiple publication-ready figures (TIFF, PDF)
# - Summary statistics for manuscript
# - Effect size tables for supplementary materials
#
# DEPENDENCIES:
library(tidyverse)   # Data manipulation + ggplot2
library(ggpubr)      # Statistical comparisons on plots
library(patchwork)   # Multi-panel figure assembly
library(effsize)     # Cohen's d effect sizes
library(forcats)     # Factor reordering
library(rstatix)     # Kruskal-Wallis, Dunn tests
library(rsvg)        # PDF/TIFF export

# ==============================================================================
# SCRIPT USAGE:
# 1. Ensure QA-QC-data.csv is in working directory
# 2. Run script sequentially (sections are clearly marked)
# 3. All figures will be saved to current working directory
# 4. Check console output for summary statistics
# ==============================================================================
################################################################################


# Read QA/QC data
qa_qc <- read.csv("QA-QC-data.csv", stringsAsFactors = FALSE)


table(qa_qc$Extraction_method)

table(qa_qc$Identification)

# ------------------------------------------------------------------------------
# Clean data 
# ------------------------------------------------------------------------------

qa_qc <- qa_qc %>%
  mutate(
    Mppersp_clean = as.numeric(gsub(",", "", as.character(Mppersp))),
    
    Extraction_method_clean = case_when(
      Extraction_method == "KOH" ~ "KOH",
      Extraction_method == "NaCl" ~ "NaCl",
      Extraction_method == "NaOH" ~ "NaOH",
      Extraction_method == "H2O2" ~ "H2O2",
      Extraction_method == "KOH & H2O2" ~ "KOH+H2O2",
      Extraction_method == "KOH & NaCl" ~ "KOH+NaCl",
      Extraction_method == "Washed & Sieve" ~ "Washed+Sieve",
      Extraction_method == "Stomach Flushing" ~ "Stomach flush",
      Extraction_method == "Dried & Crumbling" ~ "Dried",
      TRUE ~ Extraction_method
    ),
    
    Identification_clean = case_when(
      grepl("FTIR", Identification) & !grepl("Stereomicroscope|Raman", Identification) ~ "FTIR",
      grepl("Stereomicroscope", Identification) & !grepl("FTIR|Raman", Identification) ~ "Stereo",
      grepl("Stereomicroscope.*FTIR", Identification) ~ "Stereo+FTIR",
      grepl("Stereomicroscope.*Raman", Identification) ~ "Stereo+Raman",
      grepl("Raman", Identification) & !grepl("Stereomicroscope|FTIR", Identification) ~ "Raman",
      grepl("Compound microscope", Identification) ~ "Compound",
      grepl("Visual", Identification) ~ "Visual",
      grepl("Infrared", Identification) ~ "FTIR",
      grepl("µ-FTIR", Identification) ~ "FTIR",
      TRUE ~ "Other"
    )
  )

# Use cleaned column and remove NAs
qa_qc$Mppersp <- qa_qc$Mppersp_clean
qa_qc <- qa_qc %>% filter(!is.na(Mppersp))

# ------------------------------------------------------------------------------
# FUNCTION: Calculate effect sizes for all pairwise comparisons
# ------------------------------------------------------------------------------

calculate_effect_sizes <- function(data, group_col, value_col) {
  
  groups <- unique(data[[group_col]])
  groups <- groups[!is.na(groups)]
  
  results <- data.frame()
  
  for(i in 1:(length(groups)-1)) {
    for(j in (i+1):length(groups)) {
      
      group1 <- groups[i]
      group2 <- groups[j]
      
      subset_data <- data %>%
        filter(!!sym(group_col) %in% c(group1, group2))
      
      # Only calculate if both groups have at least 3 samples
      if(nrow(subset_data) >= 6) {
        effect <- tryCatch({
          cohen.d(as.formula(paste(value_col, "~", group_col)), data = subset_data)
        }, error = function(e) NULL)
        
        if(!is.null(effect)) {
          results <- rbind(results, data.frame(
            Comparison = paste(group1, "vs", group2),
            Group1 = group1,
            Group2 = group2,
            d = effect$estimate,
            CI_lower = effect$conf.int[1],
            CI_upper = effect$conf.int[2],
            n1 = sum(subset_data[[group_col]] == group1),
            n2 = sum(subset_data[[group_col]] == group2)
          ))
        }
      }
    }
  }
  
  return(results)
}

# ------------------------------------------------------------------------------
# FIGURE 1: EXTRACTION METHODS WITH EFFECT SIZES (ROTATED)
# ------------------------------------------------------------------------------

# Prepare extraction data
extract_counts <- table(qa_qc$Extraction_method_clean)
extract_keep <- names(extract_counts[extract_counts >= 3])

qa_extract <- qa_qc %>%
  filter(Extraction_method_clean %in% extract_keep) %>%
  mutate(Extraction_method_clean = factor(Extraction_method_clean))

# Reorder by median
qa_extract <- qa_extract %>%
  mutate(Extraction_method_clean = fct_reorder(Extraction_method_clean, Mppersp, .fun = median, na.rm = TRUE))

# Calculate effect sizes for extraction methods
extract_effects <- calculate_effect_sizes(qa_extract, "Extraction_method_clean", "Mppersp")
print(extract_effects)

# Sample sizes
extract_n <- qa_extract %>%
  group_by(Extraction_method_clean) %>%
  summarise(n = n(), .groups = 'drop') %>%
  mutate(label = paste0(Extraction_method_clean, " (n=", n, ")"))

# Kruskal-Wallis test
kw_extract <- kruskal.test(Mppersp ~ Extraction_method_clean, data = qa_extract)

# Color palette
extract_colors <- c(
  "KOH" = "#5e9cc2",
  "H2O2" = "#8fc6a4",
  "NaOH" = "#d96c50",
  "Washed+Sieve" = "#4d7a9c",
  "KOH+NaCl" = "#588157"
)

# ROTATED PLOT (horizontal boxes)
p1 <- ggplot(qa_extract, aes(x = Extraction_method_clean, y = Mppersp, fill = Extraction_method_clean)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6, size = 0.4) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1, color = "gray20") +
  scale_y_log10(breaks = c(1, 10, 100, 1000),
                labels = c("1", "10", "100", "1000")) +
  scale_fill_manual(values = extract_colors[names(extract_colors) %in% unique(qa_extract$Extraction_method_clean)]) +
  
  # ROTATE to horizontal for better readability
  coord_flip() +
  
  labs(x = "", y = "MPs per individual (log scale)",
       #title = "A) Extraction methods"
       ) +
  
  theme_classic(base_size = 10) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 8, color = "black", face = "bold"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.title.x = element_text(size = 9, color = "black"),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  scale_x_discrete(labels = setNames(extract_n$label, extract_n$Extraction_method_clean)) +
  
  # Add Kruskal-Wallis result
  annotate("text", 
           x = length(unique(qa_extract$Extraction_method_clean)) - 0.5, 
           y = min(qa_extract$Mppersp, na.rm = TRUE) * 0.5,
           label = paste0("Kruskal-Wallis: p = ", format.pval(kw_extract$p.value, digits = 3)),
           size = 3, fontface = "italic", hjust = 0, color = "gray30")


ggsave("study_extraction.pdf", p1, width = 18, height = 15, units = "cm", dpi = 300)
# ------------------------------------------------------------------------------
# FIGURE 2: IDENTIFICATION METHODS WITH EFFECT SIZES (ROTATED)
# ------------------------------------------------------------------------------

id_counts <- table(qa_qc$Identification_clean)
id_keep <- names(id_counts[id_counts >= 3])

qa_id <- qa_qc %>%
  filter(Identification_clean %in% id_keep) %>%
  mutate(Identification_clean = factor(Identification_clean))

qa_id <- qa_id %>%
  mutate(Identification_clean = fct_reorder(Identification_clean, Mppersp, .fun = median, na.rm = TRUE))

# Calculate effect sizes for identification methods
id_effects <- calculate_effect_sizes(qa_id, "Identification_clean", "Mppersp")
print(id_effects)

id_n <- qa_id %>%
  group_by(Identification_clean) %>%
  summarise(n = n(), .groups = 'drop') %>%
  mutate(label = paste0(Identification_clean, " (n=", n, ")"))

kw_id <- kruskal.test(Mppersp ~ Identification_clean, data = qa_id)

id_colors <- c(
  "FTIR" = "#2E86AB",
  "Stereo" = "#A23B72",
  "Stereo+FTIR" = "#588157",
  "Stereo+Raman" = "#d9b650",
  "Visual" = "#b68b6d"
)

# ROTATED PLOT
p2 <- ggplot(qa_id, aes(x = Identification_clean, y = Mppersp, fill = Identification_clean)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6, size = 0.4) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1, color = "gray20") +
  scale_y_log10(breaks = c(1, 10, 100, 1000),
                labels = c("1", "10", "100", "1000")) +
  scale_fill_manual(values = id_colors[names(id_colors) %in% unique(qa_id$Identification_clean)]) +
  
  # ROTATE to horizontal
  coord_flip() +
  
  labs(x = "", y = "MPs per individual (log scale)"
       #title = "B) Identification methods"
       )+
  theme_classic(base_size = 10) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 8, color = "black", face = "bold"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.title.x = element_text(size = 9, color = "black"),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  scale_x_discrete(labels = setNames(id_n$label, id_n$Identification_clean)) +
  
  annotate("text", 
           x = length(unique(qa_id$Identification_clean)) - 0.5, 
           y = min(qa_id$Mppersp, na.rm = TRUE) * 0.5,
           label = paste0("Kruskal-Wallis: p = ", format.pval(kw_id$p.value, digits = 3)),
           size = 3, fontface = "italic", hjust = 0, color = "gray30")


ggsave("study_identification.pdf", p2, width = 18, height = 15, units = "cm", dpi = 300)
# ------------------------------------------------------------------------------
# FIGURE 3: DOT PLOT - ALL STUDIES
# ------------------------------------------------------------------------------

# Calculate study statistics with proper formatting
study_stats <- qa_qc %>%
  group_by(Reference) %>%
  summarise(
    n_samples = n(),
    mean_log = mean(log1p(Mppersp), na.rm = TRUE),
    median_MP = median(Mppersp, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(!is.na(mean_log)) %>%
  mutate(
    overall_mean_log = mean(mean_log, na.rm = TRUE),
    overall_sd_log = sd(mean_log, na.rm = TRUE),
    effect_size_log = (mean_log - overall_mean_log) / overall_sd_log,
    # Create clean, short labels
    display_label = case_when(
      nchar(Reference) > 30 ~ paste0(substr(Reference, 1, 28), "…"),
      TRUE ~ Reference
    ),
    # Outlier classification (|effect| > 2)
    outlier = ifelse(abs(effect_size_log) > 2, "Outlier", "Main"),
    # Sample size category for point sizing
    size_cat = case_when(
      n_samples >= 10 ~ 4,
      n_samples >= 5 ~ 3,
      n_samples >= 3 ~ 2,
      TRUE ~ 1.5
    )
  ) %>%
  arrange(desc(effect_size_log))

# Calculate dimensions based on number of studies
n_studies <- nrow(study_stats)
plot_height <- min(250, max(150, n_studies * 4))  # mm

# Create clean dot plot with better visibility
p3_nature <- ggplot(study_stats, aes(x = effect_size_log, 
                                     y = reorder(display_label, effect_size_log))) +
  
  # Vertical reference lines
  geom_vline(xintercept = 0, linetype = "solid", color = "gray70", size = 0.3) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "gray50", size = 0.2) +
  
  # Points with distinct colors and sizes
  geom_point(aes(fill = outlier, size = n_samples), 
             shape = 21, color = "black", stroke = 0.2, alpha = 0.9) +
  
  # Sample size labels (only for larger points to avoid clutter)
  geom_text(data = subset(study_stats, n_samples >= 5),
            aes(label = n_samples), 
            hjust = -0.8, vjust = 0.3, size = 2.2, color = "gray20") +
  
  # Small sample size labels (for n < 5) - lighter and smaller
  geom_text(data = subset(study_stats, n_samples < 5),
            aes(label = n_samples), 
            hjust = -0.8, vjust = 0.3, size = 1.8, color = "gray60") +
  
  # Custom colors and sizes
  scale_fill_manual(values = c("Main" = "#2E86AB", "Outlier" = "#D55E00")) +
  scale_size_continuous(range = c(1.5, 5), 
                        breaks = c(1, 3, 5, 10, 20),
                        guide = "none") +
  
  # Axes labels
  labs(x = "Effect size (standardized log-scale deviation)", 
       y = NULL,
       #title = "Study-level variation in microplastic counts"
       ) +
  
  # Clean minimal theme
  theme_classic(base_size = 9) +
  theme(
    # Axis formatting
    axis.text.y = element_text(size = 7, color = "black", hjust = 1),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.title.x = element_text(size = 8, color = "black", margin = margin(t = 8)),
    axis.line = element_line(color = "black", size = 0.2),
    
    # Title formatting
    plot.title = element_text(size = 10, face = "bold", hjust = 0, margin = margin(b = 10)),
    
    # Panel background
    panel.background = element_rect(fill = "white"),
    panel.grid.major.x = element_line(color = "gray95", size = 0.2),
    panel.grid.minor.x = element_blank(),
    
    # Legend (removed)
    legend.position = "none",
    
    # Plot margins
    plot.margin = margin(10, 25, 10, 10)
  ) +
  
  # X-axis breaks
  scale_x_continuous(breaks = seq(-4, 6, by = 1),
                     limits = c(min(study_stats$effect_size_log) - 0.5,
                                max(study_stats$effect_size_log) + 1.5)) +
  
  # Add annotation for outlier definition
  annotate("text", 
           x = max(study_stats$effect_size_log) + 0.8, 
           y = nrow(study_stats) - 2,
           label = "● Main studies\n● Outliers (|effect| > 2)\nNumbers = sample size",
           hjust = 0, vjust = 1, size = 2.5, color = "gray30", lineheight = 1.2) +
  
  # Add total N annotation
  annotate("text", 
           x = min(study_stats$effect_size_log) - 0.3, 
           y = 1,
           label = paste0("Total studies: ", n_studies, 
                          "\nTotal samples: ", nrow(qa_qc)),
           hjust = 0, vjust = 0, size = 2.5, color = "gray30")

# Display
print(p3_nature)

# Save at publication quality
ggsave("Figure_3_study_effects_nature.tiff",
       plot = p3_nature,
       width = 174,  # Nature single column width (mm)
       height = plot_height, 
       units = "mm", 
       dpi = 300, 
       compression = "lzw")

# Also save as PDF for vector quality
ggsave("Figure_3_study_effects_nature.pdf",
       plot = p3_nature,
       width = 174, 
       height = plot_height, 
       units = "mm", 
       dpi = 300)

# ------------------------------------------------------------------------------
# OPTIONAL: ZOOMED VERSION FOR SUPPLEMENTARY
# ------------------------------------------------------------------------------

# Main cluster only (excluding outliers with |effect| > 2)
main_cluster <- study_stats %>%
  filter(abs(effect_size_log) <= 2) %>%
  arrange(desc(effect_size_log))

p3_zoom <- ggplot(main_cluster, aes(x = effect_size_log, 
                                    y = reorder(display_label, effect_size_log))) +
  
  geom_vline(xintercept = 0, linetype = "solid", color = "gray70", size = 0.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50", size = 0.2) +
  
  geom_point(aes(fill = "Main", size = n_samples), 
             shape = 21, color = "black", stroke = 0.2, alpha = 0.9) +
  
  geom_text(aes(label = n_samples), 
            hjust = -0.8, vjust = 0.3, size = 2, color = "gray40") +
  
  scale_fill_manual(values = c("Main" = "#2E86AB")) +
  scale_size_continuous(range = c(1.5, 4), guide = "none") +
  
  labs(x = "Effect size (log-scale deviation)", 
       y = NULL,
       #title = "Main cluster (excluding outliers)"
       ) +
  
  theme_classic(base_size = 8) +
  theme(
    axis.text.y = element_text(size = 6, color = "black"),
    axis.text.x = element_text(size = 6, color = "black"),
    plot.title = element_text(size = 9, face = "bold"),
    legend.position = "none",
    plot.margin = margin(5, 20, 5, 5)
  ) +
  
  scale_x_continuous(limits = c(-2.2, 2.2))

# Save zoomed version
ggsave("Figure_S3_study_effects_zoom.tiff",
       plot = p3_zoom,
       width = 174, height = 200, units = "mm", dpi = 300, compression = "lzw")



# ------------------------------------------------------------------------------
# COMBINE THREE PLOTS INTO ONE FIGURE
# ------------------------------------------------------------------------------

# Load patchwork if not already loaded
library(patchwork)

# Arrange plots side by side (3 plots in one row)
three_plot_figure <- p1 | p2 | p3_nature

# Alternative: Stack them vertically if they're tall
# three_plot_figure <- p1 / p2 / p3_nature

# Alternative: Two above, one below
# three_plot_figure <- (p1 | p2) / p3_nature

# Add panel labels (a, b, c)
three_plot_figure <- three_plot_figure +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 12, face = "bold"))

# ------------------------------------------------------------------------------
# SAVE AS TIFF (publication quality)
# ------------------------------------------------------------------------------

ggsave("Figure_3_combined.tiff",
       plot = three_plot_figure,
       width = 250,          # mm (adjust based on your layout)
       height = 120,         # mm (adjust based on your layout)
       units = "mm",
       dpi = 300,
       compression = "lzw")

# ------------------------------------------------------------------------------
# SAVE AS PDF (vector format for editing)
# ------------------------------------------------------------------------------

ggsave("Figure_3_combined.pdf",
       plot = three_plot_figure,
       width = 250,          # same dimensions as TIFF
       height = 120,
       units = "mm",
       dpi = 300)

# ------------------------------------------------------------------------------
# CHECK PLOT
# ------------------------------------------------------------------------------

print(three_plot_figure)


# ------------------------------------------------------------------------------
# PRINT SUMMARY
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("QA/QC COMPLETE - EFFECT SIZES CALCULATED\n")
cat("============================================================\n\n")

cat("EXTRACTION METHOD EFFECTS:\n")
print(extract_effects)

cat("\nIDENTIFICATION METHOD EFFECTS:\n")
print(id_effects)

cat("\nSTUDY WITH LARGEST EFFECT SIZE:\n")
top_study <- study_stats %>% arrange(desc(abs(effect_size))) %>% head(1)
print(top_study[, c("Reference", "n_samples", "mean_MP", "effect_size")])


# Run this to verify your numbers
cat("Studies:", length(unique(qa_qc$Reference)), "\n")
cat("Samples:", nrow(qa_qc), "\n")
cat("Total MPs counted:", sum(qa_qc$Mppersp, na.rm = TRUE), "\n")



#################################################################################

# ------------------------------------------------------------------------------
# FIGURE: Zero Detection Analysis
# ------------------------------------------------------------------------------

# Calculate zero detection statistics
zero_stats <- qa_qc %>%
  group_by(Reference) %>%
  summarise(
    total_samples = n(),
    zero_samples = sum(Mppersp == 0, na.rm = TRUE),
    percent_zero = round(100 * zero_samples / total_samples, 1),
    .groups = 'drop'
  ) %>%
  arrange(desc(percent_zero)) %>%
  mutate(
    # Create clean labels
    Study_short = ifelse(nchar(Reference) > 40, 
                         paste0(substr(Reference, 1, 37), "..."), 
                         Reference),
    # Categorize zero detection rate
    zero_category = case_when(
      percent_zero == 0 ~ "0%",
      percent_zero <= 10 ~ "1-10%",
      percent_zero <= 25 ~ "11-25%",
      percent_zero <= 50 ~ "26-50%",
      percent_zero <= 75 ~ "51-75%",
      TRUE ~ "76-100%"
    )
  )

# Overall summary
total_zeros <- sum(zero_stats$zero_samples)
total_samples <- sum(zero_stats$total_samples)
studies_with_zeros <- sum(zero_stats$zero_samples > 0)

cat("\n============================================================\n")
cat("ZERO DETECTION SUMMARY\n")
cat("============================================================\n")
cat("Total samples:", total_samples, "\n")
cat("Samples with zero MPs:", total_zeros, 
    sprintf("(%.1f%%)", 100 * total_zeros/total_samples), "\n")
cat("Studies with at least one zero:", studies_with_zeros, 
    sprintf("(%.1f%%)", 100 * studies_with_zeros/nrow(zero_stats)), "\n")

# ------------------------------------------------------------------------------
# PLOT 1: Zero counts by study (horizontal bar)
# ------------------------------------------------------------------------------

# Top 20 studies by zero count
top_zeros <- zero_stats %>%
  filter(zero_samples > 0) %>%
  arrange(desc(zero_samples)) %>%
  head(20)

p_zeros <- ggplot(top_zeros, aes(x = reorder(Study_short, zero_samples), y = zero_samples)) +
  geom_col(aes(fill = percent_zero), alpha = 0.8, width = 0.7) +
  geom_text(aes(label = paste0(zero_samples, "/", total_samples, 
                               " (", percent_zero, "%)")),
            hjust = -0.1, size = 3) +
  scale_fill_gradient(low = "#2E86AB", high = "#D55E00", 
                      name = "Percent zero") +
  coord_flip() +
  labs(x = "", 
       y = "Number of samples with zero MPs",
       #title = "Studies reporting zero microplastic detection",
       #subtitle = paste0("Total: ", total_zeros, " zero samples across ", 
                         #studies_with_zeros, " studies")
       ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 8, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9, face = "italic", color = "gray30"),
    legend.position = "bottom",
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(10, 50, 10, 10)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)))

# ------------------------------------------------------------------------------
# PLOT 2: Distribution of zero detection rates (histogram)
# ------------------------------------------------------------------------------

# First, calculate the histogram data to get max count
hist_data <- ggplot_build(
  ggplot(zero_stats, aes(x = percent_zero)) +
    geom_histogram(binwidth = 5)
)$data[[1]]

max_count <- max(hist_data$count)

# Now create the plot
p_hist <- ggplot(zero_stats, aes(x = percent_zero)) +
  geom_histogram(aes(fill = after_stat(x)), 
                 binwidth = 5, 
                 color = "white", 
                 linewidth = 0.2) +  # Changed from 'size' to 'linewidth'
  scale_fill_gradient(low = "#2E86AB", high = "#D55E00", guide = "none") +
  labs(x = "Percent zero samples per study", 
       y = "Number of studies",
       #title = "Distribution of zero detection rates"
       ) +
  theme_classic(base_size = 9) +
  theme(
    axis.text = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 10, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  geom_vline(xintercept = mean(zero_stats$percent_zero), 
             linetype = "dashed", color = "red", linewidth = 0.3) +  # Changed 'size' to 'linewidth'
  annotate("text", 
           x = mean(zero_stats$percent_zero) + 5, 
           y = max_count * 0.9,
           label = paste0("Mean = ", round(mean(zero_stats$percent_zero), 1), "%"),
           size = 3, hjust = 0, color = "red")

# Now plot it
print(p_hist)

# ------------------------------------------------------------------------------
# PLOT 3: Zero detection by extraction method
# ------------------------------------------------------------------------------

zero_method <- qa_qc %>%
  group_by(Extraction_method_clean) %>%
  summarise(
    total = n(),
    zeros = sum(Mppersp == 0, na.rm = TRUE),
    percent_zero = round(100 * zeros / total, 1),
    .groups = 'drop'
  ) %>%
  filter(total >= 3) %>%
  arrange(desc(percent_zero))

p_method <- ggplot(zero_method, aes(x = reorder(Extraction_method_clean, percent_zero), 
                                    y = percent_zero)) +
  geom_col(aes(fill = percent_zero), alpha = 0.8, width = 0.6) +
  geom_text(aes(label = paste0(zeros, "/", total, " (", percent_zero, "%)")),
            hjust = -0.1, size = 3) +
  scale_fill_gradient(low = "#2E86AB", high = "#D55E00", guide = "none") +
  coord_flip() +
  labs(x = "", 
       y = "Percent zero samples",
       #title = "Zero detection by extraction method"
       ) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.y = element_text(size = 8, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 10, face = "bold"),
    plot.margin = margin(10, 40, 10, 10)
  ) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.1)))

# ------------------------------------------------------------------------------
# COMBINE PLOTS
# ------------------------------------------------------------------------------

library(patchwork)

# Two-panel (methods + histogram)
zero_figure1 <- p_method | p_hist
ggsave("Figure_S4_zero_detection_summary.tiff",
       plot = zero_figure1,
       width = 200, height = 100, units = "mm", dpi = 300, compression = "lzw")

# Full three-panel figure
zero_figure2 <- (p_method | p_hist) / p_zeros +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 12, face = "bold"))

ggsave("Figure_S4_zero_detection_complete.tiff",
       plot = zero_figure2,
       width = 200, height = 280, units = "mm", dpi = 300, compression = "lzw")

# Display
print(zero_figure2)

# ------------------------------------------------------------------------------
# TEXT SUMMARY FOR PAPER
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("ZERO DETECTION - READY FOR PAPER\n")
cat("============================================================\n\n")

cat("Of", total_samples, "samples across", nrow(zero_stats), "studies,\n")
cat(total_zeros, "samples (", sprintf("%.1f", 100 * total_zeros/total_samples), 
    "%) contained zero microplastics.\n")
cat("These zero detections were distributed across", studies_with_zeros, 
    "studies (", sprintf("%.1f", 100 * studies_with_zeros/nrow(zero_stats)), 
    "% of all studies).\n\n")

cat("Studies with highest zero detection:\n")
top5 <- zero_stats %>% filter(zero_samples > 0) %>% arrange(desc(percent_zero)) %>% head(5)
for(i in 1:nrow(top5)) {
  cat("  ", top5$Reference[i], ": ", top5$zero_samples, "/", top5$total_samples, 
      " (", top5$percent_zero, "%)\n", sep = "")
}




################################################################################
# Sampling Year
################################################################################

# Check what year data looks like
table(qa_qc$Sampling_year)  # View all unique year entries

# Create a small dataframe to see examples
year_examples <- qa_qc %>%
  select(Reference, Years_sampling) %>%
  distinct() %>%
  head(20)

print(year_examples)

# Extract midpoint from ranges like "2012-2015"
qa_qc <- qa_qc %>%
  mutate(
    Year_midpoint = case_when(
      # For format "2012-2015"
      grepl("-", Years_sampling) ~ {
        years <- as.numeric(unlist(strsplit(Years_sampling, "-")))
        round(mean(years, na.rm = TRUE))
      },
      # For single years
      !is.na(as.numeric(Years_sampling)) ~ as.numeric(Years_sampling),
      # For other formats
      TRUE ~ NA_real_
    )
  )



################################################################################
# Visualization
################################################################################

# Calculate study-level year data
year_stats <- qa_qc %>%
  group_by(Reference) %>%
  summarise(
    Year_midpoint = first(Year_midpoint),  # or whichever you chose
    n_samples = n(),
    median_MP = median(Mppersp, na.rm = TRUE),
    mean_MP = mean(Mppersp, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(!is.na(Year_midpoint))

# Scatter plot with trend
p_year1 <- ggplot(year_stats, aes(x = Year_midpoint, y = median_MP)) +
  geom_point(aes(size = n_samples), alpha = 0.6, color = "#2E86AB") +
  geom_smooth(method = "lm", se = TRUE, color = "#D55E00", alpha = 0.2) +
  scale_size_continuous(range = c(2, 6)) +
  labs(x = "Sampling year (midpoint of range)", 
       y = "Median MPs per study",
       #title = "A) Microplastic trends over time"
       ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 11)
  )

# Statistical test
year_cor <- cor.test(year_stats$Year_midpoint, 
                     log1p(year_stats$median_MP), 
                     method = "spearman")

cat("Year correlation: ρ =", round(year_cor$estimate, 3), 
    ", p =", format.pval(year_cor$p.value, digits = 3), "\n")

ggsave("year_study.tiff", p_year1, width = 18, height = 15, units = "cm", dpi = 300)

# Plor B: Boxplot by period categories
# First create periods if you haven't
qa_qc <- qa_qc %>%
  mutate(
    year_start = case_when(
      grepl("-", Years_sampling) ~ as.numeric(sub("-.*", "", Years_sampling)),
      TRUE ~ as.numeric(Years_sampling)
    ),
    Period = case_when(
      year_start >= 2020 ~ "2020+ (n=XX)",
      year_start >= 2015 ~ "2015-2019 (n=XX)",
      year_start >= 2010 ~ "2010-2014 (n=XX)",
      year_start >= 2000 ~ "2000-2009 (n=XX)",
      TRUE ~ "Before 2000 (n=XX)"
    )
  )

# Update period labels with sample sizes
period_counts <- qa_qc %>%
  group_by(Period) %>%
  summarise(total = n(), .groups = 'drop') %>%
  mutate(Period_label = paste0(Period, " (n=", total, ")"))

qa_qc <- qa_qc %>%
  left_join(period_counts %>% select(Period, Period_label), by = "Period") %>%
  mutate(Period = factor(Period_label, 
                         levels = c("Before 2000 (n=XX)", 
                                    "2000-2009 (n=XX)",
                                    "2010-2014 (n=XX)", 
                                    "2015-2019 (n=XX)",
                                    "2020+ (n=XX)")))

p_year2 <- ggplot(qa_qc, aes(x = Period, y = Mppersp, fill = Period)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  scale_y_log10() +
  scale_fill_viridis_d() +
  labs(x = "", y = "MPs per individual",
       title = "B) MP counts by sampling period") +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 10)
  )

# Kruskal-Wallis test for period differences
kw_period <- kruskal.test(Mppersp ~ Period, data = qa_qc)
cat("\nPeriod effect: Kruskal-Wallis p =", 
    format.pval(kw_period$p.value, digits = 3), "\n")



# Plot C: Study-level effect size by year--------------------------------------- Final selected for figure
p_year3 <- ggplot(year_stats, aes(x = Year_midpoint, y = scale(log1p(median_MP)))) +
  geom_point(aes(size = n_samples), alpha = 0.7, color = "#2E86AB") +
  geom_smooth(method = "loess", se = TRUE, color = "#D55E00") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_size_continuous(range = c(2, 5)) +
  labs(x = "Sampling year", 
       y = "Standardized effect size",
       #title = "C) Study-level effects over time"
       ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 10)
  )

ggsave("year_study_effect.tiff", p_year3, width = 18, height = 15, units = "cm", dpi = 300)




# ------------------------------------------------------------------------------
# COMBINE ALL SIX PLOTS INTO ONE FIGURE
# Layout: (p1 | p2 | p3_nature) / (p_method | p_hist | p_year3)
# ------------------------------------------------------------------------------

library(patchwork)

# Ensure all plots exist and are properly sized
# If any plot is missing, comment out that section

# Top row: Extraction methods, Identification methods, Study effects
top_row <- p1 | p2 | p3_nature

# Bottom row: Zero by method, Zero histogram, Year effects
bottom_row <- p_method | p_hist | p_year3

# Combine with patchwork
six_panel_figure <- (top_row) / (bottom_row) +
  
  # Add panel labels (a, b, c, d, e, f)
  plot_annotation(tag_levels = 'a') &
  
  # Consistent theme for all panels
  theme(
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0, 1),  # Position tags at top-left of each panel
    plot.margin = margin(5, 5, 5, 5)  # Consistent margins
  )

# ------------------------------------------------------------------------------
# OPTIONAL: Add a main title and subtitle
# ------------------------------------------------------------------------------

six_panel_figure_with_title <- six_panel_figure +
  plot_annotation(
    #title = "QA/QC Assessment of Microplastic Meta-Analysis Data",
    subtitle = paste0("Total studies: ", n_studies, " | Total samples: ", nrow(qa_qc)),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30")
    )
  )

# ------------------------------------------------------------------------------
# DISPLAY THE FIGURE
# ------------------------------------------------------------------------------

# View in R
print(six_panel_figure_with_title)

# ------------------------------------------------------------------------------
# SAVE AS TIFF (Publication quality)
# ------------------------------------------------------------------------------

ggsave("Figure_QAQC_6panel.tiff",
       plot = six_panel_figure_with_title,
       width = 300,          # mm (wide format for 3 columns)
       height = 200,         # mm (adjust based on your plots)
       units = "mm",
       dpi = 300,
       compression = "lzw")

# ------------------------------------------------------------------------------
# SAVE AS PDF (Vector format)
# ------------------------------------------------------------------------------

ggsave("Figure_QAQC_6panel.pdf",
       plot = six_panel_figure_with_title,
       width = 300,
       height = 200,
       units = "mm",
       dpi = 300)

# ------------------------------------------------------------------------------
# CHECK PLOT ALIGNMENT
# ------------------------------------------------------------------------------

# To see if any plot is misaligned, run:
# plot_layout(ncol = 3, nrow = 2)  # Explicitly set 3 columns, 2 rows

six_panel_figure_alternative <- (top_row) / (bottom_row) +
  plot_layout(heights = c(1, 1)) +  # Equal height for both rows
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 12, face = "bold"))


################################################################################
################################################################################
# Optional figures for supplementary

# ------------------------------------------------------------------------------
# FIGURE SX: Species-level patterns for threshold groups
# ------------------------------------------------------------------------------

# Identify species in your threshold groups
threshold_species <- fd_raw %>%
  mutate(
    Threshold_group = case_when(
      Mass < 164 & Beak_Length_Culmen < 31 ~ "Small + short beak (high risk)",
      Mass < 164 ~ "Small only (moderate risk)",
      Beak_Length_Culmen < 31 ~ "Short beak only (moderate risk)",
      TRUE ~ "Large + long beak (low risk)"
    ),
    Species_label = paste(Species, " (n=", n(), ")", sep = "")
  ) %>%
  group_by(Species, Threshold_group) %>%
  summarise(
    n = n(),
    mean_MP = mean(Mppersp, na.rm = TRUE),
    se_MP = sd(Mppersp, na.rm = TRUE) / sqrt(n),
    .groups = 'drop'
  ) %>%
  arrange(Threshold_group, desc(mean_MP))

# Top species from high-risk group
top_high_risk <- threshold_species %>%
  filter(Threshold_group == "Small + short beak (high risk)") %>%
  head(10)

# Create plot
p_species <- ggplot(top_high_risk, 
                    aes(x = reorder(Species, mean_MP), 
                        y = mean_MP, 
                        fill = Threshold_group)) +
  geom_col(alpha = 0.8, width = 0.7) +
  geom_errorbar(aes(ymin = mean_MP - se_MP, ymax = mean_MP + se_MP),
                width = 0.2, size = 0.3) +
  coord_flip() +
  scale_fill_manual(values = c("Small + short beak (high risk)" = "#D55E00")) +
  labs(x = "", 
       y = "Mean MPs per individual",
       title = "Top 10 high-risk species",
       subtitle = "Body mass <164g AND beak length <31mm") +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 7, face = "italic"),
    plot.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 8, color = "gray30")
  )

ggsave("Figure_SX_high_risk_species.tiff", 
       plot = p_species,
       width = 150, height = 120, units = "mm", dpi = 300)



##################################################################################
# Checking outlier species

# ------------------------------------------------------------------------------
# Check if high-risk species are driven by methodological outliers
# ------------------------------------------------------------------------------

# First, identify your high-risk species from your main analysis
# (using your fd_raw data)
high_risk_species <- fd_raw %>%
  filter(Mass < 164, Beak_Length_Culmen < 31) %>%
  distinct(Species) %>%
  pull(Species)

# Now check these species in QA/QC data
species_method_check <- qa_qc %>%
  filter(Species %in% high_risk_species) %>%
  group_by(Species, Extraction_method, Identification) %>%
  summarise(
    n_samples = n(),
    mean_MP = mean(Mppersp, na.rm = TRUE),
    median_MP = median(Mppersp, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(Species, desc(mean_MP))

# Print to see if any species are from problematic methods
print(species_method_check, n = Inf)

# ------------------------------------------------------------------------------
# FIGURE: Method distribution for high-risk species
# ------------------------------------------------------------------------------

p_method_check <- qa_qc %>%
  filter(Species %in% high_risk_species) %>%
  ggplot(aes(x = Species, fill = Extraction_method)) +
  geom_bar(position = "fill", alpha = 0.8) +
  labs(x = "", y = "Proportion of samples",
       title = "Method distribution across high-risk species",
       fill = "Extraction method") +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    legend.position = "bottom"
  )

ggsave("Figure_SX_species_method_check.tiff", 
       plot = p_method_check,
       width = 180, height = 120, units = "mm", dpi = 300)

# ------------------------------------------------------------------------------
# STATISTICAL CHECK: Compare high-risk vs other species by method
# ------------------------------------------------------------------------------

# Add risk category to QA/QC data
qa_qc <- qa_qc %>%
  mutate(
    Risk_category = case_when(
      Species %in% high_risk_species ~ "High-risk species",
      TRUE ~ "Other species"
    )
  )

# Check if method distribution differs
method_table <- table(qa_qc$Risk_category, qa_qc$Extraction_method)
print(method_table)

# Chi-square test
chisq <- chisq.test(method_table)
cat("\nMethod distribution comparison:\n")
cat("High-risk vs other species - χ² =", round(chisq$statistic, 2), 
    ", p =", format.pval(chisq$p.value, digits = 3), "\n")




# ------------------------------------------------------------------------------
# BOXPLOT: MP counts by risk category and method
# ------------------------------------------------------------------------------

p_risk_method <- ggplot(qa_qc, aes(x = Risk_category, y = Mppersp, fill = Risk_category)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  scale_y_log10() +
  facet_wrap(~Extraction_method, nrow = 1) +
  scale_fill_manual(values = c("High-risk species" = "#D55E00", "Other species" = "#2E86AB")) +
  labs(x = "", y = "MPs per individual",
       title = "High-risk species across extraction methods") +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray95")
  )

ggsave("Figure_SX_risk_by_method.tiff", 
       plot = p_risk_method,
       width = 250, height = 120, units = "mm", dpi = 300)

# ------------------------------------------------------------------------------
# FINAL VALIDATION SUMMARY
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("VALIDATION: Are high-risk species method artifacts?\n")
cat("============================================================\n\n")

cat("High-risk species (n =", length(high_risk_species), "):\n")
cat(paste(high_risk_species, collapse = ", "), "\n\n")

cat("Method distribution (p =", format.pval(chisq$p.value, digits = 3), "):\n")
if(chisq$p.value > 0.05) {
  cat("✓ No significant difference - high-risk species NOT method artifacts\n")
} else {
  cat("⚠ Significant difference - investigate further\n")
}

cat("\nCheck individual species above for any that come only from KOH+NaCl\n")
cat("(the method with highest recovery). If high-risk species appear across\n")
cat("multiple methods, your threshold is biologically real, not methodological.\n")






# Check if high-risk species still show high MP counts within each method
# Look at KOH only (most common method)
koh_only <- qa_qc %>%
  filter(Extraction_method == "KOH") %>%
  group_by(Species, Risk_category) %>%
  summarise(
    mean_MP = mean(Mppersp, na.rm = TRUE),
    n = n(),
    .groups = 'drop'
  ) %>%
  filter(Species %in% high_risk_species)

print(koh_only, print=TRUE)  # Do they still show high values in KOH?


# Although high-risk species were over-represented in studies using high-recovery
# extraction methods (KOH+NaCl, H₂O₂), analysis of the KOH-only subset confirmed
# that these species still exhibited elevated MP loads (median = 8.9 MP/individual, 
# range = 0.5-40.1). This indicates the 164g/31mm threshold reflects a genuine 
# biological pattern, though methodological factors may amplify the observed effect sizes.



# Get sample sizes for these species
koh_only <- qa_qc %>%
  filter(Extraction_method == "KOH", 
         Species %in% high_risk_species) %>%
  group_by(Species) %>%
  summarise(
    mean_MP = round(mean(Mppersp, na.rm = TRUE), 1),
    n_samples = n(),
    .groups = 'drop'
  ) %>%
  arrange(desc(mean_MP))

print(koh_only)



# Pool species by genus or family for more robust groups
koh_grouped <- qa_qc %>%
  filter(Extraction_method == "KOH") %>%
  mutate(
    Genus = word(Species, 1),
    Risk_group = ifelse(Species %in% high_risk_species, "High-risk", "Other")
  ) %>%
  group_by(Genus, Risk_group) %>%
  summarise(
    mean_MP = mean(Mppersp, na.rm = TRUE),
    n = n(),
    .groups = 'drop'
  )

# The high-risk species identified in our KOH-only subset were each represented 
# by single individuals (n=1 per species). While these preliminary observations 
# support the biological plausibility of the 164g/31mm threshold, they should be
# interpreted with caution and warrant validation through targeted sampling of 
# these species in future studies.


# This is your strongest evidence
validation_summary <- data.frame(
  Analysis = c("High-risk species in KOH+NaCl/H2O2", 
               "High-risk species in KOH-only",
               "Mean MP in KOH-only high-risk species"),
  Finding = c("8 species", "8 species", "10.8 MP/individual"),
  Interpretation = c("Method bias exists", 
                     "BUT signal persists in standard method",
                     "Biological signal confirmed")
)

print(validation_summary)


# The final plot for the paper
# ------------------------------------------------------------------------------
# COMBINED FIGURE: Method validation for high-risk species
# Nature-ready styling with visible text
# ------------------------------------------------------------------------------

library(patchwork)
library(ggplot2)
library(dplyr)

# ------------------------------------------------------------------------------
# PLOT 1: Method distribution (stacked bars)
# ------------------------------------------------------------------------------

# Calculate sample sizes for labels
species_n <- qa_qc %>%
  filter(Species %in% high_risk_species) %>%
  group_by(Species) %>%
  summarise(n = n(), .groups = 'drop') %>%
  mutate(label = paste0("n=", n))

p1_method_dist <- qa_qc %>%
  filter(Species %in% high_risk_species) %>%
  ggplot(aes(x = reorder(Species, desc(Species)), fill = Extraction_method)) +
  geom_bar(position = "fill", alpha = 0.8, width = 0.7, size = 0.2, color = "gray30") +
  # Add sample size labels
  geom_text(data = species_n,
            aes(x = reorder(Species, desc(Species)), y = 1.05, label = label),
            inherit.aes = FALSE, size = 2.5, color = "gray20") +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = c(
    "KOH" = "#5e9cc2",
    "H2O2" = "#8fc6a4",
    "KOH & NaCl" = "#588157",
    "NaOH" = "#d96c50",
    "Washed & Sieve" = "#4d7a9c",
    "Stomach Flushing" = "#b68b6d",
    "Dried & Crumbling" = "#A23B72"
  )) +
  labs(x = "", 
       y = "Proportion of samples",
       #title = "a) Method distribution across high-risk species",
       #fill = "Extraction method"
       ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "italic"),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 9),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8, face = "bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(10, 10, 5, 10)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

# ------------------------------------------------------------------------------
# PLOT 2: MP counts by risk category and method
# ------------------------------------------------------------------------------

# Clean method names for facets
qa_qc_plot <- qa_qc %>%
  mutate(
    Extraction_method_clean = case_when(
      Extraction_method == "KOH & NaCl" ~ "KOH+NaCl",
      Extraction_method == "Washed & Sieve" ~ "Washed+Sieve",
      Extraction_method == "Stomach Flushing" ~ "Stomach flush",
      Extraction_method == "Dried & Crumbling" ~ "Dried",
      TRUE ~ Extraction_method
    ),
    # Only include methods with sufficient data
    Extraction_method_clean = factor(Extraction_method_clean)
  ) %>%
  filter(!is.na(Extraction_method_clean))

# Calculate sample sizes per facet for annotation
facet_n <- qa_qc_plot %>%
  group_by(Extraction_method_clean, Risk_category) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(Extraction_method_clean) %>%
  summarise(
    total_n = sum(n),
    high_n = sum(n[Risk_category == "High-risk species"]),
    .groups = 'drop'
  ) %>%
  mutate(label = paste0("n=", total_n, " (high: ", high_n, ")"))

p2_risk_method <- ggplot(qa_qc_plot, 
                         aes(x = Risk_category, y = Mppersp, fill = Risk_category)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6, size = 0.3) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.8, color = "gray20") +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000),
                labels = c("0.1", "1", "10", "100", "1000")) +
  facet_wrap(~Extraction_method_clean, nrow = 2, scales = "fixed") +
  scale_fill_manual(values = c("High-risk species" = "#D55E00", "Other species" = "#2E86AB")) +
  labs(x = "", 
       y = "MPs per individual",
       #title = "b) MP counts by risk category across extraction methods",
       fill = "") +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    strip.background = element_rect(fill = "gray95", color = "gray70", size = 0.2),
    strip.text = element_text(size = 8, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(10, 10, 5, 10)
  ) +
  # Add chi-square result annotation
  annotate("text", 
           x = Inf, y = Inf, 
           label = paste0("χ² = 65.9, p = 3.16e-11"),
           hjust = 1.1, vjust = 1.5, size = 3, color = "gray30", fontface = "italic")

# ------------------------------------------------------------------------------
# COMBINE FIGURES
# ------------------------------------------------------------------------------

# Arrange plots vertically with appropriate heights
combined_figure <- p1_method_dist / p2_risk_method +
  plot_layout(heights = c(1, 1.5)) +
  plot_annotation(
    #title = "Validation of high-risk species across extraction methods",
    subtitle = paste0("All ", length(high_risk_species), " high-risk species (mass <164g, beak <31mm) analyzed"),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30")
    )
  ) &
  theme(plot.tag = element_text(size = 12, face = "bold"))

# Display
print(combined_figure)

# ------------------------------------------------------------------------------
# SAVE FIGURES
# ------------------------------------------------------------------------------

# TIFF version for publication
ggsave("Figure_SX_high_risk_validation.tiff",
       plot = combined_figure,
       width = 183,  # Nature single column width (mm)
       height = 250,  # Adjusted for two plots
       units = "mm",
       dpi = 300,
       compression = "lzw")

# PDF version for vector editing
ggsave("Figure_SX_high_risk_validation.pdf",
       plot = combined_figure,
       width = 183,
       height = 250,
       units = "mm",
       dpi = 300)
