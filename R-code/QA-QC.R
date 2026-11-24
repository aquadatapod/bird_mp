################################################################################
# Figure S2: METHODOLOGICAL HETEROGENEITY IN MICROPLASTIC LOAD
#
# INTERPRETIVE NOTE :
# Extraction/identification method and season are each essentially constant
# WITHIN a study. Their fixed-effect tests are therefore comparisons of
# study-level group means, properly weighted via the study random intercept
# - but they CANNOT be interpreted as causal ("method X extracts more
# efficiently than method Y"), since method is confounded with every other
# unmeasured difference between the studies that used it. Report only the
# association, with this caveat stated explicitly.
################################################################################

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(patchwork)

qa_data <- read.csv(here("data", "QA-QC-data.csv"), stringsAsFactors = FALSE, na.strings = c("", "NA"))

df <- qa_data %>%
  filter(!is.na(Mppersp), Mppersp >= 0, !is.na(study_id)) %>%
  mutate(study_id = factor(study_id), logMP = log1p(Mppersp))

cat(sprintf("Analysis records: %d, studies: %d\n", nrow(df), n_distinct(df$study_id)))

# ==========================================================================
# PREPARE FACTOR — unchanged from original: all categories kept for the
# descriptive plot, only categories with >=2 studies enter the model
# ==========================================================================

prepare_factor <- function(data, variable) {
  d <- data %>%
    filter(!is.na(.data[[variable]]), .data[[variable]] != "") %>%
    mutate(factor_value = as.character(.data[[variable]]))
  
  counts <- d %>%
    group_by(factor_value) %>%
    summarise(n_records = n(), n_studies = n_distinct(study_id), .groups = "drop")
  
  valid_levels <- counts %>% filter(n_studies >= 2) %>% pull(factor_value)
  
  d <- d %>% left_join(counts, by = "factor_value") %>%
    mutate(tested = factor_value %in% valid_levels)
  
  d_model <- d %>% filter(tested) %>% mutate(factor_value = factor(factor_value))
  
  list(plot = d, model = d_model, counts = counts)
}

extraction     <- prepare_factor(df, "Extraction_method")
identification <- prepare_factor(df, "Identification")
season         <- prepare_factor(df, "Season")

# ==========================================================================
# MIXED-EFFECTS MODELS — unchanged
# ==========================================================================

m_extraction     <- lmer(logMP ~ factor_value + (1 | study_id), data = extraction$model, REML = FALSE)
m_identification <- lmer(logMP ~ factor_value + (1 | study_id), data = identification$model, REML = FALSE)
m_season         <- lmer(logMP ~ factor_value + (1 | study_id), data = season$model, REML = FALSE)

p_extraction     <- anova(m_extraction)$`Pr(>F)`[1]
p_identification <- anova(m_identification)$`Pr(>F)`[1]
p_season         <- anova(m_season)$`Pr(>F)`[1]

emm_to_df <- function(model) {
  emm <- emmeans(model, ~ factor_value)
  as.data.frame(emm) %>% mutate(MP_estimate = expm1(emmean), MP_lower = expm1(lower.CL), MP_upper = expm1(upper.CL))
}

emm_extraction_df     <- emm_to_df(m_extraction)
emm_identification_df <- emm_to_df(m_identification)
emm_season_df         <- emm_to_df(m_season)

# ==========================================================================
# Tested-vs-descriptive-only
# distinction, single muted color, no overlapping labels
# ==========================================================================

make_method_plot <- function(data, emm_data, x_title, p_value, panel_tag) {
  
  level_order <- data %>%
    group_by(factor_value) %>%
    summarise(median_MP = median(Mppersp, na.rm = TRUE), .groups = "drop") %>%
    arrange(median_MP) %>% pull(factor_value)
  
  data$factor_value <- factor(data$factor_value, levels = level_order)
  emm_data$factor_value <- factor(emm_data$factor_value, levels = level_order)
  
  labels_df <- data %>%
    distinct(factor_value, n_records, n_studies, tested) %>%
    mutate(label = sprintf("%s (n=%d, %d stud%s)%s",
                           factor_value, n_records, n_studies,
                           ifelse(n_studies == 1, "y", "ies"),
                           ifelse(tested, "", " \u2014 not tested*")))
  
  label_map <- setNames(labels_df$label, labels_df$factor_value)
  
  p_label <- ifelse(p_value < 0.001, "Mixed model: p < 0.001",
                    paste0("Mixed model: p = ", format.pval(p_value, digits = 2)))
  
  ggplot(data, aes(x = factor_value, y = Mppersp, fill = tested)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.55, linewidth = 0.3) +
    geom_jitter(width = 0.12, alpha = 0.25, size = 0.7, color = "grey30") +
    geom_point(data = emm_data, aes(x = factor_value, y = MP_estimate),
               inherit.aes = FALSE, shape = 18, size = 2.2, color = "black") +
    geom_errorbar(data = emm_data, aes(x = factor_value, ymin = MP_lower, ymax = MP_upper),
                  inherit.aes = FALSE, width = 0.12, linewidth = 0.4, color = "black") +
    scale_fill_manual(values = c("TRUE" = "#2E86AB", "FALSE" = "grey80"), guide = "none") +
    scale_x_discrete(labels = label_map) +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10),
      breaks = c(0, 1, 10, 100, 1000, 10000),
      labels = c("0", "1", "10", "100", "1,000", "10,000")
    ) +
    coord_flip() +
    labs(x = "", y = "MPs per individual", tag = panel_tag,
         subtitle = paste0(p_label, "  |  grey bars = single-study, descriptive only (*not entered in model)")) +
    theme_classic(base_size = 8.5) +
    theme(
      axis.text.y = element_text(size = 6.5, color = "black"),
      axis.text.x = element_text(size = 7, color = "black"),
      plot.tag = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 6.5, face = "italic", color = "grey30"),
      plot.margin = margin(4, 8, 4, 4)
    )
}

p_extraction_plot     <- make_method_plot(extraction$plot,     emm_extraction_df,     "Extraction method",     p_extraction,     "a")
p_identification_plot <- make_method_plot(identification$plot, emm_identification_df, "Identification method", p_identification, "b")
p_season_plot         <- make_method_plot(season$plot,         emm_season_df,         "Sampling season",       p_season,         "c")

# Stack vertically instead of side-by-side - horizontal bars need width,
# not height, so a 1-column layout reads far better than 3 columns
method_figure <- p_extraction_plot / p_identification_plot / p_season_plot

print(method_figure)

ggsave(here("outputs", "Figure_Methodological_Heterogeneity_v2.pdf"), method_figure,
       width = 16, height = 24, units = "cm", device = "pdf", useDingbats = FALSE)
ggsave(here("outputs", "Figure_Methodological_Heterogeneity_v2.tiff"), method_figure,
       width = 16, height = 24, units = "cm", dpi = 600, compression = "lzw")

cat("\n================ CORRECT INTERPRETATION (use this wording) ================\n")
cat(sprintf("Extraction method: p = %.3f - ", p_extraction))
cat(ifelse(p_extraction < 0.05,
           "associated with reported MP burden; NOTE method is confounded with\nstudy identity (nearly constant within study), so this reflects between-\nstudy differences broadly, not a controlled test of extraction efficiency.\n",
           "no significant association.\n"))
cat(sprintf("Identification method: p = %.3f - %s\n", p_identification,
            ifelse(p_identification < 0.05, "significant association.", "no significant association - detection method does not explain between-study heterogeneity.")))
cat(sprintf("Season: p = %.3f - %s\n", p_season,
            ifelse(p_season < 0.05, "significant association.", "no significant association - sampling season does not explain between-study heterogeneity.")))
cat("\nSuggested one-line Results/Supplement statement:\n")
cat(sprintf(
  "\"Reported microplastic burden differed across extraction methods (mixed\nmodel, p = %.3f), though this association cannot be separated from other\nunmeasured differences between studies using different protocols; no\nassociation was found with identification method (p = %.3f) or sampling\nseason (p = %.3f), suggesting these methodological factors do not\nsystematically bias the compiled dataset.\"\n",
  p_extraction, p_identification, p_season))



###############################################################################
# Cross checking
###############################################################################
extraction$plot %>% filter(factor_value %in% c("NaCl", "KOH & H2O2")) %>% select(factor_value, Mppersp)
identification$plot %>% filter(factor_value == "Stereomicroscope, Raman & FTIR") %>% select(factor_value, Mppersp)
season$plot %>% filter(factor_value == "Winter/Summer") %>% select(factor_value, Mppersp)





unique(emm_extraction_df$factor_value)      # should list ONLY the blue-bar categories
unique(emm_identification_df$factor_value)  # same check
unique(emm_season_df$factor_value)          # same check



# Did you get the warning, and how many records does it affect?
sum(qa_data$Mppersp == 0, na.rm = TRUE)

# Per panel - how many zero records exist in each factor's plotting data
sum(extraction$plot$Mppersp == 0, na.rm = TRUE)
sum(identification$plot$Mppersp == 0, na.rm = TRUE)
sum(season$plot$Mppersp == 0, na.rm = TRUE)