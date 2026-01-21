# =========================================================
# 03_figures.R
# Figures for manuscript (AMR classes + MDR + gene prevalence)
# =========================================================

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(scales)

# ----------------------------
# Load outputs from analysis
# ----------------------------
amr_class_pct <- read_csv("AMR_class_distribution_by_source.csv")
mdr_summary   <- read_csv("MDR_summary_by_source.csv")

# ---------------------------------------------------------
# Figure 1: 100% stacked bar of AMR classes by source
# ---------------------------------------------------------
df_class <- amr_class_pct %>%
  group_by(source) %>%
  mutate(p = n / sum(n)) %>%
  ungroup()

p_fig1 <- ggplot(df_class, aes(x = source, y = p, fill = Class)) +
  geom_col() +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_brewer(palette = "Set3") +
  labs(x = NULL, y = "Share of AMR classes (100%)", fill = "AMR class") +
  theme_minimal(base_size = 14)

ggsave("Figure1_AMRclass_100stack.png", p_fig1, width = 7, height = 4, dpi = 300)

# ---------------------------------------------------------
# Figure 3: MDR prevalence by source (side-by-side bars)
# ---------------------------------------------------------
mdr_long <- mdr_summary %>%
  transmute(
    source,
    non_MDR = non_MDR / total,
    MDR = MDR / total
  ) %>%
  pivot_longer(-source, names_to = "status", values_to = "p")

p_fig3 <- ggplot(mdr_long, aes(x = source, y = p, fill = status)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = percent(p, accuracy = 0.1)),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 4) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.05)) +
  scale_fill_brewer(palette = "Set3") +
  labs(x = NULL, y = "Proportion of isolates", fill = NULL) +
  theme_minimal(base_size = 14)

ggsave("Figure3_MDR_by_source.png", p_fig3, width = 6, height = 4, dpi = 300)

# ---------------------------------------------------------
# Supplementary (optional): gene distribution plots
# NOTE: requires gene_pct_top object saved separately if you want to reproduce
# If you want, we can add saving/loading gene_pct_top in 02_amr_analysis.R
# ---------------------------------------------------------
