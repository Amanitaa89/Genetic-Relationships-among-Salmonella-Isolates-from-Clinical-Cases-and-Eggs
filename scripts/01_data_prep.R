# =========================================================
# 01_data_prep.R
# Data preparation & metadata integration
# Project: Salmonella AMR (clinical vs eggs)
# =========================================================

library(dplyr)
library(readr)
library(tidyr)

# ---------------------------------------------------------
# Read metadata
# ---------------------------------------------------------
metadata <- read_csv("metadata_master.csv") %>%
  mutate(source = tolower(source))

# ---------------------------------------------------------
# Read AMRFinderPlus gene output (clinical isolates)
# ---------------------------------------------------------
amr_raw <- read_tsv("amrfinderplus-genes.tsv")

# ---------------------------------------------------------
# Keep high-confidence AMR hits
# ---------------------------------------------------------
amr_clean <- amr_raw %>%
  rename(isolate_id = Name) %>%
  mutate(isolate_id = as.character(isolate_id)) %>%
  filter(
    `Element type` == "AMR",
    `% Coverage of reference sequence` >= 90,
    `% Identity to reference sequence` >= 90
  )

# ---------------------------------------------------------
# Save cleaned AMR table
# ---------------------------------------------------------
write_csv(amr_clean, "amr_clean_clinical.csv")

