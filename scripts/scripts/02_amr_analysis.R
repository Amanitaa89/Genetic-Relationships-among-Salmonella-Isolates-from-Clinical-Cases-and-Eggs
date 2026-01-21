# =========================================================
# 02_amr_analysis.R
# AMR class summaries + MDR (clinical vs egg)
# =========================================================

library(dplyr)
library(readr)
library(tidyr)

# ----------------------------
# Load inputs
# ----------------------------
metadata <- read_csv("metadata_master.csv") %>%
  mutate(source = tolower(source))

amr_clean <- read_csv("amr_clean_clinical.csv")

# Egg ARGs file (gene list) + mapping table (edit name if needed)
egg_args <- read_csv("EGG_ARGs.csv")  # columns: isolate_id, Gene

# ----------------------------
# Clinical: class-level hits from AMRFinder
# ----------------------------
clin_classes <- amr_clean %>%
  distinct(isolate_id, Class) %>%
  mutate(source = "clinical")

# Clinical: gene-level hits
clin_genes <- amr_clean %>%
  distinct(isolate_id, gene = `Gene symbol`) %>%
  mutate(source = "clinical")

# ----------------------------
# Egg: gene-level hits + map to classes
# ----------------------------
egg_genes <- egg_args %>%
  transmute(isolate_id = as.character(isolate_id),
            gene = Gene,
            source = "egg") %>%
  distinct()

# Gene → Class mapping (extend if needed)
egg_gene_map <- tibble::tribble(
  ~gene,           ~Class,
  "aac(6')-Iaa",   "AMINOGLYCOSIDE",
  "fosA7",         "PHENICOL",
  "blaCMY-2",      "BETA-LACTAM",
  "sul2",          "SULFONAMIDE",
  "sul1",          "SULFONAMIDE",
  "dfrA14",        "TRIMETHOPRIM",
  "tet(A)",        "TETRACYCLINE",
  "floR_2",        "PHENICOL",
  "blaCTX-M-65",   "BETA-LACTAM",
  "qnrB19",        "QUINOLONE",
  "qnrS1",         "QUINOLONE",
  "aph(3')-Ia",    "AMINOGLYCOSIDE",
  "aph(3'')-Ib",   "AMINOGLYCOSIDE",
  "aph(4)-Ia",     "AMINOGLYCOSIDE",
  "aph(6)-Id",     "AMINOGLYCOSIDE",
  "aac(3)-Iva",    "AMINOGLYCOSIDE",
  "ant(3'')-Ia",   "AMINOGLYCOSIDE"
)

egg_classes <- egg_genes %>%
  left_join(egg_gene_map, by = "gene") %>%
  filter(!is.na(Class)) %>%
  distinct(isolate_id, Class) %>%
  mutate(source = "egg")

# ----------------------------
# Combine class-level table
# ----------------------------
amr_classes_all <- bind_rows(clin_classes, egg_classes)

# Class distribution (% within source)
amr_class_pct <- amr_classes_all %>%
  distinct(isolate_id, source, Class) %>%
  count(source, Class, name = "n") %>%
  group_by(source) %>%
  mutate(percent = round(n / sum(n) * 100, 1)) %>%
  ungroup()

write_csv(amr_class_pct, "AMR_class_distribution_by_source.csv")

# ----------------------------
# MDR (>=3 classes)
# ----------------------------
mdr_all <- amr_classes_all %>%
  distinct(isolate_id, source, Class) %>%
  count(isolate_id, source, name = "n_classes") %>%
  mutate(MDR = n_classes >= 3)

mdr_summary <- mdr_all %>%
  count(source, MDR) %>%
  tidyr::pivot_wider(names_from = MDR, values_from = n, values_fill = 0) %>%
  rename(non_MDR = `FALSE`, MDR = `TRUE`) %>%
  mutate(total = non_MDR + MDR,
         MDR_percent = round(MDR / total * 100, 1))

write_csv(mdr_summary, "MDR_summary_by_source.csv")
