# Attach requirement packages and setting WD -----------------------------------
packages_names <-
  c("tidyverse", "data.table", "readxl", "reshape2", "rstudioapi")

lapply(packages_names, require, character.only = TRUE)

setwd(dirname(getActiveDocumentContext()$path))

library(UpSetR)
library(ComplexUpset)
library(ComplexHeatmap)

rename <- dplyr::rename
select <- dplyr::select
filter <- dplyr::filter
group_by <- dplyr::group_by
mutate <- dplyr::mutate

# Load PAD, FMI and PTM data ---------------------------------------------------
FMI_samples <-
  readRDS("../../../../Sources/FMI samples assigned.rds") %>%
  as_tibble()

PTEN_data_LOF_MTH <-
  fread("../../../Common raw data/PTEN LoF annotation.csv") %>%
  filter(Assigned == "MT-H") %>%
  select(description_v2, tot_effect) %>%
  distinct() %>%
  as_tibble()

PTEN_hotspots_MTH <-
  fread("../../../Common raw data/Hotspots by MS status updated.csv") %>%
  filter(Assigned == "MT-H", htsp_status != "non") %>%
  as_tibble() %>%
  mutate(htsp_status = "htsp") %>%
  select(codon, htsp_status) %>%
  rename(codon_SV = codon)

FMI_variants_PTEN_MTH <-
  fread("../../../../Sources/FMI variants with syn.txt") %>%
  rename(did = deidentifiedSpecimenName) %>%
  left_join(select(FMI_samples, did, Assigned)) %>%
  filter(
    Assigned == "MT-H",
    gene == "PTEN",
    description_v2 != "deletion",
    variantClass_v2 != "amplification",
    variantClass_v2 != "rearrangement",
    variantType != "RE",
    codingType_SV != "synonymous"
  ) %>%
  select(did, description_v2, codon_SV, codingType_SV) %>%
  as_tibble() %>%
  mutate(codon_SV = as.numeric(codon_SV)) %>%
  mutate(
    tot_effect = ifelse(
      description_v2 %in% filter(PTEN_data_LOF_MTH, tot_effect == "lof")$description_v2,
      "lof",
      NA
    ),
    tot_effect = ifelse(
      description_v2 %in% filter(PTEN_data_LOF_MTH, tot_effect == "wt")$description_v2,
      "wt",
      tot_effect
    ),
    tot_effect = ifelse(codingType_SV == "splice", "lof", tot_effect)
  ) %>%
  left_join(PTEN_hotspots_MTH) %>%
  mutate(htsp_status = replace_na(htsp_status, "non-htsp")) %>%
  select(-codon_SV, -codingType_SV)


table(filter(FMI_variants_PTEN_MTH, tot_effect == "lof")$htsp_status)
# htsp non-htsp 
# 407      157 
table(filter(FMI_variants_PTEN_MTH, tot_effect == "wt")$htsp_status)
# htsp non-htsp 
# 2       30
table(filter(FMI_variants_PTEN_MTH, is.na(tot_effect))$htsp_status)
# htsp non-htsp 
# 6        6 
