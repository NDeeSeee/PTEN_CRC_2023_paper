# Attach requirement packages and setting WD -----------------------------------
packages_names <-
  c("tidyverse", "data.table", "readxl", "reshape2", "rstudioapi")

lapply(packages_names, require, character.only = TRUE)

setwd(dirname(getActiveDocumentContext()$path))

rename <- dplyr::rename
select <- dplyr::select
filter <- dplyr::filter
group_by <- dplyr::group_by
mutate <- dplyr::mutate

# Load FMI data ----------------------------------------------------------------
FMI_samples <-
  readRDS("../../../../Sources/FMI samples assigned.rds") %>%
  as_tibble()

FMI_variants <-
  fread("../../../../Sources/FMI variants with syn.txt") %>%
  as_tibble() %>%
  filter(codingType_SV != "synonymous") %>%
  mutate(alt = description_v2) %>%
  rename(did = deidentifiedSpecimenName) %>%
  left_join(y = select(FMI_samples, did, Assigned), by = "did")

FMI_variants_PTEN <-
  filter(FMI_variants, gene == "PTEN")

tntcont_FMI_all <- FMI_variants %>%
  select(trinucleotideAlteration_SV) %>%
  table()

tntcont_FMI_PTEN <- FMI_variants %>%
  filter(gene == "PTEN") %>%
  select(trinucleotideAlteration_SV) %>%
  table()

rbind(tntcont_FMI_PTEN, tntcont_FMI_all) %>%
  as.data.frame.matrix() %>%
  write.csv("../Charts, data, statistics/TNT context PTEN vs all FMI CRC.csv")


FMI_variants %>%
  filter(gene == "PTEN", trinucleotideAlteration_SV != "-") %>%
  select(Assigned, trinucleotideAlteration_SV) %>%
  table() %>%
  write.csv("../Charts, data, statistics/TNT context PTEN by assigned FMI CRC.csv")


