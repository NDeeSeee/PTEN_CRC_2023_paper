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
  mutate(alt = description_v2) %>%
  rename(did = deidentifiedSpecimenName) %>%
  left_join(y = select(FMI_samples, did, Assigned), by = "did")

# Custom gene set --------------------------------------------------------------
our_genes <-
  c("KRAS", "NRAS", "HRAS", "APC", "SMAD4", "PIK3CA", "TP53", "PTEN")

dNdS_summary <- tibble(
  our_genes = c("KRAS", "NRAS", "HRAS", "APC", "SMAD4", "PIK3CA", "TP53", "PTEN"),
  syn = 0,
  nontrunc = 0,
  trunc = 0
)

for (i in 1:8) {
  dNdS_summary$syn[i] <- FMI_variants %>%
    filter(gene_SV == dNdS_summary$our_genes[i],
           codingType_SV == "synonymous") %>%
    nrow()
  dNdS_summary$nontrunc[i] <- FMI_variants %>%
    filter(
      gene_SV == dNdS_summary$our_genes[i],
      variantClass_v2 == "point",
      !codingType_SV == "splice"
    ) %>%
    nrow()
  dNdS_summary$trunc[i] <- FMI_variants %>%
    filter(
      gene_SV == dNdS_summary$our_genes[i],
      variantClass_v2 == "truncation",
      !codingType_SV == "splice"
    ) %>%
    nrow()
}

dNdS_summary <- dNdS_summary %>%
  group_by(our_genes) %>%
  mutate(dNdS_nontr = round((nontrunc - syn) / syn, 3),
         dNdS_tr = round(trunc / syn, 2)) %>%
  mutate(
    dNdS_nontrL = round(prop.test(syn, nontrunc - syn)$conf.int[1], 3),
    dNdS_nontrU = round(prop.test(syn, nontrunc - syn)$conf.int[2], 3)
  ) %>%
  mutate(
    dNdS_nontrL = prop.test(syn, (nontrunc - syn))$conf.int[1],
    dNdS_nontrU = prop.test(syn, (nontrunc - syn))$conf.int[2]
  ) %>%
  ungroup()

dNdS_summary %>%
  write.csv("../Charts, data, statistics/Panel A data.csv", row.names = F)
