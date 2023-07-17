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

# Load 3D hotspots data and LPA with protein abundance assessments data --------
pten_clust <- read_excel("../Raw data/PTEN_AF_3Dhtsp_252_273.xlsx")

pten_lpa_abund <-
  read_excel("../../../Common raw data/Mighell Eng integrated LPA abundance.xlsx")

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
  filter(FMI_variants, gene == "PTEN", codingType_SV != "synonymous")


lpa_vs_abund_data <- FMI_variants_PTEN %>%
  filter(!is.na(codon_SV)) %>%
  select(codon_SV, proteinEffect_SV) %>%
  group_by(proteinEffect_SV) %>%
  summarise(count = n()) %>%
  left_join(distinct(select(
    FMI_variants_PTEN, proteinEffect_SV, codon_SV
  ))) %>%
  left_join(
    select(pten_lpa_abund, var, Fitness_score, abundance_score),
    by = c("proteinEffect_SV" = "var")
  ) %>%
  mutate(
    domain = case_when(
      codon_SV >= 190 & codon_SV <= 350 ~ "C2",
      codon_SV >= 14 & codon_SV <= 185 ~ "phos",
      codon_SV >= 1 & codon_SV <= 13 ~ "PBD",
      codon_SV > 350 ~ "tail"
    )
  ) %>%
  filter(!is.na(Fitness_score)) %>%
  arrange(desc(count))

lpa_vs_abund <- lpa_vs_abund_data %>%
  ggplot(aes(
    x = Fitness_score,
    y = abundance_score,
    size = count,
    color = domain
  ))


# Figure S3 panel A ------------------------------------------------------------
lpa_vs_abund_data %>%
  mutate(
    codon_SV = as.numeric(codon_SV),
    domain = ifelse(codon_SV <= 185 &
                      codon_SV >= 15, "phos", domain)
  ) %>%
  filter(!is.na(domain)) %>%
  ggplot(aes(
    x = Fitness_score,
    y = abundance_score,
    size = count,
    color = domain
  )) +
  geom_point(aes(fill = domain), alpha = 0.9, shape = 21) +
  scale_size(range = c(1, 24)) +
  theme_minimal() +
  scale_fill_manual(values = c(
    "PBD" = "#DFA0F5",
    "phos" = "#FBE983",
    "C2" = "#8EECFB",
    "tail" = "#8ADA97"
  )) +
  scale_color_manual(values = c(
    "PBD" = "#DFA0F5",
    "phos" = "#FBE983",
    "C2" = "#8EECFB",
    "tail" = "#8ADA97"
  ))

ggsave(
  "../Charts, data, statistics/Panel A.png",
  dpi = 400,
  height = 5,
  width = 5,
  units = "in"
)
