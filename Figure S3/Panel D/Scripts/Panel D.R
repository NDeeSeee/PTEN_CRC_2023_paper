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

PTEN_miss_hotspots <- c(15, 16, 17, 36, 93, 136)

# Load GERP score data ---------------------------------------------------------
GERP_score <-
  fread("../Raw data/PTEN_COSMIC_GERP_score.csv") %>%
  filter(!is.na(GERP.._RS)) %>%
  mutate(alt_nn = str_remove(Mutation.CDS, "c.")) %>%
  select(alt_nn, GERP.._RS) %>%
  rename(GERP_score = GERP.._RS) %>%
  distinct()

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

# MT-L only --------------------------------------------------------------------
PTEN_hotspots <-
  fread("../../../Common raw data/Hotspots by MS status updated.csv")

PTEN_variants_GERP <- FMI_variants_PTEN %>%
  rename(alt_nn = transcriptEffect_SV) %>%
  mutate(htsp_status = ifelse(
    codon_SV %in% filter(PTEN_hotspots, Assigned == "MT-L", htsp_status == "non")$codon,
    "non",
    "htsp"
  )) %>%
  mutate(htsp_status = as.factor(htsp_status)) %>%
  select(did, alt_nn, htsp_status) %>%
  left_join(GERP_score) %>%
  filter(!is.na(GERP_score))

# Figure S3 panel D ------------------------------------------------------------
PTEN_variants_GERP %>%
  ggplot(aes(x = htsp_status, y = GERP_score)) +
  geom_violin(trim = T,
              fill = "#34cfeb",
              alpha = .6) +
  geom_boxplot(width = 0.05, col = "blue") +
  theme_minimal() +
  coord_cartesian(ylim = c(4.5, 6))

ggsave(
  "../Charts, data, statistics/Panel D.png",
  units = "in",
  dpi = 500,
  width = 6.65,
  height = 5.94
)

# Statistics -------------------------------------------------------------------
t.test(
  filter(PTEN_variants_GERP, htsp_status == "htsp")$GERP_score,
  filter(PTEN_variants_GERP, htsp_status == "non")$GERP_score
)
# p-value = 3.829e-12
