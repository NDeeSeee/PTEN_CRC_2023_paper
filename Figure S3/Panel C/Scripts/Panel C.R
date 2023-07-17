# Attach requirement packages and setting WD -----------------------------------
packages_names <- c("readxl",
                    "data.table",
                    "reshape2",
                    "tidyverse",
                    "rstudioapi",
                    "readxl")

lapply(packages_names, require, character.only = TRUE)

setwd(dirname(getActiveDocumentContext()$path))

rename <- dplyr::rename
select <- dplyr::select
filter <- dplyr::filter
group_by <- dplyr::group_by
mutate <- dplyr::mutate

# Load PAD, FMI and PTM data ---------------------------------------------------
PAD_samples <-
  fread("../../../Common raw data/PAD samples.csv", sep = ",") %>%
  as_tibble()
PAD_variants <-
  fread("../../../Common raw data/PAD variants.csv", sep = ",") %>%
  as_tibble() %>%
  rename(codon_SV = codon)

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

# Load PTEN loss of function data ----------------------------------------------
# From Mighell paper
mighell_data <-
  read_xlsx("../../../Common raw data/Mighell Eng integrated LPA abundance.xlsx") %>%
  rename(
    alt = var,
    VAMP_score = abundance_score,
    VAMP_effect = abundance_score_class,
    MAVE_score = Fitness_score,
    MAVE_effect = Fitness_score_class
  ) %>%
  mutate(
    VAMP_effect = str_replace_all(VAMP_effect, "truncation-like", "lof"),
    VAMP_effect = str_replace_all(VAMP_effect, "hypomorphic", "lof"),
    VAMP_effect = str_replace_all(VAMP_effect, "wildtype-like", "wt"),
    MAVE_effect = str_replace_all(MAVE_effect, "truncation-like", "lof"),
    MAVE_effect = str_replace_all(MAVE_effect, "hypomorphic", "lof"),
    MAVE_effect = str_replace_all(MAVE_effect, "wildtype-like", "wt")
  )

# Add MAVE scores to PTEN variants data ----------------------------------------
PTEN_variants <- filter(FMI_variants, gene == "PTEN")

PTEN_variants_score <- PTEN_variants %>%
  filter(
    variantClass_v2 == "point",
    !str_detect(alt, "del"),
    !str_detect(alt, "ins"),
    !str_detect(alt, "_"),
    !str_detect(alt, "\\*")
  ) %>%
  select(alt, Assigned, codingType_SV) %>%
  left_join(select(mighell_data, alt, MAVE_score, VAMP_score)) %>%
  mutate(MAVE_score = replace_na(MAVE_score, 0),
         VAMP_score = replace_na(VAMP_score, 1))

# Load signature probabilities data --------------------------------------------
pred_MTL_prop <-
  fread("../../../Common raw data/MTL probs predicted.csv") %>%
  rename(alt = amino_acid_change, probability = rel_prob) %>%
  group_by(alt) %>%
  summarise(probability = sum(probability)) %>%
  filter(!str_detect(alt, "\\*")) %>%
  mutate(prop_pred = round(prop.table(probability) * 10000)) %>%
  select(-probability) %>%
  left_join(select(mighell_data, alt, MAVE_score, VAMP_score)) %>%
  mutate(MAVE_score = replace_na(MAVE_score, 0),
         VAMP_score = replace_na(VAMP_score, 1))

# Predicted data ---------------------------------------------------------------
VAMP_score_pred_MTL <- c()
for (i in 1:nrow(pred_MTL_prop)) {
  VAMP_score_pred_MTL <- c(VAMP_score_pred_MTL,
                           rep(pred_MTL_prop$VAMP_score[i],
                               pred_MTL_prop$prop_pred[i]))
}

pred_MTL_prop_VAMP <- tibble(VAMP_score = VAMP_score_pred_MTL,
                             data = "pred")

# Observed data ---------------------------------------------------------------
PTEN_lof_data_MTL <-
  filter(PTEN_variants_score, Assigned == "MT-L") %>%
  select(alt, MAVE_score, VAMP_score) %>%
  mutate(data = "real")

Panel_C_data <- bind_rows(pred_MTL_prop_VAMP,
                          select(PTEN_lof_data_MTL, VAMP_score, data)) %>%
  mutate(data = as.factor(data))

# Draw the chart Panel C -------------------------------------------------------
Panel_C_data %>%
  ggplot(aes(x = VAMP_score, fill = data)) +
  geom_density(alpha = .6) +
  theme_minimal() +
  scale_fill_manual(values = c("real" = "#6A3DB3", "pred" = "#EB5B51"))

# Save chart data --------------------------------------------------------------
group_by(PTEN_lof_data_MTL, alt) %>%
  summarise(
    count = n(),
    MAVE_score = unique(MAVE_score)[1],
    VAMP_score = unique(VAMP_score)[1],
    data = "observed"
  ) %>%
  bind_rows(mutate(rename(pred_MTL_prop, count = prop_pred),
                   data = "predicted")) %>%
  write.csv("../Charts, data, statistics/Panel C observed and predicted data.csv",
            row.names = F)

Panel_C_data %>%
  write.csv("../Charts, data, statistics/Panel C chart data.csv",
            row.names = F)

# Save image -------------------------------------------------------------------
ggsave(
  "../Charts, data, statistics/Panel C.png",
  dpi = 400,
  height = 5,
  width = 6.5,
  units = "in"
)
