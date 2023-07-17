# Attach requirement packages and setting WD -----------------------------------
packages_names <- c("readxl", "data.table", "reshape2", "tidyverse",
                    "rstudioapi")

lapply(packages_names, require, character.only = TRUE)

setwd(dirname(getActiveDocumentContext()$path))

rename = dplyr::rename
select = dplyr::select
filter = dplyr::filter
group_by = dplyr::group_by
mutate = dplyr::mutate

# Load PAD, FMI and PTM data ---------------------------------------------------
PAD_samples <-
  fread("../../../Common raw data/PAD samples.csv", sep = ",") %>%
  as_tibble()
PAD_variants <-
  fread("../../../Common raw data/PAD variants.csv", sep = ",") %>%
  as_tibble() %>%
  rename(codon_SV = codon) %>%
  filter(msi == "MSS")


FMI_samples <-
  readRDS("../../../../Sources/FMI samples assigned.rds") %>%
  as_tibble() %>%
  mutate(did = as.character(did)) %>%
  filter(Assigned == "MT-L") %>%
  select(-PTEN_count, -PTEN_presence)

FMI_variants <-
  fread("../../../../Sources/FMI variants with syn.csv") %>%
  as_tibble() %>%
  filter(codingType_SV != "synonymous") %>%
  mutate(alt = description_v2) %>%
  filter(Assigned == "MT-L")

PTEN_LoF_data_full <-
  fread("../../../Common raw data/LoF annotated PTEN alterations.csv") %>%
  as_tibble()

PTEN_LoF_data <-
  fread("../../../Common raw data/PTEN LoF annotation.csv") %>%
  as_tibble()

# PAD data ---------------------------------------------------------------------
PAD_variants_grouped <- PAD_variants %>%
  mutate(
    pten_lof = ifelse(
      gene == "PTEN" &
        alt %in% filter(PTEN_LoF_data_full, consensus_VM_KB == "lof")$alt,
      1,
      NA
    ),
    pten_lof = ifelse(
      gene == "PTEN" &
        alt %in% filter(PTEN_LoF_data_full, consensus_VM_KB == "wt")$alt,
      0,
      pten_lof
    ),
    pten_lof = ifelse(gene == "PTEN" & str_detect(alt, "\\*"),
                      1,
                      pten_lof)
  ) %>%
  group_by(sample_id) %>%
  summarise(PTEN_LOF = sign(sum(pten_lof == 1, na.rm = T)),
            PTEN_WT = sign(sum(pten_lof == 0, na.rm = T))) %>%
  ungroup()

PAD_co_occurrence <-
  fread("../../../Common raw data/PAD co-occurrence.csv") %>%
  left_join(select(PAD_samples, sample_id, msi), by = "sample_id") %>%
  filter(msi == "MSS") %>%
  select(
    sample_id,
    PTEN_ANY,
    PIK3R1_ANY_OUT,
    PIK3CA_ANY_OUT,
    PTEN_DEL,
    TP53_DEL,
    TP53_MUT,
    KRAS_ANY_OUT,
    SMAD4_ANY_OUT,
    BRAF_ANY_OUT,
    APC_ANY_OUT
  ) %>%
  rename(
    PIK3R1_OUT_ANY = PIK3R1_ANY_OUT,
    PIK3CA_OUT_ANY = PIK3CA_ANY_OUT,
    TP53_OUT_DEL = TP53_DEL,
    TP53_OUT_MUT = TP53_MUT,
    KRAS_OUT_ANY = KRAS_ANY_OUT,
    SMAD4_OUT_ANY = SMAD4_ANY_OUT,
    BRAF_OUT_ANY = BRAF_ANY_OUT,
    APC_OUT_ANY = APC_ANY_OUT
  ) %>%
  left_join(PAD_variants_grouped, by = "sample_id") %>%
  mutate(PTEN_LOF = replace_na(PTEN_LOF, 0),
         PTEN_WT = replace_na(PTEN_WT, 0)) %>%
  filter(!if_any(
    .cols = everything(),
    .fns = function(x) {
      is.na(x)
    }
  ))

write.csv(
  PAD_co_occurrence,
  "../Charts, data, statistics/PAD co-occurrence data.csv",
  row.names = F
)

# Chart data -------------------------------------------------------------------
PAD_co_occurrence_chart_data <- PAD_co_occurrence %>%
  pivot_longer(cols = matches("PTEN"),
               names_to = "PTEN_type",
               values_to = "PTEN_value") %>%
  pivot_longer(cols = matches("OUT"),
               names_to = "OUT_type",
               values_to = "OUT_value") %>%
  mutate(PTEN_type = as.factor(PTEN_type),
         OUT_type = as.factor(OUT_type)) %>%
  group_by(PTEN_type, OUT_type) %>%
  summarise(
    both = sum(PTEN_value == 1 & OUT_value == 1),
    PTEN_only = sum(PTEN_value == 1 & OUT_value == 0),
    OUT_only = sum(PTEN_value == 0 & OUT_value == 1),
    neither = sum(PTEN_value == 0 & OUT_value == 0)
  ) %>%
  group_by(PTEN_type, OUT_type) %>%
  mutate(fisher_test_res = paste(as.character(round(
    unlist(fisher.test(matrix(
      c(both, PTEN_only, OUT_only, neither),
      ncol = 2
    ))[c(1, 2, 3)]), 4
  )), collapse = ";")) %>%
  ungroup() %>%
  separate(fisher_test_res,
           c("p_val", "log2_lower", "log2_upper", "log2_odds_ratio"),
           sep = ";") %>%
  mutate(across(
    contains("log2"),
    .fns = function(x) {
      log2(as.numeric(x))
    }
  ),
  combination = as.factor(paste(PTEN_type, OUT_type, sep = " ")))

PAD_co_occurrence_chart_data %>% 
  write.csv("../Charts, data, statistics/Panel C PAD data.csv", row.names = F)

# Save PAD numbers version data ------------------------------------------------
PAD_co_occurrence_chart_data %>%
  select(-both, -PTEN_only, -combination, -OUT_only, -neither, -p_val) %>%
  mutate(
    OUT_type = str_remove(OUT_type, "_OUT_ANY"),
    OUT_type = str_remove(OUT_type, "_OUT"),
    PTEN_type = str_remove_all(PTEN_type, "PTEN_")
  ) %>%
  pivot_wider(
    names_from = PTEN_type,
    values_from = c(log2_odds_ratio, log2_upper, log2_lower)
  ) %>%
  pivot_longer(cols = matches("log2")) %>%
  mutate(name = str_remove(name, "log2_")) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  mutate(
    upper_LOF = upper_LOF - odds_ratio_LOF,
    lower_LOF = odds_ratio_LOF - lower_LOF,
    upper_WT = upper_WT - odds_ratio_WT,
    lower_WT = odds_ratio_WT - lower_WT,
    upper_DEL = upper_DEL - odds_ratio_DEL,
    lower_DEL = odds_ratio_DEL - lower_DEL,
    upper_ANY = upper_ANY - odds_ratio_ANY,
    lower_ANY = odds_ratio_ANY - lower_ANY
  ) %>%
  arrange(OUT_type) %>%
  mutate(
    LOF = seq(1, 23, 3),
    WT = seq(2.4, 24, 3),
    ANY = seq(0.3, 22, 3),
    DEL = seq(1.7, 23, 3)
  ) %>%
  select(
    c(
      OUT_type,
      "LOF",
      "odds_ratio_LOF",
      "upper_LOF",
      "lower_LOF",
      "WT",
      "odds_ratio_WT",
      "upper_WT",
      "lower_WT",
      "DEL",
      "odds_ratio_DEL",
      "upper_DEL",
      "lower_DEL",
      "ANY",
      "odds_ratio_ANY",
      "upper_ANY",
      "lower_ANY",
    )
  ) %>%
  mutate(across(
    .cols = where(is.numeric),
    .fns = function(x) {
      round(x, 2)
    }
  )) %>%
  write.csv("Panel C PAD chart data, Numbers edition.csv", row.names = F)


# FMI data ---------------------------------------------------------------------

FMI_co_occurrence <- FMI_variants %>%
  mutate(
    PTEN_lof = ifelse(
      gene == "PTEN" &
        alt %in% filter(PTEN_LoF_data, tot_effect == "lof")$alt,
      "lof",
      NA
    ),
    PTEN_lof = ifelse(
      gene == "PTEN" &
        alt %in% filter(PTEN_LoF_data, tot_effect == "wt")$alt,
      "wt",
      PTEN_lof
    )
  ) %>%
  group_by(did) %>%
  summarise(
    PTEN_LOF = sign(sum(gene == "PTEN" &
                          (PTEN_lof == "lof"),
                        na.rm = T)),
    PTEN_WT = sign(sum(gene == "PTEN" &
                         (PTEN_lof == "wt"),
                       na.rm = T)),
    PTEN_DEL = sign(sum(
      gene == "PTEN" & variantClass_v2 == "deletion",
      na.rm = T
    )),
    TP53_OUT_MUT = sign(sum(
      gene == "TP53" &
        (variantClass_v2 == "point" |
           variantClass_v2 == "truncation"),
      na.rm = T
    )),
    TP53_OUT_DEL = sign(sum(
      gene == "TP53" & variantClass_v2 == "deletion",
      na.rm = T
    )),
    PIK3CA_OUT_ANY = sign(sum(gene == "PIK3CA", na.rm = T)),
    PTEN_ANY = sign(sum(gene == "PTEN", na.rm = T)),
    SMAD4_OUT_ANY = sign(sum(gene == "SMAD4", na.rm = T)),
    APC_OUT_ANY = sign(sum(gene == "PIK3CA", na.rm = T)),
    KRAS_OUT_ANY = sign(sum(gene == "KRAS", na.rm = T))
  ) %>%
  left_join(x = filter(FMI_samples, Assigned == "MT-L")) %>%
  mutate(across(matches("OUT|PTEN"), function(x) {
    x <- replace_na(x, 0)
  }))

FMI_co_occurrence %>%
  write.csv("../Charts, data, statistics/Panel C FMI full data.csv",
            row.names = F)

# Chart data -------------------------------------------------------------------
FMI_co_occurrence_chart_data <-
  FMI_co_occurrence %>%
  pivot_longer(cols = matches("PTEN"),
               names_to = "PTEN_type",
               values_to = "PTEN_value") %>%
  pivot_longer(cols = matches("OUT"),
               names_to = "OUT_type",
               values_to = "OUT_value") %>%
  mutate(PTEN_type = as.factor(PTEN_type),
         OUT_type = as.factor(OUT_type)) %>%
  group_by(PTEN_type, OUT_type) %>%
  summarise(
    both = sum(PTEN_value == 1 & OUT_value == 1),
    PTEN_only = sum(PTEN_value == 1 & OUT_value == 0),
    OUT_only = sum(PTEN_value == 0 & OUT_value == 1),
    neither = sum(PTEN_value == 0 & OUT_value == 0)
  ) %>%
  group_by(PTEN_type, OUT_type) %>%
  mutate(fisher_test_res = paste(as.character(round(
    unlist(fisher.test(matrix(
      c(both, PTEN_only, OUT_only, neither),
      ncol = 2
    ))[c(1, 2, 3)]), 4
  )), collapse = ";")) %>%
  ungroup() %>%
  separate(fisher_test_res,
           c("p_val", "log2_lower", "log2_upper", "log2_odds_ratio"),
           sep = ";") %>%
  mutate(across(
    contains("log2"),
    .fns = function(x) {
      log2(as.numeric(x))
    }
  ),
  combination = as.factor(paste(PTEN_type, OUT_type, sep = " ")))

FMI_co_occurrence_chart_data %>% 
  write.csv("../Charts, data, statistics/Panel C FMI data.csv", row.names = F)

# Both data --------------------------------------------------------------------
BOTH_co_occurrence <-
  bind_rows(FMI_co_occurrence, rename(PAD_co_occurrence, did = sample_id)) %>%
  select(-sex, -site, -age, -tmb, -msi, -baitSet, -Assigned) %>%
  select(matches("OUT|PTEN")) %>%
  select(!matches("PIK3R1|BRAF|presence|count"))

BOTH_co_occurrence_chart_data <- BOTH_co_occurrence %>%
  pivot_longer(cols = matches("PTEN"),
               names_to = "PTEN_type",
               values_to = "PTEN_value") %>%
  pivot_longer(cols = matches("OUT"),
               names_to = "OUT_type",
               values_to = "OUT_value") %>%
  mutate(PTEN_type = as.factor(PTEN_type),
         OUT_type = as.factor(OUT_type)) %>%
  group_by(PTEN_type, OUT_type) %>%
  summarise(
    both = sum(PTEN_value == 1 & OUT_value == 1, na.rm = T),
    PTEN_only = sum(PTEN_value == 1 & OUT_value == 0, na.rm = T),
    OUT_only = sum(PTEN_value == 0 & OUT_value == 1, na.rm = T),
    neither = sum(PTEN_value == 0 & OUT_value == 0, na.rm = T)
  ) %>%
  group_by(PTEN_type, OUT_type) %>%
  mutate(fisher_test_res = paste(as.character(round(
    unlist(fisher.test(matrix(
      c(both, PTEN_only, OUT_only, neither),
      ncol = 2
    ))[c(1, 2, 3)]), 4
  )), collapse = ";")) %>%
  ungroup() %>%
  separate(fisher_test_res,
           c("p_val", "log2_lower", "log2_upper", "log2_odds_ratio"),
           sep = ";") %>%
  mutate(across(
    contains("log2"),
    .fns = function(x) {
      log2(as.numeric(x))
    }
  ),
  combination = as.factor(paste(PTEN_type, OUT_type, sep = " ")))

# Plotting ---------------------------------------------------------------------
ggplot(
  FMI_co_occurrence_chart_data,
  aes(x = combination, y = log2_odds_ratio, col = PTEN_type)
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = log2_lower, ymax = log2_upper),
    width = .2,
    position = position_dodge(0.05),
    linewidth = 1
  ) +
  geom_hline(yintercept = 0,
             col = "#929292",
             linewidth = 0.3) +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  xlab("") +
  ylab("Log 2 of odds ratio") +
  theme(legend.title = element_blank())

ggsave(
  "../Charts, data, statistics/Panel C FMI.png",
  units = "in",
  dpi = 500,
  width = 6.65,
  height = 5.94
)

# Table version for Numbers ----------------------------------------------------
FMI_co_occurrence_chart_data %>%
  select(-both, -PTEN_only, -combination, -OUT_only, -neither, -p_val) %>%
  mutate(
    OUT_type = str_remove(OUT_type, "_OUT_ANY"),
    OUT_type = str_remove(OUT_type, "_OUT"),
    PTEN_type = str_remove_all(PTEN_type, "PTEN_")
  ) %>%
  pivot_wider(
    names_from = PTEN_type,
    values_from = c(log2_odds_ratio, log2_upper, log2_lower)
  ) %>%
  pivot_longer(cols = matches("log2")) %>%
  mutate(name = str_remove(name, "log2_")) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  mutate(
    upper_LOF = upper_LOF - odds_ratio_LOF,
    lower_LOF = odds_ratio_LOF - lower_LOF,
    upper_WT = upper_WT - odds_ratio_WT,
    lower_WT = odds_ratio_WT - lower_WT,
    upper_DEL = upper_DEL - odds_ratio_DEL,
    lower_DEL = odds_ratio_DEL - lower_DEL,
    upper_ANY = upper_ANY - odds_ratio_ANY,
    lower_ANY = odds_ratio_ANY - lower_ANY
  ) %>%
  bind_rows(tibble(OUT_type = c("BRAF", "PIK3R1"))) %>%
  arrange(OUT_type) %>%
  mutate(
    LOF = seq(1, 23, 3),
    WT = seq(2.4, 24, 3),
    ANY = seq(0.3, 22, 3),
    DEL = seq(1.7, 23, 3)
  ) %>%
  select(
    c(
      "OUT_type",
      "LOF",
      "odds_ratio_LOF",
      "upper_LOF",
      "lower_LOF",
      "WT",
      "odds_ratio_WT",
      "upper_WT",
      "lower_WT",
      "DEL",
      "odds_ratio_DEL",
      "upper_DEL",
      "lower_DEL",
      "ANY",
      "odds_ratio_ANY",
      "upper_ANY",
      "lower_ANY",
    )
  ) %>%
  mutate(across(
    .cols = where(is.numeric),
    .fns = function(x) {
      round(x, 2)
    }
  )) %>%
  write.csv("Panel C FMI chart data, Numbers edition.csv", row.names = F)



BOTH_co_occurrence_chart_data %>%
  select(-both, -PTEN_only, -combination, -OUT_only, -neither, -p_val) %>%
  mutate(
    OUT_type = str_remove(OUT_type, "_OUT_ANY"),
    OUT_type = str_remove(OUT_type, "_OUT"),
    PTEN_type = str_remove_all(PTEN_type, "PTEN_")
  ) %>%
  pivot_wider(
    names_from = PTEN_type,
    values_from = c(log2_odds_ratio, log2_upper, log2_lower)
  ) %>%
  pivot_longer(cols = matches("log2")) %>%
  mutate(name = str_remove(name, "log2_")) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  mutate(
    upper_LOF = upper_LOF - odds_ratio_LOF,
    lower_LOF = odds_ratio_LOF - lower_LOF,
    upper_WT = upper_WT - odds_ratio_WT,
    lower_WT = odds_ratio_WT - lower_WT,
    upper_DEL = upper_DEL - odds_ratio_DEL,
    lower_DEL = odds_ratio_DEL - lower_DEL,
    upper_ANY = upper_ANY - odds_ratio_ANY,
    lower_ANY = odds_ratio_ANY - lower_ANY
  ) %>%
  bind_rows(tibble(OUT_type = c("BRAF", "PIK3R1"))) %>%
  arrange(OUT_type) %>%
  mutate(
    LOF = seq(1, 23, 3),
    WT = seq(2.4, 24, 3),
    ANY = seq(0.3, 22, 3),
    DEL = seq(1.7, 23, 3)
  ) %>%
  select(
    c(
      "OUT_type",
      "LOF",
      "odds_ratio_LOF",
      "upper_LOF",
      "lower_LOF",
      "WT",
      "odds_ratio_WT",
      "upper_WT",
      "lower_WT",
      "DEL",
      "odds_ratio_DEL",
      "upper_DEL",
      "lower_DEL",
      "ANY",
      "odds_ratio_ANY",
      "upper_ANY",
      "lower_ANY",
    )
  ) %>%
  mutate(across(
    .cols = where(is.numeric),
    .fns = function(x) {
      round(x, 2)
    }
  )) %>%
  write.csv("Panel C BOTH chart data, Numbers edition.csv", row.names = F)
