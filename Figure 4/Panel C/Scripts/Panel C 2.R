# Attach requirement packages and setting WD -----------------------------------
packages_names <- c("readxl", "data.table", "reshape2", "tidyverse",
                    "rstudioapi")

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

# Justin data ------------------------------------------------------------------
genes <-
  c("SMAD4",
    "APC",
    "BRAF",
    "KRAS",
    "PIK3CA",
    "PIK3R1",
    "TP53_mut",
    "TP53_del")

piece_of_table <- tibble()
for (i in 1:8) {
  piece_of_table <-
    read_xlsx(
      "~/Documents/GitHub/PTEN_CRC_2022/PTEN_cooccurence/FMI_FCCC_CRC_PTEN_coocc_230508/Cooccurrence_PTEN_panel.xlsx",
      sheet = i
    ) %>%
    as_tibble() %>%
    mutate(across(
      matches("both|only|neither|odds"),
      .fns = function(x) {
        replace_na(as.numeric(x), 0)
      }
    ),
    gene = genes[i]) %>%
    filter(
      str_detect(
        group,
        "MAVE_lof|MAVE_wt|VAMP_07_lof|VAMP_07_wt|DB_MAVE_VAMP_07_lof|DB_MAVE_VAMP_07_wt"
      ),
      group != "MAVE_VAMP_07_lof",
      group != "MAVE_VAMP_07_wt",
      msi_status == "MSS"
    ) %>%
    bind_rows(piece_of_table)
}

FMI_co_occurrence_chart_data <- piece_of_table %>%
  mutate(group = as.factor(group),
         gene = as.factor(gene)) %>%
  select(-msi_status, -odds_ratio, -p_value) %>%
  group_by(group, gene) %>%
  mutate(fisher_test_res = paste(as.character(round(
    unlist(fisher.test(matrix(
      c(
        n_both_pten_panel,
        n_pten_only,
        n_panel_only,
        n_neither_pten_panel
      ),
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
  )) %>%
  mutate(across(
    .cols = matches("log2"),
    .fns = function(x) {
      ifelse(is.infinite(x), 0, x)
    }
  ),
  combination = as.factor(paste(gene, group, sep = " "))) %>%
  filter(p_val != 1,
         if_all(matches("log2"), ~ .x != 0)) %>%
  mutate(
    group = case_when(
      group == "DB_MAVE_VAMP_07_lof" ~ "PTEN_LOF_consensus",
      group == "DB_MAVE_VAMP_07_wt" ~ "PTEN_WT_consensus",
      group == "MAVE_lof" ~ "PTEN_LOF_mave",
      group == "MAVE_wt" ~ "PTEN_LOF_mave",
      group == "VAMP_07_lof" ~ "PTEN_LOF_consensus",
      group == "VAMP_07_wt" ~ "PTEN_LOF_consensus",
      .default = as.character(group)
    )
  ) %>%
  mutate(
    gene = str_replace(gene, "del", "OUT_DEL"),
    gene = str_replace(gene, "mut", "OUT_MUT"),
    across(
      matches("gene"),
      .fns = function(x)
        paste(x, "_OUT_ANY", sep = "")
    ),
    gene = str_replace(gene, "TP53_OUT_DEL_OUT_ANY", "TP53_OUT_DEL"),
    gene = str_replace(gene, "TP53_OUT_MUT_OUT_ANY", "TP53_OUT_MUT")
  )


# PAD data ---------------------------------------------------------------------
PAD_variants_grouped <- PAD_variants %>%
  mutate(
    pten_lof_consensus = ifelse(
      gene == "PTEN" &
        alt %in% filter(PTEN_LoF_data_full, consensus_VM_KB == "lof")$alt,
      1,
      NA
    ),
    pten_lof_consensus = ifelse(
      gene == "PTEN" &
        alt %in% filter(PTEN_LoF_data_full, consensus_VM_KB == "wt")$alt,
      0,
      pten_lof_consensus
    ),
    pten_lof_consensus = ifelse(gene == "PTEN" &
                                  str_detect(alt, "\\*"),
                                1,
                                pten_lof_consensus),
    pten_lof_mave = ifelse(
      gene == "PTEN" &
        alt %in% filter(PTEN_LoF_data_full, MAVE_effect == "lof")$alt,
      1,
      NA
    ),
    pten_lof_mave = ifelse(
      gene == "PTEN" &
        alt %in% filter(PTEN_LoF_data_full, MAVE_effect == "wt")$alt,
      0,
      pten_lof_mave
    ),
    pten_lof_mave = ifelse(gene == "PTEN" & str_detect(alt, "\\*"),
                           1,
                           pten_lof_mave),
    pten_lof_vamp = ifelse(
      gene == "PTEN" &
        alt %in% filter(PTEN_LoF_data_full, VAMP_effect == "lof")$alt,
      1,
      NA
    ),
    pten_lof_vamp = ifelse(
      gene == "PTEN" &
        alt %in% filter(PTEN_LoF_data_full, VAMP_effect == "wt")$alt,
      0,
      pten_lof_vamp
    ),
    pten_lof_vamp = ifelse(gene == "PTEN" & str_detect(alt, "\\*"),
                           1,
                           pten_lof_vamp)
  ) %>%
  group_by(sample_id) %>%
  summarise(
    PTEN_LOF_consensus = sign(sum(pten_lof_consensus == 1, na.rm = T)),
    PTEN_WT_consensus = sign(sum(pten_lof_consensus == 0, na.rm = T)),
    PTEN_LOF_mave = sign(sum(pten_lof_mave == 1, na.rm = T)),
    PTEN_WT_mave = sign(sum(pten_lof_mave == 0, na.rm = T)),
    PTEN_LOF_vamp = sign(sum(pten_lof_vamp == 1, na.rm = T)),
    PTEN_WT_vamp = sign(sum(pten_lof_vamp == 0, na.rm = T))
  ) %>%
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
  mutate(across(
    matches("LOF|WT"),
    .fns = function(x)
      replace_na(x, 0)
  )) %>%
  filter(!if_any(
    .cols = everything(),
    .fns = function(x) {
      is.na(x)
    }
  ))


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
  combination = as.factor(paste(PTEN_type, OUT_type, sep = " "))) %>%
  rename(
    group = PTEN_type,
    gene = OUT_type,
    n_both_pten_panel = both,
    n_pten_only = PTEN_only,
    n_panel_only = OUT_only,
    n_neither_pten_panel = neither
  ) %>%
  mutate(gene = as.character(gene))


combined_df <- bind_rows(PAD_co_occurrence_chart_data,
                         FMI_co_occurrence_chart_data) %>%
  mutate(gene = as.factor(gene),
         group = as.factor(group))

combined_df_result <- combined_df %>%
  group_by(group, gene) %>%
  summarise(
    n_both_pten_panel = sum(n_both_pten_panel, na.rm = T),
    n_pten_only = sum(n_pten_only, na.rm = T),
    n_panel_only = sum(n_panel_only, na.rm = T),
    n_neither_pten_panel = sum(n_neither_pten_panel, na.rm = T)
  ) %>%
  group_by(group, gene) %>%
  mutate(fisher_test_res = paste(as.character(round(
    unlist(fisher.test(matrix(
      c(
        n_both_pten_panel,
        n_pten_only,
        n_panel_only,
        n_neither_pten_panel
      ),
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
  combination = as.factor(paste(group, gene, sep = " ")))

combined_df_result %>%
  select(group, gene, log2_lower, log2_upper, log2_odds_ratio) %>%
  mutate(log2_lower = log2_odds_ratio - log2_lower,
         log2_upper = log2_upper - log2_odds_ratio) %>%
  pivot_wider(
    names_from = group,
    values_from = c(log2_odds_ratio, log2_upper, log2_lower)
  ) %>%
  pivot_longer(cols = matches("log2")) %>%
  mutate(name = str_remove(name, "log2_"),
         name = str_remove(name, "_PTEN")) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  arrange(gene) %>%
  mutate(
    LOF_mave = seq(0.3, 22, 3),
    WT_mave = seq(0.6, 23, 3),
    LOF_vamp = seq(1.3, 23, 3),
    WT_vamp = seq(1.6, 24, 3),
    LOF_consensus = seq(2.3, 25, 3),
    WT_consensus = seq(2.6, 25, 3)
  ) %>%
  select(
    c(
      "gene",
      "LOF_mave",
      "odds_ratio_LOF_mave",
      "upper_LOF_mave",
      "lower_LOF_mave",
      "LOF_vamp",
      "odds_ratio_LOF_vamp",
      "upper_LOF_vamp",
      "lower_LOF_vamp",
      "LOF_consensus",
      "odds_ratio_LOF_consensus",
      "upper_LOF_consensus",
      "lower_LOF_consensus",
      "WT_mave",
      "odds_ratio_WT_mave",
      "upper_WT_mave",
      "lower_WT_mave",
      "WT_vamp",
      "odds_ratio_WT_vamp",
      "upper_WT_vamp",
      "lower_WT_vamp",
      "WT_consensus",
      "odds_ratio_WT_consensus",
      "upper_WT_consensus",
      "lower_WT_consensus"
    )
  ) %>%
  mutate(across(
    .cols = where(is.numeric),
    .fns = function(x) {
      round(x, 2)
    }
  )) %>%
  write.csv("Justin data 2023 co-occurrence chart data with PAD, Numbers edition.csv",
            row.names = F)
