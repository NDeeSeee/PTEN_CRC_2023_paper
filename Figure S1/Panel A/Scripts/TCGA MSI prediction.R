# R version 4.2.0 (2022-04-22)
# Encoding = UTF-8
# Packages
packages_names <-
  c("readxl", "data.table", "reshape2", "gridExtra", "tidyverse")

lapply(packages_names, require, character.only = TRUE)

# TCGA pan cancer data
setwd("~/Documents/KRAS_project/KRAS-project-update/PAD 2022 sources")

TCGA_mutations <-
  fread("coadread_tcga_pan_can_atlas_2018/data_mutations.txt") %>% as_tibble()
TCGA_samples <-
  fread("coadread_tcga_pan_can_atlas_2018/data_clinical_sample.txt",
        skip = 4) %>% as_tibble()
TCGA_clinical <-
  fread("coadread_tcga_pan_can_atlas_2018/data_clinical_patient.txt",
        skip = 4) %>% as_tibble()
TCGA_clinical_expanded <-
  fread("coadread_tcga_pan_can_atlas_2018/data_clinical_expanded.tsv") %>% as_tibble()

TCGA_samples_msi_true <-
  fread(
    "../KRAS-project/PAD_preprocessing/Oncopanels/coadread_tcga_pan_can_atlas_2018_clinical_data.tsv"
  ) %>%
  mutate(
    msi_true = ifelse(str_detect(Subtype, "MSI"), "MSI-H", Subtype),
    msi_true = ifelse(str_detect(Subtype, "POLE"), NA, msi_true),
    msi_true = ifelse(
      str_detect(Subtype, "POLE") |
        str_detect(Subtype, "MSI"),
      msi_true,
      "MSS"
    )
  ) %>%
  rename(msi_final = msi_true) %>%
  select(`Sample ID`, msi_final) %>%
  rename(sample_id = `Sample ID`)

# TCGA public
TCGA_samples_pub <-
  fread("coadread_tcga_pub/data_clinical_sample.txt", skip = 4) %>% as_tibble()

# TCGA PREDICTION MANTIS SENSOR

TCGA_msi_corr <-
  left_join(
    select(TCGA_samples_pub, SAMPLE_ID, MSI_STATUS),
    select(TCGA_samples, SAMPLE_ID, MSI_SCORE_MANTIS, MSI_SENSOR_SCORE),
    by = "SAMPLE_ID"
  ) %>%
  pivot_longer(
    cols = c(MSI_SCORE_MANTIS, MSI_SENSOR_SCORE),
    values_to = "score",
    names_to = "msi_method"
  ) %>%
  mutate(msi_method = ifelse(str_detect(msi_method, "SENSOR"), "sensor", "mantis")) %>%
  filter(!is.na(score), MSI_STATUS != "Not Evaluable") %>%
  rename(msi = MSI_STATUS, sample_id = SAMPLE_ID) %>%
  mutate(msi = str_replace(msi, "MSI-L", "MSS"))

# Algorithm to predict TCGA pancancer MSI
# default cutoffs 0.4 and 3.5

# Samples with missprediction
pivot_wider(TCGA_msi_corr, names_from = msi_method, values_from = score) %>%
  mutate(
    msi_mantis = cut(
      mantis,
      breaks = c(0, 0.4, 20),
      labels = c("MSS", "MSI-H")
    ),
    msi_sensor = cut(
      sensor,
      breaks = c(0, 3.5, 50),
      labels = c("MSS", "MSI-H")
    )
  ) %>%
  rename(msi_tcga = msi,
         mantis_score = mantis,
         sensor_score = sensor) %>%
  filter(msi_mantis != msi_sensor)

# SENSOR score distribution across MSI
TCGA_msi_corr %>%
  filter(msi_method == "sensor") %>%
  ggplot(aes(
    x = as.factor(msi),
    y = score,
    fill = as.factor(msi)
  )) +
  geom_boxplot()

filter(TCGA_msi_corr, msi_method == "sensor") %>%
  group_by(msi) %>%
  summarise(min = min(score), max = max(score))

ggsave("../KRAS-project/PAD_preprocessing/Oncopanels/sensor_msi.png")

# MANTIS score distribution across MSI
TCGA_msi_corr %>%
  filter(msi_method == "mantis") %>%
  ggplot(aes(
    x = as.factor(msi),
    y = score,
    fill = as.factor(msi)
  )) +
  geom_boxplot()

filter(TCGA_msi_corr, msi_method == "mantis") %>%
  group_by(msi) %>%
  summarise(min = min(score), max = max(score))

ggsave("../KRAS-project/PAD_preprocessing/Oncopanels/mantis_msi.png")


# Correspondence between MSI score (pancaner) and MSI status (pub)
TCGA_msi_corr <-
  left_join(
    select(TCGA_samples_pub, SAMPLE_ID, MSI_STATUS),
    select(TCGA_samples, SAMPLE_ID, MSI_SCORE_MANTIS, MSI_SENSOR_SCORE),
    by = "SAMPLE_ID"
  ) %>%
  pivot_longer(
    cols = c(MSI_SCORE_MANTIS, MSI_SENSOR_SCORE),
    values_to = "score",
    names_to = "msi_method"
  ) %>%
  mutate(msi_method = ifelse(str_detect(msi_method, "SENSOR"), "sensor", "mantis")) %>%
  filter(!is.na(score), MSI_STATUS != "Not Evaluable") %>%
  rename(msi = MSI_STATUS, sample_id = SAMPLE_ID) %>%
  mutate(msi = str_replace(msi, "MSI-L", "MSS"))

TCGA_samples %>%
  select(SAMPLE_ID, MSI_SCORE_MANTIS, MSI_SENSOR_SCORE) %>%
  mutate(
    msi_mantis = cut(
      MSI_SCORE_MANTIS,
      breaks = c(0, 0.4, 20),
      labels = c("MSS", "MSI-H")
    ),
    msi_sensor = cut(
      MSI_SENSOR_SCORE,
      breaks = c(0, 3.5, 50),
      labels = c("MSS", "MSI-H")
    )
  ) %>%
  left_join(y = select(TCGA_samples_pub, SAMPLE_ID, MSI_STATUS),
            by = "SAMPLE_ID") %>%
  rename(
    mantis_score = MSI_SCORE_MANTIS,
    sensor_score = MSI_SENSOR_SCORE,
    msi_tcga = MSI_STATUS,
    sample_id = SAMPLE_ID
  ) %>%
  mutate(
    msi_tcga = str_replace(msi_tcga, "MSI-L", "MSS"),
    msi_tcga = ifelse(msi_tcga == "Not Evaluable", NA, msi_tcga)
  ) %>%
  write.csv(
    "../KRAS-project/PAD_preprocessing/Oncopanels/TCGA_msi_prediction.csv",
    row.names = F
  )

TCGA_msi_prediction <-
  fread("../KRAS-project/PAD_preprocessing/Oncopanels/TCGA_msi_prediction.csv") %>% as_tibble()

TCGA_msi_prediction %>%
  mutate(msi_final = msi_tcga) %>%
  mutate(msi_final = ifelse(
    is.na(msi_tcga) &
      (msi_sensor == "MSS" | msi_mantis == "MSS"),
    "MSS",
    msi_tcga
  )) %>%
  mutate(msi_final = ifelse(
    is.na(msi_final) &
      (msi_sensor == "MSI-H" &
         msi_mantis == "MSI-H"),
    "MSI-H",
    msi_final
  )) %>%
  mutate(msi_final = ifelse(
    is.na(msi_final) &
      (msi_sensor == "MSI-H" |
         msi_mantis == "MSI-H"),
    "MSI-H",
    msi_final
  )) %>%
  write.csv(
    "../KRAS-project/PAD_preprocessing/Oncopanels/TCGA_msi_prediction_completed.csv",
    row.names = F
  )

TCGA_msi_full <-
  fread(
    "../KRAS-project/PAD_preprocessing/Oncopanels/TCGA_msi_prediction_completed.csv"
  ) %>%
  as_tibble() %>%
  rename(msi_pred = msi_final) %>%
  left_join(TCGA_samples_msi_true) %>%
  mutate(msi_final = ifelse(is.na(msi_final), msi_pred, msi_final))

TCGA_msi_full %>%
  write_csv("../KRAS-project/PAD_preprocessing/Oncopanels/TCGA_msi_prediction_samples.csv")

TCGA_mutations_msi <-
  filter(TCGA_mutations, Consequence != "synonymous_variant") %>%
  select(Tumor_Sample_Barcode, Hugo_Symbol, HGVSp_Short) %>%
  rename(sample_id = Tumor_Sample_Barcode,
         gene = Hugo_Symbol,
         alt = HGVSp_Short) %>%
  mutate(alt = str_remove(alt, "p.")) %>%
  left_join(
    y = select(
      TCGA_msi_full,
      sample_id,
      msi_tcga,
      mantis_score,
      sensor_score,
      msi_final
    ),
    by = "sample_id"
  ) %>%
  mutate(
    msi_MLH = ifelse(gene == "MLH1" & str_detect(alt, "[*]"), 1, 0),
    msi_RNF = ifelse(gene == "RNF43" &
                       str_detect(alt, "G659"), 1, 0)
  ) %>%
  left_join(y = select(
    rename(
      TCGA_clinical_expanded,
      sample_id = `Sample ID`,
      cna_frac = `Fraction Genome Altered`
    ),
    sample_id,
    cna_frac
  )) %>%
  select(-alt) %>%
  rename(msi_pred = msi_final) %>%
  left_join(TCGA_samples_msi_true) %>%
  mutate(msi_final = ifelse(is.na(msi_final), msi_pred, msi_final))

TCGA_mutations_msi %>%
  write_csv(
    "../KRAS-project/PAD_preprocessing/Oncopanels/TCGA_msi_prediction_completed_expanded.csv"
  )

TCGA_mutations_msi %>%
  select(-mantis_score, -sensor_score, -msi_tcga, -msi_pred, -gene) %>%
  group_by(sample_id) %>%
  summarise(
    mut_count = n(),
    msi_final = unique(msi_final)[1],
    msi_MLH = max(msi_MLH),
    msi_RNF = max(msi_RNF),
    cna_frac = max(cna_frac)
  ) %>%
  write_csv(
    "../KRAS-project/PAD_preprocessing/Oncopanels/TCGA_msi_full_for_predict_samples.csv"
  )

TCGA_mutations_msi %>%
  group_by(sample_id) %>%
  summarise(mut_count = n(), msi_final = unique(msi_final)[1]) %>%
  write_csv(
    "../KRAS-project/PAD_preprocessing/Oncopanels/TCGA_msi_prediction_mut_count.csv"
  )

# MSKCC cross-validation
MSC_mut_count <- group_by(MSC_mutations, Tumor_Sample_Barcode) %>%
  summarise(mut_count = n()) %>%
  rename(sample_id = Tumor_Sample_Barcode)

MSC_samples_msi <-
  select(MSC_samples, SAMPLE_ID, MSI_STATUS, MSI_SCORE, GENE_PANEL) %>%
  rename(sample_id = SAMPLE_ID,
         msi_mskcc = MSI_STATUS,
         assay = GENE_PANEL) %>%
  mutate(
    msi_sensor = ifelse(MSI_SCORE > 3.5, "MSI-H", MSI_SCORE),
    msi_sensor = ifelse(MSI_SCORE <= 3.5, "MSS", msi_sensor),
    assay = paste("MSK", assay, sep = "-")
  ) %>%
  left_join(y = MSC_mut_count, by = "sample_id") %>%
  mutate(
    msi_mskcc_true = as.numeric(str_replace(
      str_replace(msi_mskcc, "MSS", "0"), "MSI", "1"
    )),
    msi_sensor_true = as.numeric(str_replace(
      str_replace(msi_sensor, "MSS", "0"), "MSI-H", "1"
    ))
  )


MSC_samples_msi_pred <-
  left_join(
    x = MSC_samples_msi,
    y = select(all_cutoff_panel, assay, cutoff_under, cutoff_upper),
    by = "assay"
  ) %>%
  filter(!is.na(mut_count), !is.na(cutoff_under)) %>%
  mutate(
    msi_pred = ifelse(mut_count < cutoff_under, "MSS", NA),
    msi_pred = ifelse(mut_count > cutoff_upper, "MSI-H", msi_pred),
    msi_pred_status = as.numeric(str_replace(
      str_replace(msi_pred, "MSS", "0"), "MSI-H", "1"
    ))
  )

MSC_true_pred_1 <- MSC_samples_msi_pred %>%
  filter(!is.na(msi_pred_status), !is.na(msi_mskcc_true))

auc(MSC_true_pred_1$msi_mskcc_true,
    MSC_true_pred_1$msi_pred_status)

MSC_true_pred_2 <- MSC_samples_msi_pred %>%
  filter(!is.na(msi_pred_status), !is.na(msi_sensor_true))

auc(MSC_true_pred_2$msi_sensor_true,
    MSC_true_pred_2$msi_pred_status)
