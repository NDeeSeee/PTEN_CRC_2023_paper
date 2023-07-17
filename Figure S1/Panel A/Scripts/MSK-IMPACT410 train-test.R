#Packages
packages_names <- c("tidyverse", "readxl", "data.table", "reshape2",
                    "ggpubr", "gridExtra", "pROC", "Metrics", "party")

lapply(packages_names, require, character.only = TRUE)


# WD
setwd("~/Documents/KRAS_project/KRAS-project-update/PAD 2022 sources")

# Input data
MSC_pat_full = fread("../KRAS-project/PAD 2022 sources/MSC_after_deduplication.csv") %>% as_tibble()
MSC_samples = fread("crc_msk_2017/data_clinical_sample.txt", skip = 4) %>% as_tibble() 
TCGA_mutations_msi = fread("../KRAS-project/PAD_preprocessing/Oncopanels/TCGA_msi_prediction_completed_expanded.csv") %>% as_tibble()

# Function
get_panel_mut_count = function(panel_name) {
  panel_gene_list = t(read.table(paste0("OncoPanels/data_gene_panel_", panel_name, ".txt"), skip = 2)[-1]) %>%
    as_tibble(.name_repair = "universal") %>%
    rename(gene = ...1)
  
  panel_mut_count = filter(TCGA_mutations_msi, gene %in% panel_gene_list$gene) %>%
    group_by(sample_id) %>%
    summarise(mut_count = n(), msi_true = unique(msi_tcga)[1], mantis_score = mean(mantis_score, na.rm = T),
              sensor_score = mean(sensor_score, na.rm = T), msi_MLH = max(msi_MLH, na.rm = T), 
              msi_RNF = max(msi_RNF, na.rm = T), cna_frac = mean(cna_frac, na.rm = T),
              msi_final = unique(msi_final)[1]) %>%
    filter(!is.na(msi_final)) %>%
    mutate(mantis_score = ifelse(mantis_score >= 0.4, 1, mantis_score),
           mantis_score = ifelse(mantis_score < 0.4, 0, mantis_score),
           sensor_score = ifelse(sensor_score >= 3.5, 1, sensor_score),
           sensor_score = ifelse(sensor_score < 3.5, 0, sensor_score)) %>%
    mutate(mantis_score = as.factor(mantis_score), sensor_score = as.factor(sensor_score), 
           msi_MLH = as.factor(msi_MLH), msi_RNF = as.factor(msi_RNF),
           msi_true = as.factor(msi_true), msi_final = as.factor(msi_final),
           hypermut = ifelse(mut_count > quantile(mut_count, probs = 0.967), 1, 0)) %>%
    filter(hypermut == 0) %>%
    select(-sample_id, -msi_true, -mantis_score, -sensor_score, -hypermut)
  return(panel_mut_count)
}
  
# Validation case
val_panel_name = "MSK-IMPACT410"

train_data = get_panel_mut_count(val_panel_name)

#train_data %>%
#  write_csv("panel_410_train_data.csv")

model_410_R = ctree(msi_final ~ ., train_data, controls = ctree_control(mincriterion = 0.995))
plot(model_410_R)
#ggsave("MSK-IMPACT410 tree model.png")

# remove external duplicates from PAD_2022_rebuild.R
MSC_samples = filter(MSC_samples, SAMPLE_ID %in% MSC_pat_full$sample_id)

MSC_samples_410 = filter(MSC_samples, GENE_PANEL == "IMPACT410") %>%
  rename(sample_id = SAMPLE_ID, msi_mskcc = MSI_STATUS, msi_score = MSI_SCORE, assay = GENE_PANEL) %>%
  select(sample_id, msi_mskcc, msi_score, assay)

panel_gene_list = t(read.table(paste0("OncoPanels/data_gene_panel_", val_panel_name, ".txt"), skip = 2)[-1]) %>%
  as_tibble(.name_repair = "universal") %>%
  rename(gene = ...1)

MSK_mutations_410 = filter(MSC_mutations, Consequence != "synonymous_variant") %>%
  select(Tumor_Sample_Barcode, Hugo_Symbol, HGVSp_Short) %>%
  rename(sample_id = Tumor_Sample_Barcode, gene = Hugo_Symbol, alt = HGVSp_Short) %>%
  mutate(alt = str_remove(alt, "p.")) %>%
  left_join(y = select(MSC_samples_410, sample_id, msi_mskcc, msi_score), by ="sample_id") %>%
  mutate(msi_MLH = ifelse(gene == "MLH1" & str_detect(alt, "[*]"), 1, 0),
         msi_RNF = ifelse(gene == "RNF43" & str_detect(alt, "G659"), 1, 0)) %>%
  left_join(y = select(rename(MSC_clinical_expanded, sample_id = `Sample ID`, cna_frac = `Fraction Genome Altered`, mut_count_native = `Mutation Count`), sample_id, cna_frac, mut_count_native)) %>%
  select(-alt) 

MSK_test_410 = group_by(MSK_mutations_410, sample_id) %>%
  summarise(mut_count = n(), msi_mskcc = unique(msi_mskcc)[1], msi_score = mean(msi_score, na.rm = T),
            mut_count_native = min(mut_count_native, na.rm = T), msi_MLH = max(msi_MLH, na.rm = T), 
            msi_RNF = max(msi_RNF, na.rm = T), cna_frac = mean(cna_frac, na.rm = T)) %>%
  filter(!is.na(msi_mskcc)) %>%
  mutate(msi_final = ifelse(msi_score >= 10, "MSI-H", msi_score),
         msi_final = ifelse(msi_score < 10, "MSS", msi_final)) %>%
  mutate(msi_mskcc = ifelse(msi_mskcc == "Inconclusive", NA, msi_mskcc), 
         msi_mskcc = str_replace(msi_mskcc, "MSI", "MSI-H")) %>% 
  mutate(msi_final = as.factor(msi_final), msi_mskcc = as.factor(msi_mskcc), 
         msi_MLH = as.factor(msi_MLH), msi_RNF = as.factor(msi_RNF)) %>%
  select(-sample_id) %>%
  mutate(hypermut = ifelse(mut_count > quantile(mut_count, probs = 0.967), 1, 0)) %>%
  filter(hypermut == 0) %>%
  select(-hypermut)

MSK_test_410_origin = MSK_test_410 %>%
  select(-msi_final, -mut_count, -msi_score) %>%
  rename(msi_final = msi_mskcc, mut_count = mut_count_native)

#MSK_test_410_origin %>%
#  write_csv("~/panel_410_test_data.csv")

#Predict with R and metrics
true_values = filter(MSK_test_410_origin, !is.na(msi_final))$msi_final
predicted_values_probs = predict(model_410_R, select(filter(MSK_test_410_origin, !is.na(msi_final)), -msi_final), type = "prob")
# Here we decided to pick 0.85 probability as threshold 
predicted_values = factor(sapply(predicted_values_probs, "[[", 1) >= 0.85, levels = c(T, F), labels = c("MT-H", "MSS"))


table(true_values, predicted_values)

auc(as.numeric(true_values)-1, as.numeric(predicted_values)-1)
roc(as.numeric(true_values) ~ as.numeric(predicted_values), plot = TRUE, print.auc = TRUE)

