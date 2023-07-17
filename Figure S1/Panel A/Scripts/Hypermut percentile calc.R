# Packages
packages_names <- c(
  "tidyverse",
  "readxl",
  "data.table",
  "reshape2",
  "ggpubr",
  "gridExtra",
  "pROC",
  "Metrics",
  "party"
)

lapply(packages_names, require, character.only = TRUE)


# Counting hypermutation cutoff

test_panel_name <- "MSK-IMPACT410"

get_panel_mut_count <- function(panel_name) {
  panel_gene_list <-
    t(read.table(
      paste0("OncoPanels/data_gene_panel_", panel_name, ".txt"),
      skip = 2
    )[-1]) %>%
    as_tibble(.name_repair = "universal") %>%
    rename(gene = ...1)
  
  panel_mut_count <-
    filter(TCGA_mutations_msi, gene %in% panel_gene_list$gene) %>%
    group_by(sample_id) %>%
    summarise(
      mut_count = n(),
      msi_true = unique(msi_tcga)[1],
      mantis_score = mean(mantis_score, na.rm = T),
      sensor_score = mean(sensor_score, na.rm = T),
      msi_MLH = max(msi_MLH, na.rm = T),
      msi_RNF = max(msi_RNF, na.rm = T),
      cna_frac = mean(cna_frac, na.rm = T),
      msi_final = unique(msi_final)[1]
    ) %>%
    filter(!is.na(msi_final)) %>%
    mutate(
      mantis_score = ifelse(mantis_score >= 0.4, 1, mantis_score),
      mantis_score = ifelse(mantis_score < 0.4, 0, mantis_score),
      sensor_score = ifelse(sensor_score >= 3.5, 1, sensor_score),
      sensor_score = ifelse(sensor_score < 3.5, 0, sensor_score)
    ) %>%
    mutate(
      mantis_score = as.factor(mantis_score),
      sensor_score = as.factor(sensor_score),
      msi_MLH = as.factor(msi_MLH),
      msi_RNF = as.factor(msi_RNF),
      msi_true = as.factor(msi_true),
      msi_final = as.factor(msi_final)
    ) %>%
    select(-sample_id, -msi_true, -mantis_score, -sensor_score)
  
  return(panel_mut_count)
}

panel_mut_count <- get_panel_mut_count(test_panel_name)


panel_mut_count %>%
  mutate(hypemut = ifelse(mut_count > 100, 1, 0)) %>%
  select(hypemut) %>%
  table() %>%
  prop.table()

#         0          1
# 0.96699029 0.03300971
# == 3.3%
