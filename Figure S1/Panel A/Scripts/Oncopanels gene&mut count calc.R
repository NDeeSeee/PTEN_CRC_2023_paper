# R version 4.2.0 (2022-04-22)
# Encoding = UTF-8
# Packages
packages_names <-
  c("readxl", "data.table", "reshape2", "gridExtra", "tidyverse")

lapply(packages_names, require, character.only = TRUE)

setwd(
  "~/Documents/KRAS_project/KRAS-project-update/KRAS-project/PAD 2022 sources/OncoPanels/"
)
TCGA_mutations_msi <-
  fread(
    "../../PAD_preprocessing/Oncopanels/TCGA_msi_prediction_completed_expanded.csv"
  )

# This data was downloaded from synapse
# https://www.synapse.org/#!Synapse:syn50678294
gene_panel_files <-
  list.files()[str_detect(list.files(), "data_gene_panel")]

# Here we get only panel name and gene count
gene_panel_all <- data.frame()
for (gene_panel in gene_panel_files) {
  int_data <- read.table(gene_panel, nrows = 2, fill = T)
  int_data_df <-
    data.frame(oncopanel = int_data$V2[1], gene_count = int_data$V7[2])
  gene_panel_all <- bind_rows(gene_panel_all, int_data_df)
}
gene_panel_all <- gene_panel_all %>%
  as_tibble() %>%
  arrange(-gene_count)

# Save it as data frame
write.csv(gene_panel_all, "Oncopanels.csv", row.names = F)

# Save separately file for those panels, which targeted more than 200 genes
# Less than 200 worse for any such of prediction
gene_panel_all %>%
  filter(gene_count > 200) %>%
  write.csv("OncoPanels_200.csv", row.names = F)

onco_panels_new <- fread("OncoPanels_200.csv") %>%
  as_tibble()

# Calculate count of mut per sample for each oncopanel
# Based on TCGA mutations data
for (i in onco_panels_new$oncopanel) {
  panel_gene_list <-
    t(read.table(paste0("data_gene_panel_", i, ".txt"), skip = 2)[-1]) %>%
    as_tibble(.name_repair = "universal") %>%
    rename(gene = ...1)
  panel_mut_count <-
    filter(TCGA_mutations_msi, gene %in% panel_gene_list$gene) %>%
    group_by(sample_id) %>%
    summarise(mut_count = n(), msi = unique(msi_final)[1])
  write.csv(
    panel_mut_count,
    file = paste0("TCGA_mut_count/mut_count_panel_", i, ".csv"),
    row.names = F
  )
}
