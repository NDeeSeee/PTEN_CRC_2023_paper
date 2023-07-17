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
  filter(codingType_SV != "synonymous") %>%
  mutate(alt = description_v2) %>%
  rename(did = deidentifiedSpecimenName) %>%
  left_join(y = select(FMI_samples, did, Assigned), by = "did")

FMI_variants_PTEN <-
  filter(FMI_variants, gene == "PTEN", codingType_SV != "synonymous")

# These are ALL hotspots -------------------------------------------------------
PTEN_hotspots <-
  fread("../../../Common raw data/true_PTEN_hotspots.csv")

# Loss of function -------------------------------------------------------------
PTEN_lof_data <-
  fread("../../../Common raw data/PTEN LoF annotation.csv")


PTEN_variants <- filter(FMI_variants, gene == "PTEN") %>%
  left_join(select(PTEN_lof_data, did, alt, tot_effect), by = c("did", "alt"))

# Assign codon types -----------------------------------------------------------
PTEN_variants <- PTEN_variants %>% mutate(
  htsp_cmb = 0,
  htsp_point = 0,
  htsp_trunc = 0,
  htsp_trunc_nons = 0,
  htsp_trunc_frms = 0,
  syn_mut = 0,
  lof_mut = 0,
  wt_mut = 0
)

# All hotspots are marked as " 1 "
PTEN_variants$htsp_cmb[PTEN_variants$codon_SV %in% PTEN_hotspots$codon] <-
  1

# Then assign other categories -------------------------------------------------
PTEN_variants$htsp_point[PTEN_variants$codon_SV %in% PTEN_hotspots$codon &
                           PTEN_variants$variantClass_v2 == "point"] <- 1
PTEN_variants$htsp_trunc[PTEN_variants$codon_SV %in% PTEN_hotspots$codon &
                           PTEN_variants$variantClass_v2 == "truncation"] <- 1
PTEN_variants$htsp_trunc_nons[PTEN_variants$codon_SV %in% PTEN_hotspots$codon &
                                PTEN_variants$variantClass_v2 == "truncation" &
                                PTEN_variants$codingType_SV == "nonsense"] <- 1
PTEN_variants$htsp_trunc_frms[PTEN_variants$codon_SV %in% PTEN_hotspots$codon &
                                PTEN_variants$variantClass_v2 == "truncation" &
                                PTEN_variants$codingType_SV == "frameshift"] <- 1
PTEN_variants$syn_mut[!PTEN_variants$codon_SV %in% PTEN_hotspots$codon &
                        PTEN_variants$codingType_SV == "synonymous"] <- 1
PTEN_variants$lof_mut[PTEN_variants$tot_effect == "lof"] <- 1
PTEN_variants$wt_mut[PTEN_variants$tot_effect == "wt"] <- 1

PTEN_AF <- PTEN_variants %>%
  select(
    did,
    Assigned,
    gene,
    codon_SV,
    htsp_cmb,
    htsp_point,
    htsp_trunc,
    syn_mut,
    lof_mut,
    alleleFreq_SV
  ) %>%
  filter(!is.na(codon_SV)) %>%
  rename(PTEN_AF = alleleFreq_SV) %>%
  mutate(
    htsp_cmb = ifelse(syn_mut == 1, "syn", htsp_cmb),
    htsp_cmb = ifelse(htsp_cmb == 1, "htsp", htsp_cmb),
    htsp_cmb = ifelse(htsp_cmb == 0, "non", htsp_cmb)
  )

PTEN_lof_AF <-
  select(PTEN_AF, did, Assigned, gene, codon_SV, lof_mut, PTEN_AF) %>%
  rename(htsp_cmb = lof_mut) %>%
  mutate(htsp_cmb = ifelse(htsp_cmb == 1, "lof", "wt"))


PTEN_AF <- bind_rows(PTEN_AF, PTEN_lof_AF) %>%
  select(-htsp_point, -htsp_trunc, -syn_mut, -lof_mut)


AF_genes <- FMI_variants %>%
  filter(
    gene %in% c("APC", "KRAS", "TP53", "PIK3CA", "PTEN"),
    did %in% PTEN_AF$did,
    Assigned == "MT-L"
  ) %>%
  mutate(gene = paste(gene, "_AF", sep = "")) %>%
  group_by(did, gene) %>%
  summarise(AF = max(alleleFreq_SV)) %>%
  ungroup() %>%
  filter(!is.na(AF)) %>%
  mutate(gene = str_remove_all(gene, "_AF"),
         AF = as.numeric(AF)) %>%
  filter(!is.na(AF))

write.csv(AF_genes,
          "../Charts, data, statistics/Panel F data.csv",
          row.names = F)

# Supplementary Figure 3, Panel F chart ----------------------------------------
AF_genes %>%
  ggplot(aes(x = gene, y = AF)) +
  geom_violin(trim = T,
              fill = "#34cfeb",
              alpha = .6) +
  geom_boxplot(width = 0.1, col = "blue") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold")
  ) +
  ylab("Allele Frequency") +
  xlab("Gene")

other_genes <- unique(AF_genes$gene)[unique(AF_genes$gene) != "PTEN"]

ks_test_list <- list()
for (i in 1:length(other_genes)) {
  temp_test <- t.test(filter(AF_genes, gene == "PTEN")$AF,
                      filter(AF_genes, gene == other_genes[i])$AF)
  temp_test$data.name <- paste0(c("x: PTEN, y: ", other_genes[i]), collapse = "")
  ks_test_list[[i]] <- temp_test
}

saveRDS(ks_test_list, file = "../Charts, data, statistics/Panel F statistics ks test.RData")


ggsave(
  "../Charts, data, statistics/Panel F.png",
  units = "in",
  dpi = 500,
  width = 6.65,
  height = 5.94
)
