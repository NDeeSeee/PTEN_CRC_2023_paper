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

top_htsps <- fread("../../../Common raw data/Hotspots by MS status.csv")

top_htsps_mtl <- filter(top_htsps, 
                        htsp_status == "top_htsp", 
                        Assigned == "MT-L")$codon

top_htsps_mth <- filter(top_htsps, 
                        htsp_status == "top_htsp", 
                        Assigned == "MT-H")$codon

# Loss of function -------------------------------------------------------------
PTEN_lof_data <-
  fread("../../../Common raw data/PTEN LoF annotation.csv")

# Data processing --------------------------------------------------------------
PTEN_variants <- filter(FMI_variants, gene == "PTEN") %>%
  left_join(select(PTEN_lof_data, did, alt, tot_effect), by = c("did", "alt"))

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
PTEN_variants$htsp_cmb[PTEN_variants$codon_SV %in% PTEN_hotspots$codon] <-
  1
PTEN_variants$htsp_point[PTEN_variants$codon_SV %in% PTEN_hotspots$codon &
                           PTEN_variants$variantClass_v2 == "point"] <-
  1
PTEN_variants$htsp_trunc[PTEN_variants$codon_SV %in% PTEN_hotspots$codon &
                           PTEN_variants$variantClass_v2 == "truncation"] <-
  1
PTEN_variants$htsp_trunc_nons[PTEN_variants$codon_SV %in% PTEN_hotspots$codon &
                                PTEN_variants$variantClass_v2 == "truncation" &
                                PTEN_variants$codingType_SV == "nonsense"] <-
  1
PTEN_variants$htsp_trunc_frms[PTEN_variants$codon_SV %in% PTEN_hotspots$codon &
                                PTEN_variants$variantClass_v2 == "truncation" &
                                PTEN_variants$codingType_SV == "frameshift"] <-
  1
PTEN_variants$syn_mut[!PTEN_variants$codon_SV %in% PTEN_hotspots$codon &
                        PTEN_variants$codingType_SV == "synonymous"] <-
  1
PTEN_variants$lof_mut[PTEN_variants$tot_effect == "lof"] <- 1
PTEN_variants$wt_mut[PTEN_variants$tot_effect == "wt"] <- 1


PTEN_variants <- PTEN_variants %>%
  mutate(
    htsp_cmb = ifelse((codon_SV %in% top_htsps_mtl &
                         Assigned == "MT-L") |
                        (codon_SV %in% top_htsps_mth &
                           Assigned == "MT-H") &
                        htsp_cmb == 1,
                      codon_SV,
                      htsp_cmb
    ),
    htsp_point = ifelse((codon_SV %in% top_htsps_mtl &
                           Assigned == "MT-L") |
                          (codon_SV %in% top_htsps_mth &
                             Assigned == "MT-H") &
                          htsp_point == 1,
                        codon_SV,
                        htsp_point
    ),
    htsp_trunc = ifelse((codon_SV %in% top_htsps_mtl &
                           Assigned == "MT-L") |
                          (codon_SV %in% top_htsps_mth &
                             Assigned == "MT-H") &
                          htsp_trunc == 1,
                        codon_SV,
                        htsp_trunc
    )
  )

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
  select(-htsp_point, -htsp_trunc, -syn_mut, -lof_mut) %>% 
  mutate(PTEN_AF = as.numeric(PTEN_AF)) %>% 
  filter(!is.na(PTEN_AF))

other_genes_AF <- FMI_variants %>%
  filter(gene %in% c("APC", "KRAS", "TP53", "PIK3CA"), did %in% PTEN_AF$did) %>%
  mutate(gene = paste(gene, "_AF", sep = "")) %>%
  group_by(did, gene) %>%
  summarise(AF = max(alleleFreq_SV)) %>%
  ungroup() %>%
  pivot_wider(names_from = gene, values_from = AF) %>% 
  mutate(across(matches("AF"), .fns = function(x) as.numeric(x)))

# Here we divide PTEN AF in other gene AF --------------------------------------
PTEN_AF_comparison_VP <-
  left_join(PTEN_AF, other_genes_AF, by = "did") %>%
  mutate(
    PTEN_KRAS = PTEN_AF / KRAS_AF,
    PTEN_TP53 = PTEN_AF / TP53_AF,
    PTEN_APC = PTEN_AF / APC_AF,
    PTEN_PIK3CA = PTEN_AF / PIK3CA_AF
  ) %>%
  select(-TP53_AF,
         -APC_AF,
         -PIK3CA_AF,
         -KRAS_AF,
         -gene,
         -codon_SV,
         -PTEN_AF) %>%
  pivot_longer(
    names_to = "gene",
    cols = c("PTEN_KRAS", "PTEN_TP53", "PTEN_APC", "PTEN_PIK3CA"),
    values_to = "AF_adjusted"
  ) %>%
  filter(!is.na(AF_adjusted), !is.na(htsp_cmb))

PTEN_AF_comparison_VP_lof_wt <- PTEN_AF_comparison_VP %>% 
  filter(str_detect(htsp_cmb, "lof|wt"))

# Function for p-value statistics ----------------------------------------------
# of interaction between hotspots and other categories
get_htsp_int_table_by_htsp <-
  function(htsp_cmb_name, Assigned_selected) {
    int_data <-
      filter(PTEN_AF_comparison_VP,
             htsp_cmb == htsp_cmb_name,
             Assigned == Assigned_selected)
    genes <- sort(unique(int_data$gene))
    out_matrix <-
      matrix(data = NA,
             nrow = length(genes),
             ncol = length(genes)) %>%
      as.data.frame()
    rownames(out_matrix) <- genes
    colnames(out_matrix) <- genes
    t_test_matrix <- out_matrix
    ks_test_matrix <- out_matrix
    for (i in 1:length(genes)) {
      for (j in 1:length(genes)) {
        if (i != j &
            (nrow(filter(
              int_data, gene == colnames(out_matrix)[i]
            )) > 5) &
            (nrow(filter(
              int_data, gene == colnames(out_matrix)[j]
            )) > 5)) {
          t_test_matrix[i, j] <-
            round(t.test(
              filter(int_data, gene == colnames(out_matrix)[i])$AF_adjusted,
              filter(int_data, gene == colnames(out_matrix)[j])$AF_adjusted
            )$p.value,
            4)
          ks_test_matrix[i, j] <-
            round(ks.test(
              filter(int_data, gene == colnames(out_matrix)[i])$AF_adjusted,
              filter(int_data, gene == colnames(out_matrix)[j])$AF_adjusted
            )$p.value,
            4)
        }
      }
    }
    out_matrix[upper.tri(out_matrix)] <-
      ks_test_matrix[upper.tri(ks_test_matrix)]
    out_matrix[lower.tri(out_matrix)] <-
      t_test_matrix[lower.tri(t_test_matrix)]
    count_htsp <- table(int_data$gene)
    out_matrix <- rbind(out_matrix, count_htsp)
    rownames(out_matrix)[which(rownames(out_matrix) == tail(rownames(out_matrix), 1))] <-
      "count"
    return(out_matrix)
  }

# Function to create table with boxplot stats ----------------------------------
get_htsp_stat_table <- function(gene_name, Assigned_selected) {
  int_data <-
    filter(PTEN_AF_comparison_VP,
           gene == gene_name,
           Assigned == Assigned_selected) %>%
    mutate(htsp_cmb = as.factor(htsp_cmb))
  hotspots_AlleleFreq <-
    as.data.frame(matrix(0, nrow = 5, ncol = length(sort(
      unique(int_data$htsp_cmb)
    ))))
  colnames(hotspots_AlleleFreq) <- sort(unique(int_data$htsp_cmb))
  rownames(hotspots_AlleleFreq) <-
    c("lower_whisker",
      "lower_hinge",
      "median",
      "upper_hinge",
      "upper_whisker")
  for (i in 1:5) {
    for (j in 1:ncol(hotspots_AlleleFreq)) {
      hotspots_AlleleFreq[i, j] <-
        boxplot.stats(filter(int_data, htsp_cmb == colnames(hotspots_AlleleFreq)[j])$AF_adjusted)$stats[i]
    }
  }
  hotspots_AlleleFreq <- round(hotspots_AlleleFreq, 1)
  return(hotspots_AlleleFreq)
}

# Function to create violin+boxplot --------------------------------------------
get_htsp_int_plot_by_htsp <-
  function(htsp_cmb_name, Assigned_selected) {
    int_data <-
      filter(PTEN_AF_comparison_VP,
             htsp_cmb == htsp_cmb_name,
             Assigned == Assigned_selected) %>%
      mutate(gene = as.factor(gene), gene = str_remove_all(gene, "PTEN_"))
    out_plot <- int_data %>%
      ggplot(aes(y = AF_adjusted, x = gene)) +
      geom_violin(trim = T,
                  fill = "#34cfeb",
                  alpha = .6) +
      geom_boxplot(width = 0.1, col = "blue") +
      coord_cartesian(ylim = c(0, 3)) +
      theme_minimal() +
      theme(
        axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold")
      ) +
      ylab("AF adjusted") +
      xlab("Gene")
    return(out_plot)
  }

get_htsp_stat_table_by_htsp <-
  function(htsp_cmb_name, Assigned_selected) {
    int_data <-
      filter(PTEN_AF_comparison_VP,
             htsp_cmb == htsp_cmb_name,
             Assigned == Assigned_selected) %>%
      mutate(gene = as.factor(gene))
    if (nrow(int_data) > 10) {
      hotspots_AlleleFreq <-
        as.data.frame(matrix(0, nrow = 5, ncol = length(sort(
          unique(int_data$gene)
        ))))
      colnames(hotspots_AlleleFreq) <- sort(unique(int_data$gene))
      rownames(hotspots_AlleleFreq) <-
        c("lower_whisker",
          "lower_hinge",
          "median",
          "upper_hinge",
          "upper_whisker")
      for (i in 1:5) {
        for (j in 1:ncol(hotspots_AlleleFreq)) {
          hotspots_AlleleFreq[i, j] <-
            boxplot.stats(filter(int_data, gene == colnames(hotspots_AlleleFreq)[j])$AF_adjusted)$stats[i]
        }
      }
      hotspots_AlleleFreq <- round(hotspots_AlleleFreq, 1)
      return(hotspots_AlleleFreq)
    }
  }


# Figure 3, Panel H ------------------------------------------------------------
get_htsp_int_plot_by_htsp("lof", "MT-L")

ggsave(
  filename = paste0("../Charts, data, statistics/", "LoF", " ", "MT-L", ".png"),
  units = "in",
  dpi = 500,
  width = 7.65,
  height = 5.94
)

get_htsp_int_plot_by_htsp("wt", "MT-L")

ggsave(
  filename = paste0("../Charts, data, statistics/", "WT", " ", "MT-L", ".png"),
  units = "in",
  dpi = 500,
  width = 7.65,
  height = 5.94
)

get_htsp_int_table_by_htsp("lof", "MT-L") %>%
  write.csv(
    file = paste0("../Charts, data, statistics/", "LoF", " ", "MT-L", ".csv"),
    row.names = T
  )

get_htsp_int_table_by_htsp("wt", "MT-L") %>%
  write.csv(
    file = paste0("../Charts, data, statistics/", "WT", " ", "MT-L", ".csv"),
    row.names = T
  )

get_htsp_stat_table_by_htsp("lof", "MT-L") %>%
  write.csv(
    file = paste0(
      "../Charts, data, statistics/",
      "LoF",
      " ",
      "MT-L",
      " ",
      "boxplot",
      ".csv"
    ),
    row.names = T
  )

get_htsp_stat_table_by_htsp("wt", "MT-L") %>%
  write.csv(
    file = paste0(
      "../Charts, data, statistics/",
      "WT",
      " ",
      "MT-L",
      " ",
      "boxplot",
      ".csv"
    ),
    row.names = T
  )
