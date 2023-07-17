# Attach requirement packages and setting WD -----------------------------------
packages_names <- c("readxl", "data.table", "reshape2", "tidyverse",
                    "rstudioapi")

lapply(packages_names, require, character.only = TRUE)

setwd(dirname(getActiveDocumentContext()$path))

# Load FMI data with LoF and hotspots ------------------------------------------
FMI_samples <-
  readRDS("../../../../Sources/FMI samples assigned.rds") %>%
  as_tibble()
FMI_variants <-
  fread("../../../../Sources/FMI variants with syn.txt") %>%
  as_tibble() %>%
  mutate(alt = description_v2) %>%
  rename(did = deidentifiedSpecimenName, codon = codon_SV) %>%
  left_join(y = select(FMI_samples, did, Assigned), by = "did")

# File with PTEN hotspots grouped by microsatellite status
PTEN_hotspots_MTL <-
  fread("../../../Common raw data/Hotspots by MS status.csv") %>%
  filter(Assigned == "MT-L") %>%
  select(codon, htsp_status) %>%
  left_join(x = data.frame(codon = 1:403), by = "codon") %>%
  mutate(htsp_status = replace_na(htsp_status, "non")) %>%
  as_tibble()

# File with PTEN annotated with LoF by previous methods
PTEN_lof_data <-
  fread("../../../Common raw data/PTEN LoF annotation.csv") %>%
  rename(codon = codon_SV) %>%
  filter(Assigned == "MT-L",
         !str_detect(alt, "splice|del|_|ext|ins"))

# Define several codon classes -------------------------------------------------
PTEN_variants <- FMI_variants %>%
  filter(gene == "PTEN", codon != ".", Assigned == "MT-L") %>%
  select(did, alt, alleleFreq_SV, codingType_SV, variantClass_v2) %>%
  filter(!str_detect(alt, "splice|del|_|ext|ins")) %>%
  full_join(
    PTEN_lof_data,
    by = c("did", "alt", "codingType_SV", "variantClass_v2"),
    multiple = "all"
  ) %>%
  mutate(codon = ifelse(is.na(codon),
                        as.numeric(str_extract(alt, "[:digit:]+")),
                        codon)) %>%
  left_join(PTEN_hotspots_MTL, by = "codon") %>%
  mutate(
    tot_effect = ifelse(codingType_SV == "synonymous", "wt", tot_effect),
    syn_mut = ifelse(codingType_SV == "synonymous", 1, 0),
    lof_mut = ifelse(tot_effect == "lof", 1, 0),
    wt_mut = ifelse(tot_effect == "wt", 1, 0),
    non_htsp_mut = ifelse(htsp_status == "non", 1, 0),
    htsp_mut = ifelse(htsp_status == "other_htsp", 1, 0),
    top_htsp_mut = ifelse(htsp_status == "top_htsp", codon, 0)
  ) %>%
  pivot_longer(
    cols = c(syn_mut, lof_mut, wt_mut,
             non_htsp_mut, htsp_mut, top_htsp_mut),
    names_to = "codon_type",
    values_to = "value"
  ) %>%
  mutate(
    alleleFreq_SV = as.numeric(alleleFreq_SV),
    codon_type = ifelse(value != 0 & value != 1, value, codon_type),
    codon_type = str_remove_all(codon_type, "_mut"),
    codon_type = as.factor(codon_type)
  ) %>%
  filter(value != 0) %>%
  select(-did, -value)

# Save chart data --------------------------------------------------------------
PTEN_variants %>%
  write.csv("../Charts, data, statistics/Panel I data.csv", row.names = F)

# Draw the chart Panel C -------------------------------------------------------
PTEN_variants %>%
  ggplot(aes(x = codon_type, y = alleleFreq_SV)) +
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
  xlab("Mutation type")

# Save image -------------------------------------------------------------------
ggsave(
  "../Charts, data, statistics/Panel I.png",
  units = "in",
  dpi = 500,
  width = 11.65,
  height = 5.94
)

# Statistics -------------------------------------------------------------------

my_results_MTL <- as.data.frame(matrix(0,
                                       nrow = length(sort(
                                         unique(PTEN_variants$codon_type)
                                       )),
                                       ncol = length(sort(
                                         unique(PTEN_variants$codon_type)
                                       ))))
colnames(my_results_MTL) <- sort(unique(PTEN_variants$codon_type))
rownames(my_results_MTL) <- sort(unique(PTEN_variants$codon_type))

my_results_ks_test <- my_results_ttest <- my_results_MTL

for (i in 1:nrow(my_results_MTL)) {
  for (j in 1:nrow(my_results_MTL)) {
    my_results_ttest[i, j] <- t.test(
      filter(PTEN_variants, codon_type == colnames(my_results_MTL)[i])$alleleFreq_SV,
      filter(PTEN_variants, codon_type == colnames(my_results_MTL)[j])$alleleFreq_SV
    )$p.value
    my_results_ks_test[i, j] <- ks.test(
      filter(PTEN_variants, codon_type == colnames(my_results_MTL)[i])$alleleFreq_SV,
      filter(PTEN_variants, codon_type == colnames(my_results_MTL)[j])$alleleFreq_SV
    )$p.value
  }
}

diag(my_results_MTL) <- NA
my_results_MTL[upper.tri(my_results_MTL)] <-
  my_results_ks_test[upper.tri(my_results_MTL)]
my_results_MTL[lower.tri(my_results_MTL)] <-
  my_results_ttest[lower.tri(my_results_MTL)]

write.csv(my_results_MTL,
          "../Charts, data, statistics/Panel I statistics.csv")

hotspots_AlleleFreq <- as.data.frame(matrix(0,
                                            nrow = 5,
                                            ncol = length(sort(
                                              unique(PTEN_variants$codon_type)
                                            ))))
colnames(hotspots_AlleleFreq) <-
  sort(unique(PTEN_variants$codon_type))
rownames(hotspots_AlleleFreq) <- c("lower_whisker",
                                   "lower_hinge",
                                   "median",
                                   "upper_hinge",
                                   "upper_whisker")
for (i in 1:5) {
  for (j in 1:ncol(hotspots_AlleleFreq)) {
    hotspots_AlleleFreq[i, j] <- boxplot.stats(filter(
      PTEN_variants,
      codon_type == colnames(hotspots_AlleleFreq)[j]
    )$alleleFreq_SV)$stats[i]
  }
}
hotspots_AlleleFreq <- round(hotspots_AlleleFreq, 2)

write.csv(hotspots_AlleleFreq,
          "../Charts, data, statistics/Panel I boxplot statistics.csv")
