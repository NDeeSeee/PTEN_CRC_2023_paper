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

PTEN_miss_hotspots <- c(15, 16, 17, 36, 93, 136)

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
  fread("../../../Common raw data/true_PTEN_hotspots.csv") %>% 
  filter(!codon %in% PTEN_miss_hotspots) %>% 
  as_tibble()

# Loss of function -------------------------------------------------------------
PTEN_lof_data <-
  fread("../../../Common raw data/PTEN LoF annotation.csv")

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

# all hotspots are marked as " 1 "
PTEN_variants$htsp_cmb[PTEN_variants$codon_SV %in% PTEN_hotspots$codon] <-
  1

# Then assign other categories
PTEN_variants$htsp_point[PTEN_variants$codon_SV %in% PTEN_hotspots$codon &
                           PTEN_variants$variantClass_v2 == "point"] <- 1
PTEN_variants$htsp_trunc[PTEN_variants$codon_SV %in% PTEN_hotspots$codon &
                           PTEN_variants$variantClass_v2 == "truncation"] <- 1
PTEN_variants$syn_mut[!PTEN_variants$codon_SV %in% PTEN_hotspots$codon &
                        PTEN_variants$codingType_SV == "synonymous"] <- 1
PTEN_variants$lof_mut[PTEN_variants$tot_effect == "lof"] <- 1
PTEN_variants$wt_mut[PTEN_variants$tot_effect == "wt"] <- 1

# Save Panel G -----------------------------------------------------------------
Panel_G_data <- PTEN_variants %>%
  filter(!is.na(Assigned), alleleFreq_SV != ".") %>%
  mutate(
    htsp_point = as.factor(ifelse(htsp_point == 1, "htsp", "non")),
    alleleFreq_SV = as.numeric(alleleFreq_SV),
    Assigned = as.factor(Assigned)
  ) %>% 
  select(htsp_point, alleleFreq_SV, Assigned)

MS_types <- unique(Panel_G_data$Assigned)

t_test_list <- list()
for (i in 1:length(MS_types)) {
  temp_data <- filter(Panel_G_data, Assigned == MS_types[i])
  temp_test <- t.test(filter(temp_data, htsp_point == "non")$alleleFreq_SV,
                      filter(temp_data, htsp_point == "htsp")$alleleFreq_SV)
  temp_test$data.name <- paste0(MS_types[i], ", x: non-hotspots, y: hotspots", collapse = "")
  t_test_list[[i]] <- temp_test
}

saveRDS(t_test_list, file = "../Charts, data, statistics/Panel G statistics.RData")

Panel_G_data %>% 
  ggplot(aes(
    x = Assigned,
    y = alleleFreq_SV,
    col = htsp_point,
    fill = Assigned
  )) +
  geom_violin(trim = F,
              show.legend = T,
              linewidth = 0.7) +
  geom_boxplot(width = 0.15, position = position_dodge(.9)) +
  ylim(c(0, 1)) +
  xlab("PTEN point mutations") +
  ylab("Allele Frequency") +
  theme_minimal() +
  theme(legend.position = "top", legend.title = element_blank()) +
  theme(
    legend.title = element_blank(),
    legend.key = element_blank(),
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(
    values = c(
      "MT-L htmb" = "#FFBF15",
      "MT-H" = "#E4ACFA",
      "MT-L" = "#489FF8"
    ),
    guide = guide_legend(override.aes = list(fill = c(
      "#FFBF15",
      "#E4ACFA",
      "#489FF8"
    )))
  ) +
  scale_color_manual(values = c("htsp" = "red", "non" = "gray40"))

ggsave(
  "../Charts, data, statistics/Panels G-H/Panel G.png",
  units = "in",
  dpi = 500,
  width = 6.65,
  height = 5.94
)


# Save Panel I -----------------------------------------------------------------
Panel_I_data <- PTEN_variants %>%
  filter(!is.na(Assigned), alleleFreq_SV != ".") %>%
  mutate(
    htsp_trunc = as.factor(ifelse(htsp_trunc == 1, "htsp", "non")),
    alleleFreq_SV = as.numeric(alleleFreq_SV),
    Assigned = as.factor(Assigned)
  ) %>% 
  select(htsp_trunc, alleleFreq_SV, Assigned)

MS_types <- unique(Panel_I_data$Assigned)

t_test_list <- list()
for (i in 1:length(MS_types)) {
  temp_data <- filter(Panel_I_data, Assigned == MS_types[i])
  temp_test <- t.test(filter(temp_data, htsp_trunc == "non")$alleleFreq_SV,
                      filter(temp_data, htsp_trunc == "htsp")$alleleFreq_SV)
  temp_test$data.name <- paste0(MS_types[i], ", x: non-hotspots, y: hotspots", collapse = "")
  t_test_list[[i]] <- temp_test
}

saveRDS(t_test_list, file = "../Charts, data, statistics/Panel I statistics.RData")

Panel_I_data %>%
  ggplot(aes(
    x = Assigned,
    y = alleleFreq_SV,
    col = htsp_trunc,
    fill = Assigned
  )) +
  geom_violin(trim = F,
              show.legend = T,
              linewidth = 0.7) +
  geom_boxplot(width = 0.15, position = position_dodge(.9)) +
  ylim(c(0, 1)) +
  xlab("PTEN point mutations") +
  ylab("Allele Frequency") +
  theme_minimal() +
  theme(legend.position = "top", legend.title = element_blank()) +
  theme(
    legend.title = element_blank(),
    legend.key = element_blank(),
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(
    values = c(
      "MT-L htmb" = "#FFBF15",
      "MT-H" = "#E4ACFA",
      "MT-L" = "#489FF8"
    ),
    guide = guide_legend(override.aes = list(fill = c(
      "#FFBF15",
      "#E4ACFA",
      "#489FF8"
    )))
  ) +
  scale_color_manual(values = c("htsp" = "red", "non" = "gray40"))

ggsave(
  "../Charts, data, statistics/Panels G-H/Panel H.png",
  units = "in",
  dpi = 500,
  width = 6.65,
  height = 5.94
)

# Save Panel G-H data
select(PTEN_variants, alleleFreq_SV, htsp_trunc, htsp_point, Assigned) %>%
  write.csv("../Charts, data, statistics/Panels G-H/Panels G-H data.csv",
            row.names = F)
