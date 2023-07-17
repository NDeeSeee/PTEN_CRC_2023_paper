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

# PTEN LOF data (by 2020 PTEN project)
PTEN_lof_data <-
  fread("../../../Common raw data/PTEN_LoF_annotation_2022.csv")

all_htsps <-
  fread("../../../Common raw data/Hotspots by MS status.csv")


# Foundation medicine data
FMI_samples <-
  readRDS("../../../../Sources/FMI samples assigned.rds") %>% as_tibble()

FMI_variants <-
  fread("../../../../Sources/FMI variants with syn.txt",
        na.strings = c(".", "NA", "<NA>")) %>%
  rename(did = deidentifiedSpecimenName) %>%
  select(-specimenNumber) %>%
  left_join(y = select(FMI_samples, did, Assigned), by = "did") %>%
  mutate(
    alleleFreq_SV = as.numeric(alleleFreq_SV),
    codon_SV = as.numeric(codon_SV),
    alt = description_v2
  ) %>%
  left_join(select(PTEN_lof_data, did, alt, tot_effect), by = c("did", "alt"))

FMI_variants <-
  FMI_variants[-which(duplicated.data.frame(FMI_variants))]

PAD_samples <-
  fread("../../../Common raw data/PAD samples.csv", sep = ",") %>%
  as_tibble() %>% 
  rename(Assigned = msi, did = sample_id) %>% 
  mutate(Assigned = ifelse(Assigned == "MSS", "MT-L", Assigned))


PAD_variants <-
  fread("../../../Common raw data/PAD variants.csv", sep = ",") %>%
  as_tibble() %>%
  rename(codon_SV = codon) %>% 
  mutate(ifelse(str_detect(alt, "splice"), NA, codon_SV)) %>% 
  rename(Assigned = msi, did = sample_id) %>% 
  mutate(Assigned = ifelse(Assigned == "MSS", "MT-L", Assigned))

# Adding PTEN LoF status
FMI_variants <- FMI_variants %>%
  mutate(
    pten_lof = ifelse(
      gene == "PTEN" &
        alt %in% filter(PTEN_lof_data, tot_effect == "lof")$alt,
      "lof",
      NA
    ),
    pten_lof = ifelse(
      gene == "PTEN" &
        alt %in% filter(PTEN_lof_data, tot_effect == "wt")$alt,
      "wt",
      pten_lof
    )
  )

PAD_variants <- PAD_variants %>%
  mutate(
    pten_lof = ifelse(
      gene == "PTEN" &
        alt %in% filter(PTEN_lof_data, tot_effect == "lof")$alt,
      "lof",
      NA
    ),
    pten_lof = ifelse(
      gene == "PTEN" &
        alt %in% filter(PTEN_lof_data, tot_effect == "wt")$alt,
      "wt",
      pten_lof
    )
  )



# Execution of this code demands a lot of time ---------------------------------
# PTEN co-occurrence
# PTEN_cooccurrence <- FMI_variants %>%
#   mutate(
#     htsp_status_mss = ifelse(
#       gene == "PTEN" & Assigned == "MT-L" &
#         codon_SV %in% filter(all_htsps, Assigned == "MT-L", htsp_status == "non")$codon,
#       "non",
#       NA
#     ),
#     htsp_status_mss = ifelse(
#       gene == "PTEN" & Assigned == "MT-L" &
#         codon_SV %in% filter(all_htsps, Assigned == "MT-L", htsp_status == "other_htsp")$codon,
#       "other_htsp",
#       htsp_status_mss
#     ),
#     htsp_status_msi = ifelse(
#       gene == "PTEN" & Assigned == "MT-H" &
#         codon_SV %in% filter(all_htsps, Assigned == "MT-H", htsp_status == "non")$codon,
#       "non",
#       NA
#     ),
#     htsp_status_msi = ifelse(
#       gene == "PTEN" & Assigned == "MT-H" &
#         codon_SV %in% filter(all_htsps, Assigned == "MT-H", htsp_status == "other_htsp")$codon,
#       "other_htsp",
#       htsp_status_msi
#     )
#   ) %>%
#   group_by(did) %>%
#   summarise(
#     KRAS_presence = sign(sum(gene == "KRAS")),
#     TP53_presence = sign(sum(gene == "TP53")),
#     PIK3CA_presence = sign(sum(gene == "PIK3CA")),
#     SMAD4_presence = sign(sum(gene == "SMAD4")),
#     BRAF_presence = sign(sum(gene == "BRAF")),
#     APC_presence = sign(sum(gene == "APC")),
#     PTEN_presence = sign(sum(gene == "PTEN")),
#     PTEN_130_presence = sign(sum(gene == "PTEN" &
#                                    codon_SV == 130, na.rm = T)),
#     PTEN_136_presence = sign(sum(gene == "PTEN" &
#                                    codon_SV == 136, na.rm = T)),
#     PTEN_173_presence = sign(sum(gene == "PTEN" &
#                                    codon_SV == 173, na.rm = T)),
#     PTEN_233_presence = sign(sum(gene == "PTEN" &
#                                    codon_SV == 233, na.rm = T)),
#     PTEN_319_presence = sign(sum(gene == "PTEN" &
#                                    codon_SV == 319, na.rm = T)),
#     PTEN_323_presence = sign(sum(gene == "PTEN" &
#                                    codon_SV == 323, na.rm = T)),
#     PTEN_157_presence = sign(sum(gene == "PTEN" &
#                                    codon_SV == 157, na.rm = T)),
#     PTEN_267_presence = sign(sum(gene == "PTEN" &
#                                    codon_SV == 267, na.rm = T)),
#     PTEN_268_presence = sign(sum(gene == "PTEN" &
#                                    codon_SV == 268, na.rm = T)),
#     PTEN_mss_other_htsp_presence = sign(sum(htsp_status_mss == "other_htsp", na.rm = T)),
#     PTEN_mss_non_htsp_presence = sign(sum(htsp_status_mss == "non", na.rm = T)),
#     PTEN_msi_other_htsp_presence = sign(sum(htsp_status_msi == "other_htsp", na.rm = T)),
#     PTEN_msi_non_htsp_presence = sign(sum(htsp_status_msi == "non", na.rm = T)),
#     PTEN_synonymous_presence = sign(
#       sum(
#         gene == "PTEN" & codingType_SV == "synonymous" &
#           !codon_SV %in% all_htsps,
#         na.rm = T
#       )
#     ),
#     PTEN_lof_presence = sign(sum(tot_effect == "lof", na.rm = T)),
#     PTEN_wt_presence = sign(sum(tot_effect == "wt", na.rm = T)),
#     Assigned = unique(Assigned)[1]
#   ) %>%
#   mutate(across(contains("presence"), ~ replace_na(., 0)),
#          across(contains("presence"), as.factor))
#
#   PTEN_cooccurrence <- full_join(PTEN_cooccurrence, select(FMI_samples, did, Assigned)) %>%
#     mutate(across(.cols = matches("presence"), .fns = function(x) replace_na(x, 0)))

# write.csv(PTEN_cooccurrence, "../../../Common raw data/PTEN co-occurrence.csv", row.names = F)

PTEN_cooccurrence_FMI <-
  fread("../../../Common raw data/PTEN co-occurrence.csv")

# Execution of this code demands a lot of time ---------------------------------
# PTEN co-occurrence PAD
# PTEN_cooccurrence_PAD <- PAD_variants %>%
#   mutate(
#     htsp_status_mss = ifelse(
#       gene == "PTEN" & Assigned == "MT-L" &
#         codon_SV %in% filter(all_htsps, Assigned == "MT-L", htsp_status == "non")$codon,
#       "non",
#       NA
#     ),
#     htsp_status_mss = ifelse(
#       gene == "PTEN" & Assigned == "MT-L" &
#         codon_SV %in% filter(all_htsps, Assigned == "MT-L", htsp_status == "other_htsp")$codon,
#       "other_htsp",
#       htsp_status_mss
#     ),
#     htsp_status_msi = ifelse(
#       gene == "PTEN" & Assigned == "MT-H" &
#         codon_SV %in% filter(all_htsps, Assigned == "MT-H", htsp_status == "non")$codon,
#       "non",
#       NA
#     ),
#     htsp_status_msi = ifelse(
#       gene == "PTEN" & Assigned == "MT-H" &
#         codon_SV %in% filter(all_htsps, Assigned == "MT-H", htsp_status == "other_htsp")$codon,
#       "other_htsp",
#       htsp_status_msi
#     )
#   ) %>%
#   group_by(did) %>%
#   summarise(
#     KRAS_presence = sign(sum(gene == "KRAS")),
#     TP53_presence = sign(sum(gene == "TP53")),
#     PIK3CA_presence = sign(sum(gene == "PIK3CA")),
#     SMAD4_presence = sign(sum(gene == "SMAD4")),
#     BRAF_presence = sign(sum(gene == "BRAF")),
#     APC_presence = sign(sum(gene == "APC")),
#     PTEN_presence = sign(sum(gene == "PTEN")),
#     PTEN_130_presence = sign(sum(gene == "PTEN" &
#                                    codon_SV == 130, na.rm = T)),
#     PTEN_136_presence = sign(sum(gene == "PTEN" &
#                                    codon_SV == 136, na.rm = T)),
#     PTEN_173_presence = sign(sum(gene == "PTEN" &
#                                    codon_SV == 173, na.rm = T)),
#     PTEN_233_presence = sign(sum(gene == "PTEN" &
#                                    codon_SV == 233, na.rm = T)),
#     PTEN_319_presence = sign(sum(gene == "PTEN" &
#                                    codon_SV == 319, na.rm = T)),
#     PTEN_323_presence = sign(sum(gene == "PTEN" &
#                                    codon_SV == 323, na.rm = T)),
#     PTEN_157_presence = sign(sum(gene == "PTEN" &
#                                    codon_SV == 157, na.rm = T)),
#     PTEN_267_presence = sign(sum(gene == "PTEN" &
#                                    codon_SV == 267, na.rm = T)),
#     PTEN_268_presence = sign(sum(gene == "PTEN" &
#                                    codon_SV == 268, na.rm = T)),
#     PTEN_mss_other_htsp_presence = sign(sum(htsp_status_mss == "other_htsp", na.rm = T)),
#     PTEN_mss_non_htsp_presence = sign(sum(htsp_status_mss == "non", na.rm = T)),
#     PTEN_msi_other_htsp_presence = sign(sum(htsp_status_msi == "other_htsp", na.rm = T)),
#     PTEN_msi_non_htsp_presence = sign(sum(htsp_status_msi == "non", na.rm = T)),
#     PTEN_lof_presence = sign(sum(pten_lof == "lof", na.rm = T)),
#     PTEN_wt_presence = sign(sum(pten_lof == "wt", na.rm = T)),
#     Assigned = unique(Assigned)[1]
#   ) %>%
#   mutate(across(contains("presence"), ~ replace_na(., 0)),
#          across(contains("presence"), as.factor))
# 
# PTEN_cooccurrence_PAD <- full_join(PTEN_cooccurrence_PAD, select(PAD_samples, did, Assigned)) %>%
#   mutate(across(.cols = matches("presence"), .fns = function(x) x = as.numeric(as.character(x)))) %>% 
#   mutate(across(.cols = matches("presence"), .fns = function(x) replace_na(x, 0))) %>% 
#   mutate(across(.cols = matches("presence"), .fns = function(x) x = as.factor(x)))
# 
# write.csv(PTEN_cooccurrence_PAD, "../../../Common raw data/PTEN co-occurrence PAD.csv", row.names = F)

PTEN_cooccurrence_PAD <- fread("../../../Common raw data/PTEN co-occurrence PAD.csv") %>% 
  as_tibble()

PTEN_cooccurrence <- bind_rows(PTEN_cooccurrence_PAD, PTEN_cooccurrence_FMI)

# Glob vars --------------------------------------------------------------------
# Target genes
genes <- c("APC", "TP53", "PIK3CA", "SMAD4", "KRAS") %>%
  paste0(sep = "_presence")

htsps_mss <- c(
  filter(all_htsps, Assigned == "MT-L", htsp_status == "top_htsp")$codon,
  "mss_other_htsp",
  "mss_non_htsp",
  "lof",
  "wt",
  "synonymous"
)

htsps_msi <- c(
  filter(all_htsps, Assigned == "MT-H", htsp_status == "top_htsp")$codon,
  "msi_other_htsp",
  "msi_non_htsp",
  "lof",
  "wt",
  "synonymous"
)

htsps_mss <- paste("PTEN", htsps_mss, "presence", sep = "_")
htsps_msi <- paste("PTEN", htsps_msi, "presence", sep = "_")

# Main Function
plot_cooccurrence_vv <-
  function(hotspots,
           gene,
           df = FALSE,
           data = PTEN_cooccurrence) {
    df_for_plot_tot <- data.frame()
    for (k in MS_types) {
      odds <- bottom <- upper <- n_samples <- p_val <- c()
      for (i in hotspots) {
        if (dim(table(select(filter(
          data, Assigned == k
        ), gene, i)))[1] != 2) {
          odds <- append(odds, NA)
          upper <- append(upper, NA)
          bottom <- append(bottom, NA)
          n_samples <-
            append(n_samples, sum(table(select(
              filter(data, Assigned == k), gene, i
            ))[-1]))
          p_val <- append(p_val, NA)
        } else {
          fisher_res <-
            tryCatch(
              fisher.test(table(select(
                filter(data, Assigned == k), gene, i
              )))
            )
          odds <- round(append(odds, log2(fisher_res$estimate)), 3)
          upper <-
            round(append(upper, log2(fisher_res$conf.int[2])), 3)
          bottom <-
            round(append(bottom, log2(fisher_res$conf.int[1])), 3)
          n_samples <-
            append(n_samples, sum(table(select(
              filter(data, Assigned == k), gene, i
            ))[-1]))
          p_val <- round(append(p_val, fisher_res$p.value), 10)
        }
      }
      df_for_plot <- data.frame(
        hotspot = hotspots,
        log2_odds = odds,
        upper = upper,
        bottom = bottom,
        n_samples = n_samples,
        p_value = p_val
      ) %>%
        as_tibble() %>%
        mutate(Assigned = k)
      df_for_plot_tot <- rbind(df_for_plot_tot, df_for_plot) %>%
        mutate(hotspot = str_remove_all(str_remove_all(hotspot, "_presence"), "PTEN_"))
    }
    #  df_for_plot_tot <- df_for_plot_tot %>%
    #    mutate(gene = str_extract(gene, "[:alnum:]+"), int_type = factor(cut(sign(log2_odds), c(-1, 0, 1), include.lowest = T), labels = c("co-exc", "co-occ")),
    #           signif = factor(abs(sign(upper)-sign(bottom)), labels = c("yes", "no")))
    ans_plot <-
      ggplot(df_for_plot_tot, aes(x = hotspot, y = log2_odds, col = Assigned)) +
      geom_point(position = position_dodge(0.05)) +
      theme(axis.text.x = element_text(angle = 90)) +
      geom_errorbar(aes(ymin = bottom, ymax = upper),
                    width = .2,
                    position = position_dodge(0.05)) +
      coord_cartesian(ylim = c(-3, 4.5))
    if (df) {
      return(df_for_plot_tot)
    } else {
      return(ans_plot)
    }
  }

# MT-L -------------------------------------------------------------------------
temp_cooc_table <- data.frame()

MS_types <- "MT-L"

for (i in genes) {
  temp_cooc_table <- bind_rows(temp_cooc_table,
                               mutate(
                                 plot_cooccurrence_vv(
                                   hotspots = htsps_mss,
                                   gene = i,
                                   df = T,
                                   data = PTEN_cooccurrence
                                 ),
                                 gene = i
                               ))
}

temp_cooc_table <- as_tibble(temp_cooc_table) %>%
  mutate(
    gene = str_remove_all(gene, "_presence"),
    hotspot = as.factor(hotspot),
    gene = as.factor(gene),
    log2_odds = ifelse(str_detect(log2_odds, "Inf"), NA, log2_odds),
    log2_odds_signif = ifelse(p_value > 0.005, NA, log2_odds),
    log2_odds_signif_cutoff = min(abs(c(
      max(log2_odds_signif, na.rm = T),
      min(log2_odds_signif, na.rm = T)
    ))),
    log2_odds_signif = ifelse(
      log2_odds_signif > log2_odds_signif_cutoff,
      log2_odds_signif_cutoff,
      log2_odds_signif
    ),
    log2_odds_signif = ifelse(
      log2_odds_signif < -log2_odds_signif_cutoff,
      -log2_odds_signif_cutoff,
      log2_odds_signif
    ),
    log2_odds = ifelse(
      log2_odds > log2_odds_signif_cutoff,
      log2_odds_signif_cutoff,
      log2_odds
    ),
    log2_odds = ifelse(
      log2_odds < -log2_odds_signif_cutoff,
      -log2_odds_signif_cutoff,
      log2_odds
    ),
    p_value = ifelse(p_value == 0, p_value + 0.0000000001, p_value),
    log2_pval = log(p_value, 10),
    log2_pval_noabs = ifelse(log2_odds > 0,
                             -1 * log2_pval,
                             log2_pval)
  )

log2_odds_signif_cutoff <-
  temp_cooc_table$log2_odds_signif_cutoff[1]

write.csv(
  temp_cooc_table,
  "../Charts, data, statistics/PTEN with ALL genes co-occ MT-L.csv",
  row.names = F
)

temp_cooc_table %>%
  ggplot(aes(y = hotspot, x = gene, fill = log2_odds)) +
  geom_tile() +
  theme(
    axis.text.x = element_text(angle = 90),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    legend.title = element_blank()
  ) +
  scale_fill_gradientn(
    values = c(-0.0001, last(scales::rescale(
      c(filter(temp_cooc_table, !is.na(log2_odds))$log2_odds, 0)
    )), 1.00001),
    colours = c("salmon", "white", "royalblue4")
  ) +
  ylab("PTEN") +
  xlab("Gene")

ggsave(
  "../Charts, data, statistics/Panel B.png",
  dpi = 400,
  height = 4.2,
  width = 5.2,
  units = "in"
)


temp_cooc_table %>%
  ggplot(aes(y = hotspot, x = gene, fill = log2_odds_signif)) +
  geom_tile() +
  theme(
    axis.text.x = element_text(angle = 90),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    legend.title = element_blank()
  ) +
  scale_fill_gradientn(
    values = c(-0.0001, last(scales::rescale(
      c(filter(
        temp_cooc_table, !is.na(log2_odds_signif)
      )$log2_odds_signif, 0)
    )), 1.00001),
    colours = c("salmon", "white", "royalblue4")
  ) +
  ylab("PTEN") +
  xlab("Gene")


ggsave(
  "../Charts, data, statistics/Panel C.png",
  dpi = 400,
  height = 4.2,
  width = 5.2,
  units = "in"
)


# MT-H -------------------------------------------------------------------------
temp_cooc_table <- data.frame()

MS_types <- "MT-H"

for (i in genes) {
  temp_cooc_table <- bind_rows(temp_cooc_table,
                               mutate(
                                 plot_cooccurrence_vv(
                                   hotspots = htsps_msi,
                                   gene = i,
                                   df = T,
                                   data = PTEN_cooccurrence
                                 ),
                                 gene = i
                               ))
}

temp_cooc_table <- as_tibble(temp_cooc_table) %>%
  mutate(
    gene = str_remove_all(gene, "_presence"),
    hotspot = as.factor(hotspot),
    gene = as.factor(gene),
    # log2_odds_cutoff = min(abs(c(max(log2_odds, na.rm = T),
    #                                     min(log2_odds, na.rm = T)))),
    log2_odds = ifelse(str_detect(log2_odds, "Inf"), NA, log2_odds),
    log2_odds = ifelse(
      log2_odds > log2_odds_signif_cutoff,
      log2_odds_signif_cutoff,
      log2_odds
    ),
    log2_odds = ifelse(
      log2_odds < -log2_odds_signif_cutoff,
      -log2_odds_signif_cutoff,
      log2_odds
    )
  )

write.csv(
  temp_cooc_table,
  "../Charts, data, statistics/PTEN with ALL genes co-occ MT-H.csv",
  row.names = F
)

temp_cooc_table %>%
  ggplot(aes(y = hotspot, x = gene, fill = log2_odds)) +
  geom_tile() +
  theme(
    axis.text.x = element_text(angle = 90),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    legend.title = element_blank()
  ) +
  scale_fill_gradientn(
    values = c(-0.0001, last(scales::rescale(
      c(filter(temp_cooc_table, !is.na(log2_odds))$log2_odds, 0)
    )), 1.00001),
    colours = c("salmon", "white", "royalblue4")
  ) +
  ylab("PTEN") +
  xlab("Gene")

ggsave(
  "../Charts, data, statistics/Panel D.png",
  dpi = 400,
  height = 4.2,
  width = 5.2,
  units = "in"
)

