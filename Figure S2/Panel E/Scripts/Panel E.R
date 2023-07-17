# Attach requirement packages and setting WD -----------------------------------
packages_names <- c(
  "tidyverse", 
  "data.table", 
  "readxl", 
  "reshape2", 
  "rstudioapi",
  "Biostrings",
  "rlist",
  "seqinr",
  "ape",
  "insect",
  "biomaRt"
)

lapply(packages_names, require, character.only = TRUE)

setwd(dirname(getActiveDocumentContext()$path))

rename <- dplyr::rename
select <- dplyr::select
filter <- dplyr::filter
group_by <- dplyr::group_by
mutate <- dplyr::mutate
translate <- Biostrings::translate

# Load FMI data ----------------------------------------------------------------
FMI_samples <-
  readRDS("../../../../Sources/FMI samples assigned.rds") %>%
  as_tibble()

FMI_variants <-
  fread("../../../../Sources/FMI variants with syn.txt") %>%
  as_tibble() %>%
  mutate(alt = description_v2) %>%
  rename(did = deidentifiedSpecimenName) %>%
  left_join(y = select(FMI_samples, did, Assigned), by = "did")

PTEN_variants_syn_mtl <-
  filter(FMI_variants,
         gene == "PTEN",
         codingType_SV == "synonymous",
         Assigned == "MT-L")

PTEN_variants_syn_mth <-
  filter(FMI_variants,
         gene == "PTEN",
         codingType_SV == "synonymous",
         Assigned == "MT-H")

# PTEN open reading frame ------------------------------------------------------
PTEN_ORF <-
  DNAString(str_remove_all(read_file("../../../Common raw data/PTEN ORF.txt"), "\n"))

# Complement chain
PTEN_ORF_comp <- reverseComplement(PTEN_ORF)

# Read trinucleotide context files obtained from -------------------------------
# PTEN hotspot freq prediction.R

tricontext_prob_MTH <-
  fread("../../../Common raw data/tricontext_MTH_probs_mut_sign.csv") %>%
  as_tibble()

tricontext_prob_MTL <-
  fread("../../../Common raw data/tricontext_MTL_prob_mut_sign.csv") %>%
  as_tibble()


# Function to get predict synonymous mutations distribution --------------------
# for predict frequency distribution
get_pred_dist <- function(filename, ORF) {
  tricontext <- filename
  colnames(tricontext) <- c("tricontext", "rel_prob")
  
  ORF_comp <- reverseComplement(ORF)
  
  tricontext <- tricontext %>%
    mutate(
      tricontext_sub = str_remove_all(str_extract_all(tricontext, "\\[.+\\]"), "\\[|\\]"),
      tricontext = str_remove_all(str_remove_all(tricontext, ">.."), "\\[")
    ) %>%
    rename(trinucleotideContext_SV = tricontext,
           trinucleotideAlteration_SV = tricontext_sub) %>%
    mutate(
      trinucleotideAltContext_SV = paste0(
        str_sub(trinucleotideContext_SV, 1, 1),
        str_sub(trinucleotideAlteration_SV, 3, 3),
        str_sub(trinucleotideContext_SV, 3, 3)
      ),
      ref = str_sub(trinucleotideAlteration_SV, 1, 1),
      alt = str_sub(trinucleotideAlteration_SV, 3, 3)
    ) %>%
    as_tibble()
  
  protein_ref <- translate(ORF)
  
  # Get coordinates for forward strand
  tnt_coords_forward <- list()
  for (i in 1:nrow(tricontext)) {
    tnt_coords_forward <- list.append(
      tnt_coords_forward,
      gregexpr2(tricontext$trinucleotideContext_SV[i], toString(ORF))[[1]]
    )
  }
  
  # Locations at reverse strand
  tnt_coords_reverse <- list()
  for (i in 1:nrow(tricontext)) {
    tnt_coords_reverse <- list.append(
      tnt_coords_reverse,
      gregexpr2(tricontext$trinucleotideContext_SV[i], toString(ORF_comp))[[1]]
    )
  }
  
  # Get dN/dS ratio
  non_syn <- syn <- 0
  non_syn_prob <- syn_prob <- 0
  syn_coord <- syn_coord_prob <- c()
  non_syn_coord <- non_syn_coord_prob <- c()
  
  
  
  # Forward chain
  for (tnt_con in 1:length(tnt_coords_forward)) {
    for (tnt_coord in 1:length(tnt_coords_forward[[tnt_con]])) {
      if (tnt_coords_forward[[tnt_con]][tnt_coord][1] != -1) {
        ORF_temp <- ORF
        ORF_temp <-
          replaceLetterAt(ORF_temp, tnt_coords_forward[[tnt_con]][tnt_coord][1] + 1, tricontext$alt[tnt_con])
        ORF_temp <- translate(ORF_temp)
        if (ORF_temp == protein_ref) {
          syn_prob <- syn_prob + tricontext$rel_prob[tnt_con]
          syn <- syn + 1
          syn_coord <-
            c(syn_coord, tnt_coords_forward[[tnt_con]][tnt_coord][1] + 1)
          syn_coord_prob <-
            c(syn_coord_prob, round(tricontext$rel_prob[tnt_con], 2))
        } else {
          non_syn_prob <- non_syn_prob + tricontext$rel_prob[tnt_con]
          non_syn <- non_syn + 1
          non_syn_coord <-
            c(non_syn_coord, tnt_coords_forward[[tnt_con]][tnt_coord][1] + 1)
          non_syn_coord_prob <-
            c(non_syn_coord_prob, round(tricontext$rel_prob[tnt_con], 2))
        }
      }
    }
  }
  
  # Reverse chain
  for (tnt_con in 1:length(tnt_coords_reverse)) {
    for (tnt_coord in 1:length(tnt_coords_reverse[[tnt_con]])) {
      if (tnt_coords_reverse[[tnt_con]][tnt_coord][1] != -1) {
        ORF_temp <- ORF_comp
        ORF_temp <-
          replaceLetterAt(ORF_temp, tnt_coords_reverse[[tnt_con]][tnt_coord][1] + 1, tricontext$alt[tnt_con])
        ORF_temp <- reverseComplement(ORF_temp)
        ORF_temp <- translate(ORF_temp)
        if (ORF_temp == protein_ref) {
          syn_prob <- syn_prob + tricontext$rel_prob[tnt_con]
          syn <- syn + 1
          syn_coord <-
            c(syn_coord, length(ORF) - tnt_coords_reverse[[tnt_con]][tnt_coord][1])
          syn_coord_prob <-
            c(syn_coord_prob, round(tricontext$rel_prob[tnt_con], 2))
        } else {
          non_syn_prob <- non_syn_prob + tricontext$rel_prob[tnt_con]
          non_syn <- non_syn + 1
          non_syn_coord <-
            c(non_syn_coord, length(ORF) - tnt_coords_reverse[[tnt_con]][tnt_coord][1])
          non_syn_coord_prob <-
            c(non_syn_coord_prob, round(tricontext$rel_prob[tnt_con], 2))
        }
      }
    }
  }
  
  # Getting results
  predicted_mut_freq <-
    tibble(syn_coord = syn_coord, syn_coord_prob = syn_coord_prob) %>%
    group_by(syn_coord) %>%
    summarise(syn_coord_prob_sum = sum(syn_coord_prob),
              syn_n = n())
  
  
    return(predicted_mut_freq)
}

# Get predicted nucleotide frequency processing --------------------------------

get_pred_dist_process <- function(data) {
  processed_data <-
    get_pred_dist(filename = data,
                  ORF = PTEN_ORF) %>%
    rename(
      nucleotide_position = syn_coord,
      freq_pred = syn_coord_prob_sum,
      variants = syn_n
    ) %>%
    mutate(nucleotide_position = as.numeric(as.character(nucleotide_position))) %>%
    left_join(x = tibble(nucleotide_position = 1:length(PTEN_ORF)),
              by = "nucleotide_position") %>%
    mutate(freq_pred = replace_na(freq_pred, 0),
           variants = replace_na(variants, 0))
  return(processed_data)
}

# Separately predict values for MT-L and MT-H cohorts
syn_distr_predict_mtl <- get_pred_dist_process(tricontext_prob_MTH)

syn_distr_predict_mth <- get_pred_dist_process(tricontext_prob_MTL)

# Get real nucleotide frequency function ---------------------------------------

get_real_dist <- function(data) {
  real_dist <- select(data, transcriptEffect_SV) %>%
    as_vector() %>%
    str_extract("[:digit:]+") %>%
    as.numeric() %>%
    table() %>%
    as.data.frame() %>%
    as_tibble() %>%
    rename(nucleotide_position = ".", freq_real = Freq) %>%
    mutate(nucleotide_position = as.numeric(as.character(nucleotide_position))) %>%
    left_join(x = tibble(nucleotide_position = 1:length(PTEN_ORF)),
              by = "nucleotide_position") %>%
    mutate(freq_real = replace_na(freq_real, 0))
  return(real_dist)
}

# Separately for MT-L and MT-H -------------------------------------------------
syn_distr_real_mtl <- get_real_dist(PTEN_variants_syn_mtl)

syn_distr_real_mth <- get_real_dist(PTEN_variants_syn_mth)

# Adjust probability with samples # (1273 MT-L, 238 MT-H) and sabset proportion of overall dataset (0.95 MT-L, 0.05 MT-H)
adjust_probability <-
  function(pred_data,
           real_data,
           n_samples,
           proportion) {
    adjusted_data <-
      full_join(pred_data, real_data, by = "nucleotide_position") %>%
      mutate(freq_pred = (freq_pred / n_samples) * proportion) %>%
      select(nucleotide_position, freq_pred, freq_real)
    return(adjusted_data)
  }

compare_distr_syn_mtl <-
  adjust_probability(syn_distr_predict_mtl, syn_distr_real_mtl, 1273, 0.95) %>%
  rename(freq_pred_mtl = freq_pred, freq_real_mtl = freq_real)

compare_distr_syn_mth <-
  adjust_probability(syn_distr_predict_mth, syn_distr_real_mth, 238, 0.05) %>%
  rename(freq_pred_mth = freq_pred, freq_real_mth = freq_real)

compare_distr_syn <-
  full_join(compare_distr_syn_mtl, compare_distr_syn_mth) %>%
  mutate(
    freq_real_syn = freq_real_mtl + freq_real_mth,
    freq_pred_syn = freq_pred_mtl + freq_pred_mth,
    prop_freq_real_syn = prop.table(freq_real_syn),
    prop_freq_pred_syn = prop.table(freq_pred_syn)
  ) %>%
  select(!contains("mt"))

write_csv(compare_distr_syn,
          "../Charts, data, statistics/Panel E full data.csv")

# Save chart and data ----------------------------------------------------------
panel_E_data <- compare_distr_syn %>%
  mutate(codon = rep(1:404, each = 3)) %>%
  group_by(codon) %>%
  summarise(
    prop_freq_real_syn = sum(prop_freq_real_syn),
    prop_freq_pred_syn = sum(prop_freq_pred_syn)
  ) %>%
  # 0.019 for 1 and more real mutations
  # 0.006 for 5 and more predicted mutations
  filter(prop_freq_real_syn > 0.019 |
           prop_freq_pred_syn > 0.006) %>%
  pivot_longer(
    cols = c(prop_freq_real_syn, prop_freq_pred_syn),
    names_to = "type",
    values_to = "freq"
  ) %>%
  mutate(codon = as.factor(codon))

write_csv(panel_E_data,
          "../Charts, data, statistics/Panel E chart data.csv")

ggplot(panel_E_data, aes(x = codon, y = freq, fill = type)) +
  geom_col(position = position_dodge2()) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90), legend.title = element_blank()) +
  scale_fill_manual(
    values = c(
      "prop_freq_pred_syn" = "#EB5B51",
      "prop_freq_real_syn" = "#6A3CB2"
    ),
    labels = c("Predicted", "FMI"),
    guide = guide_legend(override.aes = list(fill = c(
      "#EB5B51", "#6A3CB2"
    )))
  ) +
  xlab("Codon position") +
  ylab("Frequency proportion")


ggsave(
  "../Charts, data, statistics/Panel E/Panel E.png",
  units = "in",
  dpi = 500,
  width = 6.65,
  height = 5.94
)
