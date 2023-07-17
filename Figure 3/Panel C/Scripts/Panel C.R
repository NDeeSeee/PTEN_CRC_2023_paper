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


# Load PAD, FMI and PTM data ---------------------------------------------------
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
# LoF data reading -------------------------------------------------------------
# MAVE and VAMP data
mighell_data <-
  read_xlsx("../../../Common raw data/Mighell Eng integrated LPA abundance.xlsx") %>%
  rename(
    alt = var,
    VAMP_score = abundance_score,
    VAMP_effect = abundance_score_class,
    MAVE_score = Fitness_score,
    MAVE_effect = Fitness_score_class
  ) %>%
  mutate(
    VAMP_effect = str_replace_all(VAMP_effect, "truncation-like", "lof"),
    VAMP_effect = str_replace_all(VAMP_effect, "hypomorphic", "lof"),
    VAMP_effect = str_replace_all(VAMP_effect, "wildtype-like", "wt"),
    MAVE_effect = str_replace_all(MAVE_effect, "truncation-like", "lof"),
    MAVE_effect = str_replace_all(MAVE_effect, "hypomorphic", "lof"),
    MAVE_effect = str_replace_all(MAVE_effect, "wildtype-like", "wt")
  ) %>%
  filter(if_any(.cols = c("MAVE_effect", "VAMP_effect"), function(x) {
    !is.na(x)
  })) %>%
  mutate(
    consensus_VM = ifelse(MAVE_effect == "lof" |
                            VAMP_effect == "lof", "lof", "wt"),
    consensus_M = ifelse(MAVE_effect == "lof", "lof", "wt"),
    consensus_V = ifelse(VAMP_effect == "lof", "lof", "wt")
  )

# Data from oncoKB
protein_effect_oncoKB <-
  read_xlsx("../../../Common raw data/PTEN CKB&oncoKB data.xlsx")[c(1, 2)] %>%
  rename(alt = Alteration, oncokb = Oncogenic) %>%
  mutate(oncokb = as.factor(ifelse(oncokb == "Oncogenic", "lof", NA)))

# Data from CKB
protein_effect_CKB <-
  read_xlsx("../../../Common raw data/PTEN CKB&oncoKB data.xlsx", sheet = 2)[c(1, 3)] %>%
  rename(alt = Variant, CKB = `Protein Effect`) %>%
  mutate(CKB = ifelse(CKB == "loss of function", "lof", NA))

# Combined oncoKB and CKB
protein_effect_KB <-
  full_join(protein_effect_oncoKB, protein_effect_CKB, "alt") %>%
  filter(
    !str_detect(alt, "fs"),
    !str_detect(alt, "dec"),
    !str_detect(alt, "wild-type"),
    !str_detect(alt, "\\*"),
    !str_detect(alt, "del"),
    !str_detect(alt, "mut"),
    !str_detect(alt, "negative"),
    !str_detect(alt, "Truncating"),
    !str_detect(alt, "loss")
  ) %>%
  filter(if_any(.cols = c("oncokb", "CKB"), function(x) {
    !is.na(x)
  })) %>%
  mutate(consensus_KB = ifelse(oncokb == "lof" |
                                 CKB == "lof", "lof", "wt"))

# Combine them all
PTEN_LoF_consensus <- full_join(
  select(mighell_data, pos, alt,
         consensus_V, consensus_M, consensus_VM),
  select(protein_effect_KB, alt, consensus_KB),
  by = "alt"
) %>%
  mutate(consensus_VM_KB = ifelse(is.na(consensus_KB),
                                  consensus_VM,
                                  consensus_KB))

# Possible values of arguments
pred_prob_filenames <-
  c("MTL probs predicted.csv",
    "MTH_predicted.csv",
    "MTL_MTH_combined_predicted.csv")
lof_types <-
  c("consensus_V",
    "consensus_M",
    "consensus_VM",
    "consensus_VM_KB")
msi_values <- c("MT-L", "MT-H")

# Main function to calculate hotspots ------------------------------------------
get_hotspots <- function(data,
                         a_length,
                         output_type = "codon",
                         p_cutoff = 0.005) {
  alt_freq <- as_tibble(as.data.frame(table(data$alt))) %>%
    rename(alt = Var1) %>%
    mutate(codon = as.numeric(str_extract(alt, "[:digit:]+")),
           alt = as.character(alt)) %>%
    filter(!is.na(codon), !str_detect(alt, "splice"), !is.na(alt)) %>%
    right_join(y = data.frame(codon = 1:a_length), by = "codon") %>%
    mutate(Freq = replace_na(Freq, 0))
  codon_freq <- as_tibble(as.data.frame(table(data$codon_SV))) %>%
    rename(codon = Var1) %>%
    mutate(codon = as.numeric(as.character(codon))) %>%
    right_join(y = data.frame(codon = 1:a_length), by = "codon") %>%
    mutate(Freq = replace_na(Freq, 0))
  htsptable <- left_join(
    x = data.frame(Var1 = 1:28),
    y = mutate(as.data.frame(table(codon_freq$Freq)),
               Var1 = as.numeric(as.character(Var1))),
    by = "Var1"
  ) %>%
    rename(number_of_mut = Var1) %>%
    mutate(
      Freq = replace_na(Freq, 0),
      cod_count = 0,
      cod_sum = 0,
      pval1 = NA
    )
  for (i in 1:28) {
    htsptable$cod_count[i] <-
      sum(codon_freq$Freq <= htsptable$number_of_mut[i])
    htsptable$cod_sum[i] <-
      sum(codon_freq$Freq[which(codon_freq$Freq <= htsptable$number_of_mut[i])])
    htsptable$pval1[i] <- sum(
      dbinom(
        x = htsptable$cod_sum[i] - htsptable$number_of_mut[i]:htsptable$cod_sum[i],
        size = htsptable$cod_sum[i] - 1,
        prob = 1 - 1 / htsptable$cod_count[i]
      )
    )
  }
  cutoff <-
    htsptable$number_of_mut[which(htsptable$pval1 < p_cutoff)][1]
  codon_freq$cutoff <- cutoff
  codon_freq$no_hotspots <-
    ifelse(codon_freq$Freq >= codon_freq$cutoff[1], 0, codon_freq$Freq)
  codon_freq$no_hotspots_2 <-
    ifelse(codon_freq$Freq >= codon_freq$cutoff[1], 0, 1)
  codon_freq <-
    mutate(codon_freq, sum_of_mut = NA, count_of_mut = NA)
  for (i in 1:(nrow(codon_freq) - 4)) {
    codon_freq$sum_of_mut[i] <- sum(codon_freq$no_hotspots[i:(i + 4)])
    codon_freq$count_of_mut[i] <-
      sum(codon_freq$no_hotspots_2[i:(i + 4)])
  }
  number_of_mut_in_those <-
    filter(htsptable, number_of_mut == codon_freq$cutoff[1])$cod_sum[1]
  codons_with_less_than_that <-
    filter(htsptable, number_of_mut == codon_freq$cutoff[1])$cod_count[1]
  codon_freq$p_value <- NA
  for (i in 1:a_length) {
    if (codon_freq$Freq[i] <= 28 & codon_freq$Freq[i] > 0) {
      codon_freq$p_value[i] <-
        filter(htsptable, number_of_mut == codon_freq$Freq[i])$pval1
    } else if (codon_freq$Freq[i] != 0) {
      codon_freq$p_value[i] <-
        filter(htsptable, number_of_mut == 28)$pval1
    } else {
      codon_freq$p_value[i] <- 1
    }
  }
  htsptable <- htsptable %>%
    mutate(
      pval2 = NA,
      pval3 = NA,
      pval4 = NA,
      pval5 = NA
    )
  for (i in 1:nrow(htsptable)) {
    htsptable$pval2[i] <-
      sum(
        dbinom(
          x = 1:(number_of_mut_in_those - htsptable$number_of_mut[i]),
          size = number_of_mut_in_those -
            1,
          prob = 1 - 2 / codons_with_less_than_that
        )
      )
    htsptable$pval3[i] <-
      sum(
        dbinom(
          x = 1:(number_of_mut_in_those - htsptable$number_of_mut[i]),
          size = number_of_mut_in_those -
            1,
          prob = 1 - 3 / codons_with_less_than_that
        )
      )
    htsptable$pval4[i] <-
      sum(
        dbinom(
          x = 1:(number_of_mut_in_those - htsptable$number_of_mut[i]),
          size = number_of_mut_in_those -
            1,
          prob = 1 - 4 / codons_with_less_than_that
        )
      )
    htsptable$pval5[i] <-
      sum(
        dbinom(
          x = 1:(number_of_mut_in_those - htsptable$number_of_mut[i]),
          size = number_of_mut_in_those -
            1,
          prob = 1 - 5 / codons_with_less_than_that
        )
      )
  }
  alt_freq$cutoff <- codon_freq$cutoff[1]
  codon_freq <- arrange(codon_freq, codon)
  if (output_type == "codon") {
    return(select(codon_freq, codon, Freq, cutoff, p_value))
  } else {
    if (output_type == "alt") {
      return(alt_freq)
    }
  }
}
# Function to calculate LoF proporion ------------------------------------------
calc_lof_prop <- function(msi, lof_type, pred_prob_filename) {
  data <- fread(pred_prob_filename)
  
  lof_data_consensus <- PTEN_LoF_consensus %>%
    mutate(consensus_effect = as.vector(select(PTEN_LoF_consensus, lof_type))[[1]])
  
  data_lof <- data %>%
    rename(alt = amino_acid_change, probability = rel_prob) %>%
    group_by(alt) %>%
    summarise(
      synonymous = max(synonymous),
      rel_prob_sum = sum(probability, na.rm = T)
    ) %>%
    left_join(select(lof_data_consensus, alt, consensus_effect), by = "alt") %>%
    filter(synonymous != 1) %>%
    mutate(consensus_effect = ifelse(str_detect(alt, "\\*"), "lof", consensus_effect)) %>%
    select(-synonymous) %>%
    mutate(codon = as.numeric(str_extract(alt, "[:digit:]+"))) %>%
    mutate(consensus_effect = as.factor(consensus_effect)) %>%
    group_by(codon, consensus_effect) %>%
    summarise(prob_sum = sum(rel_prob_sum, na.rm = T)) %>%
    ungroup() %>%
    filter(!is.na(consensus_effect)) %>%
    pivot_wider(values_from = prob_sum, names_from = consensus_effect) %>%
    mutate(across(
      .cols = c("lof", "wt"),
      .fns = function(x) {
        replace_na(x, 0)
      }
    )) %>%
    group_by(codon) %>%
    mutate(lof_prop = lof / sum(lof, wt)) %>%
    ungroup() %>%
    select(-lof, -wt)
  
  return(data_lof)
}

# Example
calc_lof_prop("MT-L",
              "consensus_VM",
              "../../../Common raw data/MTL probs predicted.csv")

# Function to calculate PTEN hotspots with LoF proporion and related probs -----
get_hotspots_with_probs <-
  function(msi,
           lof_type,
           pred_prob_filename,
           top_hotspots = 404) {
    PRED_data <- fread(pred_prob_filename) %>%
      filter(synonymous != 1) %>%
      group_by(codon_position) %>%
      summarise(rel_prob_sum = sum(rel_prob)) %>%
      mutate(rel_prob_sum_prop = prop.table(rel_prob_sum)) %>%
      rename(codon = codon_position, pred_prop = rel_prob_sum_prop) %>%
      select(-rel_prob_sum)
    
    PTEN_FMI_hotspots <-
      get_hotspots(
        filter(
          FMI_variants,
          gene == "PTEN",
          Assigned %in% msi,
          variantType2 == "sub"
        ),
        403,
        "codon"
      ) %>%
      mutate(Freq_prop = Freq / nrow(filter(FMI_variants, gene == "PTEN"))) %>%
      arrange(desc(Freq_prop)) %>%
      mutate(freq_adj = Freq / max(Freq),
             sum_freq = sum(Freq))
    
    PTEN_LoF_data <-
      calc_lof_prop(msi = msi,
                    lof_type = lof_type,
                    pred_prob_filename = pred_prob_filename)
    
    merged_data_out <-
      full_join(PTEN_FMI_hotspots, PRED_data, by = "codon") %>%
      select(codon, Freq, cutoff, sum_freq, pred_prop) %>%
      left_join(select(PTEN_LoF_data, codon, lof_prop), by = "codon") %>%
      mutate(
        pred_prop_lof = pred_prop * lof_prop,
        pred_prop_lof_prop = prop.table(pred_prop_lof),
        Freq = replace_na(Freq, 0),
        cutoff = replace_na(max(cutoff, na.rm = T)),
        sum_freq = replace_na(max(sum_freq, na.rm = T))
      )
    
    merged_data_out_top <-
      slice(merged_data_out, 1:top_hotspots) %>%
      mutate(
        pred_prop = prop.table(pred_prop),
        pred_prop_lof = pred_prop * lof_prop,
        pred_prop_lof_prop = prop.table(pred_prop_lof),
        sum_freq = sum(Freq)
      )
    
    if (top_hotspots != 404) {
      return(merged_data_out_top)
    } else {
      return(merged_data_out)
    }
  }

# Example
get_hotspots_with_probs("MT-L",
                        "consensus_VM",
                        "../../../Common raw data/MTL probs predicted.csv")

get_hotspots_with_probs("MT-L",
                        "consensus_VM",
                        "../../../Common raw data/MTL probs predicted.csv",
                        100)

# Function to calculate chisq2 for plain and LoF modified predictions ----------
get_chsq2 <- function(data) {
  codon_count <- nrow(data)
  
  prob_plain <- data$pred_prop
  prob_lof <- data$pred_prop_lof_prop
  
  freq <- data$Freq
  muts <- sum(freq)
  
  freq_limit <- 105
  
  hist1 <- rep(0, freq_limit + 1)
  for (i in 1:codon_count) {
    j <- freq[i]
    hist1[j + 1] <- hist1[j + 1] + 1
  }
  
  xpd_plain <- rep(0, freq_limit)
  for (j in 0:freq_limit) {
    xpd_plain[j + 1] <- sum(dbinom(j, muts, prob_plain))
  }
  
  xpd_lof <- rep(0, freq_limit)
  for (j in 0:freq_limit) {
    xpd_lof[j + 1] <- sum(dbinom(j, muts, prob_lof))
  }
  
  chisq2_plain <- 0
  for (j in 0:10) {
    chisq2_plain <-
      chisq2_plain + ((hist1[j + 1] - xpd_plain[j + 1]) ^ 2) / xpd_plain[j + 1]
  }
  
  chisq2_lof <- 0
  for (j in 0:10) {
    chisq2_lof <-
      chisq2_lof + ((hist1[j + 1] - xpd_lof[j + 1]) ^ 2) / xpd_lof[j + 1]
  }
  
  return(c(chisq2_plain, chisq2_lof))
}

# Example
get_hotspots_with_probs("MT-L",
                        "consensus_VM",
                        "../../../Common raw data/MTL probs predicted.csv") %>%
  get_chsq2()

# Start counting chisq2 --------------------------------------------------------
lof_type_df <- data.frame()
for (temp_lof_type in lof_types) {
  chisq2_values <-
    get_hotspots_with_probs("MT-L",
                            temp_lof_type,
                            "../../../Common raw data/MTL probs predicted.csv") %>%
    get_chsq2()
  lof_type_df <-
    bind_rows(
      lof_type_df,
      data.frame(
        mut_sign_chisq2 = chisq2_values[1],
        lof_chisq2 = chisq2_values[2],
        lof_type = temp_lof_type
      )
    )
}
lof_type_df <- as.tibble(lof_type_df)
print(lof_type_df)
# End
