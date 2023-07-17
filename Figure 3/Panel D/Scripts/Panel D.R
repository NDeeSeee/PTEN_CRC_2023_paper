# Attach requirement packages and setting WD -----------------------------------
packages_names <- c("readxl", "data.table", "reshape2", "tidyverse",
                    "rstudioapi")

lapply(packages_names, require, character.only = TRUE)

setwd(dirname(getActiveDocumentContext()$path))

rename = dplyr::rename
select = dplyr::select
filter = dplyr::filter
group_by = dplyr::group_by
mutate = dplyr::mutate

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


# Load PAD, FMI and PTM data ---------------------------------------------------
PAD_samples <-
  fread("../../../Common raw data/PAD samples.csv", sep = ",") %>%
  as_tibble()
PAD_variants <-
  fread("../../../Common raw data/PAD variants.csv", sep = ",") %>%
  as_tibble() %>%
  rename(codon_SV = codon)

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

PTEN_PTM <- fread("../../../Common raw data/PTEN PTM sites.csv") %>%
  as_tibble() %>%
  pivot_longer(
    cols = c("active_site", "PTM"),
    names_to = "impact",
    values_to = "value"
  ) %>%
  filter(value != 0) %>%
  select(-value) %>%
  mutate(impact = as.factor(impact))


# Assign PTEN domain coordinates -----------------------------------------------
PTEN_domains <- tibble(
  domain = c("PDB", "PTPase", "None", "C2", "C-tail"),
  start = c(0, 15, 185, 190, 350),
  end = c(15, 185, 190, 350, 403)
)


# Load PTEN loss of function data ----------------------------------------------
# From Mighell paper
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
  filter(if_any(.cols = c("MAVE_effect", "VAMP_effect"), function(x)
    !is.na(x))) %>%
  mutate(
    consensus_VM = ifelse(MAVE_effect == "lof" |
                            VAMP_effect == "lof", "lof", "wt"),
    consensus_M = ifelse(MAVE_effect == "lof", "lof", "wt"),
    consensus_V = ifelse(VAMP_effect == "lof", "lof", "wt")
  )

# From OncoKB  -----------------------------------------------------------------
protein_effect_oncoKB <-
  read_xlsx("../../../Common raw data/PTEN CKB&oncoKB data.xlsx")[c(1, 2)] %>%
  rename(alt = Alteration, oncokb = Oncogenic) %>%
  mutate(oncokb = as.factor(ifelse(oncokb == "Oncogenic", "lof", NA)))

# From CKB ---------------------------------------------------------------------
protein_effect_CKB <-
  read_xlsx("../../../Common raw data/PTEN CKB&oncoKB data.xlsx", sheet = 2)[c(1, 3)] %>%
  rename(alt = Variant, CKB = `Protein Effect`) %>%
  mutate(CKB = ifelse(CKB == "loss of function", "lof", NA))

# Combine OncoKB and CKB -------------------------------------------------------
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
  filter(if_any(.cols = c("oncokb", "CKB"), function(x)
    !is.na(x))) %>%
  mutate(consensus_KB = ifelse(oncokb == "lof" |
                                 CKB == "lof", "lof", "wt"))

# Define the consensus LoF status based on all the above -----------------------
PTEN_LoF_consensus <- full_join(
  select(mighell_data, pos, alt,
         consensus_V, consensus_M, consensus_VM),
  select(protein_effect_KB, alt, consensus_KB),
  by = "alt"
) %>%
  mutate(consensus_VM_KB = ifelse(is.na(consensus_KB),
                                  consensus_VM,
                                  consensus_KB))



# Function to calculate LoF proportion -----------------------------------------
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
      .fns = function(x)
        replace_na(x, 0)
    )) %>%
    group_by(codon) %>%
    mutate(lof_prop = lof / sum(lof, wt)) %>%
    ungroup() %>%
    select(-lof, -wt)
  
  return(data_lof)
}

# Possible values of arguments
pred_prob_filenames <- c(
  "MTL probs predicted.csv",
  "MTH probs predicted.csv",
  "MTL_MTH combined probs predicted.csv"
)
lof_types <-
  c("consensus_V",
    "consensus_M",
    "consensus_VM",
    "consensus_VM_KB")
msi_values <- c("MT-L", "MT-H")

# Example of usage
calc_lof_prop("MT-L",
              "consensus_VM",
              "../../../Common raw data/MTL probs predicted.csv")

# Function to calculate PTEN hotspots with LoF proportion and related probs ----
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
    
    PTEN_LoF_data <- calc_lof_prop(msi = msi,
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
    
    merged_data_out_top <- slice(merged_data_out, 1:top_hotspots) %>%
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

# Example of usage
get_hotspots_with_probs("MT-L",
                        "consensus_VM",
                        "../../../Common raw data/MTL probs predicted.csv")



# Function to get chart data ---------------------------------------------------
panel_E_chart_data <-
  function(msi,
           lof_type,
           pred_prob_filename,
           top,
           var_type) {
    FMI_PTEN_hotspot <- get_hotspots(
      filter(
        FMI_variants,
        gene == "PTEN",
        Assigned %in% msi,
        variantType2 %in% var_type
      ),
      403,
      "codon"
    ) %>%
      mutate(Freq_prop = Freq / nrow(filter(FMI_variants, gene == "PTEN"))) %>%
      arrange(desc(Freq_prop)) %>%
      mutate(freq_adj = Freq / max(Freq))
    
    
    PTEN_LoF_data <- calc_lof_prop(msi = msi,
                                   lof_type = lof_type,
                                   pred_prob_filename = pred_prob_filename)
    
    
    PRED_data <- fread(pred_prob_filename) %>%
      filter(synonymous == 0) %>%
      rename(probability = rel_prob) %>%
      group_by(codon_position) %>%
      summarise(probability = sum(probability)) %>%
      arrange(desc(probability)) %>%
      rename(codon = codon_position) %>%
      mutate(
        prob_sum_adj = probability / max(probability),
        pred_prob_sum = prop.table(probability)
      )
    
    
    final_df <- FMI_PTEN_hotspot %>%
      left_join(select(PRED_data, codon, prob_sum_adj, pred_prob_sum),
                by = "codon") %>%
      rename(FMI = freq_adj, PRED = prob_sum_adj) %>%
      filter(if_any(.cols = c("FMI", "PRED"), function(x)
        x >= sort(x, decreasing = T)[top])) %>%
      pivot_longer(
        cols = c("FMI", "PRED"),
        names_to = "data",
        values_to = "value"
      ) %>%
      mutate(data = as.factor(data)) %>%
      full_join(PTEN_PTM, by = "codon") %>%
      filter(value != 0 | is.na(value)) %>%
      left_join(PTEN_LoF_data, by = "codon") %>%
      mutate(lof_prop = ifelse(!is.na(impact), NA, lof_prop)) %>%
      mutate(value = ifelse(value > 1, 1, value)) %>%
      group_by(codon) %>%
      mutate(max_value = max(value),
             lof_prop = max_value * lof_prop) %>%
      ungroup() %>%
      # IN GRAPHICS USE ONLY
      mutate(codon_label = as.character(codon),
             codon = ifelse(codon == 233, 232, codon))
    
    return(final_df)
  }

# Possible values of arguments
var_types <- c("amplification", "deletion", "indel", "sub")

# Example of usage
panel_E_chart_data("MT-L",
                   "consensus_VM",
                   "../../../Common raw data/MTL probs predicted.csv",
                   5,
                   "sub")

# Function to plot the chart data ----------------------------------------------
panel_E_plot <- function(data,
                         domain_coeff,
                         draw_lof = FALSE,
                         text_label = FALSE) {
  panel_E_chart <- ggplot() +
    geom_rect(
      data = PTEN_domains,
      aes(
        xmin = start,
        xmax = end,
        ymin = (c(1, 1, 0.5, 1, 1) * -domain_coeff) - domain_coeff,
        ymax = (c(1, 1, 0.5, 1, 1) * domain_coeff) - domain_coeff,
        fill = as.factor(domain)
      ),
      alpha = .9
    ) +
    geom_segment(
      data = data,
      aes(
        x = codon,
        xend = codon,
        y = 0,
        yend = value
      ),
      linewidth = 1.5,
      na.rm = T,
      col = ifelse(draw_lof, "gray", "black")
    ) +
    geom_segment(
      data = data,
      aes(
        x = codon,
        xend = codon,
        y = 0,
        yend = lof_prop
      ),
      linewidth = 1.5,
      na.rm = T,
      col = "black"
    ) +
    geom_point(
      data = data,
      aes(x = codon, y = value, col = data),
      size = 5,
      na.rm = T
    ) +
    geom_point(
      data = filter(data, !is.na(impact)),
      aes(
        x = codon,
        y = -0.075,
        fill = impact,
        col = impact
      ),
      size = 5,
      shape = 24
    ) +
    theme_minimal() +
    theme(legend.position = "top") +
    geom_hline(yintercept = 0, col = "#413F42") +
    scale_fill_manual(
      values = c(
        "active_site" = "#489FF8",
        "PTM" = "NA",
        "PDB" = "#E4ACFA",
        "PTPase" = "#FDEFB2",
        "C2" = "#B9EBFD",
        "C-tail" = "#9BDCAE"
      ),
      labels = c("Active site", "",
                 "PBD", "PTPase", "C2", "C-tail"),
      guide = guide_legend(
        nrow = 2,
        override.aes = list(
          shape = c(17, NA, 0, 0, 0, 0),
          fill = c(NA, NA, "#E4ACFA", "#FDEFB2", "#B9EBFD", "#9BDCAE"),
          linetype = c(0, 0, 0, 0, 0, 0),
          col = c("#489FF8", "red", NA, NA, NA, NA)
        )
      )
    ) +
    scale_colour_manual(
      values = c(
        "FMI" = "#7139B9",
        "PRED" = "#EB5B51",
        "PTM" = "NA",
        "active_site" = "NA"
      ),
      labels = c("FMI", "Predicted", "", ""),
      guide = guide_legend(override.aes = list(
        linetype = c(0, 0, 0, 0),
        shape = c(16, 16, 0, 0),
        color = c("#7139B9", "#EB5B51", NA, NA)
      ))
    ) +
    theme(
      legend.title = element_blank(),
      axis.text = element_text(
        size = 20,
        face = "bold",
        colour = "black"
      ),
      legend.text = element_text(size = 20),
      axis.line.y = element_line(colour = "#413F42"),
      axis.title.y = element_text(
        size = 25,
        colour = "black",
        vjust = +2
      ),
      axis.title.x = element_text(
        size = 25,
        colour = "black",
        vjust = -0.15
      )
    ) +
    scale_y_continuous(
      limits = c(-0.3, 1.3),
      expand = c(0, 0),
      breaks = c(0, 1)
    ) +
    scale_x_continuous(
      limits = c(0, 407),
      expand = c(0, 0),
      breaks = c(0, 100, 200, 300, 403)
    ) +
    ylab("Relative abundance") +
    xlab("Codon position")
  
  if (text_label) {
    panel_E_chart <- panel_E_chart +
      geom_text(
        data = data,
        mapping = aes(x = codon, y = max_value,
                      label = codon_label),
        nudge_y = .15,
        size = 4,
        fontface = "bold",
        angle = 0,
        check_overlap = T
      )
  }
  
  return(panel_E_chart)
}

# Draw the chart Panel D -------------------------------------------------------
panel_E_chart_data(
  "MT-L",
  "consensus_M",
  "../../../Common raw data/MTL probs predicted.csv",
  5,
  c("sub", "indel")
) %>%
  panel_E_plot(domain_coeff = 0.15, draw_lof = T)


# Save chart data --------------------------------------------------------------
panel_E_chart_data(
  "MT-L",
  "consensus_M",
  "../../../Common raw data/MTL probs predicted.csv",
  5,
  c("sub", "indel")
) %>%
  write.csv("../Charts, data, statistics/Panel D chart data.csv", row.names = F)

# Save image -------------------------------------------------------------------
ggsave(
  "../Charts, data, statistics/Panel D.png",
  dpi = 600,
  width = 18.18,
  height = 4,
  units = "in"
)
