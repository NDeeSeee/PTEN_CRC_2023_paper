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

# PTEN ORF ---------------------------------------------------------------------
PTEN_ORF <-
  DNAString(str_remove_all(read_file("../../../Common raw data/PTEN ORF.txt"), "\n"))

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


# Function to convert vector of nucleotide coordinates -------------------------
# to codon numbers vector 
# E.g. 4, 12, 18 -> 2, 4, 6
# nucl2codon(c(4,12,18))
nucl2codon <- function(nd_coord) {
  cdn_coord <- c()
  for (i in nd_coord) {
    if (i %% 3 == 0) {
      cdn_coord <- c(cdn_coord, i %/% 3)
    } else {
      cdn_coord <- c(cdn_coord, i %/% 3 + 1)
    }
  }
  return(cdn_coord)
}


# Function to replace reference nucleotide to altered (for SNP only) -----------
# at reference codon, return altered codon
# For each element of Vector
# E.g. "ATG" altered at 3 position G > T -> function returns "ATT", etc.
  # get_alt_codon(c("ATG", "TTT", "ACT"), c(2,1,1), c("G", "A", "T"))
get_alt_codon <- function(ref_codons, alt_pos, alt_nn) {
  alt_codons <- c()
  for (i in 1:length(ref_codons)) {
    alt_codons <-
      c(alt_codons, toString(replaceLetterAt(
        DNAString(ref_codons[i]),
        at = alt_pos[i],
        letter = alt_nn[i]
      )))
  }
  return(alt_codons)
}


# Main function to plot the difference -----------------------------------------
# between codon frequency (alt and reference) chart and save .csv format data
# Input parameters
# data - dataframe with transcriptEffect column
# ORF - open reading frame as string
# codon_usage_file - file with codons and their frequency
# output_prefix - for optional saving .csv and .png
codon_usage_difference <-
  function(data,
           ORF,
           codon_usage_file,
           output_prefix) {
    # Get all codons from ORF with position number
    # Like "ATGTTC" -> "ATG":1 "TTC":2 etc.
    extracted_codons <-
      extractAt(ORF, IRanges(start = seq(1, (length(
        ORF
      ) - 2), 3), width = 3)) %>%
      as.data.frame() %>%
      as_tibble() %>%
      mutate(cdn_coord = 1:(length(ORF) / 3)) %>%
      rename(ref_codon = x)


    # Main step of this function
    alt_codons_data <- data %>%
      select(transcriptEffect_SV, Assigned) %>%
      # 1. Extract SNP coord and alt nucleotide from transcriptEffect_SV column
      # "1342G>T" -> nd_coord = 1342, alt = "T"
      mutate(
        nd_coord = as.numeric(str_extract(transcriptEffect_SV, "[:digit:]+")),
        transcriptEffect_SV = str_remove_all(transcriptEffect_SV, "[:digit:]+"),
        alt = str_sub(transcriptEffect_SV, 3, 3)
      ) %>%
      select(-transcriptEffect_SV) %>%
      arrange(nd_coord) %>%
      # 2. Associate codon numbers with SNP positions using nucl2codon function
      mutate(cdn_coord = nucl2codon(nd_coord)) %>%
      left_join(y = extracted_codons, by = "cdn_coord") %>%
      # 3. Finding substitution position at codon (out of 3) based on nucleotide position
      # E.g. 3, 4, 5 -> 3, 1, 2 etc.
      mutate(
        cdn_alt_pos = nd_coord %% 3,
        cdn_alt_pos = ifelse(cdn_alt_pos == 0, 3, cdn_alt_pos),
        # 4. Get altered codon using get_alt_codon function
        alt_codon = get_alt_codon(ref_codon, cdn_alt_pos, alt),
        # 5. Translate reference and altered codons to aminoacids
        ref_aa = as.vector(translate(DNAStringSet(ref_codon))),
        alt_aa = as.vector(translate(DNAStringSet(alt_codon)))
      ) %>%
      as_tibble()


    # Read the codon frequency table
    # Then convert RNA to DNA alphabet (U -> T)
    codon_usage_DNA <- fread(codon_usage_file) %>%
      mutate(rna_triplet = str_replace_all(rna_triplet, "U", "T")) %>%
      rename(dna_triplet = rna_triplet) %>%
      as_tibble()


    # Combine alt codons data with codon usage table
    # Then calculate difference between frequencies
    alt_codons_data_freq <- alt_codons_data %>%
      left_join(
        y = rename(
          codon_usage_DNA,
          ref_codon = dna_triplet,
          ref_cdn_freq = frequency
        ),
        by = "ref_codon"
      ) %>%
      left_join(
        y = rename(
          codon_usage_DNA,
          alt_codon = dna_triplet,
          alt_cdn_freq = frequency
        ),
        by = "alt_codon"
      ) %>%
      mutate(cdn_freq_diff = ref_cdn_freq - alt_cdn_freq)


    # Save result data at .csv format
    alt_codons_data_freq %>%
      select(
        nd_coord,
        cdn_coord,
        alt,
        cdn_alt_pos,
        ref_codon,
        alt_codon,
        ref_aa,
        alt_aa,
        ref_cdn_freq,
        alt_cdn_freq,
        cdn_freq_diff,
        Assigned
      ) %>%
      write_csv(file = paste0(output_prefix, " data.csv"))


    # Unique values of MS status
    MS_status <- names(table(alt_codons_data_freq$Assigned))


    # One-sample t-test to compare with 0
    MS_one_sample_test <- tibble()
    for (i in MS_status) {
      # Limitations that t.test requires
      if (nrow(filter(alt_codons_data_freq, Assigned == i)) > 10) {
        one_sample_test <-
          t.test(
            filter(alt_codons_data_freq, Assigned == i)$cdn_freq_diff,
            mu = 0,
            alternative = "two.sided"
          )
        one_sample_int_df <- tibble(
          MS_status = i,
          p_val = one_sample_test$p.value,
          lower_conf = one_sample_test$conf.int[1],
          upper_conf = one_sample_test$conf.int[2],
          mean_diff = one_sample_test$estimate,
          n = nrow(filter(alt_codons_data_freq, Assigned == i))
        )
        MS_one_sample_test <-
          bind_rows(MS_one_sample_test, one_sample_int_df)
      }
    }

    MS_one_sample_test %>%
      write_csv(file = paste0(output_prefix, " one sample test.csv"))


    MS_status_cmb <- combn(MS_status, 2)

    # Two-sample t-test to compare with each other
    MS_two_sample_test <- tibble()
    for (i in 1:ncol(MS_status_cmb)) {
      first_group <- MS_status_cmb[, i][1]
      second_group <- MS_status_cmb[, i][2]

      # Limitations that t.test requires
      if (nrow(filter(alt_codons_data_freq, Assigned == first_group)) > 10 &
        nrow(filter(alt_codons_data_freq, Assigned == second_group)) > 10) {
        two_sample_test <-
          t.test(
            filter(alt_codons_data_freq, Assigned == first_group)$cdn_freq_diff,
            filter(alt_codons_data_freq, Assigned == second_group)$cdn_freq_diff,
            alternative = "two.sided"
          )
        two_sample_int_df <- tibble(
          MS_status = paste0(first_group, " | ", second_group),
          p_val = two_sample_test$p.value,
          n = nrow(filter(
            alt_codons_data_freq, Assigned == first_group
          )) +
            nrow(filter(
              alt_codons_data_freq, Assigned == second_group
            ))
        )
        MS_two_sample_test <-
          bind_rows(MS_two_sample_test, two_sample_int_df)
      }
    }

    MS_two_sample_test <-
      MS_two_sample_test[!duplicated.data.frame(MS_two_sample_test), ]

    MS_two_sample_test %>%
      write_csv(file = paste0(output_prefix, " two sample test.csv"))


    # Make violin chart data
    alt_codons_data_freq_plot <- alt_codons_data_freq %>%
      filter(!is.na(Assigned)) %>%
      ggplot(aes(y = cdn_freq_diff, x = Assigned, fill = Assigned)) +
      geom_violin(alpha = .8, show.legend = F) +
      geom_boxplot(width = .1) +
      theme_minimal() +
      theme(legend.title = element_blank(), legend.key = element_blank()) +
      scale_fill_manual(
        values = c(
          "MT-H" = "#E4ACFA",
          "MT-L" = "#489FF8",
          "MT-L htmb" = "#FFBF15"
        ),
        guide = guide_legend(override.aes = list(fill = c(
          "#E4ACFA", "#489FF8", "#FFBF15"
        )))
      ) +
      xlab("") +
      ylab("Codon alteration frequency difference")

    # Save violin chart data
    ggsave(
      filename = paste0(output_prefix, ".png"),
      units = "in",
      width = 6.65,
      height = 5.94
    )

    return(alt_codons_data_freq_plot)
  }


# Synonymous mutations ---------------------------------------------------------
PTEN_variants_s <- filter(
  FMI_variants,
  gene == "PTEN",
  variantType2 == "sub",
  codingType_SV != "splice",
  codingType_SV == "synonymous"
)

# Codon usage table data from here https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606&aa=1&style=N
codon_usage_difference(
  data = PTEN_variants_s,
  ORF = PTEN_ORF,
  codon_usage_file = "../../../Common raw data/Codon frequency.csv",
  output_prefix = "../Charts, data, statistics/Panel D"
)
