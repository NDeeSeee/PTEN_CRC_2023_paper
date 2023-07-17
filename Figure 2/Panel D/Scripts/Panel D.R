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
  filter(FMI_variants, gene == "PTEN")


# now it does, and germline mutations are not there, either

mut_sig_spectra <-
  read_csv("../../../Common raw data/mut_sig_spectra.csv")
# This file was extracted from Sheet SD5 of 41467_2021_22125_MOESM4_ESM.xls,
# Supplementary Data 5. Mutational signature spectra.
# these are supplem data for Nat Comm 2021 paper
# The origins and genetic interactions of KRAS mutations are allele- and tissue-specific
# https://www.nature.com/articles/s41467-021-22125-z#Abs1

# This file was extracted from Sheet SD4 of 41467_2021_22125_MOESM4_ESM.xls,
# Supplementary Data 4. Mutational signature compositions for each tumor sample.
mutsig_by_samples_orig <-
  read_excel("41467_2021_22125_MOESM4_ESM.xlsx",
             sheet = "SD4. mutational signature level",
             skip = 1)
# skipping first row which has the table name
mutsig_by_samples <- mutsig_by_samples_orig %>%
  filter(cancer == "COAD") %>%
  select(-(1:2)) %>%
  mutate(across(starts_with("sig"), ~ as.numeric(.x))) #

mutsig_by_samples <- mutsig_by_samples %>%
  select(tumor_sample_barcode, names(which(colSums(mutsig_by_samples[-1]) > 0)))
# removed all non-COAD samples, and all signatures which do not contribute to COAD mutations

# now sorting trying to sort into MSI and POLE and none of these
# note that paper classifies MSI based on sig_6 alone?
# signature names are taken from COSMIC
# https://cancer.sanger.ac.uk/signatures/sbs/sbs6/
# other MSI associated signatures: SBS14, SBS15, SBS20, SBS21, SBS26, and SBS44.
MSIsig <-
  c("sig_6", "sig_14", "sig_15", "sig_20", "sig_21", "sig_26")
POLEsig <- c("sig_10a", "sig_10b", "sig_10c")
# assuming that the signatures are defined similarly to COSMIC
mutsig_by_samples <- mutsig_by_samples %>%
  mutate(msi_sig = rowSums(across(MSIsig)), pole_sig = rowSums(across(POLEsig))) %>%
  arrange(-msi_sig, -pole_sig)
mutsig_by_samples <- mutsig_by_samples %>%
  mutate(Assigned = case_when(msi_sig > 0.5 ~ "MT-H",
                              pole_sig > 0.5 ~ "MT-L htmb",
                              TRUE ~ "MT-L"))
# visual inspection indicated that there are no samples
# with significant contributions of several signature types
mutsig_by_samples %>% write_csv("mutsig_by_samples.csv")

mutsig_by_samples <- read_csv("mutsig_by_samples.csv")

mutsig_by_samples_MTL_summary <- mutsig_by_samples %>%
  filter(Assigned == "MT-L")
mutsig_by_samples_MTL_summary <- mutsig_by_samples_MTL_summary %>%
  select(names(which(
    colSums(mutsig_by_samples_MTL_summary[-c(1, length(mutsig_by_samples_MTL_summary))]) > 0
  )))
mutsig_by_samples_MTL_summary <- mutsig_by_samples_MTL_summary %>%
  select(starts_with("sig")) %>%
  colSums()

mutsig_by_samples_MTH_summary <- mutsig_by_samples %>%
  filter(Assigned == "MT-H")
mutsig_by_samples_MTH_summary <- mutsig_by_samples_MTH_summary %>%
  select(names(which(
    colSums(mutsig_by_samples_MTH_summary[-c(1, length(mutsig_by_samples_MTH_summary))]) > 0
  )))
mutsig_by_samples_MTH_summary <- mutsig_by_samples_MTH_summary %>%
  select(starts_with("sig")) %>%
  colSums()

# mut_sig_spectra %>% select(names(mutsig_by_samples_MTL_summary) )
# select(mut_sig_spectra,names(mutsig_by_samples_MTL_summary) ) %>% length() #13
# as.vector(mutsig_by_samples_MTL_summary)
# as.vector(mutsig_by_samples_MTL_summary) %>% length() # 13
# just tested (above, commented out) that the dimensions are the same

tricontext_MTL_prob <- bind_cols(mut_sig_spectra[1],
                                 rel_prob = rowSums(mapply(
                                   "*",
                                   select(mut_sig_spectra, names(mutsig_by_samples_MTL_summary)),
                                   as.vector(mutsig_by_samples_MTL_summary)
                                 )))
tricontext_MTL_prob %>% write_csv("tricontext_MTL_prob.csv")
tricontext_MTL_prob <- read_csv("tricontext_MTL_prob.csv")
tricontext_MTL_pred <- tricontext_MTL_prob %>%
  mutate(tnt = paste0(
    str_sub(tricontext, 1, 1),
    str_sub(tricontext, 3, 3),
    str_sub(tricontext, 7, 7)
  ),
  nt_ch = str_sub(tricontext, 3, 5))

tricontext_MTH_prob <- bind_cols(mut_sig_spectra[1],
                                 rel_prob = rowSums(mapply(
                                   "*",
                                   select(mut_sig_spectra, names(mutsig_by_samples_MTH_summary)),
                                   as.vector(mutsig_by_samples_MTH_summary)
                                 )))
tricontext_MTH_prob %>% write_csv("tricontext_MTH_prob.csv")
tricontext_MTH_prob <- read_csv("tricontext_MTH_prob.csv")
tricontext_MTH_pred <- tricontext_MTH_prob %>%
  mutate(tnt = paste0(
    str_sub(tricontext, 1, 1),
    str_sub(tricontext, 3, 3),
    str_sub(tricontext, 7, 7)
  ),
  nt_ch = str_sub(tricontext, 3, 5))

PTEN_ORF <-
  "ATGACAGCCATCATCAAAGAGATCGTTAGCAGAAACAAAAGGAGATATCAAGAGGATGGATTCGA
CTTAGACTTGACCTATATTTATCCAAACATTATTGCTATGGGATTTCCTGCAGAAAGACTTGAAGGCGTA
TACAGGAACAATATTGATGATGTAGTAAGGTTTTTGGATTCAAAGCATAAAAACCATTACAAGATATACA
ATCTTTGTGCTGAAAGACATTATGACACCGCCAAATTTAATTGCAGAGTTGCACAATATCCTTTTGAAGA
CCATAACCCACCACAGCTAGAACTTATCAAACCCTTTTGTGAAGATCTTGACCAATGGCTAAGTGAAGAT
GACAATCATGTTGCAGCAATTCACTGTAAAGCTGGAAAGGGACGAACTGGTGTAATGATATGTGCATATT
TATTACATCGGGGCAAATTTTTAAAGGCACAAGAGGCCCTAGATTTCTATGGGGAAGTAAGGACCAGAGA
CAAAAAGGGAGTAACTATTCCCAGTCAGAGGCGCTATGTGTATTATTATAGCTACCTGTTAAAGAATCAT
CTGGATTATAGACCAGTGGCACTGTTGTTTCACAAGATGATGTTTGAAACTATTCCAATGTTCAGTGGCG
GAACTTGCAATCCTCAGTTTGTGGTCTGCCAGCTAAAGGTGAAGATATATTCCTCCAATTCAGGACCCAC
ACGACGGGAAGACAAGTTCATGTACTTTGAGTTCCCTCAGCCGTTACCTGTGTGTGGTGATATCAAAGTA
GAGTTCTTCCACAAACAGAACAAGATGCTAAAAAAGGACAAAATGTTTCACTTTTGGGTAAATACATTCT
TCATACCAGGACCAGAGGAAACCTCAGAAAAAGTAGAAAATGGAAGTCTATGTGATCAAGAAATCGATAG
CATTTGCAGTATAGAGCGTGCAGATAATGACAAGGAATATCTAGTACTTACTTTAACAAAAAATGATCTT
GACAAAGCAAATAAAGACAAAGCCAACCGATACTTTTCTCCAAATTTTAAGGTGAAGCTGTACTTCACAA
AAACAGTAGAGGAGCCGTCAAATCCAGAGGCTAGCAGTTCAACTTCTGTAACACCAGATGTTAGTGACAA
TGAACCTGATCATTATAGATATTCTGACACCACTGACTCTGATCCAGAGAATGAACCTTTTGATGAAGAT
CAGCATACACAAATTACAAAAGTCTGA"

PTEN_ORF <- str_replace_all(PTEN_ORF, "\\s", "")
PTEN_ORF_comp <- PTEN_ORF %>%
  Biostrings::DNAString() %>%
  Biostrings::reverseComplement() %>%
  Biostrings::toString()
PTEN_ORF_bothways <- paste0(PTEN_ORF, PTEN_ORF_comp)

tricontext_MTL_pred <- tricontext_MTL_pred %>%
  mutate(tnt_count = str_count(PTEN_ORF_bothways, tnt),
         nt_ch_pred = rel_prob * tnt_count)
tricontext_MTL_pred %>% write_csv("tricontext_MTL_pred.csv")

tricontext_MTH_pred <- tricontext_MTH_pred %>%
  mutate(tnt_count = str_count(PTEN_ORF_bothways, tnt),
         nt_ch_pred = rel_prob * tnt_count)
tricontext_MTH_pred %>% write_csv("tricontext_MTH_pred.csv")

# note: here just reading it back
tricontext_MTL_pred <- read_csv("tricontext_MTL_pred.csv")

variants_2019_pten_by_sig <- variants_2019_pten %>%
  filter(!is.na(codon_SV),
         trinucleotideContext_SV != "indel",
         trinucleotideContext_SV != "-") %>%
  select(
    did,
    Assigned,
    codon_SV,
    proteinEffect_SV,
    trinucleotideContext_SV,
    trinucleotideAlteration_SV
  )

# creating a column which would match the mut_sig_spectra
variants_2019_pten_by_sig <- variants_2019_pten_by_sig %>%
  mutate(tricontext = paste0(
    str_sub(trinucleotideContext_SV, 1, 1),
    "[",
    trinucleotideAlteration_SV,
    "]",
    str_sub(trinucleotideContext_SV, -1, -1)
  ))
variants_2019_pten_by_sig %>% write_csv("variants_2019_pten_by_sig.csv")

variants_2019_pten_by_sig %>%
  group_by(tricontext) %>%
  summarise(mut_count = n(), codon_count = n_distinct(codon_SV))

variants_2019_pten_by_sig_with_prob <- variants_2019_pten_by_sig %>%
  filter(Assigned == "MT-L") %>%
  group_by(codon_SV, tricontext) %>%
  summarise(mut_count = n()) %>%
  left_join(tricontext_MTL_prob)

variants_2019_pten_by_sig_with_prob %>% write_csv("variants_2019_pten_by_sig_with_prob.csv")
tricontext_MTL_pred

variants_2019_pten_MTL_nt_change_counts <-
  variants_2019_pten_by_sig %>%
  droplevels() %>%
  filter(Assigned == "MT-L") %>%
  select(trinucleotideAlteration_SV) %>%
  table() %>%
  as.data.frame() %>%
  as_tibble() %>%
  rename(nt_ch = ".", true_count = Freq)

variants_2019_pten_MTL_nt_change_fractions <-
  variants_2019_pten_by_sig %>%
  droplevels() %>%
  filter(Assigned == "MT-L") %>%
  select(trinucleotideAlteration_SV) %>%
  table() %>%
  prop.table() %>%
  as.data.frame() %>%
  as_tibble() %>%
  rename(nt_ch = ".", true_fraction = Freq)

tricontext_MTL_pred_nt_change <- tricontext_MTL_pred %>%
  group_by(nt_ch) %>%
  summarise(nt_ch_pred_sum = sum(nt_ch_pred))

tricontext_MTL_pred_vs_real <- tricontext_MTL_pred %>%
  group_by(nt_ch) %>%
  summarise(nt_ch_pred_sum = sum(nt_ch_pred)) %>%
  mutate(nt_ch_prop = prop.table(nt_ch_pred_sum)) %>%
  select(-nt_ch_pred_sum) %>%
  rename(predicted_fr = nt_ch_prop) %>%
  left_join(variants_2019_pten_MTL_nt_change_fractions)
tricontext_MTL_pred_vs_real %>% write_csv("tricontext_MTL_pred_vs_real.csv")

# using chisq.test to determine whether the actual counts of nt substitutions
# are according to the proportions predicted by mutation signature patterns (which are expressed as rescaled vector of probabilities)
chisq.test(
  variants_2019_pten_MTL_nt_change_counts$true_count,
  p = tricontext_MTL_pred_nt_change$nt_ch_pred_sum,
  rescale.p = T
)
# p-value < 2.2e-16

