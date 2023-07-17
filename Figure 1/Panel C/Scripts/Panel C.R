# Attach requirement packages and setting WD -----------------------------------
packages_names <- c("readxl", "data.table", "reshape2", "tidyverse",
                    "rstudioapi")

lapply(packages_names, require, character.only = TRUE)

rename = dplyr::rename
select = dplyr::select
filter = dplyr::filter
group_by = dplyr::group_by
mutate = dplyr::mutate

setwd(dirname(getActiveDocumentContext()$path))

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

# Calculation of FMI and PAD hotspots ------------------------------------------
FMI_PTEN_hotspots <- get_hotspots(filter(FMI_variants, gene == "PTEN"),
                                  403, "codon") %>%
  mutate(Freq_prop = Freq / nrow(filter(FMI_variants, gene == "PTEN"))) %>%
  arrange(desc(Freq_prop))

PAD_PTEN_hotspots <- get_hotspots(filter(PAD_variants, gene == "PTEN"),
                                  403, "codon") %>%
  mutate(Freq_prop = Freq / nrow(filter(PAD_variants, gene == "PTEN"))) %>%
  arrange(desc(Freq_prop))

# Get hotspots from both and combine them --------------------------------------
PAD_FMI_hotspots <- sort(unique(c(
  filter(PAD_PTEN_hotspots,
         Freq >= cutoff)$codon,
  filter(FMI_PTEN_hotspots,
         Freq >= cutoff)$codon
)))

# Divide each frequency by 130 codon frequency ---------------------------------
PAD_PTEN_hotspots <- PAD_PTEN_hotspots %>%
  mutate(freq_adj = Freq / filter(PAD_PTEN_hotspots, codon == 130)$Freq)

FMI_PTEN_hotspots <- FMI_PTEN_hotspots %>%
  mutate(freq_adj = Freq / filter(FMI_PTEN_hotspots, codon == 130)$Freq)

# Add PTEN PTM coordinates to hotspots data frame ------------------------------
annotated_hotspots <-
  filter(PAD_PTEN_hotspots, codon %in% PAD_FMI_hotspots) %>%
  rename(PAD = freq_adj) %>%
  select(codon, PAD) %>%
  left_join(select(rename(
    filter(FMI_PTEN_hotspots, codon %in% PAD_FMI_hotspots),
    FMI = freq_adj
  ),
  codon, FMI),
  by = "codon") %>%
  full_join(PTEN_PTM, by = "codon") %>%
  pivot_longer(
    cols = c("PAD", "FMI"),
    names_to = "data",
    values_to = "value"
  ) %>%
  mutate(data = as.factor(data)) %>%
  ungroup()

# Just graphic-used constant ---------------------------------------------------
domain_coeff <- 0.07

# Draw the chart Panel C -------------------------------------------------------
ggplot() +
  geom_rect(
    data = PTEN_domains,
    aes(
      xmin = start,
      xmax = end,
      ymin = (c(1, 1, 0.5, 1, 1) * -domain_coeff) - domain_coeff,
      ymax = (c(1, 1, 0.5, 1, 1) * domain_coeff) - domain_coeff,
      fill = as.factor(domain)
    ),
    alpha = 0.9
  ) +
  geom_segment(
    data = annotated_hotspots,
    aes(
      x = codon,
      xend = codon,
      y = 0,
      yend = value
    ),
    linewidth = 1.5,
    na.rm = TRUE
  ) +
  geom_point(
    data = annotated_hotspots,
    aes(x = codon, y = value, col = data),
    size = 5,
    na.rm = TRUE
  ) +
  geom_point(
    data = filter(annotated_hotspots, !is.na(impact)),
    aes(
      x = codon,
      y = -0.04,
      fill = impact,
      col = impact
    ),
    size = 6.5,
    shape = 24
  ) +
  theme_minimal() +
  theme(legend.position = "top") +
  geom_hline(yintercept = 0, col = "#413F42") +
  scale_fill_manual(
    values = c(
      active_site = "#489FF8",
      PTM = "NA",
      PDB = "#E4ACFA",
      PTPase = "#FDEFB2",
      C2 = "#B9EBFD",
      `C-tail` = "#9BDCAE",
      None = "gray"
    ),
    labels = c("Active site", "", "PDB", "PTPase", "C2", "C-tail", ""),
    guide = guide_legend(
      nrow = 2,
      override.aes = list(
        shape = c(17, NA, 0, 0, 0, 0, NA),
        fill = c(NA, NA, "#E4ACFA", "#FDEFB2", "#B9EBFD", "#9BDCAE", NA),
        linetype = c(0, 0, 0, 0, 0, 0, 0),
        col = c("#489FF8", "red", NA, NA, NA, NA, NA)
      )
    )
  ) +
  scale_colour_manual(
    values = c(
      FMI = "#7139B9",
      PAD = "#FFB900",
      PTM = "NA",
      active_site = "NA"
    ),
    labels = c("FMI", "PAD", "", ""),
    guide = guide_legend(override.aes = list(
      linetype = c(0, 0, 0, 0),
      shape = c(16, 16, 0, 0),
      color = c("#7139B9", "#FFB900", NA, NA)
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
    limits = c(-0.2, 1.3),
    expand = c(0, 0),
    breaks = c(0, 1, 1.3)
  ) +
  scale_x_continuous(
    limits = c(0, 407),
    expand = c(0, 0),
    breaks = c(0, 100, 200, 300, 403)
  ) +
  ylab("Relative abundance") +
  xlab("Codon position")

# Save chart data --------------------------------------------------------------
write.csv(annotated_hotspots,
          "../Charts, data, statistics/Panel C chart data.csv",
          row.names = FALSE)

write.csv(FMI_PTEN_hotspots,
          "../Charts, data, statistics/Panel C FMI hotspots.csv",
          row.names = FALSE)

write.csv(PAD_PTEN_hotspots,
          "../Charts, data, statistics/Panel C PAD hotspots.csv",
          row.names = FALSE)

# Save image -------------------------------------------------------------------
ggsave(
  "../Charts, data, statistics/Panel C.png",
  dpi = 600,
  width = 18.18,
  height = 8,
  units = "in"
)
