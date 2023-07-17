# Attach requirement packages and setting WD -----------------------------------
packages_names <- c("readxl", "data.table", "reshape2", "tidyverse",
                    "rstudioapi")

lapply(packages_names, require, character.only = TRUE)

setwd(dirname(getActiveDocumentContext()$path))

PTEN_miss_hotspots <- c(15, 16, 17, 36, 93, 136)


# File with PTEN hotspots  -----------------------------------------------------
PTEN_hotspots_MTL <-
  read_xlsx("../../../Common raw data/PTEN hotspots MTL.xlsx") %>%
  select(htsp_cmbn_point) %>%
  mutate(htsp_status = "htsp") %>%
  left_join(x = data.frame(htsp_cmbn_point = 1:403), by = "htsp_cmbn_point") %>%
  mutate(htsp_status = replace_na(htsp_status, "non")) %>%
  as_tibble() %>%
  rename(codon = htsp_cmbn_point) %>% 
  mutate(htsp_status = ifelse(codon %in% PTEN_miss_hotspots, "non", htsp_status))

# Consurf scores obtained from -------------------------------------------------
# https://consurfdb.tau.ac.il/main_output.php?pdb_ID=1D5R&view_chain=A&unique_chain=1D5RA
consurf_scores <- fread("../Raw data/Consurf Scores.csv")

# PTEN data with LoF annotation  -----------------------------------------------
PTEN_lof_data <-
  fread("../../../Common raw data/PTEN LoF annotation.csv") %>%
  rename(codon = codon_SV) %>%
  filter(Assigned == "MT-L",
         !str_detect(alt, "splice|del|_|ext|ins")) %>%
  left_join(y = consurf_scores, multiple = "all", by = "codon") %>%
  left_join(PTEN_hotspots_MTL, by = "codon") %>%
  mutate(htsp_status = as.factor(htsp_status))

# Save chart data --------------------------------------------------------------
PTEN_lof_data %>%
  select(alt, consurf_score, htsp_status) %>%
  write.csv("../Charts, data, statistics/Panel G data.csv", row.names = F)

# Draw the chart Panel G -------------------------------------------------------
PTEN_lof_data %>% ggplot(aes(x = htsp_status, y = consurf_score)) +
  geom_violin(trim = T,
              fill = "#34cfeb",
              alpha = .6) +
  geom_boxplot(width = 0.05, col = "blue") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold")
  ) +
  ylab("Consurf score") +
  xlab("Hotspot or not")

# Save chart -------------------------------------------------------------------
ggsave(
  "../Charts, data, statistics/Panel G.png",
  units = "in",
  dpi = 500,
  width = 6.65,
  height = 5.94
)
