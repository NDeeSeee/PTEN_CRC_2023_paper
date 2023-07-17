# Attach requirement packages and setting WD -----------------------------------
packages_names <- c(
  "readxl",
  "data.table",
  "reshape2",
  "tidyverse",
  "rstudioapi",
  "survminer",
  "survival"
)

lapply(packages_names, require, character.only = TRUE)

setwd(dirname(getActiveDocumentContext()$path))

rename <- dplyr::rename
select <- dplyr::select
filter <- dplyr::filter
group_by <- dplyr::group_by
mutate <- dplyr::mutate

# Load GENIE OS data ------------------------------------------------------------
GENIE_OS_data <-
  fread("../../../Common raw data/PTEN GENIE OS data.csv")

# Fit with survival ------------------------------------------------------------
GENIE_OS_data_mtl <- filter(GENIE_OS_data)

GENIE_OS_data_mtl_wtmut <-
  filter(GENIE_OS_data,
         msi == "MSS",
         category == "PTEN_mut" | category == "PTEN_wt")

GENIE_OS_data_mtl_wtdel <-
  filter(GENIE_OS_data,
         msi == "MSS",
         category == "PTEN_del" | category == "PTEN_wt")

GENIE_OS_data_mth_wtmut <-
  filter(GENIE_OS_data,
         msi == "MSI-H",
         category == "PTEN_mut" | category == "PTEN_wt")

fit_GENIE_mtl <-
  surv_fit(Surv(time = time_months, status) ~ category, data = GENIE_OS_data_mtl)

fit_GENIE_mtl_wtmut <-
  surv_fit(Surv(time = time_months, status) ~ category, data = GENIE_OS_data_mtl_wtmut)

fit_GENIE_mtl_wtdel <-
  surv_fit(Surv(time = time_months, status) ~ category, data = GENIE_OS_data_mtl_wtdel)

fit_GENIE_mth_wtmut <-
  surv_fit(Surv(time = time_months, status) ~ category, data = GENIE_OS_data_mth_wtmut)

# Plot the charts and save -----------------------------------------------------
ggsurvplot(
  fit_GENIE_mtl,
  # survfit object with calculated statistics.
  data = GENIE_OS_data_mtl,
  # data used to fit survival curves.
  pval = FALSE,
  # show p-value of log-rank test.
  conf.int = TRUE,
  # show confidence intervals for
  # point estimates of survival curves.
  xlim = c(0, 150),
  # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in months",
  # customize X axis label.
  break.time.by = 50,
  # break X axis in time intervals by 500.
  ggtheme = theme_light(),
  legend.labs =
    c("Deletion", "Mutation", "WT") # customize plot and risk table with a theme.
  # in legend of risk table
)

ggsave(
  "../Charts, data, statistics/Panel E.png",
  width = 5,
  height = 5,
  dpi = 600,
  units = "in"
)

# PTEN MUT vs WT, MT-L ---------------------------------------------------------
ggsurvplot(
  fit_GENIE_mtl_wtmut,
  data = GENIE_OS_data_mtl_wtmut,
  pval = FALSE,
  conf.int = TRUE,
  xlim = c(0, 150),
  xlab = "Time in months",
  break.time.by = 50,
  ggtheme = theme_light(),
  legend.labs =
    c("Mutation", "WT"),
  palette =
    c("#00BA38", "#619CFF")
)

ggsave(
  "../Charts, data, statistics/Panel F.png",
  width = 5,
  height = 5,
  dpi = 600,
  units = "in"
)

# PTEN DEL vs WT, MT-L ---------------------------------------------------------
ggsurvplot(
  fit_GENIE_mtl_wtdel,
  data = GENIE_OS_data_mtl_wtdel,
  pval = FALSE,
  conf.int = TRUE,
  xlim = c(0, 150),
  xlab = "Time in months",
  break.time.by = 50,
  ggtheme = theme_light(),
  legend.labs =
    c("Deletion", "WT"),
  palette =
    c("#F8766D", "#619CFF")
)

ggsave(
  "../Charts, data, statistics/Panel G.png",
  width = 5,
  height = 5,
  dpi = 600,
  units = "in"
)

# PTEN MUT vs WT, MT-H ---------------------------------------------------------
ggsurvplot(
  fit_GENIE_mth_wtmut,
  data = GENIE_OS_data_mth_wtmut,
  pval = FALSE,
  conf.int = TRUE,
  xlim = c(0, 150),
  xlab = "Time in months",
  break.time.by = 50,
  ggtheme = theme_light(),
  legend.labs =
    c("Mutation", "WT"),
  palette =
    c("#00BA38", "#619CFF")
)

ggsave(
  "../Charts, data, statistics/Panel H.png",
  width = 5,
  height = 5,
  dpi = 600,
  units = "in"
)
