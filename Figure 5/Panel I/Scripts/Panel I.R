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

# Function Definitions ---------------------------------------------------------

# Function to Compute PTEN Mutation Presence and Prepare DataFrame for Plotting
prepare_dataframe_for_plotting <-
  function(variants_data, samples_data) {
    # Compute the PTEN mutation presence and join with sample data
    # Filter and transform the data to make it suitable for plotting
    prepared_dataframe <- variants_data %>%
      group_by(sample_id) %>%
      summarise(PTEN_mutation_presence = sign(sum(gene == "PTEN"))) %>%
      left_join(x = samples_data) %>%
      mutate(PTEN_mutation_presence = replace_na(PTEN_mutation_presence, 0)) %>%
      filter(stage != "0") %>%
      mutate(stage = case_when(stage == "I" ~ 1,
                               stage == "II" ~ 2,
                               stage == "III" ~ 3,
                               stage == "IV" ~ 4)) %>%
      group_by(stage, msi) %>%
      summarise(
        sample_count = n(),
        mutated_sample_count = sum(PTEN_mutation_presence)
      ) %>%
      mutate(
        proportion_of_mutations = mutated_sample_count / sample_count,
        standard_error = sqrt((
          proportion_of_mutations * (1 - proportion_of_mutations)
        ) / sample_count),
        confidence_interval_lower_bound = proportion_of_mutations - 1.96 * standard_error,
        confidence_interval_upper_bound = proportion_of_mutations + 1.96 * standard_error
      ) %>%
      filter(!is.na(msi))
    return(prepared_dataframe)
  }

# Function to Perform Fisher's Exact Test for Pairwise Comparisons Between Cancer Stages
perform_pairwise_fisher_tests <- function(data) {
  fisher_test_results <- data.frame()
  for (msi_level in unique(data$msi)) {
    msi_filtered_data <- subset(data, msi == msi_level)
    unique_stages <- unique(msi_filtered_data$stage)
    num_unique_stages <- length(unique_stages)
    
    if (num_unique_stages == 2) {
      # Special handling for only two stages
      stage_1 = unique_stages[1]
      stage_2 = unique_stages[2]
      stage_1_data = subset(msi_filtered_data, stage == stage_1)
      stage_2_data = subset(msi_filtered_data, stage == stage_2)
      contingency_table <- matrix(c(stage_1_data$mutated_sample_count, stage_1_data$sample_count - stage_1_data$mutated_sample_count,
                                    stage_2_data$mutated_sample_count, stage_2_data$sample_count - stage_2_data$mutated_sample_count),
                                  nrow = 2)
      fisher_test_outcome <- fisher.test(contingency_table)
      fisher_test_results <- rbind(fisher_test_results, data.frame(msi = msi_level, stage1 = stage_1, stage2 = stage_2, p_value = fisher_test_outcome$p.value))
    } else {
      # Existing logic for more than two stages
      pairwise_stage_combinations <- combn(unique_stages, 2)
      for (i in seq_along(pairwise_stage_combinations[1, ])) {
        stage_1 = pairwise_stage_combinations[1, i]
        stage_2 = pairwise_stage_combinations[2, i]
        stage_1_data = subset(msi_filtered_data, stage == stage_1)
        stage_2_data = subset(msi_filtered_data, stage == stage_2)
        contingency_table <- matrix(c(stage_1_data$mutated_sample_count, stage_1_data$sample_count - stage_1_data$mutated_sample_count,
                                      stage_2_data$mutated_sample_count, stage_2_data$sample_count - stage_2_data$mutated_sample_count),
                                    nrow = 2)
        fisher_test_outcome <- fisher.test(contingency_table)
        fisher_test_results <- rbind(fisher_test_results, data.frame(msi = msi_level, stage1 = stage_1, stage2 = stage_2, p_value = fisher_test_outcome$p.value))
      }
    }
  }
  
  fisher_test_results <- fisher_test_results %>%
    mutate(significance_annotation = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ))
  
  return(fisher_test_results)
}

# Function to Generate the Plot Using ggplot2
generate_publication_quality_plot <- function(data, test_results) {
  # Define aesthetic mappings and layers for the ggplot
  final_plot <-
    ggplot(data,
           aes(
             x = msi,
             fill = as.factor(stage),
             y = proportion_of_mutations * 100
           )) +
    geom_col(position = position_dodge2(preserve = "single")) +
    geom_errorbar(
      aes(ymin = confidence_interval_lower_bound * 100, ymax = confidence_interval_upper_bound * 100),
      width = 0.2,
      position = position_dodge(0.9)
    ) +
    scale_fill_brewer(palette = "Set1") + # Use a color palette suitable for publication
    labs(
      title = "Distribution of PTEN Mutations Across MSI and Cancer Stages",
      subtitle = "Data sourced from PAD database",
      x = "MSI Status",
      y = "Proportion of PTEN Mutations (%)",
      fill = "Cancer Stage"
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      plot.title = element_text(
        hjust = 0.5,
        size = 16,
        face = "bold"
      ),
      plot.subtitle = element_text(hjust = 0.5)
    )
  return(final_plot)
}

# Main Execution ---------------------------------------------------------------

# Read the Data (Modify these paths as per your directory structure)
PAD_sample_data <- fread("../../../Common raw data/PAD samples.csv")
PAD_variant_data <-
  fread("../../../Common raw data/PAD variants.csv")

# Prepare the DataFrame for Plotting
plotting_dataframe <-
  prepare_dataframe_for_plotting(PAD_variant_data, PAD_sample_data)

# Save the DataFrame for Chart Reproduction
write.csv(plotting_dataframe,
          "../Charts, data, statistics/Panel I data for reproduction.csv", row.names = F)

# Perform Pairwise Fisher's Exact Tests
fisher_test_results_dataframe <-
  perform_pairwise_fisher_tests(plotting_dataframe)

# Save the Fisher Test Results for Further Analysis
write.csv(
  fisher_test_results_dataframe,
  "../Charts, data, statistics/Panel I pairwise fisher test.csv", row.names = F
)

# Generate and Save the Publication-Quality Plot
publication_quality_plot <-
  generate_publication_quality_plot(plotting_dataframe, fisher_test_results_dataframe)

print(publication_quality_plot)

ggsave("../Charts, data, statistics/Panel I.png", dpi = 300, height = 4, width = 8)


# Supplement Execution ---------------------------------------------------------

# Read the Data (Modify these paths as per your directory structure)
PAD_sample_data_supp <- fread("../../../Common raw data/PAD samples.csv") %>% 
  mutate(stage = case_when(stage == "I" ~ "I",
                           stage == "II" ~ "I",
                           stage == "III" ~ "IV",
                           stage == "IV" ~ "IV"))

# Prepare the DataFrame for Plotting
plotting_dataframe_supp <-
  prepare_dataframe_for_plotting(PAD_variant_data, PAD_sample_data_supp)

# Save the DataFrame for Chart Reproduction
write.csv(plotting_dataframe_supp,
          "../Charts, data, statistics/Panel I data for reproduction supplement.csv", row.names = F)

# Perform Pairwise Fisher's Exact Tests
fisher_test_results_dataframe_supplement <-
  perform_pairwise_fisher_tests(plotting_dataframe_supp)

# Save the Fisher Test Results for Further Analysis
write.csv(
  fisher_test_results_dataframe_supplement,
  "../Charts, data, statistics/Panel I pairwise fisher test supplement.csv", row.names = F
)

# Generate and Save the Publication-Quality Plot
publication_quality_plot_supplement <-
  generate_publication_quality_plot(plotting_dataframe_supp, fisher_test_results_dataframe_supplement)

print(publication_quality_plot_supplement)

ggsave("../Charts, data, statistics/Panel I supplement.png", dpi = 300, height = 4, width = 5)
