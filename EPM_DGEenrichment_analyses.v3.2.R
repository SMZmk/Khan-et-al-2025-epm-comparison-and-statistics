library(tidyverse)
library(tidygraph)
library(dplyr)
library(tidyr)
library(igraph)
library(ggraph)
library(ggplot2)
library(ggrepel)
library(topGO)
library(org.At.tair.db)
library(viridisLite) 
library(RColorBrewer)
library(pheatmap)
#########################
#########################################################
# --- Define Dynamic Parameters ---
species_tag <- "Asag"       # e.g., "Asag", "Athaliana"
condition_tag <- "wilting" # e.g., "recovery", "stress_treatment"
output_dir <- "../studies"  # Directory to save results
p_value_threshold <- 0.001 # Significance threshold for filtering GO terms
z_threshold <- 2 # Significance treshold for filtering EPM ins DGEs enrichments
##################################################################################
# Set the working directory
setwd("/home/ibg-4/Desktop/arab_env_2024i/workflows")
# Define the file paths for Arabis nemorensis EPMs and Quadrants
file_paths_file_EPM <- c(
  "../studies/AnemAsag_EPMclu/2025dAsagSWAsagSW_E-04_2025d_gene_none-mima.csv")
file_paths_file_Q0S <- c(
  "../studies/AnemAsag_dge_quadrants_khan/Asag-D0vW-dge1_reg0.csv")
# file_paths_file_Q0W <- c(
#   "../studies/AnemAsag_dge_quadrants_khan/recovery_ctrl_quadrants.csv")
file_paths_file_ortho <- c(
  "../studies/AnemAsag_dge_quadrants_khan/arabis_sag_hlx2024h.pep.tsv")
ortho_df <- read.csv(file_paths_file_ortho, sep = "\t")
###################################################################################
epm_df <- read.csv(file_paths_file_EPM)
Q0S_df <- read.csv(file_paths_file_Q0S)
#Q0W_df <- read.csv(file_paths_file_Q0W)
###################################################################################
# Select the columns "epm" and "loc_ID" from epm_df to create a new data frame epm_df0
epm_df0 <- epm_df %>%
  dplyr::select(epm, loc_ID)
head(epm_df0)
# Merge epm_df0 and Q0S_df by loc_ID and gene_id
merged_df_Q0S <- epm_df0 %>%
  inner_join(Q0S_df, by = c("loc_ID" = "gene_id"))

# Merge epm_df0 and Q0W_df by loc_ID and gene_id
# Assuming Q0W_df has the same structure as Q0S_df
# merged_df_Q0W <- epm_df0 %>%
#   inner_join(Q0W_df, by = c("loc_ID" = "gene_id"))

# Display the first few rows of the merged data frames
head(merged_df_Q0S)
#head(merged_df_Q0W)
####################################################################################

# Calculate the total number of EPMs per loc_ID across all quadrants
total_epms_per_loc <- merged_df_Q0S %>%
  group_by(loc_ID) %>%
  summarize(total_epms = n())

head(total_epms_per_loc)

# Calculate the number of distinct EPMs and the ratio of EPMs to the number of genes containing these EPMs
epm_gene_ratio <- merged_df_Q0S %>%
  group_by(epm) %>%
  summarize(
    epm_count = n(),
    unique_genes = n_distinct(loc_ID),
    ratio = epm_count / unique_genes
  )

# Display the first few rows of the proportions data frame
head(epm_gene_ratio)

head(merged_df_Q0S)
################################################################################

# Calculate the number of distinct EPMs and the ratio of EPMs to the number of genes containing these EPMs for each quadrant
epm_gene_ratio_by_regulation <- merged_df_Q0S %>%
  group_by(epm, regulation) %>%
  summarize(
    epm_count = n(),
    unique_genes = n_distinct(loc_ID),
    ratio = epm_count / unique_genes,
    .groups = 'drop'
  )

# Display the result
print(epm_gene_ratio_by_regulation)

# Transform the data frame to the desired format
epm_ratio_wide <- epm_gene_ratio_by_regulation %>%
  dplyr::select(epm, regulation, ratio) %>%
  pivot_wider(names_from = regulation, values_from = ratio)

# Display the result
print(epm_ratio_wide)

# --- Function to Calculate Relative Z-scores and Significance ---
# This function takes the DO, NS, and UP values for a single row
# and calculates Z-scores relative to NS.
calculate_relative_z <- function(do_val, ns_val, up_val, threshold = 2) {
  
  values <- c(do_val, ns_val, up_val)
  
  # Calculate standard deviation across the 3 values for the row
  row_sd <- sd(values, na.rm = TRUE)
  
  # Initialize results
  z_do_vs_ns <- NA
  z_up_vs_ns <- NA
  sig_do <- NA # Using NA for indeterminate cases
  sig_up <- NA # Using NA for indeterminate cases
  enrichment_direction <- "Indeterminate"
  
  # Check if sd calculation is possible and meaningful
  if (!is.na(row_sd) && length(na.omit(values)) >= 2) {
    if (row_sd == 0) {
      # If sd is 0, all values are the same (or only one non-NA value)
      z_do_vs_ns <- 0
      z_up_vs_ns <- 0
      sig_do <- FALSE
      sig_up <- FALSE
      enrichment_direction <- "No_Difference"
    } else {
      # Calculate Z-scores relative to NS
      z_do_vs_ns <- (do_val - ns_val) / row_sd
      z_up_vs_ns <- (up_val - ns_val) / row_sd
      
      # Check significance
      sig_do <- abs(z_do_vs_ns) > threshold
      sig_up <- abs(z_up_vs_ns) > threshold
      
      # Determine enrichment direction based on significant deviations
      is_do_sig_higher <- sig_do && z_do_vs_ns > 0
      is_up_sig_higher <- sig_up && z_up_vs_ns > 0
      is_do_sig_lower <- sig_do && z_do_vs_ns < 0
      is_up_sig_lower <- sig_up && z_up_vs_ns < 0
      
      if (is_do_sig_higher && is_up_sig_higher) {
        enrichment_direction <- "Both_Higher_vs_NS"
      } else if (is_do_sig_higher && is_up_sig_lower) {
        enrichment_direction <- "DO_Higher_UP_Lower_vs_NS"
      } else if (is_do_sig_lower && is_up_sig_higher) {
        enrichment_direction <- "DO_Lower_UP_Higher_vs_NS"
      } else if (is_do_sig_lower && is_up_sig_lower) {
        enrichment_direction <- "Both_Lower_vs_NS"
      } else if (is_do_sig_higher) {
        enrichment_direction <- "DO_Higher_vs_NS"
      } else if (is_up_sig_higher) {
        enrichment_direction <- "UP_Higher_vs_NS"
      } else if (is_do_sig_lower) {
        enrichment_direction <- "DO_Lower_vs_NS"
      } else if (is_up_sig_lower) {
        enrichment_direction <- "UP_Lower_vs_NS"
      } else {
        enrichment_direction <- "None_Significant_vs_NS"
      }
    }
  } # else: keep initial NA/Indeterminate values if sd is NA
  
  # Return results as a tibble row (makes unnesting easier)
  return(tibble(
    Z_Score_DO_vs_NS = z_do_vs_ns,
    Z_Score_UP_vs_NS = z_up_vs_ns,
    Significant_DO_vs_NS = sig_do,
    Significant_UP_vs_NS = sig_up,
    Enrichment_vs_NS = enrichment_direction
  ))
}

# --- Apply the function row-wise ---
# Set the significance threshold
#z_threshold <- 2

# Use rowwise() to apply the function to each row
epm_results_Zscores <- epm_ratio_wide %>%
  rowwise() %>%
  # Calculate stats for the row, storing results in a list column
  mutate(stats = list(calculate_relative_z(DO, NS, UP, threshold = z_threshold))) %>%
  ungroup() %>% # Essential before unnesting
  unnest(stats) # Expand the list column into regular columns

# --- Display the results ---
print(epm_results_Zscores)

################################################################################

##################################################################################
##################################################################################
##################################################################################
# --- Start of Adjusted Code ---

# Step 1: Calculate counts and ratios per epm for EACH regulation status (UP, DO, NS)
epm_gene_ratio_by_regulation <- merged_df_Q0S %>%
  # Ensure regulation is a character and filter out potential NAs if necessary
  mutate(regulation = as.character(regulation)) %>%
  filter(regulation %in% c("UP", "DO", "NS")) %>% # Keep UP, DO, and NS genes
  group_by(epm, regulation) %>%             # Group by epm and regulation status
  summarize(
    epm_count = n(),                     # Count rows
    unique_genes = n_distinct(loc_ID),   # Count unique genes
    ratio = ifelse(unique_genes > 0, epm_count / unique_genes, 0), # Calculate ratio
    .groups = 'drop'
  )

# Step 2: Pivot to wide format
# This creates columns like epm_count_UP, unique_genes_DO, ratio_NS, etc.
epm_enrich <- epm_gene_ratio_by_regulation %>%
  pivot_wider(names_from = regulation, # Use regulation column for new names
              values_from = c(epm_count, unique_genes, ratio),
              names_sep = "_") %>%      # Creates e.g., ratio_UP, ratio_NS
  # Replace potential NAs (if an EPM lacks UP, DO, or NS genes) with 0 for counts/genes
  # Ratios will remain NA if components were missing, handled by tests below.
  mutate(
    across(starts_with("epm_count_"), ~ replace_na(.x, 0)),
    across(starts_with("unique_genes_"), ~ replace_na(.x, 0))
  )

# Note: The background calculation using epm_df0 and the subsequent merge are now REMOVED.
# 'epm_enrich' now contains UP, DO, and NS (background) information.

print("Pivoted Ratios per EPM (UP, DO, NS):")
print(epm_enrich)

##################################################################################
# Step 3: 1. Z-SCORE TEST (Comparing Rates: UP/DO vs NS)
# Compares the rate (ratio) in UP/DO vs the NS "background"
epm_enrich <- epm_enrich %>%
  mutate(
    # Z-score for UP vs NS Background
    z_UP = ifelse(unique_genes_UP > 0 & unique_genes_NS > 0 & !is.na(ratio_UP) & !is.na(ratio_NS) &
                    (ratio_UP / unique_genes_UP) >= 0 & (ratio_NS / unique_genes_NS) >= 0, # Check variances >= 0
                  # Formula: (Rate1 - Rate2) / sqrt(Var1 + Var2) where Var ~ Rate / N
                  (ratio_UP - ratio_NS) / sqrt((ratio_UP / unique_genes_UP) + (ratio_NS / unique_genes_NS)),
                  NA_real_),
    
    # Z-score for DO vs NS Background
    z_DO = ifelse(unique_genes_DO > 0 & unique_genes_NS > 0 & !is.na(ratio_DO) & !is.na(ratio_NS) &
                    (ratio_DO / unique_genes_DO) >= 0 & (ratio_NS / unique_genes_NS) >= 0, # Check variances >= 0
                  (ratio_DO - ratio_NS) / sqrt((ratio_DO / unique_genes_DO) + (ratio_NS / unique_genes_NS)),
                  NA_real_)
  )

# Flag significance using a threshold (e.g., |Z| > 2)
epm_enrich <- epm_enrich %>%
  mutate(
    sig_UP_Z = ifelse(!is.na(z_UP) & abs(z_UP) > 2, TRUE, FALSE),
    sig_DO_Z = ifelse(!is.na(z_DO) & abs(z_DO) > 2, TRUE, FALSE)
  )

print("Data with Z-Scores (vs NS background):")
print(dplyr::select(epm_enrich, epm, ratio_UP, ratio_DO, ratio_NS, z_UP, sig_UP_Z, z_DO, sig_DO_Z))

##################################################################################
# Step 4: 2. POISSON RATE RATIO TEST (Log Rate Ratio Z-score vs NS)
# Compares rates using log ratio and standard error based on counts vs NS group.
epm_enrich <- epm_enrich %>%
  mutate(
    # --- UP vs NS Background ---
    log_rr_UP = ifelse(epm_count_UP > 0 & unique_genes_UP > 0 & epm_count_NS > 0 & unique_genes_NS > 0 &
                         !is.na(ratio_UP) & !is.na(ratio_NS) & ratio_UP > 0 & ratio_NS > 0, # Check rates > 0 for log
                       log(ratio_UP / ratio_NS), # Compare UP to NS
                       NA_real_),
    se_log_rr_UP = ifelse(!is.na(log_rr_UP),
                          # Use counts from UP and NS groups for SE
                          sqrt(1 / epm_count_UP + 1 / epm_count_NS),
                          NA_real_),
    z_rate_UP = ifelse(!is.na(log_rr_UP) & !is.na(se_log_rr_UP) & se_log_rr_UP > 0, # Check SE > 0
                       log_rr_UP / se_log_rr_UP,
                       NA_real_),
    p_rate_UP = ifelse(!is.na(z_rate_UP),
                       2 * pnorm(-abs(z_rate_UP)), # Two-tailed p-value
                       NA_real_),
    
    # --- DO vs NS Background ---
    log_rr_DO = ifelse(epm_count_DO > 0 & unique_genes_DO > 0 & epm_count_NS > 0 & unique_genes_NS > 0 &
                         !is.na(ratio_DO) & !is.na(ratio_NS) & ratio_DO > 0 & ratio_NS > 0, # Check rates > 0
                       log(ratio_DO / ratio_NS), # Compare DO to NS
                       NA_real_),
    se_log_rr_DO = ifelse(!is.na(log_rr_DO),
                          # Use counts from DO and NS groups for SE
                          sqrt(1 / epm_count_DO + 1 / epm_count_NS),
                          NA_real_),
    z_rate_DO = ifelse(!is.na(log_rr_DO) & !is.na(se_log_rr_DO) & se_log_rr_DO > 0, # Check SE > 0
                       log_rr_DO / se_log_rr_DO,
                       NA_real_),
    p_rate_DO = ifelse(!is.na(z_rate_DO),
                       2 * pnorm(-abs(z_rate_DO)),
                       NA_real_)
  )

print("Data with Rate Ratio Test Results (vs NS background):")
# Check if columns were created before printing
required_cols <- c("epm", "log_rr_UP", "z_rate_UP", "p_rate_UP", "log_rr_DO", "z_rate_DO", "p_rate_DO")
if(all(required_cols %in% colnames(epm_enrich))) {
  print(dplyr::select(epm_enrich, all_of(required_cols)))
} else {
  print("Rate ratio columns were not created successfully. Check calculations and input data.")
  print("Available columns:")
  print(colnames(epm_enrich))
}


##################################################################################
# Step 5: 3. PAIRED T-TEST ACROSS EPMs (vs NS)
# Compare UP/DO rates vs NS rates across EPMs. Requires >= 2 EPMs with valid pairs.
epm_up_valid <- epm_enrich %>% filter(!is.na(ratio_UP) & !is.na(ratio_NS)) # Compare UP to NS
if (nrow(epm_up_valid) >= 2) {
  # Use tryCatch in case t.test fails (e.g., constant difference)
  ttest_UP <- tryCatch(t.test(epm_up_valid$ratio_UP, epm_up_valid$ratio_NS, paired = TRUE),
                       error = function(e) { print(paste("UP vs NS t-test failed:", e)); return(NULL) })
  if (!is.null(ttest_UP)) {
    print("Paired T-test: UP Ratios vs NS Ratios")
    print(ttest_UP)
  }
} else {
  print("Paired T-test for UP vs NS: Not enough valid paired observations.")
}

epm_do_valid <- epm_enrich %>% filter(!is.na(ratio_DO) & !is.na(ratio_NS)) # Compare DO to NS
if (nrow(epm_do_valid) >= 2) {
  ttest_DO <- tryCatch(t.test(epm_do_valid$ratio_DO, epm_do_valid$ratio_NS, paired = TRUE),
                       error = function(e) { print(paste("DO vs NS t-test failed:", e)); return(NULL) })
  if (!is.null(ttest_DO)) {
    print("Paired T-test: DO Ratios vs NS Ratios")
    print(ttest_DO)
  }
} else {
  print("Paired T-test for DO vs NS: Not enough valid paired observations.")
}

##################################################################################
# Step 6: 4. Per-EPM Exact Poisson Test (poisson.test vs NS)
# Compare counts directly (UP/DO vs NS), adjusting for exposure (unique genes).
epm_enrich <- epm_enrich %>%
  rowwise() %>%
  mutate(
    # Compare UP counts vs NS counts
    p_poisson_UP = ifelse(epm_count_UP >= 0 & unique_genes_UP > 0 & epm_count_NS >= 0 & unique_genes_NS > 0,
                          tryCatch(poisson.test(c(epm_count_UP, epm_count_NS), # Compare UP to NS counts
                                                T = c(unique_genes_UP, unique_genes_NS))$p.value,
                                   error = function(e) NA_real_),
                          NA_real_),
    # Compare DO counts vs NS counts
    p_poisson_DO = ifelse(epm_count_DO >= 0 & unique_genes_DO > 0 & epm_count_NS >= 0 & unique_genes_NS > 0,
                          tryCatch(poisson.test(c(epm_count_DO, epm_count_NS), # Compare DO to NS counts
                                                T = c(unique_genes_DO, unique_genes_NS))$p.value,
                                   error = function(e) NA_real_),
                          NA_real_)
  ) %>%
  ungroup()

print("Final Enrichment Data with All Tests (vs NS background):")
print(epm_enrich)

#################################################################################
# Step 7: Create Grouping Variable for Plotting (Based on Z-Score Significance vs NS)
epm_enrich <- epm_enrich %>%
  mutate(
    sig_UP_Z_plot = ifelse(is.na(sig_UP_Z), FALSE, sig_UP_Z),
    sig_DO_Z_plot = ifelse(is.na(sig_DO_Z), FALSE, sig_DO_Z),
    sig_combined_z = ifelse(sig_UP_Z_plot | sig_DO_Z_plot, "Significant", "Not Significant"),
    epm_group = case_when(
      sig_combined_z == "Significant" & grepl("p0", epm) ~ "p0",
      sig_combined_z == "Significant" & grepl("p1", epm) ~ "p1",
      TRUE ~ "Not Significant"
    )
  )

# Define custom colors
custom_colors <- c("p0" = "black", "p1" = "black", "Not Significant" = "grey")

#################################################################################
# Step 8: Create the Z-score Plot (vs NS)

# Replace NA Z-scores with 0 for plotting
epm_enrich_plot <- epm_enrich %>%
  mutate(z_UP_plot = replace_na(z_UP, 0),
         z_DO_plot = replace_na(z_DO, 0))

# Check if there's data to plot
if(nrow(epm_enrich_plot) > 0) {
  p <- ggplot(epm_enrich_plot, aes(x = z_UP_plot, y = z_DO_plot, color = epm_group)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text_repel(
      data = filter(epm_enrich_plot, epm_group %in% c("p0", "p1")),
      aes(label = epm),
      size = 3,
      box.padding = 0.4, point.padding = 0.6, segment.color = "grey50",
      max.overlaps = 15
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0, linetype = "solid", size = 0.3, color = "darkgrey") +
    geom_vline(xintercept = 0, linetype = "solid", size = 0.3, color = "darkgrey") +
    geom_hline(yintercept = c(-2, 2), linetype = "dotted", color = "blue", size = 0.5) +
    geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = "blue", size = 0.5) +
    scale_color_manual(values = custom_colors, name = "EPM Group (Significant vs NS)") +
    labs(title = "EPM Enrichment: UP vs DOWN Regulated Genes",
         subtitle = "Z-Scores compared to NS (Non-Significant) Ratio",
         x = "Z-Score (UP vs NS)",
         y = "Z-Score (DOWN vs NS)") +
    coord_cartesian(xlim = c(min(epm_enrich_plot$z_UP_plot, -2.5, na.rm=T) - 0.5, max(epm_enrich_plot$z_UP_plot, 2.5, na.rm=T) + 0.5),
                    ylim = c(min(epm_enrich_plot$z_DO_plot, -2.5, na.rm=T) - 0.5, max(epm_enrich_plot$z_DO_plot, 2.5, na.rm=T) + 0.5)) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
  
  # Display plot
  print(p)
  
  # Save plot
  if (!dir.exists("../studies")) {dir.create("../studies", recursive = TRUE)}
  ggsave(filename = "../studies/EPM_up-doDGE_enrichment-recovery_ZscorePlot_vs_NS.png", plot = p,
         width = 8, height = 7, dpi = 300)
} else {
  print("No data available to generate the Z-score plot.")
}


######################################################################################
# Step 9: Select and Save Enriched EPMs and Corresponding Genes (Based on vs NS significance)

epms_enriched_combined <- epm_enrich %>%
  filter(sig_combined_z == "Significant") # Based on Z-score vs NS

epms_enriched_combined_unique <- epms_enriched_combined %>%
  distinct(epm, .keep_all = TRUE)

# Add the 'direction' column back to merged_df_Q0S if it was removed, or use regulation
# Let's re-add it for clarity if needed downstream, based on regulation column
merged_df_Q0S <- merged_df_Q0S %>%
  mutate(direction = ifelse(regulation %in% c("UP", "DO"), regulation, NA_character_))

# Get the subset of original DGE results (UP/DO genes only) for significant EPMs
epm_DGEenriched_subset <- merged_df_Q0S %>%
  filter(epm %in% epms_enriched_combined_unique$epm, # EPM is significant vs NS
         !is.na(direction))                          # Only UP/DO genes

epm_DGEenriched_subset0 <- epm_DGEenriched_subset %>%
  dplyr::select(epm, loc_ID, direction) # Select columns for GO

# Print summaries
print("Significant EPMs (any direction vs NS):")
if(nrow(epms_enriched_combined_unique)>0) print(dplyr::select(epms_enriched_combined_unique, epm, z_UP, sig_UP_Z, z_DO, sig_DO_Z)) else print("None")
print("Head of genes from significant EPMs (for GO):")
if(nrow(epm_DGEenriched_subset0)>0) print(head(epm_DGEenriched_subset0)) else print("None")
print(epms_enriched_combined_unique)
# --- Save Results ---
write.csv(epm_enrich, file = "../studies/epms_DEGenriched_recovery_full_stats_vs_NS.csv", row.names = FALSE)
write.csv(epms_enriched_combined_unique, file = "../studies/epms_DEGenriched_significant_recovery_vs_NS.csv", row.names = FALSE)
write.csv(epm_DGEenriched_subset0, file = "../studies/genes_from_DEGenriched_epms_recovery_vs_NS.csv", row.names = FALSE)

#################################################################################

print("--- Analysis Chunk Completed (Using NS as Background) ---")

#################################################################################
#################################################################################
#################################################################################
################################################################################# CONTINUE HERE 02042025

head(ortho_df)
colnames(ortho_df)

new_column_names <- c("Orthogroup", "Species", "target_gene_id", "At_orthologs")
# Check if the number of new names matches the number of columns (optional but good practice)
if (length(new_column_names) == ncol(ortho_df)) {
  # Assign the new names to the data frame
  colnames(ortho_df) <- new_column_names
  
  # Print the new column names to verify the change
  print(colnames(ortho_df))
  
  # You can also view the first few rows with the new names
  print(head(ortho_df))
  
} else {
  print(paste("Error: The data frame has", ncol(ortho_df), "columns, but you provided", length(new_column_names), "names."))
}
ortho_df <- ortho_df %>%
  filter(Species == "Arabidopsis_thaliana.TAIR10.59.pep.all")
ortho_long_df <- ortho_df %>%
  # REMOVED: rename(target_gene_id = target_gene_id, At_orthologs = At_orthologs) %>%
  separate_rows(target_gene_id, sep = ",\\s*") %>% # Splits by comma potentially followed by space(s)
  separate_rows(At_orthologs, sep = ",\\s*") %>%   # Splits by comma potentially followed by space(s)
  dplyr::select(Orthogroup, Species, target_gene_id, At_orthologs) # Keep only desired columns (dplyr:: is optional)

# View the result
print(head(ortho_long_df))
# Process the data frame
ortho_cleaned_df <- ortho_long_df %>%
  # 1. Select only the target_gene_id and At_orthologs columns
  dplyr::select(target_gene_id, At_orthologs) %>%
  
  # 2. Remove version suffixes (e.g., ".1", ".2") from both columns
  #    We use mutate() with across() to apply the same function to multiple columns.
  #    str_remove() removes the first match of the pattern.
  #    The pattern "\\.\\d+$" matches a literal dot `\.` followed by
  #    one or more digits `\d+` at the very end of the string `$`.
  mutate(across(c(target_gene_id, At_orthologs), ~ str_remove(.x, "\\.\\d+$"))) %>%
  # Alternative using base R's sub():
  # mutate(across(c(target_gene_id, At_orthologs), ~ sub("\\.\\d+$", "", .x))) %>%
  
  # 3. Remove duplicate rows based on the modified columns
  distinct()

#dev.off()
##################################################################################

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# --- 1. Data Preparation ---

# Assuming ortho_long_df is already loaded and contains the expanded ortholog pairs
# Clean ortholog data (select columns, remove version suffixes, remove duplicates)
ortho_cleaned_df <- ortho_long_df %>%
  dplyr::select(target_gene_id, At_orthologs) %>%
  mutate(across(c(target_gene_id, At_orthologs), ~ str_remove(.x, "\\.\\d+$"))) %>%
  distinct()

# Assuming epm_DGEenriched_subset0 is loaded
# Merge DGE subset with cleaned orthologs
merged_data <- epm_DGEenriched_subset0 %>%
  left_join(ortho_cleaned_df,
            by = c("loc_ID" = "target_gene_id"),
            relationship = "many-to-many") %>%
  filter(!is.na(At_orthologs)) # Keep only rows with a matched Arabidopsis ortholog

cat("Head of merged data:\n")
print(head(merged_data))

# Define the universe of genes (all unique Arabidopsis orthologs in the merged dataset)
universe_genes <- unique(merged_data$At_orthologs)
cat("\nTotal unique Arabidopsis orthologs in universe:", length(universe_genes), "\n")

# Get list of unique EPM conditions
epm_list <- unique(merged_data$epm)
cat("EPM conditions found:", paste(epm_list, collapse=", "), "\n\n")


# --- 2. GO Enrichment Function ---

run_go_enrichment <- function(data, epm_conditions, universe, direction_filter,
                              ontology = "BP", node_size = 10,
                              algorithm = "elim", statistic = "fisher", top_nodes = 50) {
  
  go_results_list <- list()
  cat("Running GO Enrichment for Direction:", direction_filter, "\n")
  
  for (e in epm_conditions) {
    cat("  Processing EPM:", e, "...")
    # Subset genes of interest for this EPM and direction
    genes_interest_e <- data %>%
      filter(epm == e, direction == direction_filter) %>%
      pull(At_orthologs) %>%
      unique()
    
    # Skip if no genes found for this EPM/direction
    if (length(genes_interest_e) == 0) {
      cat(" No genes found. Skipping.\n")
      next
    }
    cat(" Found", length(genes_interest_e), "genes.")
    
    # Create named factor geneList required by topGO
    geneList <- factor(as.integer(universe %in% genes_interest_e))
    names(geneList) <- universe
    
    # Create topGO data object
    suppressMessages({ # Suppress numerous messages from topGO
      GOdata <- new("topGOdata",
                    description = paste("GO Enrichment for", e, direction_filter),
                    ontology = ontology,
                    allGenes = geneList,
                    geneSel = function(x) x == 1,
                    nodeSize = node_size,
                    mapping = "org.At.tair.db",
                    annot = annFUN.org) # Requires org.At.tair.db
    })
    
    
    # Run enrichment test
    resultFisher <- runTest(GOdata, algorithm = algorithm, statistic = statistic)
    
    # Generate results table
    goEnrichment <- GenTable(GOdata, Fisher = resultFisher, orderBy = "Fisher", topNodes = top_nodes)
    
    # Add EPM identifier and store
    if (nrow(goEnrichment) > 0) {
      goEnrichment$EPM <- e
      go_results_list[[e]] <- goEnrichment
      cat(" Found", nrow(goEnrichment), "enriched terms.\n")
    } else {
      cat (" No significant terms found.\n")
    }
    
  } # End loop over epm_conditions
  
  # Combine results from the list into a single dataframe
  combined_results <- bind_rows(go_results_list)
  cat("Finished GO Enrichment for Direction:", direction_filter, "\n\n")
  return(combined_results)
}

# --- 3. Run GO Enrichment for UP and DO ---

go_results_combined_UP <- run_go_enrichment(merged_data, epm_list, universe_genes, "UP")
go_results_combined_DO <- run_go_enrichment(merged_data, epm_list, universe_genes, "DO")

# --- 4. Save GO Results ---

# Construct dynamic filenames
go_up_csv_file <- file.path(output_dir, paste0("GO_enrichment_", species_tag, "_", condition_tag, "_UP.csv"))
go_do_csv_file <- file.path(output_dir, paste0("GO_enrichment_", species_tag, "_", condition_tag, "_DO.csv"))

# Save results
write.csv(go_results_combined_UP, file = go_up_csv_file, row.names = FALSE)
write.csv(go_results_combined_DO, file = go_do_csv_file, row.names = FALSE)
cat("Saved GO results to:\n", go_up_csv_file, "\n", go_do_csv_file, "\n\n")


# --- 5. Network Building Function (Corrected) ---

build_go_network <- function(go_results, p_val_thresh) {
  cat("Building GO network...\n")
  
  # Ensure input is a tibble for more consistent dplyr behavior
  go_results_tbl <- as_tibble(go_results)
  
  # Ensure Fisher column is numeric and filter
  filtered_go <- go_results_tbl %>%
    # Ensure Term column is character (sometimes factors can cause issues)
    mutate(Term = as.character(Term),
           Fisher = as.numeric(Fisher)) %>%
    filter(Fisher < p_val_thresh) %>%
    # Add small epsilon *before* log10 to prevent log10(0)
    mutate(minusLog10Fisher = -log10(Fisher + .Machine$double.eps))
  
  if(nrow(filtered_go) == 0) {
    cat("  No significant GO terms found below threshold", p_val_thresh, ". Cannot build network.\n")
    return(NULL) # Return NULL if no data
  }
  
  cat("  Found", nrow(filtered_go), "significant term entries.\n")
  
  # Build edge list (terms within the same EPM are connected)
  edge_list <- filtered_go %>%
    # Ensure EPM is character
    mutate(EPM = as.character(EPM)) %>%
    group_by(EPM) %>%
    filter(n() > 1) %>% # Need at least 2 terms in an EPM group to make a connection
    summarise(
      # Use combn to get all pairs of terms within the group
      pairs = list(combn(as.character(Term), 2, simplify = FALSE)),
      # Use mean -log10(p-value) of the group as edge weight
      weight = mean(minusLog10Fisher, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    unnest(pairs) %>%
    # Check if pairs is empty before trying to access elements
    filter(lengths(pairs) == 2) %>%
    mutate(
      Term1 = sapply(pairs, `[`, 1),
      Term2 = sapply(pairs, `[`, 2)
    ) %>%
    dplyr::select(from = Term1, to = Term2, EPM, weight)
  
  if(nrow(edge_list) == 0) {
    cat("  No edges could be formed (no EPM groups with >1 significant term or issue creating pairs).\n")
    # Still attempt to return a graph with nodes but no edges if nodes exist
    # Or return NULL if you prefer no graph at all in this case
    # Let's proceed to node creation and return a node-only graph or NULL
  } else {
    cat("  Built", nrow(edge_list), "edges.\n")
  }
  
  
  # Build node table (unique terms and their annotation count)
  node_df <- filtered_go %>%
    filter(!is.na(Term) & Term != "") %>% # Add check for empty strings too
    # Keep info from the row with the highest annotation count if a term appears multiple times
    group_by(Term) %>%
    arrange(desc(Annotated)) %>%
    # *** FIX APPLIED HERE: Use slice_head instead of slice ***
    slice_head(n = 1) %>%
    ungroup() %>%
    mutate(name = trimws(as.character(Term))) %>% # Ensure 'name' column exists
    dplyr::select(name, Annotated)
  
  # Handle case where node_df might be empty after filtering/slicing
  if(nrow(node_df) == 0 && nrow(edge_list) == 0) {
    cat("  No nodes could be determined. Cannot build network.\n")
    return(NULL)
  } else if (nrow(node_df) == 0 && nrow(edge_list) > 0) {
    # Create nodes from edges if node_df is empty but edges exist
    all_terms_in_edges <- unique(c(edge_list$from, edge_list$to))
    node_df <- tibble(name = all_terms_in_edges, Annotated = 1) # Assign default annotation
    cat("  Created node list from edges as initial node list was empty.\n")
  }
  
  
  # Ensure all nodes from edge list are in node list (important for graph_from_data_frame)
  # This check should happen *after* node_df is potentially created from edges
  if(nrow(edge_list) > 0) {
    all_terms_in_edges <- unique(c(edge_list$from, edge_list$to))
    missing_nodes <- setdiff(all_terms_in_edges, node_df$name)
    
    if (length(missing_nodes) > 0) {
      # Use tibble for consistency
      dummy_nodes <- tibble(name = missing_nodes, Annotated = 1) # Assign min annotation?
      node_df <- bind_rows(node_df, dummy_nodes)
      cat("  Added", length(missing_nodes), "dummy nodes present in edges but not initially in node list.\n")
    }
  }
  
  # Ensure node names are unique (should be due to group_by/slice_head, but double-check)
  node_df <- node_df %>% distinct(name, .keep_all = TRUE)
  
  cat("  Built node table with", nrow(node_df), "unique terms.\n")
  
  # Create igraph object - handle case with edges but no nodes, or nodes but no edges
  if (nrow(node_df) > 0) {
    # If there are no edges, create a graph with only vertices
    if (nrow(edge_list) == 0) {
      cat("  Creating graph with nodes but no edges.\n")
      # graph_from_data_frame requires edges, so create graph from vertices directly
      graph_obj <- make_empty_graph(n = nrow(node_df), directed = FALSE)
      V(graph_obj)$name <- node_df$name
      V(graph_obj)$Annotated <- node_df$Annotated
    } else {
      # Default case: nodes and edges exist
      graph_obj <- graph_from_data_frame(d = edge_list, vertices = node_df, directed = FALSE)
    }
    cat("Finished building GO network.\n\n")
    return(graph_obj)
  } else {
    # This case should ideally not be reached due to earlier checks, but as a fallback:
    cat("  Cannot create graph: No nodes were finalized.\n\n")
    return(NULL)
  }
}

# --- 6. Build Networks for UP and DO ---
# Assuming go_results_combined_UP, go_results_combined_DO, and p_value_threshold exist

# Make sure the input 'Term' column doesn't have factors if that might be an issue
# go_results_combined_UP$Term <- as.character(go_results_combined_UP$Term)
# go_results_combined_DO$Term <- as.character(go_results_combined_DO$Term)

cat("Processing UP regulated genes...\n")
str(go_results_combined_UP)
g_UP <- build_go_network(go_results_combined_UP, p_value_threshold)

cat("\nProcessing DOWN regulated genes...\n")
str(go_results_combined_DO) # Assuming this data frame exists and has the same structure
g_DO <- build_go_network(go_results_combined_DO, p_value_threshold)

# Now you can inspect g_UP and g_DO
# print(g_UP)
# print(g_DO)

# # --- 6. Build Networks for UP and DO ---
# str(go_results_combined_UP)
# g_UP <- build_go_network(go_results_combined_UP, p_value_threshold)
# g_DO <- build_go_network(go_results_combined_DO, p_value_threshold)

# --- 7. Network Visualization ---

# Define color palette based on all EPMs found in *significant* results
unique_epms_in_graphs <- character(0)

# Check UP graph
if (!is.null(g_UP) && "EPM" %in% edge_attr_names(g_UP)) {
  # Corrected: Use E() instead of E_()
  unique_epms_in_graphs <- union(unique_epms_in_graphs, unique(E(g_UP)$EPM))
}

# Check DO graph (apply the same correction here)
if (!is.null(g_DO) && "EPM" %in% edge_attr_names(g_DO)) {
  # Corrected: Use E() instead of E_()
  unique_epms_in_graphs <- union(unique_epms_in_graphs, unique(E(g_DO)$EPM))
}

# Create colors only if there are EPMs to color
custom_colors <- NULL
if (length(unique_epms_in_graphs) > 0) {
  # Use tryCatch for viridis in case only 1 unique EPM exists,viridis might behave differently
  tryCatch({
    custom_colors <- viridis(length(unique_epms_in_graphs))
    names(custom_colors) <- unique_epms_in_graphs
    cat("Defined colors for EPMs:", paste(names(custom_colors), collapse=", "), "\n")
  }, error = function(e) {
    cat("Error generating colors with viridis:", e$message, "\nAttempting fallback.\n")
    # Fallback for single color or issues
    if(length(unique_epms_in_graphs) == 1) {
      custom_colors <- setNames("blue", unique_epms_in_graphs) # Assign a default color
    } else {
      custom_colors <- setNames(rainbow(length(unique_epms_in_graphs)), unique_epms_in_graphs) # Use rainbow as fallback
    }
    names(custom_colors) <- unique_epms_in_graphs
    cat("Using fallback colors for EPMs:", paste(names(custom_colors), collapse=", "), "\n")
  })
  
} else {
  cat("No EPMs found in significant graph edges to define colors for.\n")
}

# Print the custom colors to verify (if generated)
if (!is.null(custom_colors)) {
  print("Custom Colors:")
  print(custom_colors)
}

# --- Redefine the function - Node size by Degree (using igraph::degree) ---

plot_go_network <- function(graph_obj, direction_label, color_palette,
                            sp_tag, cond_tag, out_dir, layout_type = "kk") {
  
  # Skip plotting if graph object is NULL
  if (is.null(graph_obj) || length(V(graph_obj)) == 0) {
    cat("Skipping network plot for", direction_label, "- No graph data.\n\n")
    return(NULL)
  }
  cat("Plotting network test (Node Size by Degree) for:", direction_label, "\n") # Updated message
  
  plot_title <- paste("Network Plot Test (Node Size by Degree) for", direction_label) # Updated title
  # Define filename for saving this version
  filename <- file.path(out_dir, paste0("Network_", sp_tag, "_", cond_tag, "_", direction_label, "_node_size_degree.png")) # Updated filename
  
  # --- Calculate node degree using explicit namespace and add as attribute ---
  if (length(V(graph_obj)) > 0) { # Ensure graph has vertices before calculating degree
    # Use igraph::degree to avoid conflicts
    node_degrees <- igraph::degree(graph_obj, mode = "all")
    graph_obj <- set_vertex_attr(graph_obj, "degree", value = node_degrees)
    cat("Node degrees calculated using igraph::degree and added as 'degree' attribute.\n")
  } else {
    cat("Graph has no vertices, cannot calculate degree.\n")
  }
  
  # --- Plot Call with Node Size mapped to Degree ---
  p <- ggraph(graph_obj, layout = layout_type) +
    # Edges (same as before)
    geom_edge_link(aes(color = EPM), width = 1.2, alpha = 0.5) +
    # Nodes: Map size to the new 'degree' attribute
    geom_node_point(aes(size = degree), color = "darkgrey") +
    # Labels (same as before)
    geom_node_text(aes(label = name), size = 3, color = "black") +
    # Add scale for node size
    scale_size_continuous(name = "Node Degree", range = c(2, 10)) + # Control size range
    theme_graph(base_family = "sans") +
    labs(title = plot_title, edge_color = "EPM") # Removed size legend title here, added in scale
  
  # Add manual color scale for edges (same as before)
  if (!is.null(color_palette) && "EPM" %in% edge_attr_names(graph_obj)) {
    p <- p + scale_edge_color_manual(values = color_palette)
  } else if ("EPM" %in% edge_attr_names(graph_obj)) {
    p <- p + scale_edge_color_viridis_d()
  }
  # --- End of plot call ---
  
  # Display the plot
  tryCatch({
    print(p)
    cat("Plot with node size by degree display attempted.\n")
  }, error = function(e) {
    cat("ERROR displaying plot with node size by degree:", e$message, "\n")
  })
  
  # --- Save this version ---
  tryCatch({
    ggsave(filename, plot = p, width = 12, height = 10, dpi = 300)
    cat("Saved network plot to:", filename, "\n\n")
  }, error = function(e) {
    cat("ERROR saving plot:", e$message, "\n")
  })
  
  return(p) # Return the plot object
}
# --- End of function redefinition ---

cat("Plot_go_network function redefined using igraph::degree.\n") # Confirmation message


dev.off()
# # --- Run this code AFTER creating g_UP ---
# # --- and BEFORE calling plot_go_network ---
# 
# cat("Inspecting g_UP object:\n")
# 
# if (!is.null(g_UP)) {
#   print(g_UP) # Check node and edge counts
#   
#   cat("\nSummary of Node 'Annotated' attribute (for size):\n")
#   if ("Annotated" %in% vertex_attr_names(g_UP)) {
#     print(summary(V(g_UP)$Annotated))
#     cat("Number of NA 'Annotated' values:", sum(is.na(V(g_UP)$Annotated)), "\n")
#   } else {
#     cat("'Annotated' attribute not found on nodes.\n")
#   }
#   
#   cat("\nSummary of Edge 'EPM' attribute (for color):\n")
#   if ("EPM" %in% edge_attr_names(g_UP)) {
#     print(summary(as.factor(E(g_UP)$EPM))) # Use as.factor for summary
#     cat("Number of NA 'EPM' values:", sum(is.na(E(g_UP)$EPM)), "\n")
#   } else {
#     cat("'EPM' attribute not found on edges.\n")
#   }
#   
#   cat("\nSummary of Node 'name' attribute (for labels):\n")
#   if ("name" %in% vertex_attr_names(g_UP)) {
#     # Only print summary if there are nodes, prevent error on empty graph summary
#     if(length(V(g_UP)$name) > 0) {
#       print(summary(V(g_UP)$name))
#       cat("Number of NA 'name' values:", sum(is.na(V(g_UP)$name)), "\n")
#     } else {
#       cat("No nodes to summarize 'name' for.\n")
#     }
#   } else {
#     cat("'name' attribute not found on nodes.\n")
#   }
#   
# } else {
#   cat("g_UP object is NULL.\n")
# }
# cat("--- Inspection finished ---\n")

# --- Now you would normally try calling plot_go_network again ---
p_network_UP <- plot_go_network(g_UP, "UP", custom_colors, species_tag, condition_tag, output_dir)
p_network_DO <- plot_go_network(g_DO, "DO", custom_colors, species_tag, condition_tag, output_dir)
####################################################

# --- 8. Heatmap Data Preparation Function ---

prepare_heatmap_data <- function(go_results, p_val_thresh) {
  cat("Preparing heatmap data...\n")
  # Ensure Fisher column is numeric and filter
  filtered_go <- go_results %>%
    mutate(Fisher = as.numeric(Fisher)) %>%
    filter(Fisher < p_val_thresh) %>%
    mutate(minusLog10Fisher = -log10(Fisher + .Machine$double.eps)) # Add epsilon
  
  if(nrow(filtered_go) == 0) {
    cat("  No significant GO terms found below threshold", p_val_thresh, ". Cannot create heatmap data.\n")
    return(NULL) # Return NULL if no data
  }
  
  # Pivot wider, then longer for ggplot heatmap format
  heatmap_data_long <- filtered_go %>%
    dplyr::select(EPM, Term, minusLog10Fisher) %>%
    # Ensure unique combinations before pivoting wider if necessary
    group_by(EPM, Term) %>%
    summarise(minusLog10Fisher = max(minusLog10Fisher), .groups = "drop") %>% # Take max p-value if term/EPM repeats
    pivot_wider(names_from = Term, values_from = minusLog10Fisher, values_fill = 0) %>%
    pivot_longer(cols = -EPM, names_to = "Term", values_to = "minusLog10Fisher")
  
  cat("Finished preparing heatmap data. Dimensions:", dim(heatmap_data_long), "\n\n")
  return(heatmap_data_long)
}

# --- 9. Prepare Heatmap Data ---

heatmap_data_long_UP <- prepare_heatmap_data(go_results_combined_UP, p_value_threshold)
heatmap_data_long_DO <- prepare_heatmap_data(go_results_combined_DO, p_value_threshold)


# --- 10. Heatmap Plotting Function ---

# --- Function to plot heatmap with dendrograms using pheatmap ---

plot_go_heatmap_pheatmap <- function(filtered_go_results, direction_label, color_vec, # Note: Takes a color vector now
                                     sp_tag, cond_tag, out_dir,
                                     cluster_rows = TRUE, cluster_cols = TRUE,
                                     dist_method = "euclidean", # pheatmap uses 'correlation', 'euclidean', etc.
                                     clust_method = "ward.D2",  # pheatmap uses 'ward.D2', 'complete', etc.
                                     fontsize = 8, angle_col = 45) { # Added appearance args
  
  # --- Basic Checks ---
  if (is.null(filtered_go_results) || nrow(filtered_go_results) == 0) {
    cat("Skipping pheatmap plot for", direction_label, "- No significant GO data provided.\n\n")
    return(NULL)
  }
  cat("Preparing pheatmap with dendrograms for:", direction_label, "\n")
  
  # --- 1. Prepare Wide Matrix (Same as before) ---
  heatmap_matrix_data <- filtered_go_results %>%
    dplyr::select(EPM, Term, minusLog10Fisher) %>%
    group_by(EPM, Term) %>%
    summarise(minusLog10Fisher = max(minusLog10Fisher, na.rm = TRUE), .groups = "drop") %>%
    mutate(minusLog10Fisher = ifelse(is.infinite(minusLog10Fisher), max(.$minusLog10Fisher[!is.infinite(.$minusLog10Fisher)], 0, na.rm=TRUE) + 1, minusLog10Fisher)) %>%
    # Replace any remaining NA/NaN with 0 before converting to matrix
    mutate(minusLog10Fisher = ifelse(is.na(minusLog10Fisher) | is.nan(minusLog10Fisher), 0, minusLog10Fisher)) %>%
    pivot_wider(names_from = Term, values_from = minusLog10Fisher, values_fill = 0)
  
  # Check for sufficient rows/cols for clustering AND plotting
  can_cluster_rows <- cluster_rows && nrow(heatmap_matrix_data) >= 2
  can_cluster_cols <- cluster_cols && (ncol(heatmap_matrix_data) -1) >= 2 # -1 for EPM column
  
  if (nrow(heatmap_matrix_data) == 0 || (ncol(heatmap_matrix_data) -1) == 0) {
    cat("  Error: Matrix has zero rows or zero term columns after preparation. Cannot create heatmap.\n")
    return(NULL)
  }
  
  # Convert to matrix with EPMs as rownames
  heatmap_matrix <- heatmap_matrix_data %>%
    column_to_rownames("EPM") %>%
    as.matrix()
  
  # Check for constant rows/columns which can cause clustering issues
  if(can_cluster_rows && any(apply(heatmap_matrix, 1, function(r) length(unique(r)) == 1))) {
    cat("  Warning: Some rows (EPMs) have constant values, might affect row clustering.\n")
  }
  if(can_cluster_cols && any(apply(heatmap_matrix, 2, function(c) length(unique(c)) == 1))) {
    cat("  Warning: Some columns (Terms) have constant values, might affect column clustering.\n")
  }
  
  
  # --- 2. Define Plot Title and Filename ---
  plot_title <- paste("Clustered Heatmap of GO Enrichment for", direction_label, "Genes")
  filename_png <- file.path(out_dir, paste0("pHeatmap_Clustered_", sp_tag, "_", cond_tag, "_", direction_label, ".png"))
  
  # --- 3. Determine dynamic parameters for pheatmap ---
  # Adjust font sizes based on matrix dimensions
  row_fontsize <- max(5, fontsize - nrow(heatmap_matrix) * 0.02)
  col_fontsize <- max(5, fontsize - ncol(heatmap_matrix) * 0.02)
  # Determine plot dimensions dynamically
  plot_width <- max(7, ncol(heatmap_matrix) * 0.15 + 4) # Base width + term width + dendrogram space
  plot_height <- max(7, nrow(heatmap_matrix) * 0.15 + 4) # Base height + epm height + dendrogram space
  
  # --- 4. Call pheatmap ---
  cat("  Generating pheatmap plot...\n")
  plot_object <- NULL # Initialize variable to store pheatmap object
  tryCatch({
    plot_object <- pheatmap(
      heatmap_matrix,
      color = color_vec, # Use the provided color vector directly
      border_color = "grey80", # Separator lines
      cluster_rows = can_cluster_rows,
      cluster_cols = can_cluster_cols,
      clustering_distance_rows = dist_method,
      clustering_distance_cols = dist_method,
      clustering_method = clust_method, # Applies to both row/col unless specified separately
      main = plot_title,
      fontsize = fontsize,
      fontsize_row = row_fontsize,
      fontsize_col = col_fontsize,
      angle_col = angle_col,
      silent = TRUE, # Suppress direct plotting if we capture object
      # Directly save to file using pheatmap's arguments
      filename = filename_png,
      width = plot_width,
      height = plot_height
    )
    cat("Saved pheatmap plot to:", filename_png, "\n\n")
    
  }, error = function(e) {
    cat("ERROR generating or saving pheatmap:", e$message, "\n")
  })
  
  # pheatmap plots directly or saves to file, but invisibly returns clustering info etc.
  # We don't explicitly print it here as it either plotted to screen or saved to file.
  # If you need to display it interactively AND save, you might remove 'filename' and wrap in png()/dev.off()
  # or capture the object and print it, then save it. Saving directly is usually sufficient.
  
  # Return NULL or potentially the plot_object if needed downstream (though it's complex)
  return(invisible(plot_object))
}

# --- End of function redefinition ---

cat("pheatmap plot function 'plot_go_heatmap_pheatmap' redefined.\n")
# --- 11. Plot Heatmaps ---

# Define color gradients for heatmaps
up_heatmap_colors <- c("white", "lightpink", "deeppink4") # Adjusted gradient
do_heatmap_colors <- c("white", "lightblue", "darkblue")   # Adjusted gradient

# Plot
p_heatmap_UP <- plot_go_heatmap_pheatmap(heatmap_data_long_UP, "UP", up_heatmap_colors, species_tag, condition_tag, output_dir)
p_heatmap_DO <- plot_go_heatmap_pheatmap(heatmap_data_long_DO, "DO", do_heatmap_colors, species_tag, condition_tag, output_dir)

cat("--- Script Finished --- \n")
