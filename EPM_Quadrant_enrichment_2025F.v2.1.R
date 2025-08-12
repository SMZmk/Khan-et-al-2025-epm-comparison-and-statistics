# Load required libraries
library(dplyr)
library(stringr)
library(tidyr)

# --- Configuration & File Loading ---

# Set the working directory to the location of your script/data
# Make sure to adjust this path if necessary.
# setwd("~/Desktop/arab_env_2024i/workflows/")

# Define input file paths
CLUSTER_FILE_PATH <- "../studies/AnemAsag_EPM_set_analyses_K5K19/EPM_cluster_assignments.csv"
file_paths_file_EPM_species_1 <- c("../studies/AnemAsag_EPMclu/2025bAnemSWAnem_AnemE-04_SW2025b_gene_none-q1q9.csv")
file_paths_file_EPM_species_2 <- c("../studies/AnemAsag_EPMclu/AsagSW_E-04_2025d_gene_none-q1q9_SZ.csv")
file_paths_file_DGE_species_1_2 <- c("../studies/AnemAsag_dge_quadrants_khan/combined_species_df_sag_nem_wilt_ctrl_quadrants.csv")
file_paths_file_ortho <- c("../studies/AnemAsag_dge_quadrants_khan/arabis_nem_helixer.v.0.peptide.tsv")
file_paths_file_k18epms <- c("../studies/AnemAsag_EPM_set_analyses_K5K9/EPM_cluster_assignments_k25.csv")

# --- Define Output File Parameters ---
OUTPUT_DIR <- "../studies/AnemAsag_EPMclu/"
OUTPUT_FILENAME <- "AnemAsag_k25_SW_comparative_enrichment_summary.csv"
FULL_OUTPUT_PATH <- file.path(OUTPUT_DIR, OUTPUT_FILENAME)


# --- File Reading & Initial Processing ---
cat("Current working directory:", getwd(), "\n")
cat("Reading input files...\n")
tryCatch({
  EPMs_species_1 <- read.csv(file_paths_file_EPM_species_1)
  EPMs_species_2 <- read.csv(file_paths_file_EPM_species_2)
  DGE_species_1_2 <- read.csv(file_paths_file_DGE_species_1_2)
  ortho_df <- read.csv(file_paths_file_ortho, sep = "\t", header = FALSE) # Assuming no header in tsv
  k18epms_df <- read.csv(file_paths_file_k18epms)
  cat("Files read successfully.\n")
}, error = function(e) {
  stop("Error reading input files: ", e$message)
})

# --- Ortholog Data Cleaning ---
new_column_names <- c("Orthogroup", "Species", "target_gene_id", "As_orthologs")
if (length(new_column_names) == ncol(ortho_df)) {
  colnames(ortho_df) <- new_column_names
} else {
  warning(paste("Ortholog file column count mismatch. Expected", length(new_column_names), "got", ncol(ortho_df)))
}

ortho_long_df <- ortho_df %>%
  filter(Species == "arabis_sag_hlx2024h.pep") %>%
  filter(!is.na(target_gene_id) & target_gene_id != "",
         !is.na(As_orthologs) & As_orthologs != "") %>%
  separate_rows(target_gene_id, sep = ",\\s*") %>%
  separate_rows(As_orthologs, sep = ",\\s*") %>%
  dplyr::select(Orthogroup, target_gene_id, As_orthologs)

ortho_cleaned_df <- ortho_long_df %>%
  dplyr::select(target_gene_id, As_orthologs) %>%
  mutate(across(c(target_gene_id, As_orthologs), ~ str_remove(.x, "\\.\\d+$"))) %>%
  filter(target_gene_id != "", As_orthologs != "") %>%
  distinct()
cat("Cleaned ortholog mapping prepared.\n")

# --- Data Merging and Preparation ---
k18epms_df$epm <- substr(k18epms_df$epm, 1, 17)
k18epms_df <- unique(k18epms_df)

DGE_species_1_2_small <- DGE_species_1_2 %>%
  dplyr::select(gene_id, quadrant)

ortho_DGEs_df <- merge(ortho_cleaned_df, DGE_species_1_2_small,
                       by.x = "target_gene_id",
                       by.y = "gene_id")

EPMs_species_1_small <- EPMs_species_1 %>% dplyr::select(loc_ID, epm)
colnames(EPMs_species_1_small) <- c("loc_ID", "epm_species_1")

EPMs_species_2_small <- EPMs_species_2 %>% dplyr::select(loc_ID, epm)
colnames(EPMs_species_2_small) <- c("loc_ID", "epm_species_2")

ortho_DGEs_epm_spec1df <- merge(ortho_DGEs_df %>% dplyr::select(target_gene_id, quadrant), EPMs_species_1_small,
                                by.x = "target_gene_id", by.y = "loc_ID")
ortho_DGEs_epm_spec2df <- merge(ortho_DGEs_df %>% dplyr::select(As_orthologs, quadrant), EPMs_species_2_small,
                                by.x = "As_orthologs", by.y = "loc_ID")

# --- Assign Clusters to Gene/EPM data ---
ortho_DGEs_epmk25_spec1df <- merge(ortho_DGEs_epm_spec1df, k18epms_df,
                                   by.x = "epm_species_1", by.y = "epm" )
ortho_DGEs_epmk25_spec2df <- merge(ortho_DGEs_epm_spec2df, k18epms_df,
                                   by.x = "epm_species_2", by.y = "epm" )

cat("Preparation of data files: done. Starting tests for enrichments.\n")

###################################################################################
### REWRITTEN SECTION: Parallel Enrichment Tests (Anemorensis - Species 1)      ###
###################################################################################
cat("\nRunning Parallel Enrichment Tests for A. nemorensis (Species 1) Clusters...\n")

df_spec1 <- ortho_DGEs_epmk25_spec1df

results_list_k25_spec1 <- list()

# Get unique clusters and quadrants
all_clusters_spec1 <- unique(df_spec1$cluster)
all_quadrants_spec1 <- unique(df_spec1$quadrant)

for (cluster_id in all_clusters_spec1) {
  for (quadrant_id in all_quadrants_spec1) {
    
    # --- Test 1: Gene-Centric Enrichment (based on unique genes) ---
    genes_in_cluster <- unique(df_spec1$target_gene_id[df_spec1$cluster == cluster_id])
    genes_in_quadrant <- unique(df_spec1$target_gene_id[df_spec1$quadrant == quadrant_id])
    total_unique_genes <- n_distinct(df_spec1$target_gene_id)
    
    a1 <- length(intersect(genes_in_cluster, genes_in_quadrant))
    b1 <- length(genes_in_cluster) - a1
    c1 <- length(genes_in_quadrant) - a1
    d1 <- total_unique_genes - (a1 + b1 + c1)
    
    fisher_gene <- fisher.test(matrix(c(a1, c1, b1, d1), nrow = 2), alternative = "greater")
    
    # --- Test 2: EPM-Centric Enrichment (based on row counts) ---
    epms_in_cluster <- sum(df_spec1$cluster == cluster_id)
    epms_in_quadrant <- sum(df_spec1$quadrant == quadrant_id)
    total_epms <- nrow(df_spec1)
    
    a2 <- sum(df_spec1$cluster == cluster_id & df_spec1$quadrant == quadrant_id)
    b2 <- epms_in_cluster - a2
    c2 <- epms_in_quadrant - a2
    d2 <- total_epms - (a2 + b2 + c2)
    
    fisher_epm <- fisher.test(matrix(c(a2, c2, b2, d2), nrow = 2), alternative = "greater")
    
    # --- Collate Results ---
    results_list_k25_spec1[[paste(cluster_id, quadrant_id, sep="_")]] <- data.frame(
      Cluster = cluster_id, Quadrant = quadrant_id,
      GeneCentric_PValue = fisher_gene$p.value,
      GeneCentric_Abundance = a1,
      Total_Unique_Genes_in_Cluster = length(genes_in_cluster),
      EPMCentric_PValue = fisher_epm$p.value,
      EPMCentric_Abundance = a2
    )
  }
}
enrichment_results_k25_spec1 <- do.call(rbind, results_list_k25_spec1)
if(!is.null(enrichment_results_k25_spec1)) rownames(enrichment_results_k25_spec1) <- NULL
cat("A. nemorensis enrichment analysis complete.\n")


###################################################################################
### REWRITTEN SECTION: Parallel Enrichment Tests (Asagitatta - Species 2)       ###
###################################################################################
cat("\nRunning Parallel Enrichment Tests for A. sagitatta (Species 2) Clusters...\n")

df_spec2 <- ortho_DGEs_epmk25_spec2df
results_list_k25_spec2 <- list()

all_clusters_spec2 <- unique(df_spec2$cluster)
all_quadrants_spec2 <- unique(df_spec2$quadrant)

for (cluster_id in all_clusters_spec2) {
  for (quadrant_id in all_quadrants_spec2) {
    
    genes_in_cluster <- unique(df_spec2$As_orthologs[df_spec2$cluster == cluster_id])
    genes_in_quadrant <- unique(df_spec2$As_orthologs[df_spec2$quadrant == quadrant_id])
    total_unique_genes <- n_distinct(df_spec2$As_orthologs)
    
    a1 <- length(intersect(genes_in_cluster, genes_in_quadrant))
    b1 <- length(genes_in_cluster) - a1
    c1 <- length(genes_in_quadrant) - a1
    d1 <- total_unique_genes - (a1 + b1 + c1)
    
    fisher_gene <- fisher.test(matrix(c(a1, c1, b1, d1), nrow = 2), alternative = "greater")
    
    epms_in_cluster <- sum(df_spec2$cluster == cluster_id)
    epms_in_quadrant <- sum(df_spec2$quadrant == quadrant_id)
    total_epms <- nrow(df_spec2)
    
    a2 <- sum(df_spec2$cluster == cluster_id & df_spec2$quadrant == quadrant_id)
    b2 <- epms_in_cluster - a2
    c2 <- epms_in_quadrant - a2
    d2 <- total_epms - (a2 + b2 + c2)
    
    fisher_epm <- fisher.test(matrix(c(a2, c2, b2, d2), nrow = 2), alternative = "greater")
    
    results_list_k25_spec2[[paste(cluster_id, quadrant_id, sep="_")]] <- data.frame(
      Cluster = cluster_id, Quadrant = quadrant_id,
      GeneCentric_PValue = fisher_gene$p.value,
      GeneCentric_Abundance = a1,
      Total_Unique_Genes_in_Cluster = length(genes_in_cluster),
      EPMCentric_PValue = fisher_epm$p.value,
      EPMCentric_Abundance = a2
    )
  }
}
enrichment_results_k25_spec2 <- do.call(rbind, results_list_k25_spec2)
if(!is.null(enrichment_results_k25_spec2)) rownames(enrichment_results_k25_spec2) <- NULL
cat("A. sagitatta enrichment analysis complete.\n")


###################################################################################
### FINAL REVISED SECTION: Create Comparative Human-Readable Summary Table      ###
###################################################################################
cat("\nGenerating comparative summary table for significantly enriched clusters...\n")

P_VALUE_THRESHOLD <- 0.05

# Function to process results for one species
process_species_results <- function(results_df, species_prefix) {
  if (is.null(results_df) || nrow(results_df) == 0) return(data.frame())
  
  # Filter by either p-value being significant
  summary_filtered <- results_df %>% 
    filter(GeneCentric_PValue < P_VALUE_THRESHOLD | EPMCentric_PValue < P_VALUE_THRESHOLD)
  
  if(nrow(summary_filtered) == 0) return(data.frame())
  
  summary <- summary_filtered %>%
    mutate(
      Abundance_in_Percent = (GeneCentric_Abundance / Total_Unique_Genes_in_Cluster) * 100,
      EPMs_per_Gene = if_else(GeneCentric_Abundance > 0, EPMCentric_Abundance / GeneCentric_Abundance, 0)
    ) %>%
    # Select and rename, formatting at the end
    dplyr::select(
      Quadrant,
      Cluster,
      !!paste0(species_prefix, "_GeneCentric_PValue") := GeneCentric_PValue,
      !!paste0(species_prefix, "_GeneCentric_Abundance") := GeneCentric_Abundance,
      !!paste0(species_prefix, "_EPMCentric_PValue") := EPMCentric_PValue,
      !!paste0(species_prefix, "_EPMCentric_Abundance") := EPMCentric_Abundance,
      !!paste0(species_prefix, "_Abundance_in_Percent") := Abundance_in_Percent,
      !!paste0(species_prefix, "_EPMs_per_Gene") := EPMs_per_Gene
    )
  return(summary)
}

# 1. Process both species' results
anem_summary <- process_species_results(enrichment_results_k25_spec1, "Anemorensis")
asag_summary <- process_species_results(enrichment_results_k25_spec2, "Asagitatta")

# 2. Check for any significant results
if (nrow(anem_summary) == 0 && nrow(asag_summary) == 0) {
  cat("\nNo significant cluster enrichments found (p < ", P_VALUE_THRESHOLD, ") to generate a summary table.\n")
} else {
  # 3. Calculate unique genes per quadrant and rename for join
  genes_per_quadrant <- DGE_species_1_2_small %>%
    group_by(quadrant) %>%
    summarise(genes_per_quadrant = n_distinct(gene_id), .groups = 'drop') %>%
    rename(Quadrant = quadrant)
  
  # 4. Merge results. This join is complex because a cluster can be significant in one species but not the other.
  # We must join by both Quadrant and Cluster to keep the rows aligned correctly.
  final_summary_table <- full_join(anem_summary, asag_summary, by = c("Quadrant", "Cluster"), relationship = "many-to-many") %>%
    left_join(genes_per_quadrant, by = "Quadrant") %>%
    dplyr::select(
      quadrant = Quadrant,
      cluster = Cluster,
      genes_per_quadrant,
      starts_with("Anemorensis_"),
      starts_with("Asagitatta_")
    ) %>%
    arrange(quadrant, as.numeric(cluster))
  
  # 5. Prepare and Print Final Table
  # Formatting numbers into strings for the final print
  final_summary_for_printing <- final_summary_table %>%
    mutate(
      across(ends_with("_PValue"), ~sprintf("%.2e", .)),
      across(ends_with("_in_Percent"), ~sprintf("%.2f%%", .)),
      across(ends_with("_EPMs_per_Gene"), ~sprintf("%.2f", .))
    ) %>%
    mutate(across(everything(), as.character)) %>%
    mutate(across(everything(), ~replace_na(., "N/A")))
  
  cat("\n\n--- Comparative Summary of Enriched EPM Clusters (p < ", P_VALUE_THRESHOLD, ") ---\n\n")
  print(as.data.frame(final_summary_for_printing), row.names = FALSE, right = FALSE)
  
  # Save the results to the specified file path
  write.csv(final_summary_for_printing, FULL_OUTPUT_PATH, row.names = FALSE)
  cat("\n\nSummary table has been saved to:", FULL_OUTPUT_PATH, "\n")
}
