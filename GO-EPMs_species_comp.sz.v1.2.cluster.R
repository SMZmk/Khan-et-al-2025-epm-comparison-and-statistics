# --- START: Modified Script for Wilting GO/Cluster Association (Inclusive) ---

# Record script start time (optional)
start_time <- Sys.time()
# Get current date for potential filename stamping
current_date_str <- format(Sys.Date(), "%Y%m%d") # Uses current date: e.g., 20250427

################################################################################
#                       Configuration                                        #
################################################################################
cat("--- Setting Configuration ---\n")

# --- Species, Conditions, and Cluster K value ---
SPEC_1 <- "Anem"
SPEC_2 <- "Asag"
COND_1 <- "wilting"
COND_2 <- "recovery" # Read for Term mapping, but logic focuses on Wilting
K_VALUE <- 25       # The 'k' value used in the cluster enrichment files (e.g., k25)

# --- Column Names (CRITICAL: Ensure these match your CSV headers) ---
GO_COL <- "GO.ID"      # GO ID column
TERM_COL <- "Term"     # GO Term description column
# IMPORTANT: "Cluster" column in NEW files replaces EPM for association
ASSOC_ID_COL <- "Cluster" # Column name for the Associated ID (now Cluster)

# --- File Paths & Working Directory ---
# Set the working directory if needed (ensure this path is correct relative to studies)
# setwd("/path/to/your/working/directory")
cat("Current working directory:", getwd(), "\n")

# Define GO enrichment file paths systematically based on new naming convention
studies_dir <- file.path("..", "studies") # Relative path to studies folder
go_file_paths <- list(
  S1WUP = file.path(studies_dir, paste0("GO_cluster_enrichment_significant_", SPEC_1, "_", COND_1, "_UP_k", K_VALUE, ".csv")),
  S1WDO = file.path(studies_dir, paste0("GO_cluster_enrichment_significant_", SPEC_1, "_", COND_1, "_DO_k", K_VALUE, ".csv")),
  S2WUP = file.path(studies_dir, paste0("GO_cluster_enrichment_significant_", SPEC_2, "_", COND_1, "_UP_k", K_VALUE, ".csv")),
  S2WDO = file.path(studies_dir, paste0("GO_cluster_enrichment_significant_", SPEC_2, "_", COND_1, "_DO_k", K_VALUE, ".csv")),
  # Recovery files read primarily for comprehensive Term mapping
  S1RUP = file.path(studies_dir, paste0("GO_cluster_enrichment_significant_", SPEC_1, "_", COND_2, "_UP_k", K_VALUE, ".csv")),
  S1RDO = file.path(studies_dir, paste0("GO_cluster_enrichment_significant_", SPEC_1, "_", COND_2, "_DO_k", K_VALUE, ".csv")),
  S2RUP = file.path(studies_dir, paste0("GO_cluster_enrichment_significant_", SPEC_2, "_", COND_2, "_UP_k", K_VALUE, ".csv")),
  S2RDO = file.path(studies_dir, paste0("GO_cluster_enrichment_significant_", SPEC_2, "_", COND_2, "_DO_k", K_VALUE, ".csv"))
)

# --- Output Directory ---
output_dir <- file.path("..", "studies") # Relative path for output
if (!dir.exists(output_dir)) {
  cat("INFO: Output directory '", output_dir, "' not found. Creating it.\n")
  dir.create(output_dir, recursive = TRUE, showWarnings = TRUE)
}
cat("Output directory set to:", output_dir, "\n")

# --- Load Required Packages ---
cat("Loading required packages (dplyr, tidyr, rlang)...\n")
packages_needed <- c("dplyr", "tidyr", "rlang")
for (pkg in packages_needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("FATAL ERROR: Package '", pkg, "' is needed. Please install it."), call. = FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
cat("Packages loaded.\n")

################################################################################
#                       File Reading                                         #
################################################################################
cat("\n--- Reading Input Files ---\n")

# --- Function to read GO files safely ---
read_go_file_safe <- function(path, file_desc) {
  cat("Reading:", file_desc, "from", basename(path), "\n")
  if (!file.exists(path)) {
    warning("File not found for ", file_desc, ": ", path, ". Returning NULL.")
    return(NULL)
  }
  tryCatch({
    # Assuming the new files are also comma-separated with header
    df <- read.csv(path, stringsAsFactors = FALSE, sep = ",", header = TRUE)
    # Basic check for expected columns (using the new names)
    expected_cols <- c(GO_COL, ASSOC_ID_COL, TERM_COL, "Annotated", "Significant", "Expected", "Fisher", "minusLog10Fisher") # Add other expected cols
    missing_cols <- setdiff(expected_cols, names(df))
    if (length(missing_cols) > 0) {
      warning("Missing one or more expected columns (", paste(missing_cols, collapse=", "), ") in ", basename(path))
    }
    # Report empty or zero-row files
    if (nrow(df) == 0) {
      warning("File ", basename(path), " is empty (zero rows).")
    }
    return(df)
  }, error = function(e) {
    warning("Error reading ", file_desc, " file '", basename(path), "': ", e$message, ". Returning NULL.")
    return(NULL)
  })
}

# Read GO files
s1wup_df <- read_go_file_safe(go_file_paths$S1WUP, paste(SPEC_1, COND_1, "UP"))
s1wdo_df <- read_go_file_safe(go_file_paths$S1WDO, paste(SPEC_1, COND_1, "DO"))
s2wup_df <- read_go_file_safe(go_file_paths$S2WUP, paste(SPEC_2, COND_1, "UP"))
s2wdo_df <- read_go_file_safe(go_file_paths$S2WDO, paste(SPEC_2, COND_1, "DO"))
# Read recovery for term mapping
s1rup_df <- read_go_file_safe(go_file_paths$S1RUP, paste(SPEC_1, COND_2, "UP"))
s1rdo_df <- read_go_file_safe(go_file_paths$S1RDO, paste(SPEC_1, COND_2, "DO"))
s2rup_df <- read_go_file_safe(go_file_paths$S2RUP, paste(SPEC_2, COND_2, "UP"))
s2rdo_df <- read_go_file_safe(go_file_paths$S2RDO, paste(SPEC_2, COND_2, "DO"))

# List of *all* DFs for Term mapping later
all_dfs_for_terms <- list(s1wup=s1wup_df, s1wdo=s1wdo_df, s2wup=s2wup_df, s2wdo=s2wdo_df,
                          s1rup=s1rup_df, s1rdo=s1rdo_df, s2rup=s2rup_df, s2rdo=s2rdo_df)

# Check essential Wilting DFs were loaded (even if empty, function returns df)
essential_wilting_dfs <- list(s1wup_df, s1wdo_df, s2wup_df, s2wdo_df)
if (any(sapply(essential_wilting_dfs, is.null))) {
  stop("FATAL ERROR: One or more essential Wilting GO dataframes failed to load entirely (returned NULL). Cannot proceed.")
}
if (all(sapply(essential_wilting_dfs, function(df) is.null(df) || nrow(df) == 0))) {
  warning("All Wilting input dataframes are empty or failed to load. Output will be empty.")
  # Allow script to continue to potentially generate an empty file or header.
}

# --- Removed: Cluster assignment file reading and processing is no longer needed ---


################################################################################
#          Identify *ALL* Significant GO IDs & Tag Occurrence                #
################################################################################
cat("\n--- Identifying ALL Significant GO IDs from Wilting and Tagging Occurrence ---\n")

# --- Helper function to get valid, unique GO IDs from a dataframe ---
get_valid_go_ids_from_df <- function(df, go_col_name) {
  if (is.null(df) || nrow(df) == 0 || !go_col_name %in% names(df)) {
    return(character(0)) # Return empty character vector if df is bad
  }
  unique(na.omit(df[[go_col_name]][df[[go_col_name]] != ""]))
}

# Get GO ID lists for each Wilting condition
go_s1wup <- get_valid_go_ids_from_df(s1wup_df, GO_COL)
go_s1wdo <- get_valid_go_ids_from_df(s1wdo_df, GO_COL)
go_s2wup <- get_valid_go_ids_from_df(s2wup_df, GO_COL)
go_s2wdo <- get_valid_go_ids_from_df(s2wdo_df, GO_COL)

# --- *** CORRECTED LOGIC: Include ALL unique GO IDs from ANY Wilting file *** ---
all_wilting_gos <- unique(c(go_s1wup, go_s1wdo, go_s2wup, go_s2wdo))
cat("Found", length(all_wilting_gos), "unique GO IDs significant in AT LEAST ONE Wilting condition.\n")

# Create the base data frame for *all* these GO IDs
if (length(all_wilting_gos) > 0) {
  go_base_df <- tibble(!!GO_COL := all_wilting_gos)
  
  # --- Tag GO IDs based on presence in each input file ---
  cat("Tagging GO IDs with their presence in input files...\n")
  go_base_df <- go_base_df %>%
    mutate(
      In_Anem_Wilt_UP = .data[[GO_COL]] %in% go_s1wup,
      In_Anem_Wilt_DO = .data[[GO_COL]] %in% go_s1wdo,
      In_Asag_Wilt_UP = .data[[GO_COL]] %in% go_s2wup,
      In_Asag_Wilt_DO = .data[[GO_COL]] %in% go_s2wdo
    )
  
  # --- Determine the Wilting Status (more comprehensive) ---
  cat("Determining Wilting status for each GO ID...\n")
  go_base_df <- go_base_df %>%
    mutate(
      Wilting_Status = case_when(
        # Start with intersection patterns
        In_Anem_Wilt_UP & In_Asag_Wilt_UP      ~ "UP_Both", # Make sure spaces here are normal
        In_Anem_Wilt_DO & In_Asag_Wilt_DO      ~ "DOWN_Both", # Make sure spaces here are normal
        In_Anem_Wilt_UP & In_Asag_Wilt_DO      ~ "UP_Anem/DOWN_Asag", # Make sure spaces here are normal
        In_Anem_Wilt_DO & In_Asag_Wilt_UP      ~ "DOWN_Anem/UP_Asag", # Make sure spaces here are normal
        # Then single species patterns (ensure no overlap with above)
        In_Anem_Wilt_UP & !In_Asag_Wilt_UP & !In_Asag_Wilt_DO ~ "UP_Anem_Only", # Make sure spaces here are normal
        In_Anem_Wilt_DO & !In_Asag_Wilt_UP & !In_Asag_Wilt_DO ~ "DOWN_Anem_Only", # Make sure spaces here are normal
        In_Asag_Wilt_UP & !In_Anem_Wilt_UP & !In_Anem_Wilt_DO ~ "UP_Asag_Only", # Make sure spaces here are normal
        In_Asag_Wilt_DO & !In_Anem_Wilt_UP & !In_Anem_Wilt_DO ~ "DOWN_Asag_Only", # Make sure spaces here are normal
        # Handle cases where one species has UP/DOWN, other has none/one
        In_Anem_Wilt_UP & In_Anem_Wilt_DO & !In_Asag_Wilt_UP & !In_Asag_Wilt_DO ~ "UP&DO_Anem_Only", # Make sure spaces here are normal
        !In_Anem_Wilt_UP & !In_Anem_Wilt_DO & In_Asag_Wilt_UP & In_Asag_Wilt_DO ~ "UP&DO_Asag_Only", # Make sure spaces here are normal
        # Add more complex overlaps if needed...
        TRUE                                   ~ "Complex/Other" # Fallback case # Make sure spaces here are normal
      )
    )
  
} else {
  warning("No GO IDs found significant in any Wilting condition. The final output will be empty.")
  go_base_df <- tibble(!!GO_COL := character()) # Create empty tibble with correct col name
}

################################################################################
#      Consolidate Associated Cluster IDs for ALL Wilting GOs                #
################################################################################
cat("\n--- Consolidating Associated Cluster IDs for ALL Wilting GO IDs ---\n")

cluster_aggregated <- NULL # Initialize

# Use the comprehensive list 'all_wilting_gos' here
if (length(all_wilting_gos) > 0) {
  
  # --- Create a single long-format table of GO-Cluster ID pairs from Wilting files ---
  select_and_filter_assoc <- function(df, go_list, go_c, assoc_c) {
    if (is.null(df) || nrow(df) == 0 || !all(c(go_c, assoc_c) %in% names(df))) return(NULL)
    df %>%
      filter(.data[[go_c]] %in% go_list) %>% # Filter using the full list
      dplyr::select(all_of(c(go_c, assoc_c))) %>%
      # Filter out NA or empty cluster IDs
      filter(!is.na(.data[[assoc_c]]) & .data[[assoc_c]] != "" & .data[[assoc_c]] != "NA") %>%
      dplyr::rename(GO.ID = !!sym(go_c), Associated_Cluster_ID = !!sym(assoc_c)) # Standardize names
  }
  
  list_of_go_assoc_dfs <- list(
    select_and_filter_assoc(s1wup_df, all_wilting_gos, GO_COL, ASSOC_ID_COL),
    select_and_filter_assoc(s1wdo_df, all_wilting_gos, GO_COL, ASSOC_ID_COL),
    select_and_filter_assoc(s2wup_df, all_wilting_gos, GO_COL, ASSOC_ID_COL),
    select_and_filter_assoc(s2wdo_df, all_wilting_gos, GO_COL, ASSOC_ID_COL)
  )
  
  # Combine all valid GO-Cluster ID pairs into one table and get unique pairs
  all_wilting_go_clusters <- bind_rows(list_of_go_assoc_dfs) %>%
    distinct(GO.ID, Associated_Cluster_ID)
  
  cat("Found", nrow(all_wilting_go_clusters), "unique GO.ID-Cluster_ID pairs across Wilting conditions for the", length(all_wilting_gos), "GOs.\n")
  
  # --- Aggregate Cluster IDs per GO ID ---
  cat("Aggregating unique Cluster IDs per GO ID...\n")
  if (nrow(all_wilting_go_clusters) > 0) {
    cluster_aggregated <- all_wilting_go_clusters %>%
      group_by(GO.ID) %>%
      summarise(
        # Consolidate unique Cluster IDs
        Associated_Clusters = paste(sort(unique(Associated_Cluster_ID)), collapse = "; "),
        .groups = 'drop' # Drop grouping after summarising
      )
    cat("Aggregation complete.\n")
  } else {
    cat("No GO.ID-Cluster_ID pairs found. Aggregation skipped.\n")
    # Create empty tibble with expected columns
    cluster_aggregated <- tibble(
      GO.ID = character(),
      Associated_Clusters = character()
    )
  }
  
  
} else {
  cat("Skipping Cluster consolidation and linking as no significant GO IDs were found in any wilting condition.\n")
  # Ensure cluster_aggregated is NULL or an empty tibble with expected columns
  cluster_aggregated <- tibble(
    GO.ID = character(),
    Associated_Clusters = character()
  ) %>% dplyr::rename(!!GO_COL := GO.ID) # Ensure correct GO col name
}


################################################################################
#                   Combine All Information & Finalize                       #
################################################################################
cat("\n--- Combining GO Tags, Status, and Aggregated Cluster Info ---\n")

# Start with the base GO table (GOs, Tags, Status)
final_wilting_table <- go_base_df

# Join the aggregated Cluster information
if (!is.null(cluster_aggregated) && nrow(cluster_aggregated) > 0) {
  # Ensure the column name matches GO_COL before joining
  if (!GO_COL %in% names(cluster_aggregated)) {
    # Attempt to rename if 'GO.ID' exists (from internal summarisation)
    if ("GO.ID" %in% names(cluster_aggregated)) {
      cluster_aggregated <- cluster_aggregated %>% rename(!!GO_COL := GO.ID)
    } else {
      stop("Cannot find join column '", GO_COL, "' or 'GO.ID' in aggregated Cluster data.")
    }
  }
  final_wilting_table <- final_wilting_table %>%
    left_join(cluster_aggregated, by = GO_COL)
} else {
  # Add empty columns if no Cluster info was generated
  cat("No Cluster information to join. Adding empty columns.\n")
  final_wilting_table <- final_wilting_table %>%
    mutate(
      Associated_Clusters = ""
    )
}

# --- Add Term descriptions ---
cat("Adding GO Term descriptions...\n")
go_term_map <- NULL
term_map_list <- list()
# Consolidate terms from *all* loaded dataframes
for(df_name in names(all_dfs_for_terms)) {
  df <- all_dfs_for_terms[[df_name]]
  # Check if the required columns exist in the dataframe before processing
  if (!is.null(df) && nrow(df) > 0 && all(c(GO_COL, TERM_COL) %in% names(df))) {
    term_map_list[[length(term_map_list) + 1]] <- df %>%
      dplyr::select(all_of(c(GO_COL, TERM_COL))) %>%
      filter(!is.na(.data[[GO_COL]]) & .data[[GO_COL]] != "" &
               !is.na(.data[[TERM_COL]]) & .data[[TERM_COL]] != "") %>%
      distinct(!!sym(GO_COL), .keep_all = TRUE) # Keep first term found per GO ID
  } else if (!is.null(df) && nrow(df) > 0) {
    warning(paste0("Term mapping skipped for dataframe '", df_name, "' due to missing GO_COL or TERM_COL."))
  }
}

if (length(term_map_list) > 0) {
  go_term_map <- bind_rows(term_map_list) %>%
    distinct(!!sym(GO_COL), .keep_all = TRUE) # Ensure overall uniqueness by GO ID
  cat("Created Term map with", nrow(go_term_map), "unique GO ID-Term pairs.\n")
  
  # Join terms to the final table
  final_wilting_table <- final_wilting_table %>%
    left_join(go_term_map, by = GO_COL) %>%
    # Bring important columns to the front
    dplyr::select(
      all_of(GO_COL),
      all_of(TERM_COL),
      Wilting_Status, # Use the corrected status column
      starts_with("In_"), # Bring tag columns next
      starts_with("Associated_"), # Then Cluster info (now only Associated_Clusters)
      everything() # Keep any other columns (like Annotated, Fisher, etc. from original input)
    )
  # Handle cases where a GO ID didn't have a term found
  final_wilting_table[[TERM_COL]] <- ifelse(is.na(final_wilting_table[[TERM_COL]]), "<Term not found>", final_wilting_table[[TERM_COL]])
  cat("Joined Term descriptions.\n")
} else {
  warning("Could not build a Term map from any input files.")
  final_wilting_table[[TERM_COL]] <- "<Term map unavailable>"
  # Reorder columns even if term map failed
  final_wilting_table <- final_wilting_table %>%
    dplyr::select(
      all_of(GO_COL),
      all_of(TERM_COL),
      Wilting_Status,
      starts_with("In_"),
      starts_with("Associated_"),
      everything()
    )
}


# --- Final Cleanup ---
# Replace NA values in aggregation columns (can happen if a GO ID had no associated cluster)
cols_to_clean_na <- c("Associated_Clusters") # Only need to clean this one now
# Ensure columns exist before cleaning NAs
cols_to_clean_na <- intersect(cols_to_clean_na, names(final_wilting_table))
if(length(cols_to_clean_na) > 0) {
  final_wilting_table <- final_wilting_table %>%
    mutate(across(all_of(cols_to_clean_na), ~replace_na(., "")))
}


# Display head of the final table
if (nrow(final_wilting_table) > 0) {
  cat("\n--- Head of the Final Wilting Table (Inclusive) ---\n")
  print(head(final_wilting_table))
} else {
  cat("\n--- Final Wilting Table is Empty ---\n")
}


################################################################################
#                     Export Final Table to CSV File                         #
################################################################################
cat("\n--- Exporting Final Wilting Table (Inclusive) to CSV File ---\n")
final_filename <- file.path(output_dir, paste0("Wilting_GO_Cluster_Associations_Inclusive_k", K_VALUE, "_", current_date_str, ".csv")) # Added k-value to filename

# Check if the table has rows before writing
if (nrow(final_wilting_table) > 0) {
  tryCatch({
    write.csv(final_wilting_table, final_filename, row.names = FALSE, na = "") # Write NAs as empty strings
    cat("Successfully saved final table to:", final_filename, "\n")
  }, error = function(e) {
    warning("ERROR saving final table file: ", e$message)
  })
} else {
  cat("Skipping export because the final table is empty.\n")
}

# --- Script End ---
end_time <- Sys.time()
cat("\n--- Analysis and Export Finished ---\n")
cat("Total script execution time:", format(end_time - start_time), "\n")

# --- END: Modified Script ---
# # Load required libraries for network visualization
# library(dplyr)
# library(tidyr)      # Make sure tidyr is loaded for separate_rows()
# library(igraph)
# library(tidygraph)
# library(ggraph)
# library(RColorBrewer) # For color palettes
# library(ggplot2)      # ggsave is part of ggplot2
# library(ggrepel)      # For geom_node_text_repel (relies on this)
# 
# cat("\n--- Generating Network Visualization with Labels ---\n")
# 
# if (nrow(final_wilting_table) > 0 && any(final_wilting_table$Associated_Clusters != "")) {
#   
#   # --- 1. Prepare the edge list (GO-Cluster pairs) ---
#   # separate_rows is used here
#   go_cluster_pairs <- final_wilting_table %>%
#     filter(Associated_Clusters != "") %>% # Only include GOs with cluster associations
#     separate_rows(Associated_Clusters, sep = "; ") %>% # separate_rows function call
#     dplyr::rename(Cluster = Associated_Clusters) %>%
#     # filter(Cluster != "cNA") %>% # Uncomment if you don't want cNA in network
#     distinct(GO.ID, Cluster, Wilting_Status, Term) # Keep relevant info
#   
#   if(nrow(go_cluster_pairs) > 0) {
#     
#     # --- 2. Prepare Node Data ---
#     go_nodes <- go_cluster_pairs %>%
#       distinct(GO.ID, Wilting_Status, Term) %>%
#       dplyr::rename(name = GO.ID) %>%
#       mutate(type = "GO")
#     
#     cluster_nodes <- go_cluster_pairs %>%
#       distinct(Cluster) %>%
#       dplyr::rename(name = Cluster) %>%
#       mutate(type = "Cluster", Wilting_Status = NA, Term = NA) # Add NA for GO-specific attributes
#     
#     all_nodes <- bind_rows(go_nodes, cluster_nodes) %>%
#       distinct(name, .keep_all = TRUE)
#     
#     # --- 3. Prepare Edge Data ---
#     edge_list <- go_cluster_pairs %>%
#       dplyr::select(from = GO.ID, to = Cluster) %>%
#       mutate(cluster_for_edge_color = to)
#     
#     
#     # --- 4. Create the tidygraph object ---
#     graph <- tbl_graph(nodes = all_nodes, edges = edge_list, directed = FALSE)
#     
#     
#     # --- 5. Visualize the network using ggraph ---
#     cat("Drawing network plot...\n")
#     
#     # Define a color palette for Wilting Status (for GO nodes)
#     actual_statuses <- unique(go_nodes$Wilting_Status)
#     status_palette_values <- RColorBrewer::brewer.pal(max(3, n_distinct(actual_statuses)), "Set3")[1:n_distinct(actual_statuses)]
#     status_colors_map <- setNames(status_palette_values, actual_statuses)
#     
#     # Define a color palette for Clusters (for edges and Cluster nodes)
#     actual_clusters <- unique(cluster_nodes$name)
#     cluster_palette_values <- scales::hue_pal()(n_distinct(actual_clusters)) # Using scales::hue_pal()
#     cluster_colors_map <- setNames(cluster_palette_values, actual_clusters)
#     
#     # Add a specific color for cNA if included and is in the actual clusters
#     if ("cNA" %in% names(cluster_colors_map)) {
#       cluster_colors_map["cNA"] <- "grey50"
#     }
#     
#     # Combine the two color maps for the overall color scale
#     combined_colors_map <- c(status_colors_map, cluster_colors_map)
#     
#     
#     # Create the plot object
#     network_plot <- ggraph(graph, layout = 'fr') +
#       geom_edge_link(aes(color = factor(cluster_for_edge_color)), alpha = 0.6, width = 0.5) +
#       
#       geom_node_point(aes(shape = type, color = if_else(type == 'GO', Wilting_Status, name)), size = 5, alpha = 0.8) +
#       
#       # Add GO Term labels ONLY for GO nodes
#       # Explicitly use ggraph:: prefix for the function
#       geom_node_text(aes(label = if_else(type == 'GO', Term, "")), # CORRECTED: Changed NULL to ""
#                      check_overlap = TRUE, # Keep this for overlap handling
#                      size = 3, # Adjust label size
#                      color = "black", # Label color
#                      nudge_y = -0.2 # Optionally adjust vertical position slightly
#                      # REMOVED: max.overlaps = 10
#       ) +
#       
#       scale_shape_manual(values = c("GO" = 19, "Cluster" = 17)) +
#       scale_color_manual(values = combined_colors_map) +
#       
#       guides(shape = guide_legend(title = "Node Type"),
#              color = guide_legend(override.aes = list(alpha = 1, size = 3)),
#       ) +
#       theme_graph() +
#       labs(title = "GO Term - Cluster Association Network",
#            subtitle = paste0("Clusters from k=", K_VALUE, " Enrichment"),
#            caption = "Nodes: GO Terms (circles), Clusters (triangles)\nGO Node Color: Wilting Status | Edge/Cluster Node Color: Cluster ID")
#     
#     # Print the plot to the graphics device (optional)
#     print(network_plot)
#     
#     
#     # --- Save the plot to a publication-quality PNG ---
#     cat("Saving network plot to PNG...\n")
#     
#     # Define filename
#     network_filename <- file.path(output_dir, paste0("Wilting_GO_Cluster_Network_k", K_VALUE, "_", current_date_str, ".png"))
#     
#     # Save using ggsave
#     tryCatch({
#       ggsave(filename = network_filename,
#              plot = network_plot,
#              device = "png",
#              dpi = 300,
#              width = 20, # Try increasing more, e.g., 20 inches
#              height = 20, # Try increasing more, e.g., 20 inches
#              units = "in")
#       cat("Network plot saved successfully to:", network_filename, "\n")
#     }, error = function(e) {
#       warning("ERROR saving network plot file: ", e$message)
#     })
#     
#     
#   } else {
#     cat("No GO-Cluster pairs found after filtering, cannot generate or save network.\n")
#   }
#   
# } else {
#   cat("Final table is empty or has no cluster associations, cannot generate or save network.\n")
# }
# 
# sessionInfo()
# exists("geom_node_text_repel", where = as.environment("package:ggraph"))
