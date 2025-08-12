setwd("/home/ibg-4/Desktop/arab_drought_jdm2024h/studies/dCRE_model_predictions/")
# Load necessary libraries
###################################################################################
file_list <- list.files(pattern = ".*gen_Atha-ref_Alyr-mod.fa_prediction.csv")
names_list <- sub("gen_Atha-ref_Alyr-mod.fa_prediction.csv", "", file_list)
###################################################################################
###################################################################################
# Function to read each CSV file and rename columns
read_and_rename <- function(file, name_prefix) {
  df <- read.csv(file, stringsAsFactors = FALSE)
  # Rename columns to include the name_prefix
  colnames(df)[-1] <- paste(colnames(df)[-1], name_prefix, sep = "_")
  return(df)
}

# Read files, rename columns, and name the list elements using the extracted names
df_list <- mapply(read_and_rename, file_list, names_list, SIMPLIFY = FALSE)

# Merge all dataframes in the list
# We assume merging on 'gene_ids' column
merged_df <- Reduce(function(x, y) merge(x, y, by = "gene_ids", all = TRUE), df_list)

# Remove rows that contain any NA values
merged_df <- na.omit(merged_df)

# Order columns alphabetically, keeping 'gene_ids' as the first column
ordered_columns <- c("gene_ids", sort(setdiff(names(merged_df), "gene_ids")))
merged_df <- merged_df[, ordered_columns]
###################################################################################
###################################################################################
# Load necessary libraries
library(dplyr)
library(tidyr)
library(corrplot)  # For visualization of correlations
library(ggplot2)   # For alternative visualization (if needed)
library(reshape2)  # For melting matrices

# Assuming `merged_df` is the cleaned and ordered dataframe

# Step 1: Calculate Pairwise Correlations
# Extract numeric columns for correlation analysis
numeric_cols <- merged_df %>% select(-gene_ids)  # Exclude 'gene_ids' column

# Calculate correlation matrix
correlation_matrix <- cor(numeric_cols, use = "pairwise.complete.obs", method = "pearson")

# Print the correlation matrix
print(correlation_matrix)

# Step 2: Visualize Correlation Matrix
# Plot the correlation matrix as a heatmap
corrplot(correlation_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         title = "Correlation Matrix of Columns", tl.cex = 0.7)

# Step 3: Compute Pairwise Similarities (Cosine Similarity)
# Function to compute cosine similarity between two vectors
cosine_similarity <- function(a, b) {
  sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
}

# Initialize a matrix to store cosine similarities
cosine_similarity_matrix <- matrix(0, ncol = ncol(numeric_cols), nrow = ncol(numeric_cols))
colnames(cosine_similarity_matrix) <- colnames(numeric_cols)
rownames(cosine_similarity_matrix) <- colnames(numeric_cols)

# Compute cosine similarities
for (i in 1:ncol(numeric_cols)) {
  for (j in i:ncol(numeric_cols)) {
    cosine_similarity_matrix[i, j] <- cosine_similarity(numeric_cols[[i]], numeric_cols[[j]])
    cosine_similarity_matrix[j, i] <- cosine_similarity_matrix[i, j]  # Symmetric matrix
  }
}

# Print the cosine similarity matrix
print(cosine_similarity_matrix)

# Step 4: Visualize Cosine Similarity Matrix
# Plot the cosine similarity matrix as a heatmap
corrplot(cosine_similarity_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         title = "Cosine Similarity Matrix of Columns", tl.cex = 0.7)

# Alternative visualization using ggplot2 and reshape2 if needed
# Melt the cosine similarity matrix for ggplot2
melted_cosine_matrix <- melt(cosine_similarity_matrix)
colnames(melted_cosine_matrix) <- c("Variable1", "Variable2", "CosineSimilarity")

# Plot using ggplot2
ggplot(melted_cosine_matrix, aes(x = Variable1, y = Variable2, fill = CosineSimilarity)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  coord_fixed() +
  labs(title = "Cosine Similarity Heatmap", x = "Variable", y = "Variable")
