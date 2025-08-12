# Function to calculate FDR
calculate_fdr <- function(ref_data, pred_data) {
  # Prepare the reference data
  ref_data <- ref_data[, c(1, ncol(ref_data))]
  colnames(ref_data) <- c("gene_id", "true_target")
  
  # Prepare the prediction data
  pred_data <- pred_data[, c(1, ncol(pred_data))]
  colnames(pred_data) <- c("gene_id", "pred_prob")
  
  # Merge reference and prediction data
  df <- merge(ref_data, pred_data, by = "gene_id")
  df <- as.data.frame(df)
  
  # Create the true_false column
  df$true_false <- ifelse((df$pred_prob > 0.5 & df$true_target == 1) |
                            (df$pred_prob < 0.5 & df$true_target == 0), TRUE, FALSE)
  
  # Count TRUE and FALSE for each true_target value
  counts <- with(df, table(true_target, true_false))
  
  # Extract counts
  false_positives_0 <- ifelse("FALSE" %in% colnames(counts), counts["0", "FALSE"], 0)
  true_positives_0 <- ifelse("TRUE" %in% colnames(counts), counts["0", "TRUE"], 0)
  
  false_positives_1 <- ifelse("FALSE" %in% colnames(counts), counts["1", "FALSE"], 0)
  true_positives_1 <- ifelse("TRUE" %in% colnames(counts), counts["1", "TRUE"], 0)
  
  # Calculate FDR
  fdr_0 <- if (false_positives_0 + true_positives_0 > 0) false_positives_0 / (false_positives_0 + true_positives_0) else NA
  fdr_1 <- if (false_positives_1 + true_positives_1 > 0) false_positives_1 / (false_positives_1 + true_positives_1) else NA
  
  # Calculate average FDR
  fdr_total <- if (!is.na(fdr_0) && !is.na(fdr_1)) (fdr_0 + fdr_1) / 2 else NA
  
  # Return FDR values
  return(c(fdr_0 = fdr_0, fdr_1 = fdr_1, fdr_total = fdr_total))
}
setwd("/home/ibg-4/Desktop/arab_drought_jdm2024h/")
# Read reference data
ref0h <- read.csv("runs/Alyrata_drought_time_series_counts/Alyr-D0h-Xn0Q0.75_expr_class.csv", header = TRUE)

# Read prediction files
pred0h <- read.csv("studies/dCRE_model_predictions/mod_Atha-D0h-Xn0Q0.75_gen_Alyrata_v.1.0.dna.toplevel.fa_prediction.csv", header = TRUE)
pred3h <- read.csv("studies/dCRE_model_predictions/mod_Atha-D3h-Xn0Q0.75_gen_Alyrata_v.1.0.dna.toplevel.fa_prediction.csv", header = TRUE)
pred6h <- read.csv("studies/dCRE_model_predictions/mod_Atha-D6h-Xn0Q0.75_gen_Alyrata_v.1.0.dna.toplevel.fa_prediction.csv", header = TRUE)
pred12h <- read.csv("studies/dCRE_model_predictions/mod_Atha-D12h-Xn0Q0.75_gen_Alyrata_v.1.0.dna.toplevel.fa_prediction.csv", header = TRUE)

# Calculate FDR for each prediction file
fdr0h <- calculate_fdr(ref0h, pred0h)
fdr3h <- calculate_fdr(ref0h, pred3h)
fdr6h <- calculate_fdr(ref0h, pred6h)
fdr12h <- calculate_fdr(ref0h, pred12h)

# Create a data frame for results
fdr_results <- data.frame(
  Prediction = c("pred0h", "pred3h", "pred6h", "pred12h"),
  FDR_0 = c(fdr0h["fdr_0"], fdr3h["fdr_0"], fdr6h["fdr_0"], fdr12h["fdr_0"]),
  FDR_1 = c(fdr0h["fdr_1"], fdr3h["fdr_1"], fdr6h["fdr_1"], fdr12h["fdr_1"]),
  FDR_Total = c(fdr0h["fdr_total"], fdr3h["fdr_total"], fdr6h["fdr_total"], fdr12h["fdr_total"])
)

# Write the results to a CSV file
write.csv(fdr_results, "studies/refAlyr12h_Atha0h3h6h12h_cross_fdr_results.csv", row.names = FALSE)

# Print the results to the console
print(fdr_results)
