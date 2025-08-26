#[2024-13-08]
# Version 2 - Now with "target" column for training
#Define Experiment
count_acro <- "Atha-D0h-Mn0Q5" 
#four letter, Species names [1] + epithet [3] + condition [2] + classification [4]  Mn0Q4: log Max Median non 0 quartiles 25% Q4
step <- "_expr_class"
#############################################
File1 <- "arabis-as-ctrl_tpm.csv"
input_dir <- "../../runs/Arabis_nem_sag_drought_mRNA_counts/"
setwd("~/Desktop/arab_env_2024i/workflows/read_count_handling")
##############################################
# Set upper percentile treshold
Q <- 0.80
#############################################
df0 <- read.csv(paste0(input_dir, File1),
               sep = ",", header = TRUE, row.names = 1)
##########################################
# Try to remove logMaxTPM and true_target columns
df0 <- tryCatch({
  subset(df0, select = -c(logMaxTPM, true_target))
}, error = function(e) {
  message("Error: 'true_target' not found. Trying to remove 'logMaxTPM' and 'target' instead.")
  tryCatch({
    # Try to remove logMaxTPM and target columns if true_target does not exist
    subset(df0, select = -c(logMaxTPM, target))
  }, error = function(e) {
    message("Error: 'target' not found. Removing only 'logMaxTPM'.")
    # Remove only logMaxTPM if both true_target and target do not exist
    subset(df0, select = -c(logMaxTPM))
  })
})

# df0 will now have the columns removed as per the successful step
#colnames(df0)[1] <- "gene_id"
df0$log10MedTPM <- log10(apply(df0, 1, median)+1)
###############################################
upper_quantile <- quantile(df0$log10MedTPM, Q)
lower_quantile <- quantile(df0$log10MedTPM, (1-Q))
count_above_upper_quantile <- sum(df0$log10MedTPM > upper_quantile)
##############################################
# Select contrast class including zeros and lowly expressed genes and balance with highly expressed genes
zero_class_count <- count_above_upper_quantile/2
subset_zero <- df0[df0$log10MedTPM == 0, ]
if (nrow(subset_zero) < zero_class_count) {
  warning("There are fewer rows with log10MedTPM equal to zero than required. Please use percentile treshold for classification.")
}
#set.seed(123)  # Optional: set seed for reproducibility
selected_rows_A <- subset_zero[sample(nrow(subset_zero), min(zero_class_count, nrow(subset_zero))), ]
subset_NONzero <- df0[df0$log10MedTPM > 0, ]
subset_NONzero_sorted <- subset_NONzero[order(subset_NONzero$log10MedTPM), ]
selected_rows_B <- head(subset_NONzero_sorted, zero_class_count)
##################################################################################
rows_above_upper_quantile <- df0[df0$log10MedTPM > upper_quantile, ]
duplicates <- intersect(rownames(rows_above_upper_quantile), rownames(selected_rows_B))
# Generate a warning if duplicates are found
if (length(duplicates) > 0) {
  warning("Duplicates found between rows with log10MedTPM greater than upper_quantile and selected_rows_B.")
}
# Display the duplicates (if any)
if (length(duplicates) > 0) {
  message("Duplicated row names: ", paste(duplicates, collapse = ", "))
}
#################################################################################
# Create true targets
df0$target <- 2
df0$target[rownames(df0) %in% rownames(rows_above_upper_quantile)] <- 1
selected_rows_combined <- rbind(selected_rows_A, selected_rows_B)
df0$target[rownames(df0) %in% rownames(selected_rows_combined)] <- 0
###############################################n
count_filtered_rows <- nrow(df0[df0$log10MedTPM > max(selected_rows_B$log10MedTPM, na.rm = TRUE) & df0$log10MedTPM < upper_quantile, ])
###############################################
print(paste0("EXPRESSION CLASS REPORT for:                                           ", count_acro))
print(paste0("Upper percentile treshold:                                             ", ((1-Q)*100),"%"))
print(paste0("Calculated upper cutoff:                                               ", upper_quantile))
print(paste0("Hypothetical lower cutoff:                                             ", lower_quantile))
print(paste0("Genes counted above cutoff:                                            ", count_above_upper_quantile))
print(paste0("Genes with a log10MedTPM equal to zero:                                ", sum(df0$target == 0)))
print(paste0("Genes with a log10MedTPM equal to zero selected and larger than zero : ", (sum(df0$target == 0))/2))
print(paste0("Lower percentile cutoff after balancing:                               ", max(selected_rows_B$log10MedTPM)))
print(paste0("Genes between lower and upper cutoff after balancing:                  ", count_filtered_rows))
##############################################
df0$gene_id <- rownames(df0)
df0 <- df0[, c("gene_id", colnames(df0)[-ncol(df0)])]
rownames(df0) <- NULL
head(df0)
###############################################
write.csv(df0, paste0(input_dir,count_acro, Q, step, ".csv"), row.names = FALSE)
###############################################
png_filename <- paste0(count_acro, Q, step, ".png")
full_path <- file.path(input_dir, png_filename)
print(paste("Saving PNG file to:", full_path))
# Open a PNG device
png(filename = full_path, width = 800, height = 600)
# Create the histogram with the adjusted color scheme
hist(df0$log10MedTPM,
     main = paste0("Distribution of log10MedTPM of ", count_acro, " ", Q),
     xlab = "log10MedTPM values",
     ylab = "Frequency",
     breaks = seq(floor(min(df0$log10MedTPM)), ceiling(max(df0$log10MedTPM)), length.out = 30),
     #     col = colors
)

# Close the graphics device
dev.off()
###############################################
#Write the report
report_filename <- paste0(count_acro, Q, step,"_report", ".txt")
full_report_path <- file.path(input_dir, report_filename)
report_conn <- file(full_report_path, open = "wt")

# Write the report content to the text file
writeLines(paste0("EXPRESSION CLASS REPORT for:                                       ", count_acro), con = report_conn)
writeLines(paste0("Upper percentile treshold:                                         ", ((1-Q)*100), "%"), con = report_conn)
writeLines(paste0("Calculated upper cutoff:                                           ", upper_quantile), con = report_conn)
writeLines(paste0("Hypothetical lower cutoff:                                         ", lower_quantile), con = report_conn)
writeLines(paste0("Genes counted above cutoff:                                        ", count_above_upper_quantile), con = report_conn)
writeLines(paste0("Genes with a log10MedTPM equal to zero:                            ", sum(df0$target == 0)), con = report_conn)
writeLines(paste0("Genes with a log10MedTPM equal to zero selected and larger than 0: ", (sum(df0$target == 0))/2), con = report_conn)
writeLines(paste0("Lower percentile cutoff after balancing:                           ", max(selected_rows_B$log10MedTPM)), con = report_conn)
writeLines(paste0("Genes between lower and upper cutoff after balancing:              ", count_filtered_rows), con = report_conn)

# Close the connection to the text file
close(report_conn)
#################################################