# Install and load required packages
# install.packages("eulerr")
# install.packages("gridExtra")
library(eulerr)
library(gridExtra)
library(ggplot2)
library(gridExtra)

# Set working directory
setwd("~/Desktop/arab_env_2024i/workflows")

# Define the directory path
directory <- "/home/ibg-4/Desktop/arab_env_2024i/studies/AnemAsag_EPMclu/"
################################################################################
pval_treshold <- 0.01
eval_treshold <- 1
################################################################################
# Read each CSV file with tab separation
epmAnemS0 <- read.csv(paste0(directory, "rdf5_epmAnemS0_cwm-motifs.jaspar_comparison_JASPAR2020.csv"), sep = "\t")
epmAsagS0 <- read.csv(paste0(directory, "rdf5_epmAsagS0_cwm-motifs.jaspar_comparison_JASPAR2020.csv"), sep = "\t")
epmAnemSS <- read.csv(paste0(directory, "rdf5_epmAnemSS_cwm-motifs.jaspar_comparison_JASPAR2020.csv"), sep = "\t")
epmAsagSS <- read.csv(paste0(directory, "rdf5_epmAsagSS_cwm-motifs.jaspar_comparison_JASPAR2020.csv"), sep = "\t")
epmAnemSW <- read.csv(paste0(directory, "rdf5_epmAnemSW_cwm-motifs.jaspar_comparison_JASPAR2020.csv"), sep = "\t")
epmAsagSW <- read.csv(paste0(directory, "rdf5_epmAsagSW_cwm-motifs.jaspar_comparison_JASPAR2020.csv"), sep = "\t")

# Sort files by the first column (subj) alphanumerically
epmAnemS0 <- epmAnemS0[order(epmAnemS0$subj), ]
epmAsagS0 <- epmAsagS0[order(epmAsagS0$subj), ]
epmAnemSS <- epmAnemSS[order(epmAnemSS$subj), ]
epmAsagSS <- epmAsagSS[order(epmAsagSS$subj), ]
epmAnemSW <- epmAnemSW[order(epmAnemSW$subj), ]
epmAsagSW <- epmAsagSW[order(epmAsagSW$subj), ]

head(epmAnemS0)
#function(df)
# Apply the filter to remove "R_" from the subj column in each file
files <- list(epmAnemS0, epmAsagS0, epmAnemSS, epmAsagSS, epmAnemSW, epmAsagSW)
files <- lapply(files, function(file) { file %>% filter(!grepl("R_", subj)) })
################################################# BOTTLENECK
filtered_files <- lapply(files, function(df) {
  df %>%
    filter(pval<pval_treshold & eval<eval_treshold)
})
files<-filtered_files
################################################# BOTTLENECK
# Extract unique "targ" values from each file
targ_list <- lapply(files, function(file) unique(file$targ))

# Get specific "targ" values for AsagS0, AsagSW, AsagSS
AsagS0 <- targ_list[[2]]
AsagSW <- targ_list[[5]]
AsagSS <- targ_list[[4]]

# Create the Euler diagram for AsagS0, AsagSW, and AsagSS
comparison_a <- euler(list(AsagS0 = AsagS0, AsagSW = AsagSW, AsagSS = AsagSS))

# Plot the Venn diagram with quantities (numbers)
a <- plot(comparison_a, main = "AsagS0 | AsagSW | AsagSS", quantities = TRUE)
print(a)
################################################################################
# Get specific "targ" values for AnemS0, AnemSW, AnemSS
AnemS0 <- targ_list[[1]]
AnemSW <- targ_list[[6]]
AnemSS <- targ_list[[3]]

# Create the Euler diagram for AnemS0, AnemSW, and AnemSS
comparison_b <- euler(list(AnemS0 = AnemS0, AnemSW = AnemSW, AnemSS = AnemSS))

# Plot the Venn diagram with quantities (numbers)
b <- plot(comparison_b, main = "AnemS0 | AnemSW | AnemSS", quantities = TRUE)
print(b)

####  #### ####   ####    ####    ####   ####  ####   ####   ####
# Create the Euler plots with custom colors for AsagS0, AsagSW, and AsagSS
comparison_a <- euler(list(AsagS0 = AsagS0, AsagSW = AsagSW, AsagSS = AsagSS))

# Create the Euler plots with custom colors for AnemS0, AnemSW, and AnemSS
comparison_b <- euler(list(AnemS0 = AnemS0, AnemSW = AnemSW, AnemSS = AnemSS))

# Create the Euler plots with custom colors for AsagS0, AsagSW, and AsagSS
comparison_c <- euler(list(AsagS0 = AsagS0, AnemS0 = AnemS0))

# Create the Euler plots with custom colors for AnemS0, AnemSW, and AnemSS
comparison_d <- euler(list(AsagSW = AsagSW, AnemSW = AnemSW))

# Create the Euler plots with custom colors for AnemS0, AnemSW, and AnemSS
comparison_e <- euler(list(AsagSS = AsagSS, AnemSS = AnemSS))

# # Set custom colors
# colors_a <- c("darkcyan", "darkcyan", "darkcyan")
# colors_b <- c("darkorange", "darkorange", "darkorange")
# 
# # Set custom border and fill colors for Asag
# fill_colors_a <- c("lightgrey", "bisque", "lightsalmon")

# Set custom border and fill colors for Asag
border_colors_a <- c("cyan3", "cyan3", "cyan3")
fill_colors_a <- c("lightgrey", "bisque", "lightblue")

# Set custom border and fill colors for Anem
border_colors_b <- c("darkcyan", "darkcyan", "darkcyan")
fill_colors_b <- c("grey", "bisque3", "lightblue4")

# Set custom border and fill colors for AnemSag
border_colors_c <- c("cyan3", "darkcyan")
fill_colors_c <- c("lightgrey", "grey")

# Set custom border and fill colors for AnemSag
border_colors_d <- c("cyan3", "darkcyan")
fill_colors_d <- c("bisque", "bisque3")

# Set custom border and fill colors for AnemSag
border_colors_e <- c("cyan3", "darkcyan")
fill_colors_e <- c("lightblue", "lightblue4")

# Use gridExtra to arrange the plots side by side
library(gridExtra)

# Create the plots with custom border and fill colors
plot_a <- plot(comparison_a, main = "AsagS0 | AsagSW | AsagSS", quantities = TRUE, 
               col = border_colors_a, fill = fill_colors_a)
plot_b <- plot(comparison_b, main = "AnemS0 | AnemSW | AnemSS", quantities = TRUE, 
               col = border_colors_b, fill = fill_colors_b)
# Create the plots with custom border and fill colors
plot_c <- plot(comparison_c, main = "AsagS0 | AnemS0", quantities = TRUE, 
               col = border_colors_a, fill = fill_colors_c)
plot_d <- plot(comparison_d, main = "AsagSW | AnemSW", quantities = TRUE, 
               col = border_colors_b, fill = fill_colors_d)
plot_e <- plot(comparison_e, main = "AsagSS | AnemSS", quantities = TRUE, 
               col = border_colors_b, fill = fill_colors_e)

# Arrange them side by side
grid.arrange(plot_a,
             plot_b,
             plot_c,
             plot_d,
             plot_e,
             ncol = 2)
###########################################################################################
# Extract unique "targ" values from each file
targ_list <- lapply(files, function(file) unique(file$targ))

# Get specific "targ" values for AsagS0, AsagSW, and AsagSS
AsagS0 <- targ_list[[2]]
AsagSW <- targ_list[[5]]
AsagSS <- targ_list[[4]]

# Create the Euler diagram for AsagS0, AsagSW, and AsagSS
comparison_a <- euler(list(AsagS0 = AsagS0, AsagSW = AsagSW, AsagSS = AsagSS))

# Get the common targ values for each intersection
common_asag_1_2_3 <- intersect(intersect(AsagS0, AsagSW), AsagSS)  # intersection of all three sets
common_asag_1_2 <- intersect(AsagS0, AsagSW)  # intersection of AsagS0 and AsagSW
common_asag_1_3 <- intersect(AsagS0, AsagSS)  # intersection of AsagS0 and AsagSS
common_asag_2_3 <- intersect(AsagSW, AsagSS)  # intersection of AsagSW and AsagSS

# Display the common "targ" values for each intersection
cat("Common between AsagS0, AsagSW, and AsagSS:\n", common_asag_1_2_3, "\n\n")
cat("Common between AsagS0 and AsagSW:\n", setdiff(common_asag_1_2, common_asag_1_2_3), "\n\n")
cat("Common between AsagS0 and AsagSS:\n", setdiff(common_asag_1_3, common_asag_1_2_3), "\n\n")
cat("Common between AsagSW and AsagSS:\n", setdiff(common_asag_2_3, common_asag_1_2_3), "\n\n")

# Get the unique targ values for each individual set
unique_asag_1 <- setdiff(AsagS0, union(AsagSW, AsagSS))  # Unique to AsagS0
unique_asag_2 <- setdiff(AsagSW, union(AsagS0, AsagSS))  # Unique to AsagSW
unique_asag_3 <- setdiff(AsagSS, union(AsagS0, AsagSW))  # Unique to AsagSS

# Display the unique "targ" values for each set
cat("Unique to AsagS0:\n", unique_asag_1, "\n\n")
cat("Unique to AsagSW:\n", unique_asag_2, "\n\n")
cat("Unique to AsagSS:\n", unique_asag_3, "\n\n")

# Repeat the same steps for the Anem files
AnemS0 <- targ_list[[1]]
AnemSW <- targ_list[[6]]
AnemSS <- targ_list[[3]]

# Create the Euler diagram for AnemS0, AnemSW, and AnemSS
comparison_b <- euler(list(AnemS0 = AnemS0, AnemSW = AnemSW, AnemSS = AnemSS))

# Get the common targ values for each intersection
common_anem_1_2_3 <- intersect(intersect(AnemS0, AnemSW), AnemSS)  # intersection of all three sets
common_anem_1_2 <- intersect(AnemS0, AnemSW)  # intersection of AnemS0 and AnemSW
common_anem_1_3 <- intersect(AnemS0, AnemSS)  # intersection of AnemS0 and AnemSS
common_anem_2_3 <- intersect(AnemSW, AnemSS)  # intersection of AnemSW and AnemSS

# Display the common "targ" values for each intersection
cat("Common between AnemS0, AnemSW, and AnemSS:\n", common_anem_1_2_3, "\n\n")
cat("Common between AnemS0 and AnemSW:\n", setdiff(common_anem_1_2, common_anem_1_2_3), "\n\n")
cat("Common between AnemS0 and AnemSS:\n", setdiff(common_anem_1_3, common_anem_1_2_3), "\n\n")
cat("Common between AnemSW and AnemSS:\n", setdiff(common_anem_2_3, common_anem_1_2_3), "\n\n")

# Get the unique targ values for each individual set
unique_anem_1 <- setdiff(AnemS0, union(AnemSW, AnemSS))  # Unique to AnemS0
unique_anem_2 <- setdiff(AnemSW, union(AnemS0, AnemSS))  # Unique to AnemSW
unique_anem_3 <- setdiff(AnemSS, union(AnemS0, AnemSW))  # Unique to AnemSS

# Display the unique "targ" values for each set
cat("Unique to AnemS0:\n", unique_anem_1, "\n\n")
cat("Unique to AnemSW:\n", unique_anem_2, "\n\n")
cat("Unique to AnemSS:\n", unique_anem_3, "\n\n")
#################################################################################
# Get the common targ values for each intersection
#common_asag_1_2_3 <- intersect(intersect(AsagS0, AsagSW), AsagSS)  # intersection of all three sets
common_asagnemS0_1_2 <- intersect(AsagS0, AnemS0)  # intersection of AsagS0 and AsagSW
common_asagnemSW_1_3 <- intersect(AsagSW, AnemSW)  # intersection of AsagS0 and AsagSS
common_asagnemSS_2_3 <- intersect(AsagSS, AnemSS)  # intersection of AsagSW and AsagSS
################################################################################
cat("Common between AnemS0 and AsagS0:\n", common_asagnemS0_1_2, "\n\n")
cat("Common between AnemSW and AsagSW:\n", common_asagnemSW_1_3, "\n\n")
cat("Common between AnemSS and AsagSS:\n", common_asagnemSS_2_3, "\n\n")
# Get the unique targ values for each individual set
unique_asaganemS0_1 <- setdiff(AsagS0, AnemS0)  # Uniq
unique_asaganemS0_2 <- setdiff(AnemS0, AsagS0)  # Unique to Anem
unique_asaganemSW_3 <- setdiff(AsagSW, AnemSW)  # 
unique_asaganemSW_4 <- setdiff(AnemSW, AsagSW)  # 
unique_asaganemSS_5 <- setdiff(AsagSS, AnemSS)  # 
unique_asaganemSS_6 <- setdiff(AnemSS, AsagSS)  # 

################################################################################
cat("Unique between AsagS0:\n", unique_asaganemS0_1, "\n\n")
cat("Unique between AnemS0:\n", unique_asaganemS0_2, "\n\n")
cat("Unique between AsagSW:\n", unique_asaganemSW_3, "\n\n")
cat("Unique between AnemSW:\n", unique_asaganemSW_4, "\n\n")
cat("Unique between AsagSS:\n", unique_asaganemSS_5, "\n\n")
cat("Unique between AnemSS:\n", unique_asaganemSS_6, "\n\n")
# Get the unique targ values for each individual set
