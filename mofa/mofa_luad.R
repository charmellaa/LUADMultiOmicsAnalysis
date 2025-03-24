# set the "LUADMultiOmicsAnalysis" folder as working directory

rm(list=ls())

dirRes <- "./mofa/MOFAResults/"

if (!dir.exists(dirRes)) {
  dir.create(dirRes)
} else {
  print(paste("The directory", dirRes, "already exists"))
}

path <- "./datasets"
dataPath <- paste0(path, "/luad")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(MOFA2)
library(ggplot2)
library(stringr)

##### INPUT FILES ######

# RNA_SEQ FILES
rna_seq_file <- paste0(dataPath, "/matrix_RNAseq_LUAD.txt")
normal_rna_file <- paste0(dataPath, "/List_Normal_LUAD_RNAseq_common.txt")
tumor_rna_file <- paste0(dataPath, "/List_Tumor_LUAD_RNAseq_common.txt")

# miRNA_SEQ FILES
mirna_seq_file <- paste0(dataPath, "/matrix_miRNAseq_LUAD.txt")
normal_mirna_file <- paste0(dataPath, "/List_Normal_LUAD_miRNAseq_common.txt")
tumor_mirna_file <- paste0(dataPath, "/List_Tumor_LUAD_miRNAseq_common.txt")

# CLINICAL DATA
clinical_file <- paste0(dataPath, "/clinical_luad.txt")

############################################################################

# LOAD RNA_SEQ DATA
rna_seq <- read.table(rna_seq_file, header = T, sep = "\t", check.names = F, row.names = 1, quote = "")
genes <- rownames(rna_seq) # vector containing gene names (+ ensembl id)

normal_rna <- scan(normal_rna_file, what = character(), quiet = TRUE) # list of normal samples in rna_seq
tumor_rna <- scan(tumor_rna_file, what = character(), quiet = TRUE) # list of tumor samples in rna_seq
common_rna <- str_extract(tumor_rna, "TCGA-\\w+-\\w+") # extract the IDs

data_normal_rna <- rna_seq[, normal_rna] # subset of the rna seq matrix with only the normal samples
data_tumor_rna <- rna_seq[, tumor_rna] # subset of the rna seq matrix with only the tumor samples

data_all_rna <- cbind(data_normal_rna, data_tumor_rna) # subset of normal and tumor rna for same patients
patients_full_rna_ids <- colnames(data_all_rna) # patient all ids (normal and tumor) rna

# filter by mean
overall_mean_rna <- apply(data_all_rna, 1, mean)
i <- which(overall_mean_rna == 0)
data_normal_rna <- data_normal_rna[-i, ]
data_tumor_rna <- data_tumor_rna[-i, ]
data_all_rna <- data_all_rna[-i, ]
genes <- genes[-i]

# transform in logarithmic base 2
data_normal_rna <- log2(data_normal_rna + 1)
data_tumor_rna <- log2(data_tumor_rna + 1)
data_all_rna <- log2(data_all_rna + 1)

#########################################################################

# LOAD MIRNA_SEQ DATA
mirna_seq <- read.table(mirna_seq_file, header = T, sep = "\t", check.names = F, row.names = 1, quote = "")
migenes <- rownames(mirna_seq) # vector

normal_mirna <- scan(normal_mirna_file, what = character(), quiet = TRUE)
tumor_mirna <- scan(tumor_mirna_file, what = character(), quiet = TRUE)
common_mirna <- str_extract(normal_mirna, "TCGA-\\w+-\\w+") # patient TCGA ids mirna

data_normal_mirna <- mirna_seq[, normal_mirna]
data_tumor_mirna <- mirna_seq[, tumor_mirna]

data_all_mirna <- cbind(data_normal_mirna, data_tumor_mirna)
patients_full_mirna_ids <- colnames(data_all_mirna) # patient all ids mirna (normal and tumor)

# filter by mean
overall_mean_mirna <- apply(data_all_mirna, 1, mean)
i <- which(overall_mean_mirna == 0)
data_normal_mirna <- data_normal_mirna[-i, ]
data_tumor_mirna <- data_tumor_mirna[-i, ]
data_all_mirna <- data_all_mirna[-i, ]
migenes <- migenes[-i]

# transform in logarithmic base 2
data_normal_mirna <- log2(data_normal_mirna + 1)
data_tumor_mirna <- log2(data_tumor_mirna + 1)
data_all_mirna <- log2(data_all_mirna + 1)

#########################################################################
# LOAD CLINICAL DATA
clinical_data <- read.table(clinical_file, header = T, sep = "\t", check.names = F, row.names = 1, quote = "")
all_patient_ids <- intersect(common_mirna, common_rna)
clinical_rna <- clinical_data[common_rna, ] # clinical data of the rna seq samples
clinical_mirna <- clinical_data[common_mirna, ] # clinical data of the miRNA seq samples
clinical_all <- clinical_data[unique(all_patient_ids), ] # clinical data of the patients for which i have RNA_seq and miRNA_seq data

##########################################################################

# PREPARE THE DATA FOR MOFA 
common_all <- intersect(common_mirna, common_rna) # patients for which i have both rna seq and mirna seq

common_normal_rna_ids <- unlist(lapply(common_all, function(x) { grep(x, normal_rna, value = TRUE) }))
common_tumor_rna_ids <- unlist(lapply(common_all, function(x) { grep(x, tumor_rna, value = TRUE) }))

common_normal_mirna_ids <- unlist(lapply(common_all, function(x) { grep(x, normal_mirna, value = TRUE) }))
common_tumor_mirna_ids <- unlist(lapply(common_all, function(x) { grep(x, tumor_mirna, value = TRUE) }))


#### MOFA2 #####

# Preprocessing
data_tumor_rna_common <- data_tumor_rna[, common_tumor_rna_ids]
data_tumor_mirna_common <- data_tumor_mirna[, common_tumor_mirna_ids]

# Modify column names for both datasets to show only the TCGA ID
colnames(data_tumor_rna_common) <- sub("^(TCGA-\\d{2}-\\d{4}-\\d{2})[A-Z]?.*", "\\1", colnames(data_tumor_rna_common))
colnames(data_tumor_mirna_common) <- sub("^(TCGA-\\d{2}-\\d{4}-\\d{2})[A-Z]?.*", "\\1", colnames(data_tumor_mirna_common))

# Find the common samples between the two datasets
common_samples <- intersect(colnames(data_tumor_rna_common), colnames(data_tumor_mirna_common))

data_tumor_rna_common <- data_tumor_rna_common[, common_samples]
data_tumor_mirna_common <- data_tumor_mirna_common[, common_samples]

data_tumor_rna_common_clean <- data_tumor_rna_common[complete.cases(data_tumor_rna_common), ]
data_tumor_mirna_common_clean <- data_tumor_mirna_common[complete.cases(data_tumor_mirna_common), ]

# variance
var_rna <- apply(data_tumor_rna_common_clean, 1, var, na.rm = TRUE)
var_mirna <- apply(data_tumor_mirna_common_clean, 1, var, na.rm = TRUE)

# Keep most variable
threshold_rna <- quantile(var_rna, 0.75, na.rm = TRUE)
threshold_mirna <- quantile(var_mirna, 0.75, na.rm = TRUE)

data_tumor_rna_common_filtered <- data_tumor_rna_common_clean[var_rna > threshold_rna, ]
data_tumor_mirna_common_filtered <- data_tumor_mirna_common_clean[var_mirna > threshold_mirna, ]

# scaling
data_tumor_rna_common_filtered <- scale(data_tumor_rna_common_filtered)
data_tumor_mirna_common_filtered <- scale(data_tumor_mirna_common_filtered)

if (anyNA(data_tumor_rna_common_filtered) | anyNA(data_tumor_mirna_common_filtered)) {
  stop("There are still missing values in the filtered data.")
}

######### FOR NORMAL ################ # same procedure as tumor #
data_normal_rna_common <- data_normal_rna[, common_normal_rna_ids]
data_normal_mirna_common <- data_normal_mirna[, common_normal_mirna_ids]

colnames(data_normal_rna_common) <- sub("^(TCGA-\\d{2}-\\d{4}-\\d{2})[A-Z]?.*", "\\1", colnames(data_normal_rna_common))
colnames(data_normal_mirna_common) <- sub("^(TCGA-\\d{2}-\\d{4}-\\d{2})[A-Z]?.*", "\\1", colnames(data_normal_mirna_common))

common_samples <- intersect(colnames(data_normal_rna_common), colnames(data_normal_mirna_common))

data_normal_rna_common <- data_normal_rna_common[, common_samples]
data_normal_mirna_common <- data_normal_mirna_common[, common_samples]

data_normal_rna_common_clean <- data_normal_rna_common[complete.cases(data_normal_rna_common), ]
data_normal_mirna_common_clean <- data_normal_mirna_common[complete.cases(data_normal_mirna_common), ]

var_rna <- apply(data_normal_rna_common_clean, 1, var, na.rm = TRUE)
var_mirna <- apply(data_normal_mirna_common_clean, 1, var, na.rm = TRUE)

threshold_rna <- quantile(var_rna, 0.75, na.rm = TRUE)
threshold_mirna <- quantile(var_mirna, 0.75, na.rm = TRUE)

data_normal_rna_common_filtered <- data_normal_rna_common_clean[var_rna > threshold_rna, ]
data_normal_mirna_common_filtered <- data_normal_mirna_common_clean[var_mirna > threshold_mirna, ]

# scaling
data_normal_rna_common_filtered <- scale(data_normal_rna_common_filtered)  
data_normal_mirna_common_filtered <- scale(data_normal_mirna_common_filtered)

if (anyNA(data_normal_rna_common_filtered) | anyNA(data_normal_mirna_common_filtered)) {
  stop("There are still missing values in the filtered data.")
}

## COMBINED 


get_patient_id <- function(sample_name) {
  return(substr(sample_name, 1, 12)) 
}

# sample type labels
sample_annotation <- data.frame(
  Sample_ID = colnames(cbind(data_normal_rna_common_filtered, data_tumor_rna_common_filtered)),
  Patient_ID = sapply(colnames(cbind(data_normal_rna_common_filtered, data_tumor_rna_common_filtered)), get_patient_id),
  Sample_Type = ifelse(grepl("-11$", colnames(cbind(data_normal_rna_common_filtered, data_tumor_rna_common_filtered))), "Normal", "Tumor")
)
# head(sample_annotation)

#combine matrices
data_combined_rna <- cbind(data_normal_rna_common_filtered, data_tumor_rna_common_filtered)
data_combined_mirna <- cbind(data_normal_mirna_common_filtered, data_tumor_mirna_common_filtered)

data_combined <- list(
  RNAseq = as.matrix(data_combined_rna),
  miRNAseq = as.matrix(data_combined_mirna)
)

dim(data_combined$RNAseq)
dim(data_combined$miRNAseq)

mofa_obj_combined <- create_mofa(data_combined)

data_opts <- get_default_data_options(mofa_obj_combined)
model_opts <- get_default_model_options(mofa_obj_combined)
train_opts <- get_default_training_options(mofa_obj_combined)
model_opts$num_factors <- 9


# Prepare MOFA
mofa_obj_combined <- prepare_mofa(
  object = mofa_obj_combined,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

# Run MOFA training
mofa_obj_combined.trained <- run_mofa(
  mofa_obj_combined, 
  use_basilisk = TRUE, 
  outfile = paste0(dirRes, "mofa_combined_model.hdf5")
)

mofa_combined_model <- load_model(paste0(dirRes, "mofa_combined_model.hdf5"))

pdf(paste0(dirRes, "mofa_combined.pdf"))

# Plot the factors colored by sample type
#plot_factors(mofa_combined_model, factors = 1:2, color_by = sample_annotation$Sample_Type)

# Extract MOFA factors
factor_values_combined <- as.data.frame(get_factors(mofa_combined_model)$group1)
factor_values_combined$Sample_ID <- rownames(factor_values_combined)

factor_values_combined <- merge(factor_values_combined, sample_annotation, by = "Sample_ID")

# Compute variance of RNA and miRNA data and plot their distributions
#rna_variance <- apply(data_combined$RNAseq, 1, var)
#mirna_variance <- apply(data_combined$miRNAseq, 1, var)
#h1 = hist(rna_variance, breaks = 50, main = "RNA-seq Variance Distribution", col = "skyblue")
#h2 = hist(mirna_variance, breaks = 50, main = "miRNA-seq Variance Distribution", col = "pink")


# PLOT THAT SHOWS WHICH FACTOR EXPLAINS THE MOST VARIANCE OVERALL
plot_variance_explained(mofa_combined_model)

# Check if factor is separating tumor vs normal samples
boxplot(Factor1 ~ Sample_Type, data = factor_values_combined, main = "Factor1 Separation", col = c("blue", "red"))
ggplot(factor_values_combined, aes(x = Factor1, fill = Sample_Type)) +
  geom_density(alpha = 0.5) +  
  scale_fill_manual(values = c("blue", "red")) +  
  labs(title = "Density Plot of Factor1 Separation", x = "Factor1 Value", y = "Density") +
  theme_minimal() +
  theme(aspect.ratio = 1/3)  


factor1_weights <- get_weights(mofa_combined_model, factors = 1)

factor1_rnaseq <- factor1_weights$RNAseq

df_rna <- data.frame(Gene = rownames(factor1_rnaseq),
                     Weight = factor1_rnaseq[, 1])

df_rna <- df_rna[order(abs(df_rna$Weight), decreasing = TRUE), ]


head(df_rna, 10)

factor1_mirna <- factor1_weights$miRNAseq
df_mirna <- data.frame(miRNA = rownames(factor1_mirna),
                       Weight = factor1_mirna[, 1])
df_mirna <- df_mirna[order(abs(df_mirna$Weight), decreasing = TRUE), ]
head(df_mirna, 10)

top_10_genes <- head(df_rna, 10)
# Extract only the gene name (everything before '|')
top_10_genes$Gene <- sub("\\|.*", "", top_10_genes$Gene)

# Plot with cleaned gene names
ggplot(top_10_genes, aes(x = reorder(Gene, abs(Weight)), y = Weight)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  ggtitle("Top 10 Genes Contributing to Factor 1") +
  xlab("Genes") +
  ylab("Weight")



top_10_mirnas <- head(df_mirna, 10)

ggplot(top_10_mirnas, aes(x = reorder(miRNA, abs(Weight)), y = Weight)) +
  geom_bar(stat = "identity", fill = "firebrick") +
  coord_flip() +
  theme_minimal() +
  ggtitle("Top 10 miRNAs Contributing to Factor 1") +
  xlab("miRNAs") +
  ylab("Weight")


dev.off()

# mapping for the stages
stage_mapping <- c("Stage IA" = 1, "Stage IB" = 2, "Stage IIA" = 3, "Stage IIB" = 4,
                   "Stage IIIA" = 5, "Stage IIIB" = 6, "Stage IV" = 7)
# Merge clinical data with MOFA factors
clinical_factors_combined <- merge(factor_values_combined, clinical_all, 
                                   by.x = "Patient_ID", 
                                   by.y = "bcr_patient_barcode")

str(clinical_factors_combined)

clinical_factors_combined$Tumor_Stage_numeric <- stage_mapping[clinical_factors_combined$ajcc_pathologic_stage]

pca <- prcomp(clinical_factors_combined[, c("Factor1", "Factor2", "Factor3", "Factor6", "Tumor_Stage_numeric")], scale. = TRUE)
summary(pca)  # To see the explained variance



# Define factors to analyze
factors_to_plot <- c("Factor6", "Factor9") # put factors chosen for analysis; here for example we explore factor 6 and factor 9

# Loop through selected factors
for (factor_name in factors_to_plot) {
  
  # Define output PDF file
  pdf(paste0(dirRes, "corrplot_spearman_", factor_name, ".pdf"), width = 8, height = 6)  # Adjust size as needed
  
  # Compute Spearman correlation and test statistical significance
  cor_test_result <- cor.test(clinical_factors_combined[[factor_name]], clinical_factors_combined$Tumor_Stage_numeric, method = "spearman")
  cor_value <- cor_test_result$estimate
  p_value <- cor_test_result$p.value
  
  # Boxplot: Factor vs AJCC Pathologic Stage (ggplot2)
  print(
    ggplot(clinical_factors_combined, aes(x = ajcc_pathologic_stage, y = .data[[factor_name]], fill = ajcc_pathologic_stage)) +
      geom_boxplot(fill = "#AFCBFF", color = "black") +
      labs(title = paste(factor_name, "vs Tumor Stage", "- Correlation:", round(cor_value, 2), "- p-value:", round(p_value, 3)),
           x = "AJCC Pathologic Stage", 
           y = factor_name) +
      theme_minimal()
  )
  
  # Regression Plot: Factor vs AJCC Pathologic Stage
  print(
    ggplot(clinical_factors_combined, aes(x = ajcc_pathologic_stage, y = .data[[factor_name]])) +  
      geom_jitter(aes(color = ajcc_pathologic_stage), width = 0.2, alpha = 0.7) +  # Jittered points to avoid overlap
      geom_smooth(aes(group = 1), method = "lm", color = "red", se = FALSE) +  # Linear trend line
      labs(title = paste(factor_name, "- Correlation:", round(cor_value, 2), "- p-value:", round(p_value, 3)),
           x = "AJCC Pathologic Stage", 
           y = factor_name) + 
      theme_minimal() +
      theme(legend.position = "none")  # Remove legend
  )
  
  # Close the PDF file
  dev.off()
}


# 
### UNCOMMENT THIS SECTION FOR FACTOR ANALYSIS FOR AVERAGED SAMPLES)
# 
# # Average factor values per patient (group by Patient_ID)
# library(dplyr)
# 
# clinical_factors_averaged <- clinical_factors_combined %>%
#   group_by(Patient_ID, ajcc_pathologic_stage, Tumor_Stage_numeric) %>%  # Group by Patient_ID
#   summarise(across(starts_with("Factor"), mean, na.rm = TRUE), .groups = "drop")  # Average Factor values
# 
# # Check new structure
# str(clinical_factors_averaged)  # Now each patient has only ONE row
# 
# # Create a PDF file
# pdf(paste0(dirRes,"factorsvsstage.pdf"), width = 10, height = 6)
# for (i in 1:9) {
#   factor_name <- paste0("Factor", i)
#   
#   plot(clinical_factors_averaged[[factor_name]], clinical_factors_averaged$Tumor_Stage_numeric,
#        xlab = factor_name, 
#        ylab = "Tumor Stage", 
#        main = paste("Scatter plot of", factor_name, "from combined samples vs Tumor Stage"))
# }
# dev.off()
# # Define factors to analyze
# factors_to_plot <- c("Factor6", "Factor9")
# 
# # Loop through selected factors using averaged data
# for (factor_name in factors_to_plot) {
#   
#   # Define output PDF file
#   pdf(paste0(dirRes, "corrplot_spearman_", factor_name, "_averaged.pdf"), width = 10, height = 6)
#   
#   # Calculate Spearman correlation on averaged data
#   cor_test_result <- cor.test(clinical_factors_averaged[[factor_name]], clinical_factors_averaged$Tumor_Stage_numeric, method = "spearman")
#   cor_value <- cor_test_result$estimate
#   p_value <- cor_test_result$p.value  
#   
#   # Boxplot: Factor vs AJCC Pathologic Stage
#   print(
#     ggplot(clinical_factors_averaged, aes(x = ajcc_pathologic_stage, y = .data[[factor_name]], fill = ajcc_pathologic_stage)) +
#       geom_boxplot(fill = "#AFCBFF", color = "black") +
#       labs(title = paste(factor_name, "vs Tumor Stage", "- Correlation:", round(cor_value, 2), "- p-value:", round(p_value, 3)),
#            x = "AJCC Pathologic Stage", 
#            y = factor_name) +
#       theme_minimal()
#   )
#   
#   # Regression Plot: Factor vs AJCC Pathologic Stage
#   print(
#     ggplot(clinical_factors_averaged, aes(x = ajcc_pathologic_stage, y = .data[[factor_name]])) +  
#       geom_jitter(aes(color = ajcc_pathologic_stage), width = 0.2, alpha = 0.7) +  # Jittered points to avoid overlap
#       geom_smooth(aes(group = 1), method = "lm", color = "red", se = FALSE) +  # Linear trend line
#       labs(title = paste(factor_name, "- Correlation:", round(cor_value, 2),  "- p-value:", round(p_value, 3)),
#            x = "AJCC Pathologic Stage", 
#            y = factor_name) + 
#       theme_minimal() +
#       theme(legend.position = "none")  # Remove legend
#   )
#   
#   # Close the PDF file
#   dev.off()
# }
