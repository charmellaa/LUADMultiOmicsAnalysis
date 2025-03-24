# set the "LUADMultiOmicsAnalysis" folder as working directory

rm(list=ls())

dirRes <- "./snf/SNFResults/"

if (!dir.exists(dirRes)) {
  dir.create(dirRes)
} else {
  print(paste("The directory", dirRes, "already exists"))
}

path <- "./datasets"
dataPath <- paste0(path, "/luad")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require(ggpubr)) install.packages("ggpubr")

library(ggpubr)
library(SNFtool)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)

rna_seq_file <- paste0(dataPath, "/matrix_RNAseq_LUAD.txt")
mirna_seq_file <- paste0(dataPath, "/matrix_miRNAseq_LUAD.txt")
clinical_file <- paste0(dataPath, "/clinical_luad.txt")

rna_seq <- read.table(rna_seq_file, header = T, sep = "\t", check.names = F, row.names = 1, quote = "")
genes <- rownames(rna_seq)

list_tumor_rna <- colnames(rna_seq)
list_tumor_rna <- list_tumor_rna[grepl("^TCGA-[^-]+-[^-]+-01", list_tumor_rna)]

rna_seq_tumor <- rna_seq[,list_tumor_rna]
colnames(rna_seq_tumor) <- str_extract(colnames(rna_seq_tumor), "TCGA-\\w+-\\w+")


mirna_seq <- read.table(mirna_seq_file, header = T, sep = "\t", check.names = F, row.names = 1, quote = "")
migenes <- rownames(mirna_seq)

list_tumor_mirna <- colnames(mirna_seq)
list_tumor_mirna <- list_tumor_mirna[grepl("^TCGA-[^-]+-[^-]+-01", list_tumor_mirna)]

mirna_seq_tumor <- mirna_seq[,list_tumor_mirna]
colnames(mirna_seq_tumor) <- str_extract(colnames(mirna_seq_tumor), "TCGA-\\w+-\\w+")


clinical_data <- read.table(clinical_file, header = T, sep = "\t", check.names = F, row.names = 1, quote = "")

common_samples <- unique(intersect(colnames(mirna_seq_tumor), colnames(rna_seq_tumor)))

common_patients <- intersect(rownames(clinical_data), common_samples)

rna_seq_tumor <- rna_seq_tumor[,common_patients]
mirna_seq_tumor <- mirna_seq_tumor[,common_patients]
clinical_tumor <- clinical_data[common_patients,]

rna_seq_tumor <- rna_seq_tumor[apply(rna_seq_tumor, 1, var) > 0, ]
mirna_seq_tumor <- mirna_seq_tumor[apply(mirna_seq_tumor, 1, var) > 0, ]

rna_norm <- standardNormalization(t(rna_seq_tumor))
mirna_norm <- standardNormalization(t(mirna_seq_tumor))

rna_dist <- 1 - cor(t(rna_norm), method = "pearson") 
mirna_dist <- 1- cor(t(mirna_norm), method = "pearson")

K <- 20
alpha <- 0.5
T <- 10

W_RNA <- affinityMatrix(rna_dist, K, alpha)
W_miRNA <- affinityMatrix(mirna_dist, K, alpha)

W_fused <- SNF(list(W_RNA, W_miRNA), K, T)

estimationResult = estimateNumberOfClustersGivenGraph(W_fused, 2:5);
# estimationResult returns 3 as best number of clusters based on eigengap
# we set it as number of clusters
num_clusters <- 3
spectral_clusters <- spectralClustering(W_fused, num_clusters)


# Store results in a dataframe
patient_clusters <- data.frame(Patient_ID = common_patients, Cluster = as.factor(spectral_clusters))
table(patient_clusters$Cluster)  # Check cluster sizes
library(igraph)

# Perform PCA on W_fused and create a plot
pca <- prcomp(W_fused, scale = TRUE)
pca_df_spectral <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Cluster = patient_clusters$Cluster)

pca_spectral <- ggplot(pca_df_spectral, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Patient Clusters",
       x = "Component 1",
       y = "Component 2")

ggsave(paste0(dirRes, "PCA_SNF.pdf"), 
       plot = pca_spectral, width = 8, height = 6)





# Merge clinical data with cluster assignments
# Ensure patient IDs match
clinical_tumor$bcr_patient_barcode <- rownames(clinical_tumor)

# Merge clinical data with clusters
merged_data <- merge(patient_clusters, clinical_tumor, 
                     by.x = "Patient_ID", 
                     by.y = "bcr_patient_barcode")
# Convert Cluster to factor
merged_data$Cluster <- as.factor(merged_data$Cluster)



merged_data$ajcc_stage_numeric <- as.numeric(factor(merged_data$ajcc_pathologic_stage, 
    levels = c("Stage I", "Stage IA", "Stage IB", "Stage II", "Stage IIA", "Stage IIB", "Stage IIIA", "Stage IIIB", "Stage IV"), ordered = TRUE))

cor.test(merged_data$ajcc_stage_numeric, as.numeric(merged_data$Cluster), method = "spearman")



# Convert columns to appropriate types
merged_data <- merged_data %>%
  mutate(
    Cluster = as.factor(Cluster),
    ajcc_pathologic_stage = as.factor(ajcc_pathologic_stage),  # Convert to factor
    age_at_diagnosis = as.numeric(age_at_diagnosis),
    cigarettes_per_day = as.numeric(cigarettes_per_day),
    prior_malignancy = factor(prior_malignancy, levels = c("no", "yes", "not reported")),
    vital_status = factor(vital_status, levels = c("Alive", "Dead"))
  )

# Handle missing data: remove rows with NA in 'age_at_diagnosis' and 'cigarettes_per_day'
merged_data_cleaned <- merged_data %>%
  filter(!is.na(age_at_diagnosis) & !is.na(cigarettes_per_day))

# Ensure 'Not Reported' is a valid level for ajcc_pathologic_stage
merged_data_cleaned$ajcc_pathologic_stage <- factor(merged_data_cleaned$ajcc_pathologic_stage, 
                                                    levels = c(levels(merged_data_cleaned$ajcc_pathologic_stage), "Not Reported"))

# Replace NA values with "Not Reported"
merged_data_cleaned$ajcc_pathologic_stage[is.na(merged_data_cleaned$ajcc_pathologic_stage)] <- "Not Reported"


# summary
summary(merged_data_cleaned$age_at_diagnosis)
summary(merged_data_cleaned$cigarettes_per_day)

# Frequency table for categorical variables
table(merged_data_cleaned$Cluster)
table(merged_data_cleaned$prior_malignancy)
table(merged_data_cleaned$ajcc_pathologic_stage)  # Summary of AJCC Pathologic Stage

# Group-based summary
print(merged_data_cleaned %>%
        group_by(Cluster, ajcc_pathologic_stage) %>%
        summarize(
          avg_age = mean(age_at_diagnosis, na.rm = TRUE),
          avg_cigarettes = mean(cigarettes_per_day, na.rm = TRUE),
          count = n()
        ), n = Inf)

# Check vital status by AJCC stage (survival across different AJCC stages)
print(merged_data_cleaned %>%
  group_by(ajcc_pathologic_stage, vital_status) %>%
  summarize(
    count = n()
  ), n = Inf)


# Cluster distribution across AJCC pathologic stage
cluster_stage_distribution <- merged_data_cleaned %>%
  group_by(Cluster, ajcc_pathologic_stage) %>%
  summarize(count = n()) %>%
  ungroup()  # To avoid grouping issues in further operations

# Display the distribution
print(cluster_stage_distribution, n = Inf)


# Vital status by cluster and AJCC stage
cluster_vital_status <- merged_data_cleaned %>%
  group_by(Cluster, ajcc_pathologic_stage, vital_status) %>%
  summarize(count = n()) %>%
  ungroup()  # Ensure grouping is not persistent for other analyses

# Display the distribution of vital status by cluster and stage
print(cluster_vital_status, n = Inf)





# Define the output PDF file
pdf(paste0(dirRes, "LUAD_Overall_Analysis.pdf"), width = 8, height = 6)

print(pca_spectral) # Scatterplot from PCA 

# Cluster Distribution Across AJCC Stages (Bar Plot)
p1 <- ggplot(cluster_stage_distribution, aes(x = ajcc_pathologic_stage, y = count, fill = Cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Cluster Distribution Across AJCC Stages",
       x = "AJCC Pathologic Stage",
       y = "Number of Patients") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p1)

# Convert Age at Diagnosis from Days to Years
merged_data_cleaned <- merged_data_cleaned %>%
  mutate(age_at_diagnosis = age_at_diagnosis / 365.25)  # Convert days to years

# Age at Diagnosis by Cluster (Boxplot)
p2 <- ggplot(merged_data_cleaned, aes(x = Cluster, y = age_at_diagnosis, fill = Cluster)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Age at Diagnosis by Cluster",
       x = "Cluster",
       y = "Age at Diagnosis (yrs)")
print(p2)

# Cigarettes per Day by Cluster (Boxplot)
p3 <- ggplot(merged_data_cleaned, aes(x = Cluster, y = cigarettes_per_day, fill = Cluster)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Cigarette Consumption by Cluster",
       x = "Cluster",
       y = "Cigarettes per Day")
print(p3)

# Vital Status Across Clusters & Stages (Bar Plot)
p4 <- ggplot(cluster_vital_status, aes(x = ajcc_pathologic_stage, y = count, fill = vital_status)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Cluster) +
  theme_minimal() +
  labs(title = "Survival Status Across Clusters and Stages",
       x = "AJCC Pathologic Stage",
       y = "Number of Patients") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("lightgreen", "darkgrey"))
print(p4)


# Survival by Cluster
ggplot(cluster_vital_status, aes(x = Cluster, y = count, fill = vital_status)) +
  geom_bar(stat = "identity", position = "dodge") +  # Proportion bar chart
  theme_minimal() +
  labs(title = "Survival by Cluster",
       x = "Cluster",
       y = "Proportion",
       fill = "Vital Status") +
  scale_fill_manual(values = c("lightgreen", "darkgrey")) +
theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Close the PDF file
dev.off()
