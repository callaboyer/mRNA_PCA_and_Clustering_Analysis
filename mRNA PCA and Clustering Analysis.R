# Install Packages
# install.packages("factoextra", "ggplot2", "cluster", "pheatmap")

# Load Libraries
library(factoextra)
library(ggplot2)
library(cluster)

# Generate Random Data
set.seed(243)

batch_data <- data.frame(
  mRNA_capping_efficiency = rnorm(30, mean = 50, sd = 5),
  mRNA_purity = rnorm(30, mean = 95, sd = 4),
  mRNA_concentration = rnorm(30, mean = 1.0, sd = 0.5),
  mRNA_residual_enzymes = runif(30, min = 50, max = 100),
  mRNA_length = rnorm(30, mean = 9000, sd = 1000)
  )

# Label Treatments for ID
groups <- c("A", "B", "C")
numbers <- 1:10

treatment_labels <- paste0("Treatment ", rep(groups, each=10), rep(numbers, times=length(groups)))
batch_data$TreatmentType <- treatment_labels

## PCA Analysis
# Perform PCA for Treatment Comparison
numeric_batch_data <- batch_data[sapply(batch_data, is.numeric)]
pca_result <- prcomp(numeric_batch_data, scale. = TRUE)
fviz_pca_ind(pca_result, geom.ind = "point", col.ind = "cos2", 
  palette = "jco", addEllipses = TRUE, title = "PCA of Batch Data")

## Hierarchical Dendrogram
# Create a Hierarchical Dendrogram 
dist_matrix <- dist(numeric_batch_data)
hc_result <- hclust(dist_matrix)

# Plot Hierarchical Dendrogram
plot(hc_result, labels=batch_data$TreatmentType, main = "PCA Hierarchical Dendrogram of mRNA Treatments") +
  theme_minimal()

## Silhouette Analysis
# Perform Silhouette Analysis on Cluster Analysis Results
silhouette_scores <- silhouette(cutree(hc_result, k=3), dist(numeric_batch_data))

# Plot Silhouette Analysis
plot(silhouette_scores, main = "PCA Silhouette Analysis of mRNA Treatments") +
  theme_minimal()

## Cluster Plot
# Create a Cluster Plot from PCA Data
num_clusters <- 3
clusters <- cutree(hc_result, k = num_clusters)
batch_data$Cluster <- as.factor(clusters)

# Prepare PCA Data for Cluster Plot
pca_data <- as.data.frame(pca_result$x)
pca_data$Cluster <- as.factor(clusters)
pca_data$TreatmentType <- batch_data$TreatmentType

# Plot Cluster Plot
print(ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster, label = TreatmentType)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, hjust = 1.5, size = 3) +
  labs(title = "PCA Cluster Plot of mRNA Treatments", x = "Principal Component 1", y = "Principal Component 2", color = "Cluster") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14))
  )

## PCA Bi-Plot
# Create a PCA Bi-Plot
if (!requireNamespace("factoextra", quietly = TRUE)) {
  install.packages("factoextra")
  }

# Plot PCA Bi-Plot
print(fviz_pca_var(pca_result, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), circle = TRUE,
  repel = TRUE, axes = c(1, 2), label = "var", labelsize = 5, geom = c("arrow","text"), addEllipses = FALSE, arrow.size = 0.5, arrow.type = "closed") +
    ggtitle("PCA Bi-Plot of mRNA Treatments") +
    xlab("Principal Component 1") +
    ylab("Principal Component 2") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14))
  )