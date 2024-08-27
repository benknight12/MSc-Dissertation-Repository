library(limma)
library(factoextra)

scort_data <- read.csv("Processed_Data_Sets/SCORT_final_data.csv",  row.names = 1)

scaled_data <- as.data.frame(scale(scort_data))
nfkb_four_genes <- c("NFKB1","RELA","NFKB2","RELB")
nfkb_four_genes2 <- c("NFkB1","RELA","NFkB2","RELB")

nfkb_scort_data <- scaled_data[,nfkb_four_genes]
head(nfkb_scort_data)

centroids <- rbind(
  c(1, 1, -1, -1),  # Upregulated NFKB1 and RELA, Downregulated NFKB2 and RELB
  c(-1, -1, 1, 1),
  c(0,0,0,0) # Downregulated NFKB1 and RELA, Upregulated NFKB2 and RELB  # Average of the messy category
)

kmeans_scort <- kmeans(nfkb_scort_data, centers=centroids, nstart = 30)
table(kmeans_scort$cluster)
save(kmeans_scort, file = "1a. SCORT_NFKB_cluster_derivation_replicating/My_NFKB_Scort_kmeans.RData")
find_means <- function(data, clusters){
  cluster_levels <- levels(as.factor(clusters))
  means <- as.data.frame(matrix(c(colMeans(data[clusters==cluster_levels[1],]),
                                  colMeans(data[clusters==cluster_levels[2],]),
                                  colMeans(data[clusters==cluster_levels[3],])),nrow=3,byrow=TRUE),row.names =c("Cluster 1","Cluster 2", "Cluster 3") )
  names(means) <- names(data)
  return(means)
}
save(find_means, file="Useful_Functions/find_means.RData")
find_means(nfkb_scort_data, kmeans_scort$cluster)


pheno <- read.csv("Processed_Data_Sets/Patient data scort.csv", row.names = 1)
find_means(pheno[,nfkb_four_genes2], pheno$NFkB.Cluster)
table(pheno$NFkB.Cluster, kmeans_scort$cluster)

## Use Adam's Clusters but use my code for geo sets
## Find COLMEANS by cluster:
find_means(pheno[,nfkb_four_genes2], pheno$NFkB.Cluster)

table(pheno$NFkB.Cluster)
## Save the clusters
write.csv(pheno$NFkB.Cluster, "1a. SCORT_NFKB_cluster_derivation_replicating/NFKB_derived_scort_clusters.csv")

## Plot the clusters
data_matrix <- as.matrix(scort_data)
cluster_object <- list(data = data_matrix, cluster = pheno$NFkB.Cluster)
plot1 <- fviz_cluster(cluster_object, data = as.matrix(scort_data),
                      geom = "point",
                      ellipse.type = "convex",
                      ggtheme = theme_bw()
)
plot1
ggsave(
  "1a. SCORT_NFKB_cluster_derivation_replicating/NFKB_derived_clustering.png",
  plot = plot1)



# Perform differential expression analysis using limma
clusters <- factor(pheno$NFkB.Cluster, levels=c("Atypical", "Canonical", "Non-Canonical"), labels = c("Atypical", "Canonical", "NonCanonical"))
clusters
scort_data$cluster <- clusters
design <- model.matrix(~0 + scort_data$cluster)
colnames(design) <- levels(scort_data$cluster)
head(data_matrix[,1:5])
fit <- lmFit(t(data_matrix), design)
head(design)
# Create contrasts to compare each cluster against the others
contrast.matrix <- makeContrasts(
  Canonical_vs_Others = Canonical - (NonCanonical + Atypical)/2,
  NonCanonical_vs_Others = NonCanonical - (Canonical + Atypical)/2,
  Atypical_vs_Others = Atypical - (Canonical + NonCanonical)/2,
  levels = design
)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract top genes for each cluster
top_genes_cluster1 <- rownames(topTable(fit2,lfc=1, coef = "Canonical_vs_Others", number = 20, adjust = "BH"))
top_genes_cluster2 <- rownames(topTable(fit2,lfc=1, coef = "NonCanonical_vs_Others", number = 20, adjust = "BH"))
top_genes_cluster3 <- rownames(topTable(fit2,lfc=1, coef = "Atypical_vs_Others", number = 20, adjust = "BH"))

top_gene_dataframe <- as.data.frame(matrix(c(top_genes_cluster1, top_genes_cluster2, top_genes_cluster3), ncol=3))
names(top_gene_dataframe) = c("Canonical","NonCanonical", "Atypical")
top_gene_dataframe

write.csv(top_gene_dataframe, "1a. SCORT_NFKB_cluster_derivation_replicating/SCORT_NFKB_Top_Gene_List.csv")

library(EnhancedVolcano)

results <- topTable(fit2, coef = "Canonical_vs_Others", number = Inf)
results_df <- as.data.frame(results)
print(summary(results_df$logFC))
print(summary(results_df$P.Value))
head(results_df)

volcano_plot <- EnhancedVolcano(
  results_df, lab = rownames(results), x = 'logFC', y = 'adj.P.Val', xlab = bquote(~Log[2]~ "fold change"),
  ylab = bquote(~-Log[10]~ "adjusted p-value"), pCutoff = 0.05, FCcutoff = 1,  title = 'Volcano plot of Differentially Expressed Features', 
  subtitle = 'Adjusted p-value < 0.05 and |log2 FC| > 1', pointSize = 1.0, labSize = 3.0, col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  colAlpha = 1, legendPosition = 'top', legendLabSize = 10, legendIconSize = 3.0, drawConnectors = TRUE, widthConnectors = 0.5,
  gridlines.major = TRUE, gridlines.minor = TRUE,max.overlaps = 50 )
volcano_plot

ggsave(
  "1a. SCORT_NFKB_cluster_derivation_replicating/NFKB_derived_volcano_plot1.png",
  plot = volc_plots$volc1)
ggsave(
  "1a. SCORT_NFKB_cluster_derivation_replicating/NFKB_derived_volcano_plot2.png",
  plot = volc_plots$volc2)

