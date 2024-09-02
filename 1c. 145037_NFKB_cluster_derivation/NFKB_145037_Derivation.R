library(EnhancedVolcano)
full_145037 <- read.csv("CB_normalized_data/normalised_145037.csv",  row.names = 1)
scaled_145037 <- as.data.frame(scale(full_145037))

nfkb_four_genes <- c("NFKB1","RELA","NFKB2","RELB")

nfkb_145037 <- scaled_145037[,nfkb_four_genes]
head(nfkb_145037)


## This divides up the samples into those with expression values that very strongly indicate cluster membership and then from these high likelihood assignments find centroid points by taking their means
label_data <- function(df) {
  labels <- rep("messy", nrow(df))
  
  # Cluster 1: Upregulated NFKB1 and RELA, Downregulated NFKB2 and RELB
  labels[df$NFKB1 > 0 & df$RELA > 0 & df$NFKB2 < 0 & df$RELB < 0] <- "up_NFKB1_RELA_down_NFKB2_RELB"
  
  # Cluster 2: Downregulated NFKB1 and RELA, Upregulated NFKB2 and RELB
  labels[df$NFKB1 < 0 & df$RELA < 0 & df$NFKB2 > 0 & df$RELB > 0] <- "down_NFKB1_RELA_up_NFKB2_RELB"
  
  return(labels)
}
centroids <- rbind(colMeans(nfkb_145037[label_data(nfkb_145037)=="up_NFKB1_RELA_down_NFKB2_RELB",]), colMeans(nfkb_145037[label_data(nfkb_145037)=="down_NFKB1_RELA_up_NFKB2_RELB",]),colMeans(nfkb_145037[label_data(nfkb_145037)=="messy",]))
centroids <- rbind(
  c(2,1,-1,-1),
  c(-2,-1,1,1),
  c(-2,0,0,0)
)

kmeans_145037_nfkb <- kmeans(nfkb_145037, centers=centroids, nstart = 30)
save(kmeans_145037_nfkb, file = "1c. 145037_NFKB_cluster_derivation/kmeans_object.RData")
load("Useful_Functions/find_means.RData")
load("1c. 145037_NFKB_cluster_derivation/kmeans_object.RData")

find_means(nfkb_145037, kmeans_145037_nfkb$cluster)
## c1 = canonical, c2= noncaonical, c2 =atypical
table(kmeans_145037$cluster)

nfkb_derived_clusters_145037 <- factor(kmeans_145037$cluster, levels=c(1,2,3), labels = c("Canonical", "NonCanonical", "Atypical"))

## Save the clusters
write.csv(nfkb_derived_clusters_145037, "1c. 145037_NFKB_cluster_derivation/NFKB_derived_scort_clusters.csv")


## Plot the clusters
data_matrix <- as.matrix(full_145037)
cluster_object <- list(data = data_matrix, cluster = nfkb_derived_clusters_145037)
plot1 <- fviz_cluster(cluster_object, data = data_matrix, geom = "point", ellipse.type = "convex", ggtheme = theme_bw())

plot1
ggsave(
  "1c. 145037_NFKB_cluster_derivation/NFKB_derived_clustering.png",
  plot = plot1)

load("Useful_Functions/Diff_expression_and_volc_plot_from_clusters_and_data.RData")
diff_gene_and_volc_plots(full_145037, nfkb_derived_clusters_145037, folder_name = "1c. 145037_NFKB_cluster_derivation/", 1, 0.05)
