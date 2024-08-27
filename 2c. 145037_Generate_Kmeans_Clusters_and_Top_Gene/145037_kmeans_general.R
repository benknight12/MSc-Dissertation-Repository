full_145037 <- read.csv("CB_normalized_data/normalised_145037.csv",row.names = 1)

scaled_145037 <- as.data.frame(scale(full_145037))

cluster_kmeans<-kmeans(scaled_145037, centers = 3, nstart=20)
save(cluster_kmeans, file="2c. 145037_Generate_Kmeans_Clusters_and_Top_Gene/kmeans_object.RData")
table(cluster_kmeans$cluster)

load("Useful_Functions/find_means.RData")
nfkb_four_genes <- c("NFKB1", "RELA",  "NFKB2", "RELB")
find_means(scaled_145037[,nfkb_four_genes], clusters = cluster_kmeans$cluster)

## Means indicate Cluster 1 is NonCanonical, cluster 2 is Atypical and cluster 3 is Canonical
clusters <- factor(cluster_kmeans$cluster, levels = c(1,2,3), labels = c("NonCanonical", "Atypical", "Canonical"))
write.csv(clusters, "2c. 145037_Generate_Kmeans_Clusters_and_Top_Gene/Cluster_list.csv")
library(factoextra)
plot1 <- fviz_cluster(cluster_kmeans, data = as.matrix(scaled_145037),
                      geom = "point",
                      ellipse.type = "convex",
                      ggtheme = theme_bw()
)
plot1
ggsave(
  "2c. 145037_Generate_Kmeans_Clusters_and_Top_Gene/NFKB_derived_clustering.png",
  plot = plot1)

load("Useful_Functions/Diff_expression_and_volc_plot_from_clusters_and_data.RData")
## Run volcano plot and gene lists on unscaled data - doesn't have much of an effect
diff_gene_and_volc_plots(full_145037, clusters, "2c. 145037_Generate_Kmeans_Clusters_and_Top_Gene/")

