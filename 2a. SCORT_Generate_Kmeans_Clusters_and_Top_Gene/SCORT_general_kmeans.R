scort_data <- read.csv("CB_normalized_data/scort_normalized.csv",row.names = 1)
dim(scort_data)
scaled_scort <- as.data.frame(scale(scort_data))
pheno <- read.csv("Processed_Data_Sets/Patient data scort.csv", row.names = 1)
scort_kmeans<-kmeans(scaled_scort, centers = 3, nstart=20)
save(scort_kmeans, file="2a. SCORT_Generate_Kmeans_Clusters_and_Top_Gene/kmeans_object.RData")

table(scort_kmeans$cluster)

load("Useful_Functions/find_means.RData")
nfkb_four_genes <- c("NFKB1", "RELA",  "NFKB2", "RELB")
find_means(scaled_scort[,nfkb_four_genes], clusters = scort_kmeans$cluster)
## Means indicate Cluster 3 is aytpical, cluster 2 is canonical and cluster 1 is noncanonical
clusters <- factor(scort_kmeans$cluster, levels = c(1,2,3), labels = c("Canonical", "Atypical", "NonCanonical"))
write.csv(clusters, file = "2a. SCORT_Generate_Kmeans_Clusters_and_Top_Gene/cluster_list.csv")
table(clusters)
data_matrix <- as.matrix(scort_data)
cluster_object <- list(data = data_matrix, cluster = as.factor(clusters))
scort_kmeans$cluster <- clusters
plot1 <- fviz_cluster(scort_kmeans, data = as.matrix(scaled_scort),                  
                      geom = "point",
                      ellipse.type = "convex",
                      ggtheme = theme_minimal(),
                      legend.title = "Cluster Groups",
                      legend.labs =c("Canonical", "Atypical", "NonCanonical")
                      )
plot1
ggsave(
  "2a. SCORT_Generate_Kmeans_Clusters_and_Top_Gene/NFKB_derived_clustering.png",
  plot = plot1)

load("Useful_Functions/Diff_expression_and_volc_plot_from_clusters_and_data.RData")
diff_gene_and_volc_plots(scort_data, clusters, "2a. SCORT_Generate_Kmeans_Clusters_and_Top_Gene/")
dim(scort_data)
