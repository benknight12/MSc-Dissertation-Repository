full_109057 <- read.csv("Processed_Data_Sets/109057_set(unscaled).csv", row.names = 1)
scaled_109057 <- as.data.frame(scale(full_109057))
load("2a. SCORT_Generate_Kmeans_Clusters_and_Top_Gene/kmeans_object.RData")
load("1a. SCORT_NFKB_cluster_derivation_replicating/My_NFKB_Scort_kmeans.RData")
my_nfkb_kmeans_scort <- kmeans_scort


my_nfkb_kmeans_scort$centers

nfkb_four_genes <- c("NFKB1","RELA","NFKB2","RELB")


predict_kmeans <- function(newdata, kmeans_model) {
  library(fdm2id)
  # Calculate the distance between each point in newdata and the cluster centers from kmeans_model
  dist_matrix <- as.matrix(dist(rbind(kmeans_model$centers, newdata)))[-(1:nrow(kmeans_model$centers)), 1:nrow(kmeans_model$centers)]
  
  # Assign each point in newdata to the nearest cluster center
  cluster_assignments <- apply(dist_matrix, 1, which.min)
  
  return(cluster_assignments)
}
save(predict_kmeans, file = "Useful_Functions/Predict_Kmeans.RData")
new_clus <- predict_kmeans(scaled_109057[,nfkb_four_genes], my_nfkb_kmeans_scort)
find_means(data = scaled_109057[,nfkb_four_genes], new_clus)
## 1= canonical, 2=noncanonical, 3=atypical
clusters <- factor(new_clus, levels = c(1,2,3), labels = c("Canonical", "NonCanonical", "Atypical"))


load("Useful_Functions/Diff_expression_and_volc_plot_from_clusters_and_data.RData")
head(full_109057[,1:5])

diff_gene_and_volc_plots(full_109057, clusters, "1d. 109057_predict_fromscort_nfkb/")

