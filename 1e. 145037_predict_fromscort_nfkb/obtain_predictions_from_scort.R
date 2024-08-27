full_145037 <- read.csv("Processed_Data_Sets/145037_set(unscaled).csv", row.names = 1)
scaled_145037 <- as.data.frame(scale(full_145037))
load("2a. SCORT_Generate_Kmeans_Clusters_and_Top_Gene/kmeans_object.RData")
load("1a. SCORT_NFKB_cluster_derivation_replicating/My_NFKB_Scort_kmeans.RData")
my_nfkb_kmeans_scort <- kmeans_scort


my_nfkb_kmeans_scort$centers
nfkb_four_genes <- c("NFKB1","RELA","NFKB2","RELB")


load("Useful_Functions/Predict_Kmeans.RData")
new_clus <- predict_kmeans(scaled_145037[,nfkb_four_genes], my_nfkb_kmeans_scort)
find_means(data = scaled_145037[,nfkb_four_genes], new_clus)

## 1= canonical, 2=noncanonical, 3=atypical
clusters <- factor(new_clus, levels = c(1,2,3), labels = c("Canonical", "NonCanonical", "Atypical"))


load("Useful_Functions/Diff_expression_and_volc_plot_from_clusters_and_data.RData")


diff_gene_and_volc_plots(full_145037, clusters, "1e. 145037_predict_fromscort_nfkb/")

