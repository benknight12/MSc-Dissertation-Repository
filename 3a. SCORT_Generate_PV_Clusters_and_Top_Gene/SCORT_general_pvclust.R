scort_data <- read.csv("CB_normalized_data/scort_normalized.csv",row.names = 1)
dim(scort_data)
scaled_scort <- as.data.frame(scale(scort_data))
pheno <- read.csv("Processed_Data_Sets/Patient data scort.csv", row.names = 1)

library(pvclust)

# custom_dist <- function(x) {
#   # Calculate the correlation matrix
#   cor_matrix <- cor(x, method = "pearson")
#   
#   # Calculate the distance matrix using 1 - r^2
#   dist_matrix <- as.dist(1 - cor_matrix^2)
#   
#   return(dist_matrix)
# }
# save(custom_dist, file = "Useful_Functions/Custom_Distance.RData")
rownames(scaled_scort) <- 1:nrow(scaled_scort)
# Fit the pvclust model using the custom distance function
pv_result <- pvclust(t(scaled_scort), method.hclust = "ward.D2", nboot = 1, method.dist = "euclidean")
save(pv_result, file="3a. SCORT_Generate_PV_Clusters_and_Top_Gene/pv_object.RData")

# Plot the result
plot(pv_result)

table(cutree(pv_result$hclust, k=3))
clusters <- cutree(pv_result$hclust, k=3)
table(clusters)
load("Useful_Functions/find_means.RData")

nfkb_four_genes <- c("NFKB1", "RELA",  "NFKB2", "RELB")

find_means(scaled_scort[,nfkb_four_genes], clusters = clusters)
## Means indicate Cluster 1 is Atpyical, cluster 2 is noncanonical and cluster 3 is canonical
clusters <- factor(clusters, levels = c(1,2,3), labels = c("Atypical","NonCanonical", "Canonical"))
write.csv(clusters, "3a. SCORT_Generate_PV_Clusters_and_Top_Gene/pc_clust_clusters.csv")
table(pheno$NFkB.Cluster, clusters)

data_matrix <- as.matrix(scort_data)
plot_object <- list(data = data_matrix, cluster = clusters)

plot1 <- fviz_cluster(plot_object, data = data_matrix,
                      geom = "point",
                      ellipse.type = "convex",
                      ggtheme = theme_bw()
)
plot1
ggsave(
  "3a. SCORT_Generate_PV_Clusters_and_Top_Gene/PV_SCORT_clustering.png",
  plot = plot1)

load("Useful_Functions/Diff_expression_and_volc_plot_from_clusters_and_data.RData")
diff_gene_and_volc_plots(scort_data, clusters, "3a. SCORT_Generate_PV_Clusters_and_Top_Gene/")
dim(scort_data)

