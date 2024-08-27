full_145037 <- read.csv("CB_normalized_data/normalised_145037.csv",row.names = 1)

scaled_145037 <- as.data.frame(scale(full_145037))

load("Useful_Functions/Custom_Distance.RData")
rownames(scaled_145037) <- 1:nrow(scaled_145037)
# Fit the pvclust model using the custom distance function
pv_result <- pvclust(t(scaled_145037), method.hclust = "ward.D2", nboot = 1, method.dist = "euclidean")
save(pv_result, file="3c. 145037_Generate_PV_Clusters_and_Top_Gene/pv_object.RData")
# Plot the result
plot(pv_result)
heatmap(as.matrix(scaled_145037), scale = "column")


table(cutree(pv_result$hclust, k=3))
clusters <- cutree(pv_result$hclust, k=3)


load("Useful_Functions/find_means.RData")
nfkb_four_genes <- c("NFKB1", "RELA",  "NFKB2", "RELB")
find_means(scaled_145037[,nfkb_four_genes], clusters = clusters)

## Means indicate Cluster 1 is Canonical, cluster 2 is Atypical and cluster 3 is NonCanonical
clusters <- factor(clusters, levels = c(1,2,3), labels = c("Canonical", "Atypical", "NonCanonical"))
write.csv(clusters, "3c. 145037_Generate_PV_Clusters_and_Top_Gene/cluster_list.csv")

data_matrix <- as.matrix(full_145037)
plot_object <- list(data = data_matrix, cluster = clusters)
plot1 <- fviz_cluster(plot_object, data = data_matrix,
                      geom = "point",
                      ellipse.type = "convex",
                      ggtheme = theme_bw()
)
plot1
ggsave(
  "3c. 145037_Generate_PV_Clusters_and_Top_Gene/NFKB_derived_clustering.png",
  plot = plot1)

load("Useful_Functions/Diff_expression_and_volc_plot_from_clusters_and_data.RData")
## Run volcano plot and gene lists on unscaled data - doesn't have much of an effect
### Always does for canonical vs others
diff_gene_and_volc_plots(full_145037, clusters, "3c. 145037_Generate_PV_Clusters_and_Top_Gene/")


