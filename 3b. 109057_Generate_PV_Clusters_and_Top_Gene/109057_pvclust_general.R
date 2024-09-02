full_109057 <- read.csv("CB_normalized_data/normalised_109057.csv",row.names = 1)

scaled_109057 <- as.data.frame(scale(full_109057))

load("Useful_Functions/Custom_Distance.RData")
rownames(scaled_109057) <- 1:nrow(scaled_109057)
# Fit the pvclust model using the custom distance function
pv_result <- pvclust(t(scaled_109057), method.hclust = "ward.D2", nboot = 1, method.dist = "euclidean")
save(pv_result, file="3b. 109057_Generate_PV_Clusters_and_Top_Gene/pv_object.RData")


# Plot the result
plot(pv_result)

table(cutree(pv_result$hclust, k=3))
clusters <- cutree(pv_result$hclust, k=3)

load("3b. 109057_Generate_PV_Clusters_and_Top_Gene/pv_object.RData")

load("Useful_Functions/find_means.RData")
nfkb_four_genes <- c("NFKB1", "RELA",  "NFKB2", "RELB")
find_means(scaled_109057[,nfkb_four_genes], clusters = clusters)

## Means indicate Cluster 1 is Atypical, cluster 2 is Canonical and cluster 3 is NonCanonical
clusters <- factor(clusters, levels = c(1,2,3), labels = c("Atypical", "Canonical", "NonCanonical"))
write.csv(clusters, "3b. 109057_Generate_PV_Clusters_and_Top_Gene/cluster_list.csv")

data_matrix <- as.matrix(full_109057)
plot_object <- list(data = data_matrix, cluster = clusters)
plot1 <- fviz_cluster(plot_object, data = data_matrix,
                      geom = "point",
                      ellipse.type = "convex",
                      ggtheme = theme_bw()
)
plot1
ggsave(
  "3b. 109057_Generate_PV_Clusters_and_Top_Gene/NFKB_derived_clustering.png",
  plot = plot1)

load("Useful_Functions/Diff_expression_and_volc_plot_from_clusters_and_data.RData")
## Run volcano plot and gene lists on unscaled data - doesn't have much of an effect
diff_gene_and_volc_plots(full_109057, clusters, "3b. 109057_Generate_PV_Clusters_and_Top_Gene/")


