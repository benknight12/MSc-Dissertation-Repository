library(EnhancedVolcano)
full_109057 <- read.csv("CB_normalized_data/normalised_109057.csv",  row.names = 1)
scaled_109057 <- as.data.frame(scale(full_109057))

nfkb_four_genes <- c("NFKB1","RELA","NFKB2","RELB")

nfkb_109057 <- scaled_109057[,nfkb_four_genes]
head(nfkb_109057)
summary(nfkb_109057)

## This divides up the samples into those with expression values that very strongly indicate cluster membership and then from these high likelihood assignments find centroid points by taking their means
label_data <- function(df) {
  labels <- rep("messy", nrow(df))
  
  # Cluster 1: Upregulated NFKB1 and RELA, Downregulated NFKB2 and RELB
  labels[df$NFKB1 > 0 & df$RELA > 0 & df$NFKB2 < 0 & df$RELB < 0] <- "up_NFKB1_RELA_down_NFKB2_RELB"
  
  # Cluster 2: Downregulated NFKB1 and RELA, Upregulated NFKB2 and RELB
  labels[df$NFKB1 < 0 & df$RELA < 0 & df$NFKB2 > 0 & df$RELB > 0] <- "down_NFKB1_RELA_up_NFKB2_RELB"
  
  return(labels)
}
centroids <- rbind(colMeans(nfkb_109057[label_data(nfkb_109057)=="up_NFKB1_RELA_down_NFKB2_RELB",]), colMeans(nfkb_109057[label_data(nfkb_109057)=="down_NFKB1_RELA_up_NFKB2_RELB",]),colMeans(nfkb_109057[label_data(nfkb_109057)=="messy",]))


kmeans_109057_nfkb <- kmeans(nfkb_109057, centers=centroids, nstart = 30)
save(kmeans_109057_nfkb, file = "1b. 109057_NFKB_cluster_derivation/kmeans_object.RData")
load("Useful_Functions/find_means.RData")
find_means(nfkb_109057, kmeans_109057_nfkb$cluster)
table(kmeans_109057$cluster)

### The means suggest cluster 1 is Canonical, cluster 2 is NonCanonical and cluster 3 is Atypical
nfkb_derived_clusters_109057 <- factor(kmeans_109057$cluster, levels=c(1,2,3), labels = c("Canonical", "NonCanonical", "Atypical"))

## Save the clusters
write.csv(nfkb_derived_clusters_109057, "1b. 109057_NFKB_cluster_derivation/NFKB_derived_scort_clusters.csv")





## Plot the clusters
data_matrix <- as.matrix(full_109057)
cluster_object <- list(data = data_matrix, cluster = nfkb_derived_clusters_109057)
plot1 <- fviz_cluster(cluster_object, data = data_matrix,
                      geom = "point",
                      ellipse.type = "convex",
                      ggtheme = theme_bw()
)
plot1
ggsave(
  "1b. 109057_NFKB_cluster_derivation/NFKB_derived_clustering.png",
  plot = plot1)

## Extract differential gene lists

## Define a function for taking clusters and full data then performing differentical gene analysis and rpoduce volcano plots
## (clusters must be Canonical, NonCanonical ant Atypical - in this order)
diff_gene_and_volc_plots <- function(data_frame, clusters, folder_name, lfc_boundary = log(1.5, base=2), pval_boundary=0.05, special_labels =NULL, overide = FALSE){
  data_matrix <- as.matrix(data_frame)
  # Perform differential expression analysis using limma
  design <- model.matrix(~0 + clusters)
  colnames(design) <- levels(clusters)
  
  fit <- lmFit(t(data_matrix), design)
  # Create contrasts to compare each cluster against the others
  contrast.matrix <- makeContrasts(
    Canonical_vs_Others = Canonical - (NonCanonical + Atypical)/2,
    NonCanonical_vs_Others = NonCanonical - (Canonical + Atypical)/2,
    Atypical_vs_Others = Atypical - (Canonical + NonCanonical)/2,
    levels = design
  )
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  print(fit2$coefficients[nfkb_four_genes,"Canonical_vs_Others"]/log(2,exp(1)))
  # Extract top genes for each cluster
  top_genes_cluster1 <- (topTable(fit2,lfc = lfc_boundary,p.value = 0.05, coef = "Canonical_vs_Others", number = Inf, adjust = "BH"))
  top_genes_cluster2 <- (topTable(fit2,lfc = lfc_boundary,p.value = 0.05, coef = "NonCanonical_vs_Others", number = Inf, adjust = "BH"))
  top_genes_cluster3 <- (topTable(fit2,lfc = lfc_boundary,p.value = 0.05, coef = "Atypical_vs_Others", number = Inf, adjust = "BH"))

  max_length <- max(length(rownames(top_genes_cluster1)), length(rownames(top_genes_cluster2)), length(rownames(top_genes_cluster3)))

  top_gene_dataframe <- data.frame(
    Canonical = c(rownames(top_genes_cluster1), rep(NA, max_length - length(rownames(top_genes_cluster1)))),
    NonCanonical = c(rownames(top_genes_cluster2), rep(NA, max_length - length(rownames(top_genes_cluster2)))),
    Atypical = c(rownames(top_genes_cluster3), rep(NA, max_length - length(rownames(top_genes_cluster3))))
  )
  print(top_gene_dataframe)
  
  write.csv(top_gene_dataframe, paste0(folder_name,"/Top_Gene_List.csv"))
  
  
  results <- topTable(fit2, coef = "Canonical_vs_Others", number = Inf)
  
  # Ensure the results are in a data frame format
  results_df <- as.data.frame(results)
  
  if (is.null(special_labels)){
   volcano_plot <- EnhancedVolcano(
    results_df, lab = rownames(results), x = 'logFC', y = 'adj.P.Val', xlab = bquote(~Log[2]~ "fold change"),
    ylab = bquote(~-Log[10]~ "adjusted p-value"), pCutoff = 0.05, FCcutoff = lfc_boundary,  title = 'Volcano plot of Differentially Expressed Features', 
    subtitle = 'Adjusted p-value < 0.05 and |log2 FC| > 0.585', pointSize = 1.0, labSize = 3.0, col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
    colAlpha = 1, legendPosition = 'top', legendLabSize = 10, legendIconSize = 3.0, drawConnectors = TRUE, widthConnectors = 0.5,
    gridlines.major = TRUE, gridlines.minor = TRUE,max.overlaps = 50)
  } else {
    if (overide == FALSE){
    special_labels<-intersect(rownames(results_df[results_df$logFC>=1 & results_df$adj.P.Val<0.05,]),special_labels)
    volcano_plot <- EnhancedVolcano(
      results_df, lab = rownames(results), x = 'logFC', y = 'adj.P.Val', xlab = bquote(~Log[2]~ "fold change"),
      ylab = bquote(~-Log[10]~ "adjusted p-value"), pCutoff = 0.05, FCcutoff = lfc_boundary,  title = 'Volcano plot of Differentially Expressed Features', 
      subtitle = 'Adjusted p-value < 0.05 and |log2 FC| > 0.585', pointSize = 1.0, labSize = 3.0, col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
      colAlpha = 1, legendPosition = 'top', legendLabSize = 10, legendIconSize = 3.0, drawConnectors = TRUE, widthConnectors = 0.5,
      gridlines.major = TRUE, gridlines.minor = TRUE,max.overlaps = 50, selectLab = special_labels)
    print(special_labels)
    } else {
      volcano_plot <- EnhancedVolcano(
        results_df, lab = rownames(results), x = 'logFC', y = 'adj.P.Val', xlab = bquote(~Log[2]~ "fold change"),
        ylab = bquote(~-Log[10]~ "adjusted p-value"), pCutoff = 0.05, FCcutoff = lfc_boundary,  title = 'Volcano plot of Differentially Expressed Features', 
        subtitle = 'Adjusted p-value < 0.05 and |log2 FC| > 0.585', pointSize = 1.0, labSize = 3.0, col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
        colAlpha = 1, legendPosition = 'top', legendLabSize = 10, legendIconSize = 3.0, drawConnectors = TRUE, widthConnectors = 0.5,
        gridlines.major = TRUE, gridlines.minor = TRUE,max.overlaps = 400, selectLab = special_labels)
    }
  }
  print(volcano_plot)
  
  
  print(fit2$coefficients[nfkb_four_genes,"Canonical_vs_Others"]/log(2,exp(1)))
  print(results[nfkb_four_genes, "logFC"])
  
  ggsave(
    paste0(folder_name,"/volcano_plot1.png"),
    plot = volcano_plot)
}

diff_gene_and_volc_plots(full_109057, nfkb_derived_clusters_109057, "1b. 109057_NFKB_cluster_derivation")

save(diff_gene_and_volc_plots, file = "Useful_Functions/Diff_expression_and_volc_plot_from_clusters_and_data.RData")
