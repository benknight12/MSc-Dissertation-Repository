diff_gene_and_volc_plots2 <- function(data_frame, clusters, folder_name, lfc_boundary = log(1.5, base=2), pval_boundary=0.05, special_labels =NULL, overide = FALSE){
  data_matrix <- as.matrix(data_frame)
  # Perform differential expression analysis using limma
  design <- model.matrix(~0 + clusters)
  colnames(design) <- levels(clusters)
  
  fit <- lmFit(t(data_matrix), design)
  # Create contrasts to compare each cluster against the others
  contrast.matrix <- makeContrasts(
    Canonical_vs_NonCanonical = Canonical - NonCanonical,
    NonCanonical_vs_Atypical = NonCanonical - Atypical,
    Atypical_vs_Canonical = Atypical - Canonical,
    levels = design
  )
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  # Extract top genes for each cluster
  top_genes_cluster1 <- (topTable(fit2,lfc = lfc_boundary,p.value = 0.05, coef = "Canonical_vs_NonCanonical", number = Inf, adjust = "BH"))
  top_genes_cluster2 <- (topTable(fit2,lfc = lfc_boundary,p.value = 0.05, coef = "NonCanonical_vs_Atypical", number = Inf, adjust = "BH"))
  top_genes_cluster3 <- (topTable(fit2,lfc = lfc_boundary,p.value = 0.05, coef = "Atypical_vs_Canonical", number = Inf, adjust = "BH"))
  
  max_length <- max(length(rownames(top_genes_cluster1)), length(rownames(top_genes_cluster2)), length(rownames(top_genes_cluster3)))
  
  top_gene_dataframe <- data.frame(
    CanonicalvsNonCanonical = c(rownames(top_genes_cluster1), rep(NA, max_length - length(rownames(top_genes_cluster1)))),
    NonCanonicalvsAtypical = c(rownames(top_genes_cluster2), rep(NA, max_length - length(rownames(top_genes_cluster2)))),
    AtypicalvsCanonical = c(rownames(top_genes_cluster3), rep(NA, max_length - length(rownames(top_genes_cluster3))))
  )
  print(top_gene_dataframe)
  
  write.csv(top_gene_dataframe, paste0(folder_name,"/Top_Gene_List.csv"))
  
  
  results <- topTable(fit2, coef = "Canonical_vs_NonCanonical", number = Inf)
  results2 <- topTable(fit2, coef = "Atypical_vs_Canonical", number = Inf)
  
  # Ensure the results are in a data frame format
  results_df <- as.data.frame(results)
  results_df2 <- as.data.frame(results2)
  
  if (is.null(special_labels)){
    volcano_plot <- EnhancedVolcano(
      results_df, lab = rownames(results), x = 'logFC', y = 'adj.P.Val', xlab = bquote(~Log[2]~ "fold change"),
      ylab = bquote(~-Log[10]~ "adjusted p-value"), pCutoff = 0.05, FCcutoff = lfc_boundary,  title = 'Volcano plot of Differentially Expressed Features', 
      subtitle = 'Adjusted p-value < 0.05 and |log2 FC| > 0.585', pointSize = 1.0, labSize = 3.0, col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
      colAlpha = 1, legendPosition = 'top', legendLabSize = 10, legendIconSize = 3.0, drawConnectors = TRUE, widthConnectors = 0.5,
      gridlines.major = TRUE, gridlines.minor = TRUE,max.overlaps = 50)
    volcano_plot2 <- EnhancedVolcano(
      results_df2, lab = rownames(results2), x = 'logFC', y = 'adj.P.Val', xlab = bquote(~Log[2]~ "fold change"),
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
      special_labels2<-intersect(rownames(results_df2[results_df2$logFC>=1 & results_df2$adj.P.Val<0.05,]),special_labels)
      volcano_plot <- EnhancedVolcano(
        results_df2, lab = rownames(results2), x = 'logFC', y = 'adj.P.Val', xlab = bquote(~Log[2]~ "fold change"),
        ylab = bquote(~-Log[10]~ "adjusted p-value"), pCutoff = 0.05, FCcutoff = lfc_boundary,  title = 'Volcano plot of Differentially Expressed Features', 
        subtitle = 'Adjusted p-value < 0.05 and |log2 FC| > 0.585', pointSize = 1.0, labSize = 3.0, col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
        colAlpha = 1, legendPosition = 'top', legendLabSize = 10, legendIconSize = 3.0, drawConnectors = TRUE, widthConnectors = 0.5,
        gridlines.major = TRUE, gridlines.minor = TRUE,max.overlaps = 50, selectLab = special_labels)
      print(special_labels2)
    } else {
      volcano_plot <- EnhancedVolcano(
        results_df, lab = rownames(results), x = 'logFC', y = 'adj.P.Val', xlab = bquote(~Log[2]~ "fold change"),
        ylab = bquote(~-Log[10]~ "adjusted p-value"), pCutoff = 0.05, FCcutoff = lfc_boundary,  title = 'Volcano plot of Differentially Expressed Features', 
        subtitle = 'Adjusted p-value < 0.05 and |log2 FC| > 0.585', pointSize = 1.0, labSize = 3.0, col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
        colAlpha = 1, legendPosition = 'top', legendLabSize = 10, legendIconSize = 3.0, drawConnectors = TRUE, widthConnectors = 0.5,
        gridlines.major = TRUE, gridlines.minor = TRUE,max.overlaps = 400, selectLab = special_labels)
      volcano_plot2 <- EnhancedVolcano(
        results_df2, lab = rownames(results2), x = 'logFC', y = 'adj.P.Val', xlab = bquote(~Log[2]~ "fold change"),
        ylab = bquote(~-Log[10]~ "adjusted p-value"), pCutoff = 0.05, FCcutoff = lfc_boundary,  title = 'Volcano plot of Differentially Expressed Features', 
        subtitle = 'Adjusted p-value < 0.05 and |log2 FC| > 0.585', pointSize = 1.0, labSize = 3.0, col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
        colAlpha = 1, legendPosition = 'top', legendLabSize = 10, legendIconSize = 3.0, drawConnectors = TRUE, widthConnectors = 0.5,
        gridlines.major = TRUE, gridlines.minor = TRUE,max.overlaps = 400, selectLab = special_labels)
    }
  }
  print(volcano_plot)
  print(volcano_plot2)
  print(fit2$coefficients[nfkb_four_genes,"Canonical_vs_NonCanonical"]/log(2,exp(1)))
  print(results[nfkb_four_genes, "logFC"])
  ggsave(
    paste0(folder_name,"/volcano_plot1.png"),
    plot = volcano_plot)
  ggsave(
    paste0(folder_name,"/volcano_plot2.png"),
    plot = volcano_plot2)
}

find_means(scaled_scort[,nfkb_four_genes], pred_from_perf)
clusters <- as.factor(final_clusters)
clusters <- factor(pred_from_perf, levels = c(1,2,3), c("Canonical","NonCanonical","Atypical"))
diff_gene_and_volc_plots(scort_data, clusters, "12. gene list comp/109057/", lfc_boundary = log(2,base=2),special_labels = nfkb_four_genes, overide = T)

save(diff_gene_and_volc_plots2, file = "Useful_Functions/Diff_expression_and_volc_plot_from_clusters_and_data2.RData")
