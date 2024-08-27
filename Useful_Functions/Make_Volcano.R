make_volcano <- function(processed_data, clusters, cluster_names){
  data_matrix <- as.matrix(processed_data)
  
  ## 'Renormalise' dataframe with scale
  processed_data <- as.data.frame(scale(processed_data))
  
  print(dim(processed_data))
  ## Add column to dataframe with the cluster labels
  processed_data$cluster <- factor(clusters, labels = cluster_names)
  
  library(limma)
  
  design <- model.matrix(~0 + processed_data$cluster)
  colnames(design) <- levels(processed_data$cluster)
  fit <- lmFit(t(data_matrix), design)
  
  contrast.matrix <- makeContrasts(
    Canonical_vs_Others = Canonical - (NonCanonical + Atypical)/2,
    NonCanonical_vs_Others = NonCanonical - (Canonical + Atypical)/2,
    Atypical_vs_Others = Atypical - (Canonical + NonCanonical)/2,
    levels = design
  )
  
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  results <- topTable(fit2, coef = "Canonical_vs_Others", number = Inf)
  
  print(head(results))
  
  # Add a column for significance
  results$Significance <- ifelse(results$adj.P.Val < 0.05 & abs(results$logFC) > 1, "Significant", "Not Significant")
  print(results[results$Significance=="Significant",])
  library(EnhancedVolcano)
  volc1 <- EnhancedVolcano(
    results,
    lab = rownames(results),
    x = 'logFC',
    y = 'adj.P.Val',
    xlab = bquote(~Log[2]~ "fold change"),
    ylab = bquote(~-Log[10]~ "adjusted p-value"),
    pCutoff = 0.05,
    FCcutoff = 0.9,
    title = 'Volcano plot of Differentially Expressed Features',
    subtitle = 'Adjusted p-value < 0.05 and |log2 FC| > 1',
    pointSize = 1.0,
    labSize = 3.0,
    col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
    colAlpha = 1,
    legendPosition = 'top',
    legendLabSize = 10,
    legendIconSize = 3.0,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    gridlines.major = TRUE,
    gridlines.minor = TRUE
  )
  plot(volc1)
  
  Canon_NonCanon <- topTable(fit2, coef = "Canonical_vs_Others", number = Inf)
  C_NC_lfc <- topTable(fit2, coef = "Canonical_vs_Others", number = Inf, lfc = 1)
  
  library(tidyr)
  library(dplyr)
  library(ggrepel)
  
  volc2 <- Canon_NonCanon %>%
    mutate(DE = abs(logFC) > 1 & adj.P.Val < 0.05) %>%
    ggplot(aes(x = logFC, y = -log10(P.Value))) +
    geom_point(aes(colour = DE), alpha = 1, size = 1) +
    scale_colour_manual(values = c("black", "blue")) +
    geom_vline(xintercept = c(1, -1), colour = "blue", linetype = 2) +
    geom_hline(yintercept = 1.3, colour = "blue", linetype = 2) +
    theme_bw()
  
  plot(volc2)
  
  return(list(volc1 = volc1, volc2 = volc2, significiant = results[results$Significance=="Significant",]))
}
save(make_volcano, file = "Useful_Functions/make_volcano.RData")
