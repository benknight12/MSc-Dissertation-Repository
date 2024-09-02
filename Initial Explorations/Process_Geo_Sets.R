library(oligo)
gse109057_celdata <- read.celfiles(list.celfiles('../9_07_Folder_Reset/RAW_Files/GSE109057/GSE109057_RAW/',full.names=TRUE,listGzipped = TRUE))
gse145037_celdata <- read.celfiles(list.celfiles('../9_07_Folder_Reset/RAW_Files/GSE145037/GSE145037_RAW/',full.names=TRUE,listGzipped = TRUE))                                  


make_gene_data_frame <- function(normalized,package){
  
  if (class(normalized)=="ExpressionSet"){
    gene_data <- exprs(normalized)
  }else{
    gene_data <- normalized
  }
  
  probe_ids <- rownames(gene_data)
  gene_symbols <-mapIds(package, keys = probe_ids, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
  rownames(gene_data) <- gene_symbols
  
  
  unique_rowID <- make.unique(rownames(gene_data))
  
  # Check for missing data in rowGroup and remove rows with missing data
  valid_rows <- !is.na(rownames(gene_data))
  gene_data <- gene_data[valid_rows, ]
  unique_rowID <- unique_rowID[valid_rows]
  library(WGCNA)
  collapsed <- collapseRows(
    datET = gene_data,
    rowGroup = rownames(gene_data),
    rowID = unique_rowID,
    method = "Average"
  )
  
  collapsed_matrix <- collapsed$datETcollapsed
  
  
  print(dim(collapsed_matrix))
  
  return(as.data.frame(t(collapsed_matrix)))
}
save(make_gene_data_frame, file = "Useful_Functions/process_data_function.RData")

library(primeview.db)
normalized145037 <- rma(gse145037_celdata)
processed145037 <-  make_gene_data_frame(normalized145037, primeview.db)

normalized109057 <- rma(gse109057_celdata)
processed109057 <-  make_gene_data_frame(normalized109057, primeview.db)

write.csv(processed109057, "Processed_Data_Sets/109057_set(unscaled).csv")
write.csv(processed145037, "Processed_Data_Sets/145037_set(unscaled).csv")