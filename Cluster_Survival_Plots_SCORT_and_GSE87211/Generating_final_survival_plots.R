make_extreme_binary <- function(df) {
  for (col in colnames(df)) {
    median_value <- median(df[,col], na.rm = TRUE)
    print(col)
    binary_col <- ifelse(df[,col] > 0.5, "Upregulated", ifelse(df[,col] < -0.5, "Downregulated", "Neither"))
    df[[paste0(col, '_binary')]] <- binary_col
  }
  return(df)
}
x <- make_extreme_binary(data.frame(PCSK2 = scaled_scort[,"PCSK2"]))
x
fit_survival_plots <- function(gene, dataset, set = "scort"){
  dataset$gene_factor <- as.factor(as.data.frame(dataset[,gene])[,1])
  dataset <- dataset[dataset$gene_factor != "Neither", ]
  dfs_data <- as.data.frame(dataset[dataset$DFS<40,])
  folder_name <- paste0("8b. ScortSurvival_for_top_genes/", gene)
  if (set != "scort"){
    folder_name <- paste0("8b. ScortSurvival_for_top_genes/", set,"/", gene)
    dfs_data <- as.data.frame(dataset)
  }
  
  dir.create(folder_name, showWarnings = FALSE)
  mod1 <- survfit(Surv(OS, OS.Status=='DECEASED') ~ gene_factor, data = dataset)
  surv_diff <- survdiff(Surv(OS, OS.Status=='DECEASED') ~ gene_factor, data = dataset)
  p_value1 <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
  
  mod2 <- survfit(Surv(DFS, DFS.Status=='Recurred') ~ gene_factor, data = dfs_data)
  surv_diff <- survdiff(Surv(DFS, DFS.Status=='Recurred') ~ gene_factor, data = dfs_data)
  p_value2 <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
  
  out <- list(c())
  if(p_value1 < 1){
    OS_plot <- ggsurvplot(mod1, data = dataset, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = paste0("Kaplan-Meier Curve for Overall Survival ",gene), legend.title = "Gene Expression", legend.labs = c("Downregulated", "Upregulated"))
    print(paste0(folder_name,"/OS_plot.png"))
    png(paste0(folder_name,"/OS_plot.png"))
    print(OS_plot, newpage = FALSE)
    dev.off()
    out <- list(os = OS_plot)
  }
  if(p_value2 < 1){
    DFS_plot <- ggsurvplot(mod2, data = dfs_data, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Disease-Free Survival Probability", title = paste0("Kaplan-Meier Curve for Disease-Free Survival ",gene), legend.title = "Gene Expression",legend.labs = c("Downregulated", "Upregulated"))
    png(paste0(folder_name,"/DFS_plot.png"))
    print(DFS_plot, newpage = FALSE)
    dev.off() 
    out <- list(out, dfs=DFS_plot)
  }
  return(out)
}
expected_gene_expression <- rep("downregulated", length(gene_list))
names(expected_gene_expression) <- gene_list
expected_gene_expression[c("CKS2","RPS17")]<-"upregulated"
make_extreme_binary(scaled_scort[gene_list])
length(gene_list)

df = scaled_scort[,gene_list]
binary_genes <- make_extreme_binary(df = scaled_scort[,gene_list])[(length(gene_list)+1):(2*length(gene_list))]
colnames(binary_genes) <- gene_list
head(binary_genes)
pheno_scort <-read.csv("Processed_Data_Sets/scort.csv", row.names = 1)
survival_scort <- cbind(binary_genes, pheno_scort[,c("OS","OS.Status","DFS","DFS.Status","TRG")])
head(survival_scort)
library(survival)
library(survminer)
levels(as.factor(survival_scort[,"ATP1A2"]))
fit_survival_plots("ATP1A2", survival_scort)
gene_list
for(i in gene_list){
  f <- fit_survival_plots(i, survival_scort)
  print(i)
  print(f$os)
  print(f$dfs)
}

relevantcols <- pheno_87211[,c("disease.free.time..month..ch1", "cancer.recurrance.after.surgery.ch1", "survival.time..month..ch1", "death.due.to.tumor.ch1", "depth.of.invasion.after.rct.ch1")]
df = scaled_87211[,gene_list]

binary_genes_87211 <- make_extreme_binary(df = scaled_87211[complete.cases(relevantcols),gene_list])[,(length(gene_list)+1):(2*(length(gene_list)))]
head(binary_genes_87211)
head(scaled_87211[,"PCSK2"])
survival_gse87211 <- cbind(binary_genes_87211, relevantcols[complete.cases(relevantcols),])
names(survival_gse87211) <- c(gene_list, c("DFS", "DFS.Status", "OS", "OS.Status", "TRG"))

survival_gse87211$DFS.Status <- factor(survival_gse87211$DFS.Status, c(0,1), c("DiseaseFree","Recurred"))
survival_gse87211$OS.Status <- factor(survival_gse87211$OS.Status, c(0,1), c("LIVING","DECEASED"))





for(i in gene_list){
  g <- fit_survival_plots(i, survival_scort)
  f <- fit_survival_plots(i, survival_gse87211, set="87211")
  if (!is.null(g$os) &!is.null(f$os) ){
    print(g$os)
    print(f$os)
    print(i)
  }
  if (!is.null(g$dfs) &!is.null(f$dfs) ){
    print(g$dfs)
    print(f$dfs)
    print(i)
  }
}




