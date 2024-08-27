scaled_scort <- read.csv("CB_normalized_data/scaled_scort.csv", row.names = 1)
pheno_scort <- read.csv("Processed_Data_Sets/scort.csv", row.names = 1)

scaled_87211 <- read.csv("CB_normalized_data/scaled_87211.csv", row.names = 1) 
pheno_87211 <- read.csv("Processed_Data_Sets/pheno87211.csv", row.names = 1)

fit_survival_plots <- function(gene, dataset, set = "scort"){
  dataset$gene_factor <- as.factor(as.data.frame(dataset[,gene])[,1])
  dataset <- dataset[dataset$gene_factor != "neither", ]
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
    OS_plot <- ggsurvplot(mod1, data = dataset, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = paste0("Kaplan-Meier Curve for Overall Survival ",gene), legend.title = "Gene Expression", legend.labs = c("Canonical", "Not Canonical"))
    print(paste0(folder_name,"/OS_plot.png"))
    png(paste0(folder_name,"/OS_plot.png"))
    print(OS_plot, newpage = FALSE)
    dev.off()
    out <- list(os = OS_plot)
  }
  if(p_value2 < 1){
    DFS_plot <- ggsurvplot(mod2, data = dfs_data, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Disease-Free Survival Probability", title = paste0("Kaplan-Meier Curve for Disease-Free Survival ",gene), legend.title = "Gene Expression",legend.labs = c("Canonical", "Not Canonical"))
    png(paste0(folder_name,"/DFS_plot.png"))
    print(DFS_plot, newpage = FALSE)
    dev.off() 
    out <- list(out, dfs=DFS_plot)
  }
  return(out)
}

load("Final_gene_list.Rdata")
gene_list <- similar1

gene_list <-gene_list[!is.na(gene_list)]
load("Useful_Functions/create_binaries.RData")

binary_genes <- create_binary_columns(df = scaled_scort[,gene_list])[,(length(gene_list)+1):(2*(length(gene_list)))]
colnames(binary_genes) <- gene_list
survival_scort <- cbind(binary_genes, pheno_scort[,c("OS","OS.Status","DFS","DFS.Status","TRG")])
head(survival_scort)
dim(survival_scort)
fit_survival_plots("ATP1A2", survival_scort)

survival_scort[survival_scort$OS.Status=="DECEASED",]
pheno_scort[pheno_scort$OS.Status=="DECEASED","OS"]

for(i in gene_list){
  f <- fit_survival_plots(i, as.data.frame(survival_scort))
  print(i)
  print(f$os)
  print(f$dfs)
}
#gene_list <- c( gene_list, "MAD2L1")
relevantcols <- pheno87211[,c("disease.free.time..month..ch1", "cancer.recurrance.after.surgery.ch1", "survival.time..month..ch1", "death.due.to.tumor.ch1", "depth.of.invasion.after.rct.ch1")]
binary_genes_87211 <- create_binary_columns(df = scaled_87211[,gene_list])[,(length(gene_list)+1):(2*(length(gene_list)))]
survival_gse87211 <- cbind(binary_genes_87211, relevantcols)
survival_gse87211 <- survival_gse87211[complete.cases(survival_gse87211),]
names(survival_gse87211) <- c(gene_list, c("DFS", "DFS.Status", "OS", "OS.Status", "TRG"))
names(binary_genes_87211)
binary_genes_87211$MAD2L1_binary
survival_gse87211$DFS.Status <- factor(survival_gse87211$DFS.Status, c(0,1), c("DiseaseFree","Recurred"))
survival_gse87211$OS.Status <- factor(survival_gse87211$OS.Status, c(0,1), c("LIVING","DECEASED"))

x <-fit_survival_plots("ATP1A2", survival_gse87211, set = "87211")


for(i in gene_list){
  f <- fit_survival_plots(i, as.data.frame(survival_gse87211), set = "87211")
  print(i)
  print(f$os)
  print(f$dfs)
}


### Table of TRG proportions
emprty <- as.data.frame(matrix(rep(0,8*length(gene_list)),ncol=4,nrow = 2*length(gene_list)))
rownames(emprty) <- c(paste0(gene_list, "_up"),paste0(gene_list, "_down"))
emprty <- emprty[order(rownames(emprty)),]

for (i in gene_list){
  upset <- cbind(as.data.frame(survival_scort[survival_scort[,i]=="upregulated",i]), as.data.frame(survival_scort[survival_scort[,i]=="upregulated","TRG"]))
  downset <- cbind(as.data.frame(survival_scort[survival_scort[,i]=="downregulated",i]), as.data.frame(survival_scort[survival_scort[,i]=="downregulated","TRG"]))
  names(upset) <- c("cat","TRG")
  names(downset) <- c("cat","TRG")
  upproportions <- c(sum(upset$TRG=="pCR")/length(upset$TRG), sum(upset$TRG=="Good partial response")/length(upset$TRG),sum(upset$TRG=="Partial response")/length(upset$TRG),sum(upset$TRG=="Minimal or no response")/length(upset$TRG))
  downproportions <- c(sum(downset$TRG=="pCR")/length(upset$TRG), sum(downset$TRG=="Good partial response")/length(upset$TRG),sum(downset$TRG=="Partial response")/length(upset$TRG),sum(downset$TRG=="Minimal or no response")/length(upset$TRG))
  emprty[rownames(emprty) == paste0(i,"_up"),]<-upproportions
  emprty[rownames(emprty) == paste0(i,"_down"),]<-downproportions
 }
colnames(emprty) <- c("pCR", "Good partial reponse", "Partial response", "Minimal or no response")
emprty

survival_gse87211$TRG
emprty87211 <- as.data.frame(matrix(rep(0,10*length(gene_list)),ncol=5,nrow = 2*length(gene_list)))
rownames(emprty87211) <- c(paste0(gene_list, "_up"),paste0(gene_list, "_down"))
emprty87211 <- emprty87211[order(rownames(emprty87211)),]
head(survival_gse87211)
survival_gse87211$TRG <- factor(survival_gse87211$TRG, c(0,1,2,3,4), c("pCR","Good partial response", "Partial repsonse", "Minimal response", "No Response"))
gene_list<- gene_list[!is.na(gene_list)]
for (i in gene_list){
  upset <- cbind(as.data.frame(survival_gse87211[survival_gse87211[,i]=="upregulated",i]), as.data.frame(survival_gse87211[survival_gse87211[,i]=="upregulated","TRG"]))
  downset <- cbind(as.data.frame(survival_gse87211[survival_gse87211[,i]=="downregulated",i]), as.data.frame(survival_gse87211[survival_gse87211[,i]=="downregulated","TRG"]))
  names(upset) <- c("cat","TRG")
  names(downset) <- c("cat","TRG")
  upproportions <- c(sum(upset$TRG=="pCR")/length(upset$TRG), sum(upset$TRG=="Good partial response")/length(upset$TRG),sum(upset$TRG=="Partial repsonse")/length(upset$TRG),sum(upset$TRG=="Minimal response")/length(upset$TRG),sum(upset$TRG=="No Response")/length(upset$TRG))
  downproportions <- c(sum(downset$TRG=="pCR")/length(downset$TRG), sum(downset$TRG=="Good partial response")/length(downset$TRG),sum(downset$TRG=="Partial repsonse")/length(downset$TRG),sum(downset$TRG=="Minimal response")/length(downset$TRG),sum(downset$TRG=="No Response")/length(downset$TRG))
  emprty87211[rownames(emprty87211) == paste0(i,"_up"),]<-upproportions
  emprty87211[rownames(emprty87211) == paste0(i,"_down"),]<-downproportions
}
colnames(emprty87211) <- c("pCR", "Good partial reponse", "Partial response", "Minimal response", "No Response")
sum(emprty87211[2,])
emprty87211
write.csv(emprty87211,"~/Downloads/trg_props2.csv")

canset <- as.data.frame(survival_scort[final_clusters=="Canonical",])
nonset <- as.data.frame(survival_scort[final_clusters=="NonCanonical",])
atyset <- as.data.frame(survival_scort[final_clusters=="Atypical",])
head(canset)


### all genes predicting survival
survival_scort
gene_list
## genes all should be down except CKS2 and RPS17 are down (for canonical)
expected_gene_expression <- rep("downregulated", length(gene_list))
names(expected_gene_expression) <- gene_list
expected_gene_expression[c("CKS2","RPS17")]<-"upregulated"
expected_gene_expression
expected <- as.data.frame(matrix(rep(FALSE,length(gene_list)*nrow(survival_scort)), ncol = length(gene_list), nrow = nrow(survival_scort)))
rownames(expected) <- rownames(survival_scort)
colnames(expected) <- gene_list
expected
for (i in gene_list){
  expected[survival_scort[,i]==expected_gene_expression[i],i] <- TRUE
}
survival_scort$indicator <- ifelse(rowSums(as.data.frame(lapply(expected,as.numeric)))>=10, TRUE, FALSE)
OS_plot <- ggsurvplot(survfit(Surv(OS, OS.Status=='DECEASED') ~ indicator, data = survival_scort), data = survival_scort, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival", legend.title = "Gene Expression", legend.labs = c("Less than 10","More Than 9 'Canonical' expressed genes"))
print(OS_plot)
png("15. gene lists from adam/random_os_plot.png")
print(OS_plot)
dev.off()
DFS_plot <- ggsurvplot(survfit(Surv(DFS, DFS.Status=='Recurred') ~ indicator, data = survival_scort), data = survival_scort, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival", legend.title = "Gene Expression", legend.labs = c("Less than 10","More Than 9 'Canonical' expressed genes"))
png("15. gene lists from adam/random_dfs_plot.png")
print(DFS_plot)
dev.off()
### 
### all genes predicting binary therapy response
mod_set <- scaled_scort[,gene_list]
mod_set$binary_response <- as.factor(as.numeric(as.factor(pheno_scort$Binary_response))-1)
formula <- as.formula(paste0("binary_response ~ ", paste(gene_list, collapse = "+")))
bin_mod <- glm(formula, data=mod_set, family = binomial )
summary(bin_mod)
exp(cbind(coef(bin_mod),confint(bin_mod)))
sjPlot::tab_model(bin_mod)

mod_set$clusters <- pheno_scort$NFkB.Cluster=="Canonical"
clus_mod <- bin_mod <- glm(binary_response ~ clusters, data=mod_set, family = binomial )
summary(clus_mod)
exp(cbind(coef(clus_mod),confint(clus_mod)))

mod_set87211 <- scaled_87211[rownames(survival_gse87211),gene_list]
finclus_87211<-read.csv("13. Final_cluster_derivation/final_87211_clusters.csv",row.names = 1)
mod_set87211$clusters <- as.factor(finclus_87211[rownames(survival_gse87211),])
mod_set87211$canonical <- as.factor(finclus_87211[rownames(survival_gse87211),]=="Canonical")
mod_set87211$binary_response <- as.factor(as.numeric(as.factor(pheno_87211[rownames(survival_gse87211),"depth.of.invasion.after.rct.ch1"]<=1))-1)
mod_set87211$binary_response     
clus_mod2 <- glm(binary_response ~ clusters, data=mod_set87211, family = binomial )
summary(clus_mod2)
exp(cbind(coef(clus_mod2),confint(clus_mod2)))
sjPlot::tab_model(clus_mod2)
clus_mod3 <- glm(binary_response ~ canonical, data=mod_set87211, family = binomial )
sjPlot::tab_model(clus_mod3)
names(pheno_87211)


table(final_clusters)
find_means(scaled_scort[,nfkb_four_genes], final_clusters$x)
sum(canset$TRG=="pCR")
length(canset$TRG)
trg_means <- rbind(c(sum(canset$TRG=="pCR")/length(canset$TRG), sum(canset$TRG=="Good partial response")/length(canset$TRG),sum(canset$TRG=="Partial response")/length(canset$TRG),sum(canset$TRG=="Minimal or no response")/length(canset$TRG)),
                   c(sum(nonset$TRG=="pCR")/length(nonset$TRG), sum(nonset$TRG=="Good partial response")/length(nonset$TRG),sum(nonset$TRG=="Partial response")/length(nonset$TRG),sum(nonset$TRG=="Minimal or no response")/length(nonset$TRG)),
                   c(sum(atyset$TRG=="pCR")/length(atyset$TRG), sum(atyset$TRG=="Good partial response")/length(atyset$TRG),sum(atyset$TRG=="Partial response")/length(atyset$TRG),sum(atyset$TRG=="Minimal or no response")/length(atyset$TRG))
                   )
trg_means


load("Useful_Functions/find_means.RData")
find_means(scaled_scort[c(nfkb_four_genes, "CKS2")], final_clusters$x)
write.csv(emprty,"~/Downloads/trg_props.csv")


find_means(scaled_scort[,c(nfkb_four_genes, gene_list[order(gene_list)])], final_clusters$x)


### try produce again but htis time we only do up and down regulation for ><+-0.5
make_extreme_binary <- function(df) {
  for (col in colnames(df)) {
    options <-c("upregulated","downregulated")
    median_value <- median(df[[col]], na.rm = TRUE)
    print(expected_gene_expression[col])
    if (expected_gene_expression[col]=="downregulated"){
      binary_col <- ifelse(df[[col]] > 0.5, "Canonical", "NotCanonical")
    }else if(expected_gene_expression[col]=="upregulated"){
      binary_col <- ifelse(df[[col]] < -0.5, "Canonical", "NotCanonical")
    }
    df[[paste0(col, '_binary')]] <- binary_col
  }
  return(df)
}
make_extreme_binary(scaled_scort[gene_list])
length(gene_list)
gene_list <- c(gene_list, "ELAVL3")
expected_gene_expression
binary_genes <- make_extreme_binary(df = scaled_scort[,gene_list])[(length(gene_list)+1):(2*length(gene_list))]
colnames(binary_genes) <- c(gene_list)
head(binary_genes)
survival_scort <- cbind(binary_genes, pheno_scort[,c("OS","OS.Status","DFS","DFS.Status","TRG")])
head(survival_scort)

levels(as.factor(survival_scort[,"ELAVL3"]))
fit_survival_plots("ELAVL3", survival_scort)
gene_list
for(i in gene_list){
  f <- fit_survival_plots(i, survival_scort)
  print(i)
  print(f$os)
  print(f$dfs)
}

final_clusters <- (read.csv("13. Final_cluster_derivation/final_scort_clusters.csv", row.names = 1))$x
final_clusters <- pheno_scort$NFkB.Cluster
data <- data.frame(
  value = c(scaled_scort[final_clusters=="Canonical","RELA"], scaled_scort[final_clusters=="NonCanonical","RELA"], scaled_scort[final_clusters=="Atypical","RELA"]),
  group = factor(c(rep("Canonical", length(scaled_scort[final_clusters=="Canonical","RELA"])),
                   rep("NonCanonical", length(scaled_scort[final_clusters=="NonCanonical","RELA"])),
                   rep("Atypical", length(scaled_scort[final_clusters=="Atypical","RELA"])) ))
)
ggplot(data, aes(x = group, y = value)) +
  geom_violin(trim = FALSE, fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(title = "Violin Plots of Three Vectors",
       x = "Group",
       y = "Value")



relevantcols <- pheno87211[,c("disease.free.time..month..ch1", "cancer.recurrance.after.surgery.ch1", "survival.time..month..ch1", "death.due.to.tumor.ch1", "depth.of.invasion.after.rct.ch1")]
binary_genes_87211 <- make_extreme_binary(df = scaled_87211[,gene_list])[,(length(gene_list)+1):(2*(length(gene_list)))]
survival_gse87211 <- cbind(binary_genes_87211, relevantcols)
survival_gse87211 <- survival_gse87211[complete.cases(survival_gse87211),]
names(survival_gse87211) <- c(gene_list, c("DFS", "DFS.Status", "OS", "OS.Status", "TRG"))
head(survival_gse87211)

survival_gse87211$DFS.Status <- factor(survival_gse87211$DFS.Status, c(0,1), c("DiseaseFree","Recurred"))
survival_gse87211$OS.Status <- factor(survival_gse87211$OS.Status, c(0,1), c("LIVING","DECEASED"))

load("Final_gene_list.Rdata")
gene_list <- similar1
gene_list <- gene_list[!is.na(gene_list)]

for(i in gene_list){
  f <- fit_survival_plots(i, survival_gse87211, set="87211")
  print(i)
  print(f$os)
  print(f$dfs)
}
gene_list

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

