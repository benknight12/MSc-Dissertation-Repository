## 3. Each of top-top genes conduct cox proportional hazards survival between up and downregulated for each of these genes
scaled_scort <- as.data.frame(scale(read.csv("Processed_Data_Sets/SCORT_final_data.csv",row.names = 1)))
pheno <-read.csv("Processed_Data_Sets/Patient data scort.csv")

## SCORT Canonical
load("12. gene list comp/Final_Gene_Lists.RData")

print(unique(c(final_lists_all_sets$tier1, final_lists_all_sets$tier2, final_lists_all_sets$tier3, final_lists_all_sets$tier4)))
gene_list <- final_lists_all_sets$tier1
gene_list <- gene_list[!is.na(gene_list)]

survival_gene_df <- scaled_scort[,c(gene_list,"ATP1A2")]
head(survival_gene_df)
colMeans(survival_gene_df) # means basically 0 as scaled so use 0 as split point for high/low
library(robustbase)
colMedians(as.matrix(survival_gene_df)) # medians fairly close to 0 as well - close enough maybe just use 0

# Create binarys
create_binary_columns <- function(df) {
  for (col in colnames(df)) {
    median_value <- median(df[[col]], na.rm = TRUE)
    binary_col <- ifelse(df[[col]] > median_value, 'upregulated', 'downregulated')
    df[[paste0(col, '_binary')]] <- binary_col
  }
  return(df)
}
save(create_binary_columns, file = "Useful_Functions/create_binaries.RData")

survival_gene_df <- cbind(create_binary_columns(survival_gene_df), OS = pheno$OS, OS.Status = pheno$OS.Status, DFS = pheno$DFS, DFS.Status = pheno$DFS.Status, TRG = pheno$TRG, binary_response = pheno$Binary_response)
head(survival_gene_df)
library(survival)
library(survminer)
## See that all NA's are 36 and all 36 are either NA or Living so we can impute NA's as living - must be end of collection period for most individuals
survival_gene_df$OS[is.na(survival_gene_df$OS.Status)]
survival_gene_df$OS.Status[survival_gene_df$OS == 36]
survival_gene_df$OS[is.na(survival_gene_df$DFS.Status)]

survival_gene_df$OS.Status[is.na(survival_gene_df$OS.Status)] <- "LIVING"



dfs_data <-survival_gene_df[survival_gene_df$DFS<40,]
levels(as.factor(survival_gene_df_gse$DFS.Status))
fit_survival_plots <- function(gene, dataset, set = "scort"){
  gene_column <- as.data.frame(dataset[,gene])
  dataset$gene_factor <- as.factor(gene_column[,1])
  dfs_data <- as.data.frame(dataset[dataset$DFS<40,])
  folder_name <- paste0("8b. ScortSurvival_for_top_genes/", gene)
  if (set != "scort"){
    folder_name <- paste0("8b. ScortSurvival_for_top_genes/", set,"/", gene)
    dfs_data <- as.data.frame(dataset)
  }
  dir.create(folder_name, showWarnings = FALSE)
  OS_plot <- ggsurvplot(survfit(Surv(OS, OS.Status=='DECEASED') ~ gene_factor, data = dataset), data = dataset, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = paste0("Kaplan-Meier Curve for Overall Survival ",gene), legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
  DFS_plot <- ggsurvplot(survfit(Surv(DFS, DFS.Status=='Recurred') ~ gene_factor, data = dfs_data), data = dfs_data, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = paste0("Kaplan-Meier Curve for Disease-Free Survival ",gene), legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
  png(paste0(folder_name,"/OS_plot.png"))
  print(OS_plot, newpage = FALSE)
  dev.off()
  png(paste0(folder_name,"/DFS_plot.png"))
  print(DFS_plot, newpage = FALSE)
  dev.off()  
  return(list(os = OS_plot, dfs = DFS_plot))
  }

load("Final_gene_list.Rdata")
gene_list <- similar1

scaled_scort <- read.csv("CB_normalized_data/scaled_scort.csv",row.names = 1)
pheno_scort <- read.csv("Processed_Data_Sets/scort.csv",row.names = 1)
names(pheno_scort)
pheno_scort$OS[pheno_scort$OS.Status=="DECEASED"]
pheno_scort$OS.Status[is.na(pheno_scort$OS.Status)] <- "LIVING"
#complete.cases(pheno_scort[,c("DFS", "DFS.Status", "OS", "OS.Status", "TRG")])
pheno_scort$DFS.Status
small <- scaled_scort[,gene_list[!is.na(gene_list)]]
interim <- cbind(create_binary_columns(small), pheno_scort[,c("DFS", "DFS.Status", "OS", "OS.Status", "TRG")])
survival_gene_df <- interim[, (ncol(small)+1):ncol(interim)]
names(survival_gene_df) <- c(names(small),c("DFS", "DFS.Status", "OS", "OS.Status", "TRG"))
head(interim)

fit_survival_plots("ATP1A2", as.data.frame(survival_gene_df))

for(i in gene_list[!is.na(gene_list)]){
  f <- fit_survival_plots(i, as.data.frame(survival_gene_df))
  print(i)
  print(f$os)
  print(f$dfs)
}


### 87211

scaled_87211 <- read.csv("CB_normalized_data/scaled_87211.csv",row.names = 1)
pheno_87211 <- read.csv("Processed_Data_Sets/pheno87211.csv", row.names = 1)
names(pheno_87211)

table(is.na(pheno_87211$cancer.recurrance.after.surgery.ch1))
pheno_87211$cancer.recurrance.after.surgery.ch1
pheno_87211$disease.free.time..month..ch1

small <- scaled_87211[,gene_list[!is.na(gene_list)]]
interim <- cbind(create_binary_columns(small), pheno_87211[,c("disease.free.time..month..ch1", "cancer.recurrance.after.surgery.ch1", "survival.time..month..ch1", "death.due.to.tumor.ch1", "depth.of.invasion.after.rct.ch1")])
survival_gene_df_gse <- interim[, (ncol(small)+1):ncol(interim)]
names(survival_gene_df_gse) <- c(names(small),c("DFS", "DFS.Status", "OS", "OS.Status", "TRG"))
survival_gene_df_gse
survival_gene_df_gse <- survival_gene_df_gse[complete.cases(survival_gene_df_gse),]


for(i in gene_list[!is.na(gene_list)]){
  f <- fit_survival_plots(i, as.data.frame(survival_gene_df_gse), set="87211")
  print(i)
  print(f$os)
  print(f$dfs)
}


fit_survival_plots("CKS2", survival_gene_df_gse, set="87211")


survival_gene_df_gse$TRG
boxplot(list(survival_gene_df_gse$TRG[survival_gene_df_gse$PKNOX2=="upregulated"],
             survival_gene_df_gse$TRG[survival_gene_df_gse$PKNOX2=="downregulated"]))

df <- data.frame(
  value = survival_gene_df_gse$TRG,
  group = survival_gene_df_gse$PSD)
df
df <- df %>%
  group_by(value, group)%>%
  summarise(count = n()) %>%
  ungroup()

ggplot(df, aes(x=value, y=count, group = group, color=as.factor(group)))+
  geom_line(position = "dodge")



### create dataframe of mean survival for up and down regulation of genes
up <- c()
down <-c()
for (j in similar[!is.na(similar)]){
  up <- c(up, mean(survival_gene_df$OS[survival_gene_df[,j]=="upregulated"]))
  down <- c(down, mean(survival_gene_df$OS[survival_gene_df[,j]=="downregulated"]))
}


means <- as.data.frame(t(as.matrix(data.frame(up = up,down=down))))
colnames(means) <- similar[!is.na(similar)]

means
final_clusters <- read.csv("13. Final_cluster_derivation/final_scort_clusters.csv",row.names = 1)
c(nfkb_four_genes,!is.na(similar))
find_means(scaled_scort[,c(nfkb_four_genes,similar[!is.na(similar)])], final_clusters$x)
## this tells us if up or down reg for canonical (all down except 18&19):
upordown <- rep("D", 22)
upordown[c(18,19)] <- "U"
rbind(means, upordown)

pheno$TRG[final_clusters =="Canonical"]
mean(as.numeric(as.factor(pheno$TRG[final_clusters =="Canonical"])))
mean(as.numeric(as.factor(pheno$TRG[final_clusters =="NonCanonical"])))
mean(as.numeric(as.factor(pheno$TRG[final_clusters =="Atypical"])))

pheno_87211 <-read.csv("Processed_Data_Sets/pheno87211.csv", row.names = 1)

finclus_87211 <-finclus_87211[!is.na(pheno_87211$depth.of.invasion.after.rct.ch1)]
pheno_87211 <- pheno_87211[!is.na(pheno_87211$depth.of.invasion.after.rct.ch1),]

mean(as.numeric(as.factor(pheno_87211$depth.of.invasion.after.rct.ch1[finclus_87211 =="Canonical"])))
mean(as.numeric(as.factor(pheno_87211$depth.of.invasion.after.rct.ch1[finclus_87211 =="NonCanonical"])))
mean(as.numeric(as.factor(pheno_87211$depth.of.invasion.after.rct.ch1[finclus_87211 =="Atypical"])))

pheno_data <- data.frame(
  depth_of_invasion = pheno_87211$depth.of.invasion.after.rct.ch1,
  group = finclus_87211
)
pheno_data_filtered <- pheno_data %>%
  filter(group %in% c("Canonical", "NonCanonical", "Atypical"))
ggplot(pheno_data_filtered, aes(x = group, y = depth_of_invasion, fill = group)) +
  geom_violin(scale = "width", position = position_identity(), alpha = 0.6, color = "black") +
  geom_boxplot(width = 0.1, position = position_identity(), color = "black", alpha = 0.7) +  # Optional: add a boxplot overlay
  geom_jitter(shape = 16, position = position_jitter(0.1), alpha = 0.4) +  # Optional: add individual data points
  labs(title = "Depth of Invasion after RCT by Group",
       x = "Group",
       y = "Depth of Invasion") +
  theme_minimal() +
  theme(legend.position = "none")

pheno_scort_data <- data.frame(
  depth_of_invasion = as.numeric(factor(pheno$TRG, levels(as.factor(pheno$TRG)), c(1,3,2,0))),
  group = final_clusters$x
)
pheno_scort_data_filtered <- pheno_scort_data %>%
  filter(group %in% c("Canonical", "NonCanonical", "Atypical"))
ggplot(pheno_scort_data_filtered, aes(x = group, y = depth_of_invasion, fill = group)) +
  geom_violin(scale = "width", position = position_identity(), alpha = 0.6, color = "black") +
  geom_boxplot(width = 0.1, position = position_identity(), color = "black", alpha = 0.7) +  # Optional: add a boxplot overlay
  geom_jitter(shape = 16, position = position_jitter(0.1), alpha = 0.4) +  # Optional: add individual data points
  labs(title = "Depth of Invasion after RCT by Group",
       x = "Group",
       y = "Depth of Invasion") +
  theme_minimal() +
  theme(legend.position = "none")


final_clusters

new_mini <-cbind(clusters = as.factor(final_clusters=="Canonical"), pheno[c("OS","OS.Status", "DFS", "DFS.Status")], scaled_scort[,nfkb_four_genes])
survival_gene_df$clusters <- as.factor(final_clusters=="Canonical")
#survival_gene_df$clusters <- as.factor(pheno$NFkB.Cluster)
table(final_clusters$x)
dfsnew_mini <- new_mini[new_mini$OS<50,]


dfsnew_mini[complete.cases(new_mini),]
os <-ggsurvplot(survfit(Surv(OS, OS.Status=='DECEASED') ~ clusters, data = new_mini[complete.cases(new_mini),]), data = new_mini[complete.cases(new_mini),], pval = TRUE, xlab = "Time (months)",  risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival For Clusters", legend.title = "Cluster", legend.labs = c("Not Canonical", "Canonical"))
dfs <-ggsurvplot(survfit(Surv(DFS, DFS.Status=='Recurred') ~ clusters, data = dfsnew_mini[complete.cases(new_mini),]), data = dfsnew_mini[complete.cases(new_mini),], pval = TRUE, xlab = "Time (months)",  risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease Free Survival For Clusters", legend.title = "Cluster", legend.labs = c("Not Canonical", "Canonical"))
save(similar,file="Final_gene_list.Rdata")

png("os_scort_clus.png")
print(os, newpage = FALSE)
dev.off()
png("dfs_scort_clus.png")
print(dfs, newpage = FALSE)
dev.off()

new_mini[(new_mini$OS>30&new_mini$OS<40),]
new_mini$OS.Status == "DECEASED"

pheno_87211 <-read.csv("Processed_Data_Sets/pheno87211.csv", row.names = 1)

finclus_87211 <- c(read.csv("13. Final_cluster_derivation/final_87211_clusters.csv", row.names = 1)$x)
table(finclus_87211)


new_mini87211 <-cbind(clusters = as.factor(finclus_87211=="Canonical"), pheno_87211[,c("survival.time..month..ch1","death.due.to.tumor.ch1","disease.free.time..month..ch1","cancer.recurrance.after.surgery.ch1")], scaled_87211[,nfkb_four_genes])
names(new_mini87211) <- c("clusters", "OS","OS.Status", "DFS", "DFS.Status", nfkb_four_genes)
head(new_mini87211)
dfsnew_mini87211 <- new_mini87211[new_mini87211$OS<50,]
ggsurvplot(survfit(Surv(OS, OS.Status==1) ~ clusters, data = new_mini87211[complete.cases(new_mini87211),]), data = new_mini87211[complete.cases(new_mini87211),], pval = TRUE, xlab = "Time (months)",  risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival For Clusters, in GSE87211", legend.title = "Cluster", legend.labs = c("Canonical", "Not Canonical"))
ggsurvplot(survfit(Surv(DFS, DFS.Status==1) ~ clusters, data = new_mini87211[complete.cases(new_mini87211),]), data = new_mini87211[complete.cases(new_mini87211),], pval = TRUE, xlab = "Time (months)",  risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease Free Survival For Clusters, in GSE87211", legend.title = "Cluster", legend.labs = c("Canonical", "Not Canonical"))





## in scort
## MAD2L1
OS_plot_MAD2L1 <- ggsurvplot(survfit(Surv(OS, OS.Status=='DECEASED') ~ MAD2L1_binary, data = survival_gene_df), data = survival_gene_df, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival MAD2L1", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
DFS_plot_MAD2L1 <- ggsurvplot(survfit(Surv(DFS, DFS.Status=='Recurred') ~ MAD2L1_binary, data = dfs_data), data = dfs_data, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival MAD2L1", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/MAD2L1/OS_plot.png")
print(OS_plot_MAD2L1, newpage = FALSE)
dev.off()
png("8b. ScortSurvival_for_top_genes/MAD2L1/DFS_plot.png")
print(DFS_plot_MAD2L1, newpage = FALSE)
dev.off()

## HIBADH
OS_plot_HIBADH <- ggsurvplot(survfit(Surv(OS, OS.Status=='DECEASED') ~ HIBADH_binary, data = survival_gene_df), data = survival_gene_df, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival HIBADH", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/HIBADH/OS_plot.png")
print(OS_plot_HIBADH, newpage = FALSE)
dev.off()
DFS_plot_HIBADH <- ggsurvplot(survfit(Surv(DFS, DFS.Status=='Recurred') ~ HIBADH_binary, data = dfs_data), data = dfs_data, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival HIBADH", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/HIBADH/DFS_plot.png")
print(DFS_plot_HIBADH, newpage = FALSE)
dev.off()

## ZBTB8OS
OS_plot_ZBTB8OS <- ggsurvplot(survfit(Surv(OS, OS.Status=='DECEASED') ~ ZBTB8OS_binary, data = survival_gene_df), data = survival_gene_df, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival ZBTB8OS", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/ZBTB8OS/OS_plot.png")
print(OS_plot_ZBTB8OS, newpage = FALSE)
dev.off()
DFS_plot_ZBTB8OS <- ggsurvplot(survfit(Surv(DFS, DFS.Status=='Recurred') ~ ZBTB8OS_binary, data = dfs_data), data = dfs_data, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival ZBTB8OS", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/ZBTB8OS/DFS_plot.png")
print(DFS_plot_ZBTB8OS, newpage = FALSE)
dev.off()

## KCTD18
OS_plot_KCTD18 <- ggsurvplot(survfit(Surv(OS, OS.Status=='DECEASED') ~ KCTD18_binary, data = survival_gene_df), data = survival_gene_df, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival KCTD18", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/KCTD18/OS_plot.png")
print(OS_plot_KCTD18, newpage = FALSE)
dev.off()
DFS_plot_KCTD18 <- ggsurvplot(survfit(Surv(DFS, DFS.Status=='Recurred') ~ KCTD18_binary, data = dfs_data), data = dfs_data, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival KCTD18", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/KCTD18/DFS_plot.png")
print(DFS_plot_KCTD18, newpage = FALSE)
dev.off()

## ARL6IP1
OS_plot_ARL6IP1 <- ggsurvplot(survfit(Surv(OS, OS.Status=='DECEASED') ~ ARL6IP1_binary, data = survival_gene_df), data = survival_gene_df, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival ARL6IP1", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/ARL6IP1/OS_plot.png")
print(OS_plot_ARL6IP1, newpage = FALSE)
dev.off()
DFS_plot_ARL6IP1 <- ggsurvplot(survfit(Surv(DFS, DFS.Status=='Recurred') ~ ARL6IP1_binary, data = dfs_data), data = dfs_data, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival ARL6IP1", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/ARL6IP1/DFS_plot.png")
print(DFS_plot_ARL6IP1, newpage = FALSE)
dev.off()

## NIPSNAP2
OS_plot_NIPSNAP2 <- ggsurvplot(survfit(Surv(OS, OS.Status=='DECEASED') ~ NIPSNAP2_binary, data = survival_gene_df), data = survival_gene_df, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival NIPSNAP2", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/NIPSNAP2/OS_plot.png")
print(OS_plot_NIPSNAP2, newpage = FALSE)
dev.off()
DFS_plot_NIPSNAP2 <- ggsurvplot(survfit(Surv(DFS, DFS.Status=='Recurred') ~ NIPSNAP2_binary, data = dfs_data), data = dfs_data, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival NIPSNAP2", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/NIPSNAP2/DFS_plot.png")
print(DFS_plot_NIPSNAP2, newpage = FALSE)
dev.off()

## MCM6
OS_plot_MCM6 <- ggsurvplot(survfit(Surv(OS, OS.Status=='DECEASED') ~ MCM6_binary, data = survival_gene_df), data = survival_gene_df, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival MCM6", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/MCM6/OS_plot.png")
print(OS_plot_MCM6, newpage = FALSE)
dev.off()
DFS_plot_MCM6 <- ggsurvplot(survfit(Surv(DFS, DFS.Status=='Recurred') ~ MCM6_binary, data = dfs_data), data = dfs_data, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival MCM6", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/MCM6/DFS_plot.png")
print(DFS_plot_MCM6, newpage = FALSE)
dev.off()

## HACD3
OS_plot_HACD3 <- ggsurvplot(survfit(Surv(OS, OS.Status=='DECEASED') ~ HACD3_binary, data = survival_gene_df), data = survival_gene_df, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival HACD3", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/HACD3/OS_plot.png")
print(OS_plot_HACD3, newpage = FALSE)
dev.off()
DFS_plot_HACD3 <- ggsurvplot(survfit(Surv(DFS, DFS.Status=='Recurred') ~ HACD3_binary, data = dfs_data), data = dfs_data, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival HACD3", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/HACD3/DFS_plot.png")
print(DFS_plot_HACD3, newpage = FALSE)
dev.off()

## TTF2
OS_plot_TTF2 <- ggsurvplot(survfit(Surv(OS, OS.Status=='DECEASED') ~ TTF2_binary, data = survival_gene_df), data = survival_gene_df, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival TTF2", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/TTF2/OS_plot.png")
print(OS_plot_TTF2, newpage = FALSE)
dev.off()
DFS_plot_TTF2 <- ggsurvplot(survfit(Surv(DFS, DFS.Status=='Recurred') ~ TTF2_binary, data = dfs_data), data = dfs_data, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival TTF2", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/TTF2/DFS_plot.png")
print(DFS_plot_TTF2, newpage = FALSE)
dev.off()

## CENPK
OS_plot_CENPK <- ggsurvplot(survfit(Surv(OS, OS.Status=='DECEASED') ~ CENPK_binary, data = survival_gene_df), data = survival_gene_df, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival CENPK", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/CENPK/OS_plot.png")
print(OS_plot_CENPK, newpage = FALSE)
dev.off()
DFS_plot_CENPK <- ggsurvplot(survfit(Surv(DFS, DFS.Status=='Recurred') ~ CENPK_binary, data = dfs_data), data = dfs_data, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival CENPK", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/CENPK/DFS_plot.png")
print(DFS_plot_CENPK, newpage = FALSE)
dev.off()

## ANLN
OS_plot_ANLN <- ggsurvplot(survfit(Surv(OS, OS.Status=='DECEASED') ~ ANLN_binary, data = survival_gene_df), data = survival_gene_df, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival ANLN", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/ANLN/OS_plot.png")
print(OS_plot_ANLN, newpage = FALSE)
dev.off()
DFS_plot_ANLN <- ggsurvplot(survfit(Surv(DFS, DFS.Status=='Recurred') ~ ANLN_binary, data = dfs_data), data = dfs_data, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival ANLN", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/ANLN/DFS_plot.png")
print(DFS_plot_ANLN, newpage = FALSE)
dev.off()

## SUCLA2
OS_plot_SUCLA2 <- ggsurvplot(survfit(Surv(OS, OS.Status=='DECEASED') ~ SUCLA2_binary, data = survival_gene_df), data = survival_gene_df, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival SUCLA2", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/SUCLA2/OS_plot.png")
print(OS_plot_SUCLA2, newpage = FALSE)
dev.off()
DFS_plot_SUCLA2 <- ggsurvplot(survfit(Surv(DFS, DFS.Status=='Recurred') ~ SUCLA2_binary, data = dfs_data), data = dfs_data, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival SUCLA2", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/SUCLA2/DFS_plot.png")
print(DFS_plot_SUCLA2, newpage = FALSE)
dev.off()

## CXCL3
OS_plot_CXCL3 <- ggsurvplot(survfit(Surv(OS, OS.Status=='DECEASED') ~ CXCL3_binary, data = survival_gene_df), data = survival_gene_df, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival CXCL3", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/CXCL3/OS_plot.png")
print(OS_plot_CXCL3, newpage = FALSE)
dev.off()
DFS_plot_CXCL3 <- ggsurvplot(survfit(Surv(DFS, DFS.Status=='Recurred') ~ CXCL3_binary, data = dfs_data), data = dfs_data, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival CXCL3", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/CXCL3/DFS_plot.png")
print(DFS_plot_CXCL3, newpage = FALSE)
dev.off()

## MSMO1
OS_plot_MMP12 <- ggsurvplot(survfit(Surv(OS, OS.Status=='DECEASED') ~ MMP12_binary, data = survival_gene_df), data = survival_gene_df, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival MMP12", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/MMP12/OS_plot.png")
print(OS_plot_MMP12, newpage = FALSE)
dev.off()
DFS_plot_MMP12 <- ggsurvplot(survfit(Surv(DFS, DFS.Status=='Recurred') ~ MMP12_binary, data = dfs_data), data = dfs_data, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival MMP12", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/MMP12/DFS_plot.png")
print(DFS_plot_MMP12, newpage = FALSE)
dev.off()

## EREG
OS_plot_EREG <- ggsurvplot(survfit(Surv(OS, OS.Status=='DECEASED') ~ EREG_binary, data = survival_gene_df), data = survival_gene_df, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival EREG", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/EREG/OS_plot.png")
print(OS_plot_EREG, newpage = FALSE)
dev.off()
DFS_plot_EREG <- ggsurvplot(survfit(Surv(DFS, DFS.Status=='Recurred') ~ EREG_binary, data = dfs_data), data = dfs_data, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival EREG", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/EREG/DFS_plot.png")
print(DFS_plot_EREG, newpage = FALSE)
dev.off()

## ST6GALNAC1
OS_plot_ST6GALNAC1 <- ggsurvplot(survfit(Surv(OS, OS.Status=='DECEASED') ~ ST6GALNAC1_binary, data = survival_gene_df), data = survival_gene_df, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival ST6GALNAC1", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/ST6GALNAC1/OS_plot.png")
print(OS_plot_ST6GALNAC1, newpage = FALSE)
dev.off()
DFS_plot_ST6GALNAC1 <- ggsurvplot(survfit(Surv(DFS, DFS.Status=='Recurred') ~ ST6GALNAC1_binary, data = dfs_data), data = dfs_data, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival ST6GALNAC1", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/ST6GALNAC1/DFS_plot.png")
print(DFS_plot_ST6GALNAC1, newpage = FALSE)
dev.off()

## MMP3
OS_plot_MMP3 <- ggsurvplot(survfit(Surv(OS, OS.Status=='DECEASED') ~ MMP3_binary, data = survival_gene_df), data = survival_gene_df, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival MMP3", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/MMP3/OS_plot.png")
print(OS_plot_MMP3, newpage = FALSE)
dev.off()
DFS_plot_MMP3 <- ggsurvplot(survfit(Surv(DFS, DFS.Status=='Recurred') ~ MMP3_binary, data = dfs_data), data = dfs_data, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival MMP3", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8b. ScortSurvival_for_top_genes/MMP3/DFS_plot.png")
print(DFS_plot_MMP3, newpage = FALSE)
dev.off()







pheno_87211 <- read.csv("Processed_Data_Sets/pheno87211.csv", row.names = 1)
processed87211 <-read.csv("CB_normalized_data/normalized_87211.csv",row.names = 1)
scaled_87211 <- as.data.frame(scale(processed87211))
## in gse87211
head(pheno_87211)
survival_gene_df_87211 <- scaled_87211[,c("ATP1A2","HIBADH", "ZBTB8OS", "KCTD18", "ARL6IP1", "NIPSNAP2", "MCM6", "HACD3", "TTF2", "CENPK", "ANLN", "MSMO1", "SUCLA2", "MAD2L1", "CXCL3", "MMP12", "EREG", "ST6GALNAC1", "MMP3")]
survival_gene_df_87211 <- cbind(create_binary_columns(survival_gene_df_87211), OS = pheno_87211$survival.time..month..ch1)
survival_gene_df_87211 <- cbind(survival_gene_df_87211, OS.Status = as.factor(pheno_87211$death.due.to.tumor.ch1), DFS = pheno_87211$disease.free.time..month..ch1, DFS.Status = as.factor(pheno_87211$cancer.recurrance.after.surgery.ch1))
head(survival_gene_df_87211)
dfs_data_87211 <-survival_gene_df_87211[survival_gene_df_87211$DFS<40,]

## MAD2L1
OS_plot_MAD2L1 <- ggsurvplot(survfit(Surv(OS, OS.Status==1) ~ MAD2L1_binary, data = survival_gene_df_87211), data = survival_gene_df_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival MAD2L1", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/MAD2L1/OS_plot.png")
print(OS_plot_MAD2L1, newpage = FALSE)
dev.off()
DFS_plot_MAD2L1 <- ggsurvplot(survfit(Surv(DFS, DFS.Status==1) ~ MAD2L1_binary, data = dfs_data_87211), data = dfs_data_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival MAD2L1", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/MAD2L1/DFS_plot.png")
print(DFS_plot_MAD2L1, newpage = FALSE)
dev.off()

## HIBADH
OS_plot_HIBADH <- ggsurvplot(survfit(Surv(OS, OS.Status==1) ~ HIBADH_binary, data = survival_gene_df_87211), data = survival_gene_df_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival HIBADH", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/HIBADH/OS_plot.png")
print(OS_plot_HIBADH, newpage = FALSE)
dev.off()
DFS_plot_HIBADH <- ggsurvplot(survfit(Surv(DFS, DFS.Status==1) ~ HIBADH_binary, data = dfs_data_87211), data = dfs_data_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival HIBADH", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/HIBADH/DFS_plot.png")
print(DFS_plot_HIBADH, newpage = FALSE)
dev.off()

## ZBTB8OS
OS_plot_ZBTB8OS <- ggsurvplot(survfit(Surv(OS, OS.Status==1) ~ ZBTB8OS_binary, data = survival_gene_df_87211), data = survival_gene_df_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival ZBTB8OS", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/ZBTB8OS/OS_plot.png")
print(OS_plot_ZBTB8OS, newpage = FALSE)
dev.off()
DFS_plot_ZBTB8OS <- ggsurvplot(survfit(Surv(DFS, DFS.Status==1) ~ ZBTB8OS_binary, data = dfs_data_87211), data = dfs_data_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival ZBTB8OS", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/ZBTB8OS/DFS_plot.png")
print(DFS_plot_ZBTB8OS, newpage = FALSE)
dev.off()

## KCTD18
OS_plot_KCTD18 <- ggsurvplot(survfit(Surv(OS, OS.Status==1) ~ KCTD18_binary, data = survival_gene_df_87211), data = survival_gene_df_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival KCTD18", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/KCTD18/OS_plot.png")
print(OS_plot_KCTD18, newpage = FALSE)
dev.off()
DFS_plot_KCTD18 <- ggsurvplot(survfit(Surv(DFS, DFS.Status==1) ~ KCTD18_binary, data = dfs_data_87211), data = dfs_data_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival KCTD18", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/KCTD18/DFS_plot.png")
print(DFS_plot_KCTD18, newpage = FALSE)
dev.off()

## ARL6IP1
OS_plot_ARL6IP1 <- ggsurvplot(survfit(Surv(OS, OS.Status==1) ~ ARL6IP1_binary, data = survival_gene_df_87211), data = survival_gene_df_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival ARL6IP1", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/ARL6IP1/OS_plot.png")
print(OS_plot_ARL6IP1, newpage = FALSE)
dev.off()
DFS_plot_ARL6IP1 <- ggsurvplot(survfit(Surv(DFS, DFS.Status==1) ~ ARL6IP1_binary, data = dfs_data_87211), data = dfs_data_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival ARL6IP1", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/ARL6IP1/DFS_plot.png")
print(DFS_plot_ARL6IP1, newpage = FALSE)
dev.off()

## NIPSNAP2
OS_plot_NIPSNAP2 <- ggsurvplot(survfit(Surv(OS, OS.Status==1) ~ NIPSNAP2_binary, data = survival_gene_df_87211), data = survival_gene_df_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival NIPSNAP2", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/NIPSNAP2/OS_plot.png")
print(OS_plot_NIPSNAP2, newpage = FALSE)
dev.off()
DFS_plot_NIPSNAP2 <- ggsurvplot(survfit(Surv(DFS, DFS.Status==1) ~ NIPSNAP2_binary, data = dfs_data_87211), data = dfs_data_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival NIPSNAP2", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/NIPSNAP2/DFS_plot.png")
print(DFS_plot_NIPSNAP2, newpage = FALSE)
dev.off()

## MCM6
OS_plot_MCM6 <- ggsurvplot(survfit(Surv(OS, OS.Status==1) ~ MCM6_binary, data = survival_gene_df_87211), data = survival_gene_df_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival MCM6", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/MCM6/OS_plot.png")
print(OS_plot_MCM6, newpage = FALSE)
dev.off()
DFS_plot_MCM6 <- ggsurvplot(survfit(Surv(DFS, DFS.Status==1) ~ MCM6_binary, data = dfs_data_87211), data = dfs_data_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival MCM6", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/MCM6/DFS_plot.png")
print(DFS_plot_MCM6, newpage = FALSE)
dev.off()

## HACD3
OS_plot_HACD3 <- ggsurvplot(survfit(Surv(OS, OS.Status==1) ~ HACD3_binary, data = survival_gene_df_87211), data = survival_gene_df_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival HACD3", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/HACD3/OS_plot.png")
print(OS_plot_HACD3, newpage = FALSE)
dev.off()
DFS_plot_HACD3 <- ggsurvplot(survfit(Surv(DFS, DFS.Status==1) ~ HACD3_binary, data = dfs_data_87211), data = dfs_data_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival HACD3", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/HACD3/DFS_plot.png")
print(DFS_plot_HACD3, newpage = FALSE)
dev.off()

## TTF2
OS_plot_TTF2 <- ggsurvplot(survfit(Surv(OS, OS.Status==1) ~ TTF2_binary, data = survival_gene_df_87211), data = survival_gene_df_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival TTF2", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/TTF2/OS_plot.png")
print(OS_plot_TTF2, newpage = FALSE)
dev.off()
DFS_plot_TTF2 <- ggsurvplot(survfit(Surv(DFS, DFS.Status==1) ~ TTF2_binary, data = dfs_data_87211), data = dfs_data_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival TTF2", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/TTF2/DFS_plot.png")
print(DFS_plot_TTF2, newpage = FALSE)
dev.off()

## CENPK
OS_plot_CENPK <- ggsurvplot(survfit(Surv(OS, OS.Status==1) ~ CENPK_binary, data = survival_gene_df_87211), data = survival_gene_df_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival CENPK", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/CENPK/OS_plot.png")
print(OS_plot_CENPK, newpage = FALSE)
dev.off()
DFS_plot_CENPK <- ggsurvplot(survfit(Surv(DFS, DFS.Status==1) ~ CENPK_binary, data = dfs_data_87211), data = dfs_data_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival CENPK", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/CENPK/DFS_plot.png")
print(DFS_plot_CENPK, newpage = FALSE)
dev.off()

## ANLN
OS_plot_ANLN <- ggsurvplot(survfit(Surv(OS, OS.Status==1) ~ ANLN_binary, data = survival_gene_df_87211), data = survival_gene_df_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival ANLN", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/ANLN/OS_plot.png")
print(OS_plot_ANLN, newpage = FALSE)
dev.off()
DFS_plot_ANLN <- ggsurvplot(survfit(Surv(DFS, DFS.Status==1) ~ ANLN_binary, data = dfs_data_87211), data = dfs_data_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival ANLN", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/ANLN/DFS_plot.png")
print(DFS_plot_ANLN, newpage = FALSE)
dev.off()

## MSMO1
OS_plot_MSMO1 <- ggsurvplot(survfit(Surv(OS, OS.Status==1) ~ MSMO1_binary, data = survival_gene_df_87211), data = survival_gene_df_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival MSMO1", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/MSMO1/OS_plot.png")
print(OS_plot_MSMO1, newpage = FALSE)
dev.off()
DFS_plot_MSMO1 <- ggsurvplot(survfit(Surv(DFS, DFS.Status==1) ~ MSMO1_binary, data = dfs_data_87211), data = dfs_data_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival MSMO1", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/MSMO1/DFS_plot.png")
print(DFS_plot_MSMO1, newpage = FALSE)
dev.off()

## SUCLA2
OS_plot_SUCLA2 <- ggsurvplot(survfit(Surv(OS, OS.Status==1) ~ SUCLA2_binary, data = survival_gene_df_87211), data = survival_gene_df_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival SUCLA2", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/SUCLA2/OS_plot.png")
print(OS_plot_SUCLA2, newpage = FALSE)
dev.off()
DFS_plot_SUCLA2 <- ggsurvplot(survfit(Surv(DFS, DFS.Status==1) ~ SUCLA2_binary, data = dfs_data_87211), data = dfs_data_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival SUCLA2", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/SUCLA2/DFS_plot.png")
print(DFS_plot_SUCLA2, newpage = FALSE)
dev.off()

## CXCL3
OS_plot_CXCL3 <- ggsurvplot(survfit(Surv(OS, OS.Status==1) ~ CXCL3_binary, data = survival_gene_df_87211), data = survival_gene_df_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival CXCL3", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/CXCL3/OS_plot.png")
print(OS_plot_CXCL3, newpage = FALSE)
dev.off()
DFS_plot_CXCL3 <- ggsurvplot(survfit(Surv(DFS, DFS.Status==1) ~ CXCL3_binary, data = dfs_data_87211), data = dfs_data_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival CXCL3", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/CXCL3/DFS_plot.png")
print(DFS_plot_CXCL3, newpage = FALSE)
dev.off()

## MMP12
OS_plot_MMP12 <- ggsurvplot(survfit(Surv(OS, OS.Status==1) ~ MMP12_binary, data = survival_gene_df_87211), data = survival_gene_df_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival MMP12", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/MMP12/OS_plot.png")
print(OS_plot_MMP12, newpage = FALSE)
dev.off()
DFS_plot_MMP12 <- ggsurvplot(survfit(Surv(DFS, DFS.Status==1) ~ MMP12_binary, data = dfs_data_87211), data = dfs_data_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival MMP12", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/MMP12/DFS_plot.png")
print(DFS_plot_MMP12, newpage = FALSE)
dev.off()

## EREG
OS_plot_EREG <- ggsurvplot(survfit(Surv(OS, OS.Status==1) ~ EREG_binary, data = survival_gene_df_87211), data = survival_gene_df_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival EREG", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/EREG/OS_plot.png")
print(OS_plot_EREG, newpage = FALSE)
dev.off()
DFS_plot_EREG <- ggsurvplot(survfit(Surv(DFS, DFS.Status==1) ~ EREG_binary, data = dfs_data_87211), data = dfs_data_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival EREG", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/EREG/DFS_plot.png")
print(DFS_plot_EREG, newpage = FALSE)
dev.off()

## ST6GALNAC1
OS_plot_ST6GALNAC1 <- ggsurvplot(survfit(Surv(OS, OS.Status==1) ~ ST6GALNAC1_binary, data = survival_gene_df_87211), data = survival_gene_df_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival ST6GALNAC1", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/ST6GALNAC1/OS_plot.png")
print(OS_plot_ST6GALNAC1, newpage = FALSE)
dev.off()
DFS_plot_ST6GALNAC1 <- ggsurvplot(survfit(Surv(DFS, DFS.Status==1) ~ ST6GALNAC1_binary, data = dfs_data_87211), data = dfs_data_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival ST6GALNAC1", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
png("8. Survival_for_top_genes/ST6GALNAC1/DFS_plot.png")
print(DFS_plot_ST6GALNAC1, newpage = FALSE)
dev.off()

## ATP1A2
OS_plot_ATP1A2 <- ggsurvplot(survfit(Surv(OS, OS.Status==1) ~ ATP1A2_binary, data = survival_gene_df_87211), data = survival_gene_df_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival ATP1A2", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
OS_plot_ATP1A2
png("8. Survival_for_top_genes/MMP3/OS_plot.png")
print(OS_plot_ATP1A2, newpage = FALSE)
dev.off()
DFS_plot_ATP1A2 <- ggsurvplot(survfit(Surv(DFS, DFS.Status==1) ~ ATP1A2_binary, data = dfs_data_87211), data = dfs_data_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Disease-Free Survival ATP1A2", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
DFS_plot_ATP1A2
png("8. Survival_for_top_genes/MMP3/DFS_plot.png")
print(DFS_plot_ATP1A2, newpage = FALSE)
dev.off()




scort_survival_lots <- list(OS_plot_NXNL1 = OS_plot_NXNL1, DFS_plot_NXNL1 = DFS_plot_NXNL1, OS_plot_SYNGAP1 = OS_plot_SYNGAP1, DFS_plot_SYNGAP1 = DFS_plot_SYNGAP1, OS_plot_LDHAL6A = OS_plot_LDHAL6A, DFS_plot_LDHAL6A = DFS_plot_LDHAL6A, OS_plot_RAC3 = OS_plot_RAC3, DFS_plot_RAC3 = DFS_plot_RAC3, OS_plot_OGDH = OS_plot_OGDH, DFS_plot_OGDH = DFS_plot_OGDH, OS_plot_EML4 = OS_plot_EML4, DFS_plot_EML4 = DFS_plot_EML4, OS_plot_LOC105375116 = OS_plot_LOC105375116, DFS_plot_LOC105375116 = DFS_plot_LOC105375116, OS_plot_LOC101929400 = OS_plot_LOC101929400, DFS_plot_LOC101929400 = DFS_plot_LOC101929400, OS_plot_IGK = OS_plot_IGK, DFS_plot_IGK = DFS_plot_IGK)
save(scort_survival_lots, file = "8. Survival_for_top_genes/scort_survival_plots.RData")

###### No survival/death data for geo sets so can't verify in these.