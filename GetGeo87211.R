### try gse87211
library(GEOquery)
library(oligo)
library(WGCNA)
library(AnnotationForge)
library(AnnotationDbi)
library(BiocManager)
BiocManager::install("HsAgilentDesign026652.db", force=T)    
library(HsAgilentDesign026652.db) 
library(affycoretools)   

# load series and platform data from GEO
z <- getGEO("GSE87211", AnnotGPL = TRUE)
z<-z[[1]]
load("Useful_Functions/annotate_function.RData")
class(z)
processed87211 <- make_gene_data_frame(z, HsAgilentDesign026652.db)
class(processed87211)

intersect(names(processed87211), names(processed109057))
## Remove nfkb2 probe row, annotate normally then add it back
nfkb2_row <- exprs(z)[grep("A_23_P202156", rownames(z)),]
exprs(z) <- rbind(exprs(z)[-grep("A_23_P202156", rownames(z)),], A_23_P202156 = rep(0,length(nfkb2_row)))
processed87211 <- make_gene_data_frame(z, HsAgilentDesign026652.db)
processed87211 <- cbind(processed87211, NFKB2 = nfkb2_row)
summary(processed87211$NFKB2)


scaled_87211 <-as.data.frame(scale(processed87211))

nfkb_four_genes <- c("NFKB1","RELA","NFKB2","RELB")

nfkb_87211 <- scaled_87211[,nfkb_four_genes]
head(nfkb_87211)
## This divides up the samples into those with expression values that very strongly indicate cluster membership and then from these high likelihood assignments find centroid points by taking their means
label_data <- function(df) {
  labels <- rep("messy", nrow(df))
  
  # Cluster 1: Upregulated NFKB1 and RELA, Downregulated NFKB2 and RELB
  labels[df$NFKB1 > 0 & df$RELA > 0 & df$NFKB2 < 0 & df$RELB < 0] <- "up_NFKB1_RELA_down_NFKB2_RELB"
  
  # Cluster 2: Downregulated NFKB1 and RELA, Upregulated NFKB2 and RELB
  labels[df$NFKB1 < 0 & df$RELA < 0 & df$NFKB2 > 0 & df$RELB > 0] <- "down_NFKB1_RELA_up_NFKB2_RELB"
  
  return(labels)
}
centroids <- rbind(colMeans(nfkb_87211[label_data(nfkb_87211)=="up_NFKB1_RELA_down_NFKB2_RELB",]), colMeans(nfkb_87211[label_data(nfkb_87211)=="down_NFKB1_RELA_up_NFKB2_RELB",]),colMeans(nfkb_87211[label_data(nfkb_87211)=="messy",]))


kmeans_87211_nfkb <- kmeans(nfkb_87211, centers=rbind(c(3,3,-3,-3),
                                                      c(-4,-4,4,4),
                                                      c(0,0,0,0)), nstart = 30)
save(kmeans_87211_nfkb, file = "11. GSE87211/NFκB_derived//nfkb_derived_kmeans_object.RData")
load("11. GSE87211/NFκB_derived/nfkb_derived_kmeans_object.RData")
load("Useful_Functions/find_means.RData")

find_means(nfkb_87211, kmeans_87211_nfkb$cluster)
table(kmeans_87211_nfkb$cluster)

### The means suggest cluster 1 is Canonical, cluster 2 is NonCanonical and cluster 3 is Atypical
nfkb_derived_clusters_87211 <- factor(kmeans_87211_nfkb$cluster, levels=c(1,2,3), labels = c("Canonical", "NonCanonical", "Atypical"))

## Save the clusters
write.csv(nfkb_derived_clusters_87211, "11. GSE87211/NFκB_derived//NFKB_derived_scort_clusters.csv")

library(factoextra)
library(ggplot2)
## Plot the clusters
data_matrix <- as.matrix(processed87211)
cluster_object <- list(data = data_matrix, cluster = nfkb_derived_clusters_87211)
plot1 <- fviz_cluster(cluster_object, data = data_matrix,
                      geom = "point",
                      ellipse.type = "convex",
                      ggtheme = theme_bw()
)
plot1
ggsave(
  "11. GSE87211/NFKB_derived_clustering.png",
  plot = plot1)

## run analysis for nfkb clusters
library(EnhancedVolcano)
load("Useful_Functions/Diff_expression_and_volc_plot_from_clusters_and_data.RData")
diff_gene_and_volc_plots(processed87211, nfkb_derived_clusters_87211, "11. GSE87211/NFκB_derived/")



## run kmeans clustering and analysis
kmeans_87211 <- kmeans(scaled_87211, centers = 3, nstart = 20)
clusters<-kmeans_87211$cluster
find_means(scaled_87211[,nfkb_four_genes], clusters)
## canonical = 2, atypical =3, noncanonical =1
clusters <- factor(clusters, levels= c(1,2,3), labels = c("NonCanonical", "Canonical","Atypical"))
diff_gene_and_volc_plots(processed87211, clusters, "11. GSE87211/kmeans/")


## run pvclust and analysis
library(pvclust)
pv_result_87211 <- pvclust(t(scaled_87211), method.hclust = "ward.D2", nboot = 1, method.dist = "euclidean")
clusters<-cutree(pv_result_87211$hclust, k=3)
find_means(scaled_87211[,nfkb_four_genes], clusters)
## canonical = 2, atypical =1, noncanonical =3
clusters <- factor(clusters, levels= c(1,2,3), labels = c("Atypical", "Canonical", "NonCanonical"))
diff_gene_and_volc_plots(processed87211, clusters, "11. GSE87211/pvclust/")








## Investigate genes already identified (new script)
## predict from centroids
scort_data <- read.csv("Processed_Data_Sets/SCORT_final_data.csv",row.names = 1)
scaled_scort <- as.data.frame(scale(scort_data))

red_set2 <- scaled_scort[,common_features]
# just cluster on nfkb
load("1a. SCORT_NFKB_cluster_derivation_replicating/My_NFKB_Scort_kmeans.RData")
library(flexclust)
new_clus_87211 <- predict(kmeans_scort, scaled_87211[,nfkb_four_genes])
table(new_clus_87211)
find_means(scaled_87211[,nfkb_four_genes], new_clus_87211)
## Canonical - 1, NonCanonical -2, Atypical -3
clusters <- factor(new_clus_87211, levels = c(1,2,3), labels = c("Atypical", "NonCanonical", "Canonical"))

diff_gene_and_volc_plots(processed87211, clusters, "11. GSE87211/scort_clusters/")

load("4. Combinations_of_Full_Data_Gene_Lists/Combined_Gene_Lists.RData")
final_list$canonical_all
top_genes_scort_nfkb <- read.csv("1a. SCORT_NFKB_cluster_derivation_replicating/Top_Gene_List.csv",row.names = 1)
dim(top_genes_scort_nfkb)
top_genes_scort_nfkb$Canonical
intersect(top_genes_scort_nfkb$Canonical, final_list$canonical_all)
diff_gene_and_volc_plots(processed87211, clusters, "11. GSE87211/scort_clusters/special_volcano/", special_labels = top_genes_scort_nfkb$Canonical, overide=T)


final_list$canonical_special


## just from nfkb kmeans
find_means(scaled_87211[,nfkb_four_genes], kmeans_87211_nfkb$cluster)
# 1 - canonical, 2 - noncannoical, 3 - atypical
clusters<- factor(kmeans_87211_nfkb$cluster, levels = c(1,2,3), labels = c("Canonical","NonCanonical", "Atypical"))
diff_gene_and_volc_plots(processed87211, clusters, "11. GSE87211/scort_clusters/special_volcano/", special_labels = top_genes_scort_nfkb$Canonical, overide=F)



### Create clusterings based on fixed centroids: create cluster object with centroids where i want them
centroid_c1 <- rep(0,ncol(processed87211))
centroid_c2 <- rep(0,ncol(processed87211))
centroid_c3 <- rep(0,ncol(processed87211))
centroid_c1[grep("NFKB1", names(processed87211))] <- 1
centroid_c1[grep("NFKB2", names(processed87211))] <- -1
centroid_c1[grep("RELA", names(processed87211))] <- 1
centroid_c1[grep("RELB", names(processed87211))] <- -1

centroid_c2[grep("NFKB1", names(processed87211))] <- -1
centroid_c2[grep("NFKB2", names(processed87211))] <- 1
centroid_c2[grep("RELA", names(processed87211))] <- -1
centroid_c2[grep("RELB", names(processed87211))] <- 1

custom_object <-as.data.frame(rbind(centroid_c1, centroid_c2, centroid_c3, centroid_c3))
names(custom_object) <- names(scaled_87211)
trial_kmeans <- kmeans(custom_object, centers=3)
trial_clus_87211 <- predict(trial_kmeans, scaled_87211)
find_means(scaled_87211[,nfkb_four_genes], trial_clus_87211)
clusters <- factor(trial_clus_87211, levels = c(1,2,3), labels = c("NonCanonical", "Atypical", "Canonical"))
diff_gene_and_volc_plots(processed87211, clusters, "11. GSE87211/scort_clusters/special_volcano/", special_labels = top_genes_scort_nfkb$Canonical, overide = T)


write.csv(scaled_87211, "Processed_Data_Sets/scaled_87211.csv")
write.csv(processed87211, "Processed_Data_Sets/full_87211.csv")

z <- getGEO("GSE87211", AnnotGPL = TRUE)
pheno87211<-pData(z[[1]])
head(pheno87211)
write.csv(pheno87211, "Processed_Data_Sets/pheno87211.csv")
survival_gene_df <- data.frame(OS = as.numeric(pheno$`survival time (month):ch1`), OS.Status=pheno$`death due to tumor:ch1`, UXS1 = scaled_87211[,"UXS1"])
survival_gene_df <- survival_gene_df[!is.na(survival_gene_df$OS),]



load("Useful_Functions/create_binaries.RData")
survival_gene_df <- cbind(survival_gene_df, create_binary_columns(as.data.frame(survival_gene_df$UXS1)))
head(survival_gene_df)
table(is.na(survival_gene_df$OS))
names(survival_gene_df) <- c("OS", "OS.Status", "UXS1","dup" , "UXS1_binary")
library(survminer)
OS_plot_UXS1 <- ggsurvplot(survfit(Surv(OS, OS.Status==1) ~ survival_gene_df$UXS1_binary, data = survival_gene_df), data = survival_gene_df, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival IGK", legend.title = "Gene Expression", legend.labs = c("upregulated", "downregulated"))
OS_plot_UXS1




