### Try to get conclusive comparisons between gene lists etc
# scort
nkfb_scort <- read.csv("1a. SCORT_NFKB_cluster_derivation_replicating/Top_Gene_List.csv", row.names = 1)
kmeans_scort <- read.csv("2a. SCORT_Generate_Kmeans_Clusters_and_Top_Gene/Top_Gene_List.csv", row.names = 1)
pvclust_scort <- read.csv("3a. SCORT_Generate_PV_Clusters_and_Top_Gene//Top_Gene_List.csv", row.names = 1)

intersect(intersect(nkfb_scort$Canonical, kmeans_scort$Canonical),pvclust_scort$Canonical)
intersect(pvclust_scort$Canonical, kmeans_scort$Canonical)
intersect(nkfb_scort$Canonical, kmeans_scort$Canonical)

load("1a. SCORT_NFKB_cluster_derivation_replicating/My_NFKB_Scort_kmeans.RData")
nfkb_my_clusters <- kmeans_scort$cluster
nfkb_my_clusters
find_means(scaled_scort[,nfkb_four_genes], nfkb_my_clusters)
nfkb_my_clusters <- factor(nfkb_my_clusters, levels = c(1,2,3), labels = c("Canonical", "NonCanonical", "Atypical"))
nkfb_clusters_scort <- read.csv("1a. SCORT_NFKB_cluster_derivation_replicating/NFKB_derived_scort_clusters.csv", row.names = 1)
kmeans_clusters_scort <- read.csv("2a. SCORT_Generate_Kmeans_Clusters_and_Top_Gene/cluster_list.csv",row.names = 1)
pvclust_clusters_scort <- read.csv("3a. SCORT_Generate_PV_Clusters_and_Top_Gene/pc_clust_clusters.csv",row.names = 1)


table(nkfb_clusters_scort$x, nfkb_my_clusters)
table(nkfb_clusters_scort$x, kmeans_clusters_scort$x)
table(pvclust_clusters_scort$x, kmeans_clusters_scort$x)
table(pvclust_clusters_scort$x, nkfb_clusters_scort$x)

find_means(scaled_scort[,nfkb_four_genes], nkfb_clusters_scort$x)
find_means(scaled_scort[,nfkb_four_genes], nfkb_my_clusters)
find_means(scaled_scort[,nfkb_four_genes], kmeans_clusters_scort$x)
find_means(scaled_scort[,nfkb_four_genes], pvclust_clusters_scort$x)

pheno <- read.csv("Processed_Data_Sets/scort.csv", row.names = 1)
find_means(scaled_scort[nfkb_four_genes], pheno$NFkB.Cluster)
centroids2 <- rbind(c(1,-1,-1), c(-1,1,1), c(0,0,0))
centroids2
custom_kmeans <- list(
  centers = centroids2,
  cluster = scaled_scort$cluster, totss = NULL, withinss = NULL, tot.withinss = NULL, betweenss = NULL, size = as.numeric(table(scaled_scort$cluster)),iter = 1, ifault = 0
)
class(custom_kmeans) <- "kmeans"
pheno <- read.csv("Processed_Data_Sets/scort.csv", row.names = 1)
pred_from_perf <- predict(custom_kmeans, scaled_scort[,nfkb_four_genes])
table(pred_from_perf)
find_means(scaled_scort[,nfkb_four_genes], pred_from_perf) 
new <- kmeans(scaled_scort[,nfkb_four_genes],3, iter.max=20)
find_means(scaled_scort[,nfkb_four_genes], new$cluster) 
table(factor(new$cluster, c(2,3,1), c("Atypical", "Canonical", "NonCanonical")), pheno$NFkB.Cluster)

perf <- factor(pred_from_perf, c(3,1,2), c("Atypical", "Canonical", "NonCanonical"))
table(factor(pred_from_perf, c(3,1,2), c("Atypical", "Canonical", "NonCanonical")), pheno$NFkB.Cluster)
for(i in 1:length(pheno$NFkB.Cluster)){
  if (pheno$NFkB.Cluster[i] == perf[i]){
    cluster <- perf[i]
  }else{
    print(scaled_scort[i,nfkb_four_genes])
  }
}
library(pvclust)
pv_result <- pvclust::pvclust(t(scaled_scort), method.hclust = "ward.D2", nboot = 1000, method.dist = "euclidean")
plot(pv_result)
pvrect(pv_result, 0.95)
full_kmeans_scort <- kmeans(scaled_scort, centers = 3, nstart = 20)
find_means(scaled_scort[,nfkb_four_genes], full_kmeans_scort$cluster)

find_means(scaled_scort[,nfkb_four_genes], pheno$NFkB.Cluster)

full_kmeans_scortclus <- factor(full_kmeans_scort$cluster, c(3,1,2), c("Atypical", "Canonical","NonCanonical"))
table(full_kmeans_scortclus, pheno$NFkB.Cluster)
table(pheno$NFkB.Cluster)
library(pvclust)
# full_pv_scort <- pvclust(scaled_scort[,-ncol(scaled_scort)], method.hclust = "ward.D2", nboot = 1, method.dist = "euclidean")
dist_mat <- dist(scaled_scort[, -ncol(scaled_scort)], method = 'euclidean')
full_pv_scort <- hclust(dist_mat, method = "ward.D2", members = NULL)
find_means(scaled_scort[,nfkb_four_genes], cutree(pv_result$hclust, k=3))
plot(full_pv_scort)
full_pv_scortclus <- factor(cutree(pv_result$hclust, k=3), c(1,3,2), c("Atypical","Canonical", "Non-Canonical"))
table(full_pv_scortclus, pheno$NFkB.Cluster)
summary(scaled_scort[,"NFKB2"])

boxplot(x = scaled_scort[,nfkb_four_genes])
with3_kmeans <-  kmeans(scaled_scort[,c("NFKB1","NFKB2","RELB")], centers = centroids2, nstart = 20)
table(with3_kmeans$cluster)
find_means(scaled_scort[,nfkb_four_genes], with3_kmeans$cluster)
clusters <- factor(with3_kmeans$cluster, c(1,2,3), c("Canonical", "NonCanonical", "Atypical"))
veccan <- c(scaled_scort[clusters=="Canonical", "RELA"])
vecnon <- c(scaled_scort[clusters=="NonCanonical", "RELA"])
vecaty <- c(scaled_scort[clusters=="Atypical", "RELA"])
df <- data.frame(
  value = c(veccan, vecnon, vecaty),
  group = factor(rep(c("Canonical", "NonCanonical", "Atypical"),
                     times = c(length(veccan), length(vecnon), length(vecaty))))
)
df
ggplot(df, aes(x = group, y = value, fill = group)) +
  geom_violin(trim = FALSE) +
  labs(title = "Violin Plots of Different Length Vectors",
       x = "Vectors",
       y = "Values") +
  theme_minimal()

table(final_clusters)

scaled_scort[clusters=="Canonical" & scaled_scort$RELA < -0.5, nfkb_four_genes]
clusters[clusters=="Canonical" & scaled_scort$RELA < -0.5] <- "Atypical"
scaled_scort[(clusters=="NonCanonical" & scaled_scort$RELA > 0.5), nfkb_four_genes] 
clusters[clusters=="NonCanonical" & scaled_scort$RELA > 0.5] <- "Atypical"
clusters[c(62,86,93,170)] <- "Canonical"
clusters[c(9)] <- "NonCanonical"
table(clusters)

find_means(scaled_scort[,nfkb_four_genes], clusters)

table((clusters=="Canonical") & (scaled_scort$RELA < 0))

table(full_kmeans_scort$cluster,cutree(full_pv_scort, k=3))
table(cutree(full_pv_scort, k=3))
table(factor(new$cluster, c(2,1,3), c("Canonical", "NonCanonical","Atypical")), factor(cutree(full_pv_scort, k=3),c(3,2,1),c("Canonical", "NonCanonical","Atypical")))
table(factor(new$cluster, c(2,1,3), c("Canonical", "NonCanonical","Atypical")), factor(full_kmeans_scort$cluster,c(1,2,3),c("Canonical", "NonCanonical","Atypical")))
table(new$cluster)
table(full_kmeans_scort$cluster)

final_clusters<- clusters
final_clusters
find_means(scaled_scort[,nfkb_four_genes], final_clusters)
# final_clusters<-c()
# for(i in 1:nrow(scaled_scort)){
#   if(nkfb_clusters_scort[i,1]==kmeans_clusters_scort[i,1] || nkfb_clusters_scort[i,1]==pvclust_clusters_scort[i,1]){
#     cluster <- nkfb_clusters_scort[i,1]
#   }else{
#     if(pvclust_clusters_scort[i,1]==kmeans_clusters_scort[i,1]){
#       cluster <- kmeans_clusters_scort[i,1]
#     }else{
#       cluster<-paste0(nkfb_clusters_scort[i,1],"_x")
#     }
#   }
#   final_clusters<-c(final_clusters,cluster)
# }
# final_clusters
# 
# find_means(scaled_scort[, nfkb_four_genes], final_clusters)
# ## see all noncanonical_* so easy to extract indexes
# indices <- grep("Non-Canonical_x",final_clusters)
# indices
# expressions <- scaled_scort[indices, nfkb_four_genes]
# expressions
# ## viewing these show: 66, 15 atypical, 65, 113, 120 noncanonical, 
# final_clusters[c(65, 113, 120)] <- "NonCanonical"
# final_clusters[c(66, 15)] <- "Atypical"

table(final_clusters)
diff_gene_and_volc_plots(scort_data, as.factor(final_clusters), "12. gene list comp/")

write.csv(final_clusters, file = "12. gene list comp/Final_scort_Cluster_List.csv")
conclusive_scort_genes <- read.csv("12. gene list comp/Top_Gene_List.csv",row.names = 1)
conclusive_scort_genes$NonCanonical
conclusive_scort_genes$Canonical
  
## Function to predict from `final cluster` centroids  
library(dplyr)
scort_data<- read.csv("CB_normalized_data/scort_normalized.csv",row.names = 1)
scaled_scort<- as.data.frame(scale(scort_data))
scaled_scort$cluster <- final_clusters
centroids <- scaled_scort %>%
  group_by(cluster) %>%
  summarise(across(1:18092, mean))
dim(centroids)
fixed_centroids <- as.matrix(centroids[, -1])
rownames(fixed_centroids) <- centroids$cluster


custom_kmeans <- list(
  centers = fixed_centroids,
  cluster = scaled_scort$cluster, totss = NULL, withinss = NULL, tot.withinss = NULL, betweenss = NULL, size = as.numeric(table(scaled_scort$cluster)),iter = 1, ifault = 0
)
class(custom_kmeans) <- "kmeans"

## Check in geo109057
scaled_109057 <- read.csv("CB_normalized_data/scaled_109057.csv", row.names = 1)
scaled_145037 <- read.csv("CB_normalized_data/scaled_145037.csv", row.names = 1)
scaled_87211 <- read.csv("CB_normalized_data/scaled_87211.csv", row.names = 1)


library(fdm2id)
new_geo_109057_clus <- predict(custom_kmeans, scaled_109057)
table(new_geo_109057_clus)
find_means(scaled_109057[,nfkb_four_genes], new_geo_109057_clus)
## Suggessts noncanonical is c3, canoniacl is c2 and atpical c1
clusters <- factor(new_geo_109057_clus, levels = c(1,2,3), c("Atypical","Canonical","NonCanonical"))

diff_gene_and_volc_plots(processed109057, clusters, "12. gene list comp/109057/")
signif_109057 <- read.csv("12. gene list comp/109057/Top_Gene_List.csv",row.names = 1)
confirmed_by109 <-intersect(signif_109057$Canonical, conclusive_scort_genes$Canonical)


## Check in geo145037
new_geo_145037_clus <- predict(custom_kmeans, scaled_145037[,common_features_scort_primeview])
table(new_geo_145037_clus)
find_means(scaled_145037[,nfkb_four_genes], new_geo_145037_clus)
## Suggessts noncanonical is c3, canonical is c2 and atpical c1
clusters <- factor(new_geo_145037_clus, levels = c(1,2,3), c("Atypical","Canonical","NonCanonical"))

diff_gene_and_volc_plots(processed145037, clusters, "12. gene list comp/145037//")
signif_145037 <- read.csv("12. gene list comp/145037//Top_Gene_List.csv",row.names = 1)
confirmed_by145 <- intersect(signif_145037$Canonical, conclusive_scort_genes$Canonical)
diff_gene_and_volc_plots(processed145037, clusters, "12. gene list comp/145037/", special_labels = confirmed_by_one, overide = T)

confirmed_by_109and145 <- intersect(intersect(signif_145037$Canonical, conclusive_scort_genes$Canonical), signif_109057$Canonical)


### 87211
common_features_scort_agilent <- intersect(colnames(scaled_scort),colnames(scaled_87211))
custom_kmeans_87211 <- list(
  centers = fixed_centroids[,common_features_scort_agilent],
  cluster = scaled_scort$cluster, totss = NULL, withinss = NULL, tot.withinss = NULL, betweenss = NULL, size = as.numeric(table(scaled_scort$cluster)),iter = 1, ifault = 0
)
class(custom_kmeans_87211) <- "kmeans"
new_geo_87211_clus <- predict(custom_kmeans, scaled_87211[,common_features_scort_agilent])
table(new_geo_87211_clus)
find_means(scaled_87211[,nfkb_four_genes], new_geo_87211_clus)
## Suggessts noncanonical is c3, canoniacl is c2 and atpical c1
clusters <- factor(new_geo_87211_clus, levels = c(1,2,3), c("Atypical","Canonical","NonCanonical"))
write.csv(clusters,"12. gene list comp/Final_87211_Clusters.csv")

diff_gene_and_volc_plots(processed87211, clusters, "12. gene list comp/87211//")
signif_87211 <- read.csv("12. gene list comp/87211/Top_Gene_List.csv",row.names = 1)
confirmed_by872 <-intersect(signif_87211$Canonical, conclusive_scort_genes$Canonical)

confirmed_by_all <- intersect(confirmed_by_109and145, confirmed_by872)

confirmed_by_two <- unique(c(confirmed_by_109and145, intersect(confirmed_by109, confirmed_by872), intersect(confirmed_by145, confirmed_by872)))
confirmed_by_one <- unique(c(confirmed_by109, confirmed_by145, confirmed_by872))

confirmed_by_one
confirmed_by_two
confirmed_by_all

centroids <- rbind(c(1, 1, -1, -1),  c(-1, -1, 1, 1),  c(0,0,0,0))
obj <- kmeans(nfkb_87211, centers=centroids)
nfkb_clusters_87211 <- obj$cluster


## lfc >2.5 times gives:: 
find_means(scaled_109057[,nfkb_four_genes], new_geo_109057_clus)
clusters <- factor(new_geo_109057_clus, levels = c(1,2,3), c("Atypical","Canonical","NonCanonical"))
diff_gene_and_volc_plots(full_109057, clusters, "12. gene list comp/109057/", lfc_boundary = log(2.5,base=2),special_labels = nfkb_four_genes, overide = T)

find_means(scaled_145037[,nfkb_four_genes], new_geo_145037_clus)
clusters <- factor(new_geo_145037_clus, levels = c(1,2,3), c("Atypical","NonCanonical","Canonical"))
diff_gene_and_volc_plots(full_145037, clusters, "12. gene list comp/145037/", lfc_boundary = log(2.5,base=2),special_labels = nfkb_four_genes, overide = T)

find_means(scaled_87211[,nfkb_four_genes], new_geo_87211_clus)
clusters <- factor(new_geo_87211_clus, levels = c(1,2,3), c("Atypical","Canonical","NonCanonical"))
diff_gene_and_volc_plots(full_87211, clusters, "12. gene list comp/87211/", lfc_boundary = log(2.5,base=2),special_labels = nfkb_four_genes, overide = T)

find_means(scaled_scort[,nfkb_four_genes], final_clusters)
clusters <- factor(final_clusters, levels=c(1,2,3), labels = c("Atypical", "Canonical", "NonCanonical"))
clusters <- as.factor(final_clusters)
diff_gene_and_volc_plots(scort_data, clusters, "12. gene list comp/SCORT/", lfc_boundary = log(2.5,base=2),special_labels = nfkb_four_genes, overide = T)


signif_109057 <- read.csv("12. gene list comp/109057/Top_Gene_List.csv",row.names = 1)
signif_87211 <- read.csv("12. gene list comp/87211/Top_Gene_List.csv",row.names = 1)
conclusive_scort_genes <- read.csv("12. gene list comp/SCORT/Top_Gene_List.csv",row.names = 1)
signif_145037 <- read.csv("12. gene list comp/145037//Top_Gene_List.csv",row.names = 1)

similar <- intersect(intersect(signif_109057$Canonical, signif_87211$Canonical), intersect(conclusive_scort_genes$Canonical, signif_145037$Canonical))
similar
all_shared_genes <- intersect(common_features_scort_primeview, common_features_scort_agilent)
prob_gene_shared <- sum(!is.na(signif_145037$Canonical))/length(all_shared_genes)*
sum(!is.na(signif_109057$Canonical))/length(all_shared_genes)*
sum(!is.na(signif_87211$Canonical))/length(all_shared_genes)*
sum(!is.na(conclusive_scort_genes$Canonical))/length(all_shared_genes)

## We'd expect this many deferentially expressed genes across all datasets so fairly confident this is significant
prob_gene_shared*length(all_shared_genes)
prob_gene_shared
length(all_shared_genes)
## Perform binomial calculation with 18080 trials and p=0.001564942 but acc get 47:
pvalue_47genes = stats::binom.test(13, n=length(all_shared_genes), p=prob_gene_shared)

## Pvalue:
pvalue_47genes$p.value



### lfc1 tests:
all_shared_genes <- intersect(common_features_scort_primeview, common_features_scort_agilent)
prob_gene_shared_lfc1 <- sum(!is.na(scort_lfc1$Canonical))/length(all_shared_genes)*
  sum(!is.na(geo145037_lfc1$Canonical))/length(all_shared_genes)*
  sum(!is.na(signif_87211$Canonical))/length(all_shared_genes)*
  sum(!is.na(geo109057_lfc1$Canonical))/length(all_shared_genes)

prob_gene_shared_lfc1*length(all_shared_genes)

pheno_scort$OS.Status[is.na(pheno_scort$OS.Status)] <- "LIVING"
table((pheno_scort$OS.Status))
table(pheno_87211$`death due to tumor:ch1`)
pheno_87211[pheno_87211$`death due to tumor:ch1`=="NA",]


heatmap(as.matrix(scaled_87211[,intersect(names(scaled_87211), confirmed_by_all)]), scale = "column")
heatmap(as.matrix(scaled_109057[,intersect(names(scaled_109057), confirmed_by_all)]), scale = "column")
heatmap(as.matrix(scaled_145037[,intersect(names(scaled_145037), confirmed_by_all)]), scale = "column")
heatmap(as.matrix(scaled_scort[,intersect(names(scaled_scort), confirmed_by_all)]), scale = "column")




#### load final clusters
finclus_scort <- c(read.csv("13. Final_cluster_derivation/final_scort_clusters.csv", row.names = 1)$x)
finclus_109057 <- c(read.csv("13. Final_cluster_derivation/final_109057_clusters.csv", row.names = 1)$x)
finclus_145037 <- c(read.csv("13. Final_cluster_derivation/final_145037_clusters.csv", row.names = 1)$x)
finclus_87211 <- c(read.csv("13. Final_cluster_derivation/final_87211_clusters.csv", row.names = 1)$x)

find_means(scaled_scort[,nfkb_four_genes], finclus_scort)
find_means(scaled_109057[,nfkb_four_genes], finclus_109057)
find_means(scaled_145037[,nfkb_four_genes], finclus_145037)
find_means(scaled_87211[,nfkb_four_genes], finclus_87211)

diff_gene_and_volc_plots(full_109057, as.factor(finclus_109057), "14. gene lists from 13/109057/", lfc_boundary = log(1.75,base=2),special_labels = nfkb_four_genes, overide = T)
diff_gene_and_volc_plots(full_145037, as.factor(finclus_145037), "14. gene lists from 13/145037/", lfc_boundary = log(1.75,base=2),special_labels = nfkb_four_genes, overide = T)
diff_gene_and_volc_plots(full_87211, as.factor(finclus_87211), "14. gene lists from 13/87211/", lfc_boundary = log(1.75,base=2),special_labels = nfkb_four_genes, overide = T)
diff_gene_and_volc_plots(scort_data, as.factor(finclus_scort), "14. gene lists from 13/SCORT/", lfc_boundary = log(1.75,base=2),special_labels = nfkb_four_genes, overide = T)

signif_109057 <- read.csv("14. gene lists from 13/109057/Top_Gene_List.csv",row.names = 1)
signif_87211 <- read.csv("14. gene lists from 13/87211/Top_Gene_List.csv",row.names = 1)
conclusive_scort_genes <- read.csv("14. gene lists from 13/SCORT/Top_Gene_List.csv",row.names = 1)
signif_145037 <- read.csv("14. gene lists from 13/145037//Top_Gene_List.csv",row.names = 1)

head(conclusive_scort_genes)
sum(!is.na(conclusive_scort_genes$Canonical))
sum(!is.na(signif_109057$Canonical))
length(intersect(conclusive_scort_genes$Canonical, signif_109057$Canonical))
intersect(conclusive_scort_genes$Canonical, signif_109057$Canonical)
sum(!is.na(signif_145037$Canonical))
length(intersect(conclusive_scort_genes$Canonical, signif_145037$Canonical))
length(intersect(intersect(conclusive_scort_genes$Canonical, signif_145037$Canonical),signif_109057$Canonical))
sum(!is.na(signif_87211$Canonical))
length(intersect(conclusive_scort_genes$Canonical, signif_87211$Canonical))

similar1 <- intersect(intersect(signif_109057$Canonical, signif_87211$Canonical), intersect(conclusive_scort_genes$Canonical, signif_145037$Canonical))
as.data.frame(similar1)
length(similar1)
similar[!(similar %in% intersect(similar1, similar))]
similar1[!(similar1 %in% intersect(similar1, similar))]
load("Final_gene_list.Rdata")

save(similar1, file = "Final_gene_list.Rdata")
similar[!(similar %in% intersect(similar, similar1))]
similar1
## 1.5 gives 106
## 1.75 gives 22
## 2 gives 4
similar1


gene_list <-similar1[!is.na(similar1)]
mod_set <- scaled_scort[,gene_list]
mod_set$binary_response <- as.factor(as.numeric(as.factor(pheno_scort$Binary_response))-1)
formula <- as.formula(paste0("binary_response ~ ", paste(gene_list, collapse = "+")))
bin_mod <- glm(formula, data=mod_set, family = binomial )
summary(bin_mod)
exp(cbind(coef(bin_mod),confint(bin_mod)))

names(scort_positions)[!names(scort_positions)%in%gene_list]

rownames(conclusive_scort_genes)[conclusive_scort_genes$Canonical%in%gene_list]
scort_positions <- unlist(sapply(gene_list, function(x){mean(grep(x, conclusive_scort_genes$Canonical))}))
positions_109057 <- unlist(sapply(gene_list, function(x){mean(grep(x, signif_109057$Canonical))}))
positions_145037 <- unlist(sapply(gene_list, function(x){mean(grep(x, signif_145037$Canonical))}))
positions_87211 <- unlist(sapply(gene_list, function(x){mean(grep(x, signif_87211$Canonical))}))
totals <- as.data.frame(cbind(scort_positions, positions_109057, positions_145037, positions_87211))
table(rownames(signif_109057)[signif_109057$Canonical=="UTS2B"])

totals
totals$total <- rowSums(totals[,c("scort_positions","positions_109057", "positions_145037")])
ordered_gene_list <-rownames(totals)[order(totals$total)]


mod_set$canonical <- pheno_scort$NFkB.Cluster=="Canonical"
clus_mod <- bin_mod <- glm(binary_response ~ canonical, data=mod_set, family = binomial )
summary(clus_mod)
sjPlot::tab_model(clus_mod)
sjPlot::plot_model(clus_mod)
exp(cbind(coef(clus_mod),confint(clus_mod)))

mod_set87211 <- scaled_87211[rownames(survival_gse87211),gene_list]
finclus_87211<-read.csv("13. Final_cluster_derivation/final_87211_clusters.csv",row.names = 1)
mod_set87211$clusters <- as.factor(finclus_87211[rownames(survival_gse87211),]=="Canonical")
clus_mod2 <- glm(binary_response ~ clusters, data=mod_set, family = binomial )
summary(clus_mod2)
exp(cbind(coef(clus_mod2),confint(clus_mod2)))








all_shared_genes <- intersect(common_features_scort_primeview, common_features_scort_agilent)
prob_gene_shared <- sum(!is.na(signif_145037$Canonical))/length(all_shared_genes)*
  sum(!is.na(signif_109057$Canonical))/length(all_shared_genes)*
  sum(!is.na(signif_87211$Canonical))/length(all_shared_genes)*
  sum(!is.na(conclusive_scort_genes$Canonical))/length(all_shared_genes)

## We'd expect this many deferentially expressed genes across all datasets so fairly confident this is significant
prob_gene_shared*length(all_shared_genes)
prob_gene_shared
length(all_shared_genes)
pvalue_47genes <- stats::binom.test(20, n=length(all_shared_genes), p=prob_gene_shared)
pvalue_47genes$p.value
## Perform binomial calculation with 18080 trials and p=0.001564942 but acc get 47:
pvalue_47genes3 <- stats::binom.test(3, n=length(all_shared_genes), p=prob_gene_shared)
pvalue_47genes2 <- stats::binom.test(2, n=length(all_shared_genes), p=prob_gene_shared)
pvalue_47genes1 <- stats::binom.test(1, n=length(all_shared_genes), p=prob_gene_shared)
pvalue_47genes0 <- stats::binom.test(0, n=length(all_shared_genes), p=prob_gene_shared)
1-(pvalue_47genes3$p.value + pvalue_47genes2$p.value + pvalue_47genes1$p.value +pvalue_47genes0$p.value )
pvalue_47genes0$p.value
log(1.75,2)

df <- scaled_scort[,gene_list]
names(df)
expected_gene_expression
obtain_goodorbad <- function(df){
  new_df <- data.frame(matrix(rep(0,nrow(df)*ncol(df)), ncol=ncol(df), nrow=nrow(df)), row.names = row.names(df))
  colnames(new_df) <- names(df)
  for (col in names(df)){
    if(expected_gene_expression[col]=="downregulated"){
      new_df[,col] <- ifelse(df[,col]>0,"Not-Canonical","Canonical")
    }
    if(expected_gene_expression[col]=="upregulated"){
      new_df[,col] <- ifelse(df[,col]<0,"Not-Canonical","Canonical")
    }
  }
  return(new_df)
}
expected_gene_expression
x <- obtain_goodorbad(df)
x
## 1 marks not canlnocal
x_new <- as.data.frame(sapply(x, function(y){as.numeric(as.factor(y))-1}))
x_new
x_new$total <- as.factor(ifelse(rowSums(x_new)<7, "HighCorrelationCanonical", "Not"))
table(x_new$total)
names(survival_scort)
x_new <- cbind(x_new, survival_scort[,c("DFS", "DFS.Status", "OS", "OS.Status",  "TRG")])
head(x_new)

dfs_x_new<-x_new[x_new$DFS<40,]
ggsurvplot(survfit(Surv(OS, OS.Status=='DECEASED') ~ total, data = x_new), data = x_new, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = paste0("Kaplan-Meier Curve for Overall Survival "), legend.title = "Gene Expression", legend.labs = c("High Correlation with Canonical", "Not"))
ggsurvplot(survfit(Surv(DFS, DFS.Status=='Recurred') ~ total, data = dfs_x_new),data = dfs_x_new, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = paste0("Kaplan-Meier Curve for Overall Survival "), legend.title = "Gene Expression", legend.labs = c("High Correlation with Canonical", "Not"))

df_87211 <- scaled_87211[,gene_list]
x_87211 <- obtain_goodorbad(df_87211)
## 1 marks not canlnocal
x_new_87211 <- as.data.frame(sapply(x_87211, function(y){as.numeric(as.factor(y))-1}))
x_new_87211$total <- as.factor(ifelse(rowSums(x_new_87211)<5, "HighCorrelationCanonical", "Not"))
x_new_87211 <- cbind(x_new_87211[complete.cases(relevantcols),], survival_gse87211[,c("DFS", "DFS.Status", "OS", "OS.Status",  "TRG")])

dfs_x_new_87211<- x_new_87211
ggsurvplot(survfit(Surv(OS, OS.Status=='DECEASED') ~ total, data = x_new_87211), data = x_new_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = paste0("Kaplan-Meier Curve for Overall Survival "), legend.title = "Gene Expression", legend.labs = c("High Correlation with Canonical", "Not"))
ggsurvplot(survfit(Surv(DFS, DFS.Status=='Recurred') ~ total, data = dfs_x_new_87211),data = dfs_x_new_87211, pval = TRUE, xlab = "Time (months)", conf.int = TRUE, risk.table = TRUE, ylab = "Overall Survival Probability", title = paste0("Kaplan-Meier Curve for Overall Survival "), legend.title = "Gene Expression", legend.labs = c("High Correlation with Canonical", "Not"))

x_new
x_new$binary_response <- as.factor(as.numeric(as.factor(ifelse(survival_scort$TRG%in%c("pCR", "Good partial response"), "good", "bad")))-1)
head(x_new_new)
as.formula(paste0("binary_response ~", paste0(gene_list, collapse= "+")))
x_new_new <- scaled_scort[,gene_list]
x_new_new$binary_response <- as.factor(as.numeric(as.factor(ifelse(survival_scort$TRG%in%c("pCR", "Good partial response"), "good", "bad")))-1)


test1 <- glm(as.formula(paste0("binary_response ~", paste0(gene_list, collapse= "+"))), data = x_new_new, family = "binomial")
sjPlot::tab_model(test1)

x_new_87211 <- scaled_87211[complete.cases(relevantcols),gene_list]
x_new_87211$binary_response <- as.factor(as.numeric(as.factor(ifelse(survival_gse87211$TRG%in%c(0, 1), "good", "bad")))-1)
head(cbind(x_new_87211, survival_gse87211$TRG))
test2 <- glm(as.formula(paste0("binary_response ~", paste0(gene_list, collapse= "+"))), data = x_new_87211, family = "binomial")
sjPlot::tab_model(test2)

library(MKmisc)
library(ResourceSelection)
HLgof.test(fitted(test1), obs= x_new_new$binary_response)
hoslem.test(x_new_new$binary_response, fitted(test1), g=10)


sum(x_new_new[,i] == expected_gene_expression[i])

temp_frame <- x_new_new[,-ncol(x_new_new):-(ncol(x_new_new)-1)]
count <- sum(temp_frame[,i] == expected_gene_expression[i])

x_new_new$expected_reg <- rep(0,nrow(x_new_new))

sum(temp_frame[,i] == expected_gene_expression[i])

temp_frame[i,]
expected_gene_expression
sum(temp_frame[i,] ==expected_gene_expression)


temp_frame <- (create_binary_columns(x_new_new[,-ncol(x_new_new):-(ncol(x_new_new)-1)]))[,(1+length(expected_gene_expression)):(2*length(expected_gene_expression))]
for(i in 1:(nrow(x_new_new))){
  count <- sum(temp_frame[i,] == expected_gene_expression)
  x_new_new$expected_reg[i] <- count
  print(count)
}
x_new_new$expected_reg
x_new_new
names(x_new_new)[grep("expected_reg", names(x_new_new))]<-"Number Canonical Matches"
x_new_new$`Number Canonical Matches`
mod1 <- glm(binary_response~`Number Canonical Matches`, family = binomial, data=x_new_new)
sjPlot::tab_model(mod1)


temp_87211 <- (create_binary_columns(x_new_87211[,-ncol(x_new_87211)]))[,(1+length(expected_gene_expression)):(2*length(expected_gene_expression))]
for(i in 1:(nrow(x_new_87211))){
  count <- sum(temp_87211[i,] == expected_gene_expression)
  x_new_87211$expected_reg[i] <- count
  print(count)
}
x_new_87211$expected_reg
x_new_87211$binary_response
pheno_87211$gender.ch1
x_new_87211$bin <- ifelse(x_new_87211$expected_reg>13, "Many Matches", "Few Matches")
x_new_87211$bin
names(x_new_87211)[grep("expected_reg", names(x_new_87211))]<-"Number Canonical Matches"
mod2 <- glm(binary_response~`Number Canonical Matches`, family = binomial, data=x_new_87211)
sjPlot::tab_model(mod2)
