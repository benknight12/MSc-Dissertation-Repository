### Compare the clusterings for each dataset
scort_nfkb_clusters <- read.csv("Processed_Data_Sets/Patient data scort.csv", row.names = 1)[,"NFkB.Cluster"]
scort_kmeans_clusters <- read.csv("2a. SCORT_Generate_Kmeans_Clusters_and_Top_Gene/cluster_list.csv",row.names = 1)[,1]
scort_pvclust_clusters <- read.csv("3a. SCORT_Generate_PV_Clusters_and_Top_Gene/pc_clust_clusters.csv",row.names = 1)[,1]

table(scort_nfkb_clusters, scort_kmeans_clusters)
table(scort_nfkb_clusters, scort_pvclust_clusters)
table(scort_kmeans_clusters, scort_pvclust_clusters)
sum(diag(table(scort_pvclust_clusters, scort_kmeans_clusters)))/sum(table(scort_pvclust_clusters, scort_kmeans_clusters))
#####

nfkb_clusters_109057 <- read.csv("1b. 109057_NFKB_cluster_derivation/NFKB_derived_scort_clusters.csv",row.names = 1)[,1]
kmeans_clusters_109057 <- read.csv("2b. 109057_Generate_Kmeans_Clusters_and_Top_Gene/clusters_lists.csv",row.names = 1)[,1]
pvclust_clusters_109057 <- read.csv("3b. 109057_Generate_PV_Clusters_and_Top_Gene/cluster_list.csv",row.names = 1)[,1]

table(nfkb_clusters_109057, kmeans_clusters_109057)
table(nfkb_clusters_109057, pvclust_clusters_109057)
table(kmeans_clusters_109057, pvclust_clusters_109057)

######

nfkb_clusters_145037 <- read.csv("1c. 145037_NFKB_cluster_derivation/NFKB_derived_scort_clusters.csv",row.names = 1)[,1]
kmeans_clusters_145037 <- read.csv("2c. 145037_Generate_Kmeans_Clusters_and_Top_Gene/Cluster_list.csv",row.names = 1)[,1]
pvclust_clusters_145037 <- read.csv("3c. 145037_Generate_PV_Clusters_and_Top_Gene/cluster_list.csv",row.names = 1)[,1]

table(nfkb_clusters_145037, kmeans_clusters_145037)
table(nfkb_clusters_145037, pvclust_clusters_145037)
table(kmeans_clusters_145037, pvclust_clusters_145037)


##( Not much agreement between kmeans and pvclust in geo ) ##

##### Now compare gene lists - using scoring system -  first just scort to find genelist for this - then add others
top_gene_dataframe <- read.csv("1a. SCORT_NFKB_cluster_derivation_replicating/Top_Gene_List.csv", row.names = 1)
kmeans_clust_list <- read.csv("2a. SCORT_Generate_Kmeans_Clusters_and_Top_Gene/Top_Gene_List.csv", row.names = 1)
pv_clust_list <- read.csv("3a. SCORT_Generate_PV_Clusters_and_Top_Gene/Top_Gene_List.csv", row.names = 1)

intersect(top_gene_dataframe$Canonical, pv_clust_list$Canonical)
intersect(kmeans_clust_list$Canonical, pv_clust_list$Canonical)


canonical <- as.data.frame(matrix(c(intersect(pv_clust_list$Canonical, kmeans_clust_list$Canonical),rep(0, length(intersect(pv_clust_list$Canonical, kmeans_clust_list$Canonical)))),ncol = 2))
noncanonical <- as.data.frame(matrix(c(intersect(pv_clust_list$NonCanonical, kmeans_clust_list$NonCanonical),rep(0, length(intersect(pv_clust_list$NonCanonical, kmeans_clust_list$NonCanonical)))),ncol = 2))
atypical <- as.data.frame(matrix(c(intersect(pv_clust_list$Atypical, kmeans_clust_list$Atypical),rep(0, length(intersect(pv_clust_list$Atypical, kmeans_clust_list$Atypical)))),ncol = 2))

intersect(pv_clust_list$Canonical, kmeans_clust_list$Canonical)
for (i in intersect(pv_clust_list$Canonical, kmeans_clust_list$Canonical)){
  if (!is.na(i)){
    value1 <- grep(i, pv_clust_list$Canonical,  fixed = TRUE)[1]
    value2 <- grep(i, kmeans_clust_list$Canonical,  fixed = TRUE)[1]
    canonical[grep(i, intersect(pv_clust_list$Canonical, kmeans_clust_list$Canonical)),2] <- value1 + value2
  }
}

canonical$V2<- as.numeric(canonical$V2)
canonical[order(canonical$V2),]


intersect(pv_clust_list$NonCanonical, kmeans_clust_list$NonCanonical)
for (i in intersect(pv_clust_list$NonCanonical, kmeans_clust_list$NonCanonical)){
  if (!is.na(i)){
    value1 <- grep(i, pv_clust_list$NonCanonical,  fixed = TRUE)[1]
    value2 <- grep(i, kmeans_clust_list$NonCanonical,  fixed = TRUE)[1]
    noncanonical[grep(i, intersect(pv_clust_list$NonCanonical, kmeans_clust_list$NonCanonical)),2] <- value1 + value2
  }
}
noncanonical$V2<- as.numeric(noncanonical$V2)
noncanonical[order(noncanonical$V2),]

for (i in intersect(pv_clust_list$Atypical, kmeans_clust_list$Atypical)){
  if (!is.na(i)){
    value1 <- grep(i, pv_clust_list$Atypical,  fixed = TRUE)[1]
    value2 <- grep(i, kmeans_clust_list$Atypical,  fixed = TRUE)[1]
    atypical[grep(i, intersect(pv_clust_list$Atypical, kmeans_clust_list$Atypical)),2] <- value1 + value2
  }
}
atypical$V2<- as.numeric(atypical$V2)
atypical[order(atypical$V2),]


top_gene_dataframe
## Compare with nfkb derived
focus_canonical <- as.data.frame(matrix(c(intersect(top_gene_dataframe$Canonical, canonical$V1), canonical$V2[canonical$V1 %in% intersect(top_gene_dataframe$Canonical, canonical$V1)]), byrow=T, nrow=2))
focus_noncanonical <- as.data.frame(matrix(c(intersect(top_gene_dataframe$NonCanonical, noncanonical$V1), noncanonical$V2[noncanonical$V1 %in% intersect(top_gene_dataframe$NonCanonical, noncanonical$V1)]), byrow=T, nrow=2))
focus_atypical <- as.data.frame(matrix(c(intersect(top_gene_dataframe$Atypical, atypical$V1), atypical$V2[atypical$V1 %in% intersect(top_gene_dataframe$Atypical, atypical$V1)]), byrow=T, nrow=2))
dim(focus_canonical)
focus_canonical
focus_noncanonical
focus_atypical

## Extract a final list of those which are common between kmeans and pvclust, then also those that are also in nfkb as 'special'
final_list <- list(canonical_all = canonical$V1, canonical_scores = canonical$V2, canonical_special = focus_canonical[1,],
                   noncanonical_all = noncanonical$V1, noncanonical_scores = noncanonical$V2, noncanonical_special = focus_noncanonical[1,],
                   atypical_all = atypical$V1, atypical_scores = atypical$V2, atypical_special = focus_atypical[1,])
as.data.frame((t(focus_canonical))[order(as.numeric(as.data.frame(t(focus_canonical))[,2])),])
as.data.frame(t(focus_noncanonical)[order(as.numeric(as.data.frame(t(focus_noncanonical))[,2])),])
as.data.frame(t(focus_atypical)[order(as.numeric(as.data.frame(t(focus_atypical))[,2])),])

final_list$canonical_all
final_list$noncanonical_all
final_list$atypical_all

final_list$canonical_special
final_list$noncanonical_special
save(final_list, file = "4. Combinations_of_Full_Data_Gene_Lists/Combined_Gene_Lists.RData")



#### GEO -- no aggreement basically so motivates the predictions from scort centroids to obtain alternate clusters
top_gene_dataframe_109 <- read.csv("1b. 109057_NFKB_cluster_derivation/Top_Gene_List.csv", row.names = 1)
kmeans_clust_list_109 <- read.csv("2b. 109057_Generate_Kmeans_Clusters_and_Top_Gene/Top_Gene_List.csv", row.names = 1)
pv_clust_list_109 <- read.csv("3b. 109057_Generate_PV_Clusters_and_Top_Gene/Top_Gene_List.csv", row.names = 1)
intersect(top_gene_dataframe_109$Canonical, pv_clust_list_109$Canonical)
intersect(top_gene_dataframe_109$Canonical,  canonical$V1)
top_gene_dataframe_109$Canonical

top_gene_dataframe_145 <- read.csv("1c. 145037_NFKB_cluster_derivation/Top_Gene_List.csv", row.names = 1)
kmeans_clust_list_145 <- read.csv("2c. 145037_Generate_Kmeans_Clusters_and_Top_Gene/Top_Gene_List.csv", row.names = 1)
pv_clust_list_145 <- read.csv("3c. 145037_Generate_PV_Clusters_and_Top_Gene/Top_Gene_List.csv", row.names = 1)

intersect(kmeans_clust_list_145$Canonical, canonical$V1)
canonical$V1
