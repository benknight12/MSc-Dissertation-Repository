## 1. Mean survival and recurrence/TRG for the clusters - each dataset (just scort for now)

labelled_scort <- read.csv("Supervised_Sets/scort_supervised.csv",row.names = 1)

head(labelled_scort[,1:20])

scort_overallsurvival_means_adam <- data.frame(Canonical = mean(labelled_scort$OS[labelled_scort$NFkB.Cluster == "Canonical"]),
                                          NonCanonical = mean(labelled_scort$OS[labelled_scort$NFkB.Cluster == "Non-Canonical"]),
                                          Atypical = mean(labelled_scort$OS[labelled_scort$NFkB.Cluster == "Atypical"]))
scort_overallsurvival_means_adam

kmeans_whole_clusters <- read.csv("2a. SCORT_Generate_Kmeans_Clusters_and_Top_Gene/cluster_list.csv", row.names = 1)
scort_overallsurvival_means_kmeans <- data.frame(Canonical = mean(labelled_scort$OS[kmeans_whole_clusters == "Canonical"]),
                                               NonCanonical = mean(labelled_scort$OS[kmeans_whole_clusters == "NonCanonical"]),
                                               Atypical = mean(labelled_scort$OS[kmeans_whole_clusters == "Atypical"]))
scort_overallsurvival_means_kmeans

pvclust_whole_clusters <- read.csv("3a. SCORT_Generate_PV_Clusters_and_Top_Gene/pc_clust_clusters.csv", row.names = 1)
scort_overallsurvival_means_pvclust <- data.frame(Canonical = mean(labelled_scort$OS[pvclust_whole_clusters == "Canonical"]),
                                                 NonCanonical = mean(labelled_scort$OS[pvclust_whole_clusters == "NonCanonical"]),
                                                 Atypical = mean(labelled_scort$OS[pvclust_whole_clusters == "Atypical"]))
scort_overallsurvival_means_pvclust

all_means <- round(as.data.frame(rbind(scort_overallsurvival_means_adam, scort_overallsurvival_means_kmeans, scort_overallsurvival_means_pvclust, scort_finalclusters_means), row.names = c("NFkB Derived", "Kmeans", "PVClust", "Final")),3)
kable(all_means, caption = "Survival Means between clusters") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F) %>%
  column_spec(1, bold = TRUE)

clusters <- labelled_scort$NFkB.Cluster
options <- levels(as.factor(labelled_scort$TRG))
scort_trg_table_adam <- as.data.frame(rbind(c(table(factor(labelled_scort$TRG[clusters == "Canonical"], levels = options)), proportion_pCR = table(factor(labelled_scort$TRG[clusters == "Canonical"], levels = options))[4]/sum(table(factor(labelled_scort$TRG[clusters == "Canonical"], levels = options))), proportion_good_response = table(labelled_scort$Binary_response[clusters == "Canonical"])[2]/sum(table(labelled_scort$Binary_response[clusters == "Canonical"]))),
                                            c(table(factor(labelled_scort$TRG[clusters == "NonCanonical"], levels = options)), proportion_pCR = table(factor(labelled_scort$TRG[clusters == "NonCanonical"], levels = options))[4]/sum(table(factor(labelled_scort$TRG[clusters == "NonCanonical"], levels = options))), proportion_good_response = table(labelled_scort$Binary_response[clusters == "NonCanonical"])[2]/sum(table(labelled_scort$Binary_response[clusters == "NonCanonical"]))),
                                            c(table(factor(labelled_scort$TRG[clusters == "Atypical"], levels = options)), proportion_pCR = table(factor(labelled_scort$TRG[clusters == "Atypical"], levels = options))[4]/sum(table(factor(labelled_scort$TRG[clusters == "Atypical"], levels = options))), proportion_good_response = table(labelled_scort$Binary_response[clusters == "Atypical"])[2]/sum(table(labelled_scort$Binary_response[clusters == "Atypical"])))),
                                      row.names=c("Canonical", "NonCanonical", "Atypical"))
scort_trg_table_adam

clusters <- kmeans_whole_clusters
scort_trg_table_kmeans <- as.data.frame(rbind(c(table(factor(labelled_scort$TRG[clusters == "Canonical"], levels = options)), proportion_pCR = table(factor(labelled_scort$TRG[clusters == "Canonical"], levels = options))[4]/sum(table(factor(labelled_scort$TRG[clusters == "Canonical"], levels = options))), proportion_good_response = table(labelled_scort$Binary_response[clusters == "Canonical"])[2]/sum(table(labelled_scort$Binary_response[clusters == "Canonical"]))),
                                              c(table(factor(labelled_scort$TRG[clusters == "NonCanonical"], levels = options)), proportion_pCR = table(factor(labelled_scort$TRG[clusters == "NonCanonical"], levels = options))[4]/sum(table(factor(labelled_scort$TRG[clusters == "NonCanonical"], levels = options))), proportion_good_response = table(labelled_scort$Binary_response[clusters == "NonCanonical"])[2]/sum(table(labelled_scort$Binary_response[clusters == "NonCanonical"]))),
                                              c(table(factor(labelled_scort$TRG[clusters == "Atypical"], levels = options)), proportion_pCR = table(factor(labelled_scort$TRG[clusters == "Atypical"], levels = options))[4]/sum(table(factor(labelled_scort$TRG[clusters == "Atypical"], levels = options))), proportion_good_response = table(labelled_scort$Binary_response[clusters == "Atypical"])[2]/sum(table(labelled_scort$Binary_response[clusters == "Atypical"])))),
                                        row.names=c("Canonical", "NonCanonical", "Atypical"))
scort_trg_table_kmeans

clusters <- pvclust_whole_clusters
scort_trg_table_pvclust <- as.data.frame(rbind(c(table(factor(labelled_scort$TRG[clusters == "Canonical"], levels = options)), proportion_pCR = table(factor(labelled_scort$TRG[clusters == "Canonical"], levels = options))[4]/sum(table(factor(labelled_scort$TRG[clusters == "Canonical"], levels = options))), proportion_good_response = table(labelled_scort$Binary_response[clusters == "Canonical"])[2]/sum(table(labelled_scort$Binary_response[clusters == "Canonical"]))),
                                               c(table(factor(labelled_scort$TRG[clusters == "NonCanonical"], levels = options)), proportion_pCR = table(factor(labelled_scort$TRG[clusters == "NonCanonical"], levels = options))[4]/sum(table(factor(labelled_scort$TRG[clusters == "NonCanonical"], levels = options))), proportion_good_response = table(labelled_scort$Binary_response[clusters == "NonCanonical"])[2]/sum(table(labelled_scort$Binary_response[clusters == "NonCanonical"]))),
                                               c(table(factor(labelled_scort$TRG[clusters == "Atypical"], levels = options)), proportion_pCR = table(factor(labelled_scort$TRG[clusters == "Atypical"], levels = options))[4]/sum(table(factor(labelled_scort$TRG[clusters == "Atypical"], levels = options))), proportion_good_response = table(labelled_scort$Binary_response[clusters == "Atypical"])[2]/sum(table(labelled_scort$Binary_response[clusters == "Atypical"])))),
                                         row.names=c("Canonical", "NonCanonical", "Atypical"))
scort_trg_table_pvclust



##### Apply to final clusterings determined in 12
final_clusters <- read.csv("12. gene list comp/Final_scort_Cluster_List.csv",row.names = 1)
scort_finalclusters_means <- data.frame(Canonical = mean(labelled_scort$OS[final_clusters$x == "Canonical"]),
                                        NonCanonical = mean(labelled_scort$OS[final_clusters$x == "NonCanonical"]),
                                        Atypical = mean(labelled_scort$OS[final_clusters$x == "Atypical"]))
scort_finalclusters_means

scort_finalclusters_trg_proportions <- scort_trg_table_pvclust <- as.data.frame(rbind(c(table(factor(labelled_scort$TRG[final_clusters$x == "Canonical"], levels = options)), proportion_pCR = table(factor(labelled_scort$TRG[final_clusters$x == "Canonical"], levels = options))[4]/sum(table(factor(labelled_scort$TRG[final_clusters$x == "Canonical"], levels = options))), proportion_good_response = table(labelled_scort$Binary_response[final_clusters$x == "Canonical"])[2]/sum(table(labelled_scort$Binary_response[final_clusters$x == "Canonical"]))),
                                                                                      c(table(factor(labelled_scort$TRG[final_clusters$x == "NonCanonical"], levels = options)), proportion_pCR = table(factor(labelled_scort$TRG[final_clusters$x == "NonCanonical"], levels = options))[4]/sum(table(factor(labelled_scort$TRG[final_clusters$x == "NonCanonical"], levels = options))), proportion_good_response = table(labelled_scort$Binary_response[final_clusters$x == "NonCanonical"])[2]/sum(table(labelled_scort$Binary_response[final_clusters$x == "NonCanonical"]))),
                                                                                      c(table(factor(labelled_scort$TRG[final_clusters$x == "Atypical"], levels = options)), proportion_pCR = table(factor(labelled_scort$TRG[final_clusters$x == "Atypical"], levels = options))[4]/sum(table(factor(labelled_scort$TRG[final_clusters$x == "Atypical"], levels = options))), proportion_good_response = table(labelled_scort$Binary_response[final_clusters$x == "Atypical"])[2]/sum(table(labelled_scort$Binary_response[final_clusters$x == "Atypical"])))),
                                                                                row.names=c("Canonical", "NonCanonical", "Atypical"))
scort_finalclusters_trg_proportions

all_means_and_tables <- list(scort_overallsurvival_means_adam = scort_overallsurvival_means_adam, scort_overallsurvival_means_kmeans = scort_overallsurvival_means_kmeans, scort_overallsurvival_means_pvclust = scort_overallsurvival_means_pvclust, final_scort_means = scort_finalclusters_means,
                             scort_trg_table_adam = scort_trg_table_adam, scort_trg_table_kmeans = scort_trg_table_kmeans, scort_trg_table_pvclust=scort_trg_table_pvclust, final_trg_table = scort_finalclusters_trg_proportions)
scort_trg_table_adam$proportion_good_response.response

trg_proportions <- round(as.data.frame(rbind(scort_trg_table_adam$proportion_good_response.response, scort_trg_table_kmeans$proportion_good_response.response, scort_trg_table_pvclust$proportion_good_response.response, scort_finalclusters_trg_proportions$proportion_good_response.response), row.names = c("NFkB Derived", "Kmeans", "PVClust", "Final")),3)
colnames(trg_proportions) <- c("Canonical", "Non-Canonical", "Atypical")
kable(trg_proportions, caption = "TRG proportions between clusters") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F) %>%
  column_spec(1, bold = TRUE)
save(all_means_and_tables, file="6. Survival_and_TRG_cluster_analysis/Means_and_Tables.RData")



### For GSE87211


full_87211 <- read.csv("Processed_Data_Sets/full_87211.csv",row.names = 1)
pheno87211 <- read.csv("Processed_Data_Sets/pheno87211.csv",row.names = 1)

head(pheno87211)
relevant_columns <- c("disease.free.time..month..ch1", "cancer.recurrance.after.surgery.ch1","death.due.to.tumor.ch1", "survival.time..month..ch1")

labelled_87211 <- as.data.frame(cbind(pheno87211[,relevant_columns], full_87211))
names(labelled_87211)[1:4] <- c("DFS", "DFS.Status", "OS.Status", "OS")

exclude <-unlist(lapply(labelled_87211$OS, is.na))
exclude
labelled_87211 <- labelled_87211[!exclude, ]
table(is.na(labelled_87211$OS))
dim(labelled_87211[,1:20])

clusters <- read.csv("12. gene list comp/Final_87211_Clusters.csv", row.names = 1)
clusters <- as.factor(clusters$x[!exclude])
levels(clusters)
overallsurvival_means_87211 <- data.frame(Canonical = mean(labelled_87211$OS[clusters == "Canonical"]),
                                               NonCanonical = mean(labelled_87211$OS[clusters == "NonCanonical"]),
                                               Atypical = mean(labelled_87211$OS[clusters == "Atypical"]))
table(is.na(labelled_87211$OS))
overallsurvival_means_87211
labelled_87211$OS[clusters == "Canonical"]

mean(labelled_87211$OS[clusters == "NonCanonical"])
mean(labelled_87211$OS[clusters == "Canonical"])
mean(labelled_87211$OS[clusters == "Atypical"])
dim(labelled_87211)
write.csv(labelled_87211, "Supervised_Sets/87211.csv")
write.csv(clusters, "12. gene list comp/87211_reduced_clusters.csv")
