processed109057 <- read.csv("Processed_Data_Sets/109057_set(unscaled).csv",row.names = 1)
processed145037 <- read.csv("Processed_Data_Sets/145037_set(unscaled).csv",row.names = 1)
scort_data <- read.csv("CB_normalized_data/scort_normalized.csv",row.names = 1)
scaled_scort <- as.data.frame(scale(scort_data))
scaled_109057 <- as.data.frame(scale(processed109057))
scaled_145037 <- as.data.frame(scale(processed145037))
common_features<-intersect(colnames(scort_data),colnames(scaled_109057))

load("4. Combinations_of_Full_Data_Gene_Lists/Combined_Gene_Lists.RData")

shortestlist <- unique(unlist(c(final_list$canonical_special, final_list$atypical_special, final_list$noncanonical_special))[-1])
shortlist <- unique(c(final_list$canonical_all, final_list$noncanonical_all, final_list$atypical_all))

mini_109057 <- scaled_109057[,intersect(shortlist,names(mini_109057))]

reduced_kmeans <- kmeans(mini_109057,centers=3)
clusters <- reduced_kmeans$cluster


load("Useful_Functions/find_means.RData")
find_means(scaled_109057[,nfkb_four_genes], clusters)
## Doesn't produce great clusters


### Alternative approach using prediction from kmeans
load("Useful_Functions/Predict_Kmeans.RData")
load("2a. SCORT_Generate_Kmeans_Clusters_and_Top_Gene/kmeans_object.RData")
pheno <- read.csv("Processed_Data_Sets/Patient data scort.csv")
## refit kmeans using reduced feature set
red_set <- scaled_scort[,common_features]
dim(red_set)
red_kmeans <- kmeans(red_set, centers=3, nstart=20)
###compare to other scort kmeans
table(scort_kmeans$cluster,red_kmeans$cluster) # pretty similar
library(knitr)
library(kableExtra)
og_clus <- factor(scort_kmeans$cluster, levels = c(1,2,3), labels = c("Atypical", "Non-Canonical", "Canonical"))
red_clus <- factor(red_kmeans$cluster, levels = c(1,2,3), labels = c("Atypical", "Non-Canonical", "Canonical"))
cross_tab_df <- as.data.frame.matrix(table(og_clus, red_clus))
colnames(cross_tab_df) <- paste("Reduced Set Cluster", colnames(cross_tab_df))
rownames(cross_tab_df) <- paste("Original Cluster", rownames(cross_tab_df))
kable(cross_tab_df, caption = "Cross Tabulation of Original Clusters vs. Reduced Set Clusters") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F) %>%
  add_header_above(c(" " = 1, "Reduced Set Clusters" = ncol(cross_tab_df))) %>%
  column_spec(1, bold = TRUE)

table(red_kmeans$cluster, pheno$NFkB.Cluster) # not awful - can say 1 = Atypical, 2 = noncanonical and 3 = Canonical
table(scort_kmeans$cluster, pheno$NFkB.Cluster) # 1 = atypical, 2 = noncan, 3 = can



new_clus_109057 <- predict_kmeans(scaled_109057[,common_features], red_kmeans)
table(new_clus_109057)

load("Useful_Functions/find_means.RData")
find_means(scaled_109057[,nfkb_four_genes], new_clus_109057)
### Produces reasonable clusters with cluster 1= Canoical, 2 = Atypical, 3=Noncanonical
clusters <- factor(new_clus_109057, levels = c(1,2,3), labels = c("Canonical", "Atypical", "NonCanonical"))
library(EnhancedVolcano)
## find differential genes for thes clusters:
load("Useful_Functions/Diff_expression_and_volc_plot_from_clusters_and_data.RData")
diff_gene_and_volc_plots(full_109057, clusters, "5. GEO_exploration_of_scort_genes/109057")
dim(full_109057)



### View expression means for the genes identified in SCORT: SYNGAP1 LDHAL6A LOC105375116 IGK LOC101929400 RAC3 BRME1 OGDH NXNL1 EML4 SPMIP6
find_means(scaled_109057[,intersect(unlist(final_list$canonical_special), names(scaled_109057))], clusters)

## Although not identified as the most deferentially expressed
##plot these now:
diff_gene_and_volc_plots(full_109057, clusters, "5. GEO_exploration_of_scort_genes/109057/Special Volcano/", special_labels = unlist(final_list$canonical_all))
final_list$canonical_all
diff_gene_and_volc_plots(full_109057, clusters, "5. GEO_exploration_of_scort_genes/109057/Special Volcano/", special_labels = "LYZ", overide = T)


########## Replicate this analysis for 145037
new_clus_145037 <- predict_kmeans(scaled_145037[,common_features], red_kmeans)
table(new_clus_145037)


find_means(scaled_145037[,nfkb_four_genes], new_clus_145037)
### Produces reasonable clusters with cluster 1= Canonical, 2 = atypical, 3=NonCanonical
clusters <- factor(new_clus_145037, levels = c(1,2,3), labels = c("Canonical", "Atypical", "NonCanonical"))

## find differential genes for thes clusters:
diff_gene_and_volc_plots(full_145037, clusters, "5. GEO_exploration_of_scort_genes/145037")

### View expression means for the genes identified in SCORT: SYNGAP1 LDHAL6A LOC105375116 IGK LOC101929400 RAC3 BRME1 OGDH NXNL1 EML4 SPMIP6
find_means(scaled_145037[,intersect(unlist(final_list$canonical_special), names(scaled_145037))], clusters)
## Although not identified as the most deferentially expressed
##plot these now:
diff_gene_and_volc_plots(full_145037, clusters, "5. GEO_exploration_of_scort_genes/145037/Special Volcano/", special_labels = unlist(final_list$canonical_all))

diff_gene_and_volc_plots(full_145037, clusters, "5. GEO_exploration_of_scort_genes/145037/Special Volcano/", special_labels = c("MAT2B", "CIAO2A", "DR1", "EEIG2", "AREG"), overide = TRUE)

###### Compare these gene lists 
scort_pred_109057_gene_list <- read.csv("5. GEO_exploration_of_scort_genes/109057/Top_Gene_List.csv",row.names = 1)
scort_pred_145037_gene_list <- read.csv("5. GEO_exploration_of_scort_genes/145037/Top_Gene_List.csv",row.names = 1)
scort_pred_109057_gene_list
scort_pred_145037_gene_list
print(intersect(scort_pred_109057_gene_list$Canonical, scort_pred_145037_gene_list$Canonical))
print(intersect(final_list$canonical_all, scort_pred_145037_gene_list$Canonical))
print(intersect(final_list$canonical_all, scort_pred_109057_gene_list$Canonical))


print(intersect(scort_pred_109057_gene_list$NonCanonical, scort_pred_145037_gene_list$NonCanonical))
print(intersect(scort_pred_109057_gene_list$Atypical, scort_pred_145037_gene_list$Atypical))



load("4. Combinations_of_Full_Data_Gene_Lists/Combined_Gene_Lists.RData")
trial1 <- read.csv("2a. SCORT_Generate_Kmeans_Clusters_and_Top_Gene/Top_Gene_List.csv",row.names = 1)
trial2 <- read.csv("3a. SCORT_Generate_PV_Clusters_and_Top_Gene/Top_Gene_List.csv",row.names = 1)
trial<- rbind(trial1,trial2)

print(intersect(scort_pred_109057_gene_list$Canonical, trial$Canonical))
print(intersect(scort_pred_145037_gene_list$Canonical, trial$Canonical))

final_list$canonical_special



#### Using the special volcano plots we identify the following as verified with by either 109057 or 145037
verified_by_109057 <- c("RAC3", "OGDH", "LDHAL6A", "NXNL1")
verified_by_145037 <- c("SYNGAP1", "NXNL1", "EML4")
both_verified <- "NXNL1"

## proceed with: SYNGAP1, LDHAL6A, RAC3, OGDH, NXNL1 & EML4 
## proceed also in scort (not in geo sets) LOC105375116, IGK & LOC101929400
canonical_lists <- list(GEO_Verified = c("SYNGAP1", "LDHAL6A", "RAC3", "OGDH", "NXNL", "EML4"),
              SCORT_GEO_Verified = c("SYNGAP1", "LDHAL6A", "RAC3", "OGDH", "NXNL1", "EML4", "LOC105375116", "LOC101929400", "IGK"),
              SCORT_List_Special = c(final_list$canonical_special),
              SCORT_List_ALL = c(final_list$canonical_all))

save(canonical_lists, file = "Canonical_Final_Gene_List.RData")
