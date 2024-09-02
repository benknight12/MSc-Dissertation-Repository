## Run with ADAM clusters
pheno_scort <- read.csv("Processed_Data_Sets/scort.csv",row.names = 1)
scaled_109057 <- read.csv("CB_normalized_data/scaled_109057.csv", row.names = 1)
scaled_145037 <- read.csv("CB_normalized_data/scaled_145037.csv", row.names = 1)
scaled_87211 <- read.csv("CB_normalized_data/scaled_87211.csv", row.names = 1)
scort_data <- read.csv("CB_normalized_data/scort_normalized.csv", row.names = 1)
full_109057 <- read.csv("CB_normalized_data/normalised_109057.csv", row.names = 1)
full_145037 <- read.csv("CB_normalized_data/normalised_145037.csv", row.names = 1)
full_87211 <- read.csv("CB_normalized_data/normalized_87211.csv", row.names = 1)

clusters <- pheno_scort$NFkB.Cluster
pheno_scort$NFkB.Cluster
scort_adam_clusters <- factor(pheno_scort$NFkB.Cluster, c("Atypical", "Canonical", "Non-Canonical"), c("Atypical", "Canonical", "NonCanonical"))
clusters

fixed_centroids <- as.data.frame(find_means(scaled_scort, clusters = clusters))
custom_kmeans <- list(
  centers = fixed_centroids,
  cluster = clusters, totss = NULL, withinss = NULL, tot.withinss = NULL, betweenss = NULL, size = as.numeric(table(scaled_scort$cluster)),iter = 1, ifault = 0
)
class(custom_kmeans) <- "kmeans"
library(fdm2id)
dim(scaled_109057)

new <- predict(custom_kmeans, scaled_109057)
find_means(scaled_109057[,nfkb_four_genes], new) 
## View clusters pretty good for kmeans on whole set - only RELA mean is dodgy for noncanonical and canonical
clusters <- factor(new, c(1,2,3), c("Atypical", "Canonical","NonCanonical"))
scaled_109057[clusters=="Canonical" & scaled_109057$RELA < -0.5, nfkb_four_genes]
scaled_109057[clusters=="NonCanonical" & scaled_109057$RELA > 0.5, nfkb_four_genes]
clusters[clusters=="Canonical" & scaled_109057$RELA < -0.5] <- "Atypical"
clusters[clusters=="NonCanonical" & scaled_109057$RELA > 0.5] <- "Atypical"
find_means(scaled_109057[,nfkb_four_genes], clusters)
## Better
final_109057_clusters <- clusters


new <- predict(custom_kmeans, scaled_145037)
find_means(scaled_145037[,nfkb_four_genes], new) 
## View clusters pretty good for kmeans on whole set - only RELA mean is dodgy for noncanonical and canonical
clusters <- factor(new, c(1,2,3), c("Atypical", "Canonical","NonCanonical"))
scaled_145037[clusters=="Canonical" & scaled_145037$NFKB1 < -0.5, nfkb_four_genes]
clusters[clusters=="Canonical" & scaled_145037$NFKB1 < -0.5] <- "Atypical"
find_means(scaled_145037[,nfkb_four_genes], clusters)
## Better
final_145037_clusters <- clusters

new <- predict(custom_kmeans, scaled_87211)
find_means(scaled_87211[,nfkb_four_genes], new) 
## View clusters pretty good for kmeans on whole set - only RELA mean is dodgy for noncanonical and canonical
clusters <- factor(new, c(1,2,3), c("Atypical", "Canonical","NonCanonical"))
clusters[clusters=="Canonical" & scaled_87211$NFKB1 < -0.5] <- "Atypical"
clusters[clusters=="NonCanonical" & scaled_87211$RELB < -0.5] <- "Atypical"
clusters[clusters=="Canonical" & scaled_87211$RELA < -0.5] <- "Atypical"
clusters[clusters=="NonCanonical" & scaled_87211$RELB < -0.5] <- "Atypical"
find_means(scaled_87211[,nfkb_four_genes], clusters)
## Better
final_87211_clusters <- clusters

library(EnhancedVolcano)
load("Useful_Functions/Diff_expression_and_volc_plot_from_clusters_and_data.RData")
diff_gene_and_volc_plots(full_109057, as.factor(final_109057_clusters), "15. gene lists from adam/109057/", lfc_boundary = log(1.85,base=2),special_labels = nfkb_four_genes, overide = T)
diff_gene_and_volc_plots(full_145037, as.factor(final_145037_clusters), "15. gene lists from adam/145037/", lfc_boundary = log(1.85,base=2),special_labels = nfkb_four_genes, overide = T)
diff_gene_and_volc_plots(full_87211, as.factor(final_87211_clusters), "15. gene lists from adam/87211/", lfc_boundary = log(1.85,base=2),special_labels = nfkb_four_genes, overide = T)
diff_gene_and_volc_plots(scort_data, as.factor(scort_adam_clusters), "15. gene lists from adam/SCORT/", lfc_boundary = log(1.85,base=2),special_labels = nfkb_four_genes, overide = T)

signif_109057 <- read.csv("15. gene lists from adam/109057/Top_Gene_List.csv",row.names = 1)
signif_87211 <- read.csv("15. gene lists from adam/87211/Top_Gene_List.csv",row.names = 1)
conclusive_scort_genes <- read.csv("15. gene lists from adam/SCORT/Top_Gene_List.csv",row.names = 1)
signif_145037 <- read.csv("15. gene lists from adam/145037//Top_Gene_List.csv",row.names = 1)

sum(!is.na(conclusive_scort_genes$Canonical))
sum(!is.na(signif_109057$Canonical))
length(intersect(conclusive_scort_genes$Canonical, signif_109057$Canonical))
sum(!is.na(signif_145037$Canonical))
length(intersect(conclusive_scort_genes$Canonical, signif_145037$Canonical))
length(intersect(intersect(conclusive_scort_genes$Canonical, signif_145037$Canonical),signif_109057$Canonical))
sum(!is.na(signif_87211$Canonical))
length(intersect(conclusive_scort_genes$Canonical, signif_87211$Canonical))

similar1 <- intersect(intersect(signif_109057$Canonical, signif_87211$Canonical), intersect(conclusive_scort_genes$Canonical, signif_145037$Canonical))
similar1


head(full_87211[,(ncol(full_87211)-1:ncol(full_87211))])

head(labelled_87211[,1:5])
dim(labelled_87211)


lapply(reduced_pheno,function(x) {table(is.na(x))})
names(pheno_87211)
reduced_pheno <- pheno_87211[,c("depth.of.invasion.before.rct.ch1","cancer.recurrance.after.surgery.ch1","death.due.to.tumor.ch1","disease.free.time..month..ch1", "survival.time..month..ch1")]
names(reduced_pheno) <- c("TRG", "DFS.Status", "OS.Status", "DFS", "OS")
head(reduced_pheno)
reduced_pheno$cluster <- final_87211_clusters


table(reduced_pheno[is.na(reduced_pheno)])
sum(is.na(reduced_pheno)

    