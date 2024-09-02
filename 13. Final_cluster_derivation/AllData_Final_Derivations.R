scort_data <- read.csv("CB_normalized_data/scort_normalized.csv", row.names = 1)
full_109057 <- read.csv("CB_normalized_data/normalised_109057.csv", row.names = 1)
full_145037 <- read.csv("CB_normalized_data/normalised_145037.csv", row.names = 1)
full_87211 <- read.csv("CB_normalized_data/normalized_87211.csv", row.names = 1)
scaled_scort <- as.data.frame(scale(scort_data))
scaled_109057 <- as.data.frame(scale(full_109057))
scaled_145037 <- as.data.frame(scale(full_145037))
scaled_87211 <- as.data.frame(scale(full_87211))

centroids2 <- rbind(c(1,-1,-1), c(-1,1,1), c(0,0,0))
centroids <- rbind(c(1,1,-1,-1),c(-1,-1,1,1), c(0,0,0,0))

## SCORT 
# Check cluster issues for RELA
new <- kmeans(scaled_scort[,nfkb_four_genes],3, iter.max=20)
full_kmeans_scort <- kmeans(scaled_scort, centers = 3, nstart = 20)
dist_mat <- dist(scaled_scort, method = 'euclidean')
full_pv_scort <- hclust(dist_mat, method = "ward.D2", members = NULL)

find_means(scaled_scort[,nfkb_four_genes], new$cluster) 
find_means(scaled_scort[,nfkb_four_genes], full_kmeans_scort$cluster)
find_means(scaled_scort[,nfkb_four_genes], cutree(full_pv_scort, k=3))


## View clusters pretty good for other 3 genes but persistant issue with RELA
# Cluster on other 3
with3_kmeans <-  kmeans(scaled_scort[,c("NFKB1","NFKB2","RELB")], centers = centroids2, nstart = 20)
find_means(scaled_scort[,nfkb_four_genes], with3_kmeans$cluster)
clusters <- factor(with3_kmeans$cluster, c(1,2,3), c("Canonical", "NonCanonical", "Atypical"))
veccan <- c(scaled_scort[clusters=="Canonical", "RELA"])
vecnon <- c(scaled_scort[clusters=="NonCanonical", "RELA"])
vecaty <- c(scaled_scort[clusters=="Atypical", "RELA"])
df <- data.frame(value = c(veccan, vecnon, vecaty),  group = factor(rep(c("Canonical", "NonCanonical", "Atypical"),                 times = c(length(veccan), length(vecnon), length(vecaty)))))

ggplot(df, aes(x = group, y = value, fill = group)) +
  geom_violin(trim = FALSE) +
  labs(title = "Violin Plots of Different Length Vectors",
       x = "Vectors",
       y = "Values") +
  theme_minimal()
## See still lots of issues with RELA distribution so reclassify the Canonical and Non-Canonical for those with 

# Redefine wrong RELA expression with a defined tolerance (0.5)
scaled_scort[clusters=="Canonical" & scaled_scort$RELA < -0.5, nfkb_four_genes]
clusters[clusters=="Canonical" & scaled_scort$RELA < -0.5] <- "Atypical"
scaled_scort[(clusters=="NonCanonical" & scaled_scort$RELA > 0.5), nfkb_four_genes] 
clusters[clusters=="NonCanonical" & scaled_scort$RELA > 0.5] <- "Atypical"
clusters[c(62,86,93,170,135,34,70)] <- "Canonical"
clusters[c(9)] <- "NonCanonical"
find_means(scaled_scort[,nfkb_four_genes], clusters)

veccan <- c(scaled_scort[clusters=="Canonical", "RELA"])
vecnon <- c(scaled_scort[clusters=="NonCanonical", "RELA"])
vecaty <- c(scaled_scort[clusters=="Atypical", "RELA"])
df <- data.frame(value = c(veccan, vecnon, vecaty),  group = factor(rep(c("Canonical", "NonCanonical", "Atypical"),                 times = c(length(veccan), length(vecnon), length(vecaty)))))

ggplot(df, aes(x = group, y = value, fill = group)) +
  geom_violin(trim = FALSE) +
  labs(title = "Violin Plots of Different Length Vectors",
       x = "Vectors",
       y = "Values") +
  theme_minimal()
## Better
final_scort_clusters <- clusters
write.csv(final_scort_clusters, "13. Final_cluster_derivation/final_scort_clusters.csv")



table(final_scort_clusters)
fixed_centroids <- as.data.frame(find_means(scaled_scort, clusters = final_scort_clusters))
head(fixed_centroids[,1:5])
custom_kmeans <- list(
  centers = fixed_centroids,
  cluster = final_scort_clusters, totss = NULL, withinss = NULL, tot.withinss = NULL, betweenss = NULL, size = as.numeric(table(scaled_scort$cluster)),iter = 1, ifault = 0
)
class(custom_kmeans) <- "kmeans"


# GEO 109057
new <- predict(custom_kmeans, scaled_109057)

find_means(scaled_109057[,nfkb_four_genes], new) 
## View clusters pretty good for kmeans on whole set - only RELA mean is dodgy for noncanonical and canonical
clusters <- factor(new, c(1,2,3), c("Canonical", "NonCanonical","Atypical"))
table(clusters)
# tolerance 0.2
scaled_109057[clusters=="Canonical" & scaled_109057$RELA < -0.5, nfkb_four_genes]
scaled_109057[clusters=="NonCanonical" & scaled_109057$RELA > 0.5, nfkb_four_genes]
clusters[clusters=="Canonical" & scaled_109057$RELA < -0.5] <- "Atypical"
clusters[clusters=="NonCanonical" & scaled_109057$RELA > 0.5] <- "Atypical"

find_means(scaled_109057[,nfkb_four_genes], clusters)
table(clusters)

## Better
final_109057_clusters <- clusters
write.csv(final_109057_clusters, "13. Final_cluster_derivation/final_109057_clusters.csv")





# GEO 145037
new <- predict(custom_kmeans, scaled_145037)

find_means(scaled_145037[,nfkb_four_genes], new) 
## View clusters pretty good for kmeans on whole set - only RELA mean is dodgy for Canonical so look at these
clusters <- factor(new, c(1,2,3), c("Canonical", "NonCanonical","Atypical"))
table(clusters)
# tolerance 0.2
scaled_145037[clusters=="Canonical" & scaled_145037$NFKB1 < -0.5, nfkb_four_genes]
clusters[clusters=="Canonical" & scaled_145037$NFKB1 < -0.5] <- "Atypical"

find_means(scaled_145037[,nfkb_four_genes], clusters)
table(clusters)

## Better
final_145037_clusters <- clusters
write.csv(final_145037_clusters, "13. Final_cluster_derivation/final_145037_clusters.csv")




# GEO 87211
new <- predict(custom_kmeans, scaled_87211)
find_means(scaled_87211[nfkb_four_genes], new)
new <- kmeans(scaled_87211[,nfkb_four_genes], centers = rbind(c(1,1,-1,-1),c(-1,-1,1,1),c(0,0,0,0)), nstart = 20)
find_means(scaled_87211[,nfkb_four_genes], new$cluster) 

## Clusters need work - RELA and RELB are bad
clusters <- factor(new$cluster, c(1,2,3), c("Canonical", "NonCanonical","Atypical"))
table(clusters)
# tolerance 0.5 - issues with canonical nfkb 1 and noncanonical relb
scaled_87211[clusters=="Canonical" & scaled_87211$NFKB1 < -0.5, nfkb_four_genes]
clusters[clusters=="Canonical" & scaled_87211$NFKB1 < -0.5] <- "Atypical"

scaled_87211[clusters=="NonCanonical" & scaled_87211$RELB < -0.5, nfkb_four_genes]
clusters[clusters=="NonCanonical" & scaled_87211$RELB < -0.5] <- "Atypical"

scaled_87211[clusters=="Canonical" & scaled_87211$RELA < -0.5, nfkb_four_genes]
clusters[clusters=="Canonical" & scaled_87211$RELA < -0.5] <- "Atypical"

scaled_87211[clusters=="NonCanonical" & scaled_87211$NFKB2 < -0.5, nfkb_four_genes]
clusters[clusters=="NonCanonical" & scaled_87211$RELB < -0.5] <- "Atypical"
## some manual adjustment
clusters[c("GSM2325301", "GSM2325504", "GSM2325297", "GSM2325418", "GSM2325549", "GSM2325551", "GSM2325538")] <- "Canonical"



find_means(as.data.frame(scaled_87211[,similar[!is.na(similar)]]), clusters)

find_means(scaled_87211[,nfkb_four_genes], clusters)
table(clusters)

## Better
final_87211_clusters <- clusters
length(final_87211_clusters)
write.csv(final_87211_clusters, "13. Final_cluster_derivation/final_87211_clusters.csv")

install.packages("citation")
library(citation)
x <- citation("pvclust")
toBibtex(citation("pvclust"))
toBibtex(citation("oligo"))
toBibtex(citation("AnnotationForge"))


pheno_87211 <- read.csv("Processed_Data_Sets/pheno87211.csv",row.names = 1)
mean(pheno_87211[complete.cases(pheno_87211)&clusters == "Canonical","depth.of.invasion.after.rct.ch1"])
mean(pheno_87211[complete.cases(pheno_87211)&clusters == "NonCanonical","depth.of.invasion.after.rct.ch1"])
mean(pheno_87211[complete.cases(pheno_87211)&clusters == "Atypical","depth.of.invasion.after.rct.ch1"])
pheno$TRG
complete.cases(pheno$TRG)
trgs <-5- as.numeric(as.factor(pheno$TRG))
mean(trgs[!is.na(trgs)&final_clusters == "Canonical"])
mean(trgs[!is.na(trgs)&final_clusters == "NonCanonical"])
mean(trgs[!is.na(trgs)&final_clusters == "Atypical"])

mean(pheno$OS[!is.na(pheno$OS)&final_clusters == "Canonical"])
mean(pheno$OS[!is.na(pheno$OS)&final_clusters == "NonCanonical"])
mean(pheno$OS[!is.na(pheno$OS)&final_clusters == "Atypical"])



