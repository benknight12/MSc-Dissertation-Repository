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


kmeans_87211_nfkb <- kmeans(nfkb_87211, centers=rbind(c(1,1,-1,-1),
                                                      c(-1,-1,1,1),
                                                      c(0,0,0,0)), nstart = 30)
save(kmeans_87211_nfkb, file = "11. GSE87211/nfkb_derived_kmeans_object.RData")
load("Useful_Functions/find_means.RData")
find_means(nfkb_87211, kmeans_87211_nfkb$cluster)
table(kmeans_87211_nfkb$cluster)

### The means suggest cluster 1 is Canonical, cluster 2 is NonCanonical and cluster 3 is Atypical
nfkb_derived_clusters_87211 <- factor(kmeans_87211_nfkb$cluster, levels=c(1,2,3), labels = c("Canonical", "NonCanonical", "Atypical"))

## Save the clusters
write.csv(nfkb_derived_clusters_87211, "11. GSE87211//NFKB_derived_scort_clusters.csv")

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
load("Useful_Functions/Diff_expression_and_volc_plot_from_clusters_and_data.RData")




## run kmeans clustering and analysis



## run pvclust and analysis









## Investigate genes already identified (new script)
