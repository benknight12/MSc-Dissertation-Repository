processed109057 <- read.csv("Processed_Data_Sets/109057_set(unscaled).csv",row.names = 1)
processed145037 <- read.csv("Processed_Data_Sets/145037_set(unscaled).csv",row.names = 1)
scort_data <- read.csv("Processed_Data_Sets/SCORT_final_data.csv",row.names = 1)
scaled_scort <- as.data.frame(scale(scort_data))
scaled_109057 <- as.data.frame(scale(processed109057))
scaled_145037 <- as.data.frame(scale(processed145037))

pheno <- read.csv("Processed_Data_Sets/Patient data scort.csv")


## obtain geo pheno data
library(GEOquery)
gset <- getGEO("GSE109057", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL15207", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
pheno_109057 <- pData(gset)
relevant_colnames_109057 <- c("characteristics_ch1.1", "characteristics_ch1.2", "contact_city", "contact_country", "age:ch1", "Sex:ch1", "tissue:ch1")
pheno_109057 <- pheno_109057[,relevant_colnames_109057]

gset <- getGEO("GSE145037", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL15207", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
pheno_145037 <- pData(gset)
relevant_colnames_145037 <- c("age:ch1", "clincal n positive:ch1", "clincal t stage:ch1", "response to the crt:ch1", "Sex:ch1")
pheno_145037 <- pheno_145037[,relevant_colnames_145037]

scort_relevant_columns <- c("NFkB.Cluster", "Gender", "Age", "DFS", "DFS.Status", "OS", "OS.Status", "Hypoxia.Score", "Meth.Global.Mean", "PretreatmentMStage", "PretreatmentNStage", "PretreatmentTStage", "RNAProliferationScore", "TRG", "Binary_response")
pheno_scort <- pheno[,scort_relevant_columns]


## Merge to form 3 supervised learning sets
labelled_scort <- cbind(pheno_scort, scaled_scort)
labelled_109057 <- cbind(pheno_109057, scaled_109057)
labelled_145037 <- cbind(pheno_145037, scaled_145037) 

write.csv(labelled_scort, "Supervised_Sets/scort_supervised.csv")
write.csv(labelled_109057, "Supervised_Sets/109057_supervised.csv")
write.csv(labelled_145037, "Supervised_Sets/145037_supervised.csv")






  