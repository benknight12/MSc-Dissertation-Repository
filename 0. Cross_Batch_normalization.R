### Combine the data sets reduce to common features then perform cross batch normalization to give more comparable datasets
scort_data <- read.csv("Processed_Data_Sets/Scort_Data.csv", row.names = 1)
full_87211 <- read.csv("Processed_Data_Sets/full_87211.csv", row.names = 1)
full_109057 <- read.csv("Processed_Data_Sets/109057_set(unscaled).csv", row.names = 1)
full_145037 <- read.csv("Processed_Data_Sets/145037_set(unscaled).csv", row.names = 1)

shared_genes <- intersect(intersect(names(scort_data), names(full_87211)), intersect(names(full_109057), names(full_145037)))

shared_genes

data_combined <- as.data.frame(rbind(cbind(scort_data[,shared_genes],batch_label = rep("scort_batch", nrow(scort_data))),
                       cbind(full_87211[,shared_genes],batch_label = rep("87211_batch", nrow(full_87211))),
                       cbind(full_109057[,shared_genes],batch_label = rep("109057_batch", nrow(full_109057))),
                       cbind(full_145037[,shared_genes],batch_label = rep("145037_batch", nrow(full_145037)))))
dim(data_combined)
head(data_combined[,18000:ncol(data_combined)])

library(sva)
batch_info <- data_combined$batch_label
combined_data_corrected <- ComBat(dat=t(data_combined[,-ncol(data_combined)]), batch=batch_info, mod=NULL)
combined_data_corrected <- t(combined_data_corrected)

## separate and then scale the data
normalised_scort <- as.data.frame(combined_data_corrected[data_combined$batch_label == "scort_batch",])
write.csv(normalised_scort, "CB_normalized_data/scort_normalized.csv")
scaled_scort <- as.data.frame(scale(normalised_scort))
write.csv(scaled_scort, "CB_normalized_data/scaled_scort.csv")

normalised_87211 <- as.data.frame(combined_data_corrected[data_combined$batch_label == "87211_batch",])
write.csv(normalised_87211, "CB_normalized_data/normalized_87211.csv")
scaled_87211 <- as.data.frame(scale(normalised_87211))
write.csv(scaled_87211, "CB_normalized_data/scaled_87211.csv")

normalised_109057 <- as.data.frame(combined_data_corrected[data_combined$batch_label == "109057_batch",])
write.csv(normalised_109057, "CB_normalized_data/normalised_109057.csv")
scaled_109057 <- as.data.frame(scale(normalised_109057))
write.csv(scaled_109057, "CB_normalized_data/scaled_109057.csv")

normalised_145037 <- as.data.frame(combined_data_corrected[data_combined$batch_label == "145037_batch",])
write.csv(normalised_145037, "CB_normalized_data/normalised_145037.csv")
scaled_145037 <- as.data.frame(scale(normalised_145037))
write.csv(scaled_145037, "CB_normalized_data/scaled_145037.csv")

scaled_data[,shared_genes]
names(scaled_data)[(shared_genes %in% names(scaled_data))==FALSE]

test <- scaled_data[,names(scaled_data) %in% shared_genes]
dim(test)

## see large difference following crossbatch normalization so should rerun
test$A1BG - scaled_scort$A1BG
