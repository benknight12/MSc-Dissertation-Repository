## 5. Neural Network between scort genes (all) and therapy response with MIFS feature selection and obtain which genes are weighted most heavily

library(caret)
library(keras)
library(tidyverse)  
library(glmnet)
library(xgboost)

# Data prep
labelled_scort <- read.csv("CB_normalized_data/scort_normalized.csv", row.names = 1)
pheno <- read.csv("Processed_Data_Sets/scort.csv", row.names = 1)

X <- labelled_scort[, 16:ncol(labelled_scort)]
y <- pheno$TRG


y <- as.factor(y)
y_numeric <- as.numeric(y) - 1
elastic_net_model <- cv.glmnet(as.matrix(X), y_numeric, alpha = 0.5, family = "multinomial")
coef_list <- coef(elastic_net_model, s = "lambda.min")
selected_features_elastic_net <- c()
for (i in 1:length(coef_list)) {
  non_zero_coef <- rownames(coef_list[[i]])[which(coef_list[[i]] != 0)]
  non_zero_coef <- non_zero_coef[non_zero_coef != "(Intercept)"]
  selected_features_elastic_net <- unique(c(selected_features_elastic_net, non_zero_coef))
}
X_selected_elastic_net <- X[, selected_features_elastic_net]
trainIndex <- createDataPartition(y, p = .8, list = FALSE, times = 1)
X_train_selected_enet <- X_selected_elastic_net[trainIndex,]
X_test_selected_enet <- X_selected_elastic_net[-trainIndex,]
y_train <- y[trainIndex]
y_test <- y[-trainIndex]
y_train_cat <- to_categorical(as.numeric(y_train) - 1)
y_test_cat <- to_categorical(as.numeric(y_test) - 1)




### Train all models
## enet
model <- keras_model_sequential() %>%
  layer_dense(units = 64, activation = 'relu', input_shape = ncol(X_train_selected_enet)) %>%
  layer_dropout(rate = 0.5) %>%
  layer_dense(units = 32, activation = 'relu') %>%
  layer_dropout(rate = 0.5) %>%
  layer_dense(units = length(levels(y_train)), activation = 'softmax')
model %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = tf$keras$optimizers$legacy$Adam(),
  metrics = c('accuracy')
)
history_enet <- model %>% fit(
  as.matrix(X_train_selected_enet), y_train_cat,
  epochs = 100,
  batch_size = 32,
  validation_split = 0.2
)

## XGBoost
y_xgb <- as.numeric(y) - 1
dtrain <- xgb.DMatrix(data = as.matrix(X), label = y_xgb)
params <- list(
  objective = "multi:softmax",
  num_class = length(unique(y)),
  eval_metric = "mlogloss"
)
xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 100
)
importance_matrix <- xgb.importance(model = xgb_model)
selected_features_xgb <- importance_matrix$Feature
X_selected_xgb <- X[, selected_features_xgb]
X_train_selected_xgb <- X_selected_xgb[trainIndex,]
X_test_selected_xgb <- X_selected_xgb[-trainIndex,]
model_xgb <- keras_model_sequential() %>%
  layer_dense(units = 64, activation = 'relu', input_shape = ncol(X_train_selected_xgb)) %>%
  layer_dropout(rate = 0.5) %>%
  layer_dense(units = 32, activation = 'relu') %>%
  layer_dropout(rate = 0.5) %>%
  layer_dense(units = length(levels(y_train)), activation = 'softmax')
model_xgb %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = tf$keras$optimizers$legacy$Adam(),
  metrics = c('accuracy')
)
history_xgb <- model_xgb %>% fit(
  as.matrix(X_train_selected_xgb), y_train_cat,
  epochs = 100,
  batch_size = 32,
  validation_split = 0.2
)

## Neural Net
basic_model <- keras_model_sequential() %>%
  layer_dense(units = 64, activation = 'relu', input_shape = ncol(X_train)) %>%
  layer_dropout(rate = 0.5) %>%
  layer_dense(units = 32, activation = 'relu') %>%
  layer_dropout(rate = 0.5) %>%
  layer_dense(units = length(levels(y_train)), activation = 'softmax')
basic_model %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = tf$keras$optimizers$legacy$Adam(),
  metrics = c('accuracy')
)
basic_history <- basic_model %>% fit(
  as.matrix(X_train), y_train_cat,
  epochs = 100,
  batch_size = 32,
  validation_split = 0.2
)


# Evaluate the model on the test data
score_xgb <- model_xgb %>% evaluate(as.matrix(X_test_selected_xgb), y_test_cat)
cat('Test accuracy (XGBoost):',  score_xgb[["accuracy"]], "\n")


score_enet <- model %>% evaluate(as.matrix(X_test_selected_enet), y_test_cat)
cat('Test accuracy (Elastic Net):', score_enet[["accuracy"]], "\n")

basic_score <- basic_model %>% evaluate(as.matrix(X_test), y_test_cat)
cat('Test accuracy (Basic Neural Network):', basic_score[["accuracy"]], "\n")


## Get important features
important_genes_enet <- data.frame(
  gene = rownames(coef_list[[1]]),
  coefficient = coef_list[[1]][,1]
) %>%
  filter(coefficient != 0) %>%
  arrange(desc(abs(coefficient))) %>%
  head(100)  # Get top 20 most important genes
print(important_genes_enet)

important_genes_xgb <- xgb.importance(model = xgb_model) %>%
  arrange(desc(Gain)) %>%
  head(100)  # Get top 20 most important genes

print(important_genes_xgb)

weights_nn <- get_weights(basic_model)[[1]]
important_genes_nn <- data.frame(
  gene = colnames(X_train),
  importance = apply(abs(weights_nn), 1, sum)
) %>%
  arrange(desc(importance)) %>%
  head(100)  # Get top 20 most important genes

print(cbind(important_genes_enet$gene, important_genes_xgb$Feature, important_genes_nn$gene))

intersect(important_genes_enet$gene, important_genes_xgb$Feature)
intersect(important_genes_enet$gene, important_genes_nn$Feature)
intersect(important_genes_xgb$gene, important_genes_nn$Feature)

og_gene_list <- c("YME1L1", "GGH", "PSMD14", "TFB2M", "MAD2L1", "CIAO2A", "MAT2B", "SLBP", "BPNT2", "AREG", "FABP5", "LAPTM4B", "CCL20", "GPR160", "LYZ", "CALD1", "MMP12", "SCML1", "C1S", "BUB1", "MELK", "HKDC1")
og_gene_list <- c("ATP1A2","PKNOX2","CNTFR", "FYB2", "GALP","CA7", "CIMAP1D", "PSD", "PCSK2", "ANO5",  "KCNA2","PYY2",  "GP2", "PLP1", "CHRM1","HAPLN2", "ATP13A4", "RPS17","CKS2", "UTS2B","NLGN1","PPP6R1")
intersect(og_gene_list, important_genes_enet$gene)
intersect(og_gene_list, important_genes_xgb$gene)
intersect(og_gene_list, important_genes_nn$gene)


## Try fitting a new model with cluster membership to obtain genes as verification - just xgboost
X <- labelled_scort
final_cluster <- read.csv("13. Final_cluster_derivation/final_scort_clusters.csv",row.names=1)
y <- final_cluster$x

table(as.numeric(as.factor(y)))
y_xgb <- as.numeric(as.factor(y)) - 1
dtrain <- xgb.DMatrix(data = as.matrix(X), label = y_xgb)
params <- list(
  objective = "multi:softmax",
  num_class = length(unique(y)),
  eval_metric = "mlogloss"
)
xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 100
)
importance_matrix <- xgb.importance(model = xgb_model)
selected_features_xgb <- importance_matrix$Feature
X_selected_xgb <- X[, selected_features_xgb]
X_train_selected_xgb <- X_selected_xgb[trainIndex,]
X_test_selected_xgb <- X_selected_xgb[-trainIndex,]
y_test <- y[-trainIndex]
model_xgb <- keras_model_sequential() %>%
  layer_dense(units = 64, activation = 'relu', input_shape = ncol(X_train_selected_xgb)) %>%
  layer_dropout(rate = 0.5) %>%
  layer_dense(units = 32, activation = 'relu') %>%
  layer_dropout(rate = 0.5) %>%
  layer_dense(units = length(levels(y_train)), activation = 'softmax')
model_xgb %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = tf$keras$optimizers$legacy$Adam(),
  metrics = c('accuracy')
)
history_xgb <- model_xgb %>% fit(
  as.matrix(X_train_selected_xgb), y_train_cat,
  epochs = 100,
  batch_size = 32,
  validation_split = 0.2
)

xgb_probabilities <- predict(model_xgb, as.matrix(X_test_selected_xgb))
xgb_probabilities[,1:3]
if (is.matrix(xgb_probabilities[,1:3])) {
  xgb_predicted_classes <- apply(xgb_probabilities[,1:3], 1, which.max) - 1
} else {
  xgb_predicted_classes <- xgb_probabilities[,1:3]
}
xgb_predicted_classes
xgb_predicted_classes <- factor(xgb_predicted_classes, levels = c(1,2,0), labels = levels(y_test))

table(xgb_predicted_classes)
xgb_predicted_classes
y_test
confusion_matrix_xgb <- table(xgb_predicted_classes, y_test)
confusionMatrix(xgb_predicted_classes, y_test)
print(confusion_matrix_xgb)

weights_nn <- get_weights(model_xgb)[[1]]
important_genes_nn <- data.frame(
  gene = colnames(X_train_selected_xgb),
  importance = apply(abs(weights_nn), 1, sum)
) %>%
  arrange(desc(importance)) %>%
  head(100)  
print(important_genes_nn)

intersect(important_genes_nn$gene, og_gene_list)
