labelled_scort <- read.csv("CB_normalized_data/scort_normalized.csv", row.names = 1)
pheno <- read.csv("Processed_Data_Sets/scort.csv", row.names = 1)
Xtrain <- scale(labelled_scort)
ytrain <- as.factor(pheno$TRG)

full_87211 <- read.csv("../scaled_87211.csv", row.names = 1)
pheno_87211 <- read.csv("Processed_Data_Sets/pheno87211.csv", row.names = 1)
Xtest <- scale(full_87211)
ytest <- as.factor(pheno_87211$depth.of.invasion.after.rct.ch1)
ytest[ytest==1] <- 2
ytest <- factor(ytest, c(0,2,3,4), c("pCR","Good partial response", "Partial response","Minimal or no response"))
table(ytest)/length(ytest)
table(ytrain)/length(ytrain)
Xtrain <- as.matrix(Xtrain)
Xtest <- as.matrix(Xtest)
table(ytest)
ytrain_binary <- pheno$Binary_response
ytest_binary <- rep("Positive_Response",length(ytest))
ytest_binary[ytest%in% c("Partial response", "Minimal or no response")] <- "Negative_Response"

variances <- apply(Xtrain, 2, var)
top_features <- order(variances, decreasing = TRUE)[1:2000]
Xtrain_reduced <- Xtrain[, top_features]
Xtest_reduced <- Xtest[, top_features]





cv_fit <- cv.glmnet(Xtrain, ytrain_binary, family = "multinomial", alpha = alpha_value, type.measure = "class")
best_lambda <- cv_fit$lambda.min
elastic_net_model <- glmnet(Xtrain, ytrain_binary, family = "multinomial", alpha = alpha_value, lambda = best_lambda)
table(ytrain)

predictions <- predict(elastic_net_model, newx = Xtest, s = "lambda.min", type = "class")
table(predictions)
table(ytest_binary)
# Confusion matrix
confusion_matrix <- table(predictions, ytest_binary)
print("Confusion Matrix:")
print(confusion_matrix)

accuracy <- (confusion_matrix[1,1]+ confusion_matrix[2,2])/sum(confusion_matrix)
accuracy
intersect(gene_list, names(labelled_scort)[top_features])







### just use genelist

Xtrain_genelist <- Xtrain[, gene_list]
Xtest_genelist <- Xtest[, gene_list]

alpha_grid <- seq(0, 1, by = 0.1)
cv_results <- list()
best_lambda_per_alpha <- numeric(length(alpha_grid))
alpha_best <- NULL
min_cv_error <- Inf
for (alpha_value in alpha_grid) {
  # Train Elastic Net model using cross-validation
  cv_fit <- cv.glmnet(Xtrain_genelist, ytrain_binary, family = "multinomial", alpha = alpha_value, type.measure = "class")
  
  # Store the best lambda for this alpha value
  best_lambda_per_alpha[which(alpha_grid == alpha_value)] <- cv_fit$lambda.min
  
  # Check if the current model has the lowest CV error
  if (cv_fit$cvm[cv_fit$lambda == cv_fit$lambda.min] < min_cv_error) {
    min_cv_error <- cv_fit$cvm[cv_fit$lambda == cv_fit$lambda.min]
    alpha_best <- alpha_value
  }
  
  # Store results
  cv_results[[as.character(alpha_value)]] <- cv_fit
}
best_cv_fit <- cv_results[[as.character(alpha_best)]]
best_lambda <- best_cv_fit$lambda.min
alpha_best
# Train the final Elastic Net model with the best alpha and lambda
elastic_net_model <- glmnet(Xtrain_genelist, ytrain_binary, family = "multinomial", alpha = alpha_best, lambda = best_lambda)

# cv_fit <- cv.glmnet(Xtrain_genelist, ytrain_binary, family = "multinomial", alpha = alpha_value, type.measure = "class")
# best_lambda <- cv_fit$lambda.min
# elastic_net_model <- glmnet(Xtrain_genelist, ytrain_binary, family = "multinomial", alpha = alpha_value, lambda = best_lambda)
# table(ytrain)

predictions <- predict(elastic_net_model, newx = Xtest_genelist, s = "lambda.min", type = "class")
table(predictions)
table(ytest_binary)
# Confusion matrix
confusion_matrix <- table(predictions, ytest_binary)
print("Confusion Matrix:")
print(confusion_matrix)
accuracy <- (confusion_matrix[1,1]+ confusion_matrix[2,2])/sum(confusion_matrix)
accuracy
intersect(gene_list, names(labelled_scort)[top_features])
