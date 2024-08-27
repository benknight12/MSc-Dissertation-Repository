## 2. Logistic regression between features (top_genes) and cluster membership - run with different top gene lists
library(nnet)

labelled_scort <- read.csv("Supervised_Sets/scort_supervised.csv",row.names = 1)
head(labelled_scort[,1:20])
load("12. gene list comp/Final_Gene_Lists.RData")
conclusive_scort_genes <- read.csv("12. gene list comp/Top_Gene_List.csv",row.names = 1)
final_clusters <- read.csv("12. gene list comp/Final_scort_Cluster_List.csv", row.names = 1)
conclusive_scort_genes

labelled_scort$binary_canonical <-as.factor(final_clusters=="Canonical")

feature_names <- conclusive_scort_genes$Canonical
feature_names <- feature_names[!is.na(feature_names)]
formula <- paste0("binary_canonical ~ ", paste(feature_names, collapse = " + "))

# Log reg for binary canonical vs not

table(labelled_scort$binary_canonical)  

log_reg_model <- glm(formula, data=labelled_scort, family = binomial)
summary(log_reg_model)
exp(cbind(coef(log_reg_model),confint(log_reg_model)))

library(DescTools)
pseudo_r2b <- PseudoR2(log_reg_model)
print(pseudo_r2b)

 #### scrap this probs