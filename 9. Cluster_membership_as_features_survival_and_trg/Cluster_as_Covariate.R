## 4. a - Fit linear model/anova for survival as outcome and each cluster membership as covariate
##    b - also chi squared/log reg between cluster membership and TRG - assess performance with R2 measure/pvalues
library(car) 
library(survminer)
library(sjPlot)
library(nnet)
## final scort clusters
labelled_scort <- read.csv("Supervised_Sets/scort_supervised.csv", row.names = 1)
final_scort_clusters <- read.csv("12. gene list comp/Final_scort_Cluster_List.csv",row.names = 1)

labelled_scort <- cbind(final_scort_clusters, labelled_scort)
names(labelled_scort)[1] <- "cluster"
binary_canonical <- final_scort_clusters=="Canonical"
table(final_scort_clusters)
labelled_scort <- cbind(binary_canonical, labelled_scort)
names(labelled_scort)[1] <- "binary_canonical"
labelled_scort$OS.Status[is.na(labelled_scort$OS.Status)] <- "LIVING"
labelled_scort$OS.Status <- as.factor(labelled_scort$OS.Status)

survival_cluster_scortmodel <- lm(OS ~ binary_canonical, data = labelled_scort)
cbind(coef(survival_cluster_scortmodel),confint(survival_cluster_scortmodel))
summary(survival_cluster_scortmodel)



########## 

# anova_result <- aov(OS ~ NFkB.Cluster, data = labelled_scort)
# summary(anova_result)
# ## Check assumptions:
# plot(anova_result, which = 1)
# qqnorm(residuals(anova_result))
# qqline(residuals(anova_result))
# leveneTest(OS ~ NFkB.Cluster, data = labelled_scort)
# TukeyHSD(anova_result)
# ### Check these some time - i think they indicate no large issues - maybe qq not very normal

########## Anova acc not great bc censored data - cox ph far better really:

library(survival)
cox_model <- coxph(Surv(OS, OS.Status=="DECEASED") ~ binary_canonical, data = labelled_scort)
summary(cox_model)
exp(cbind(coef(cox_model),confint(cox_model)))
# Plot survival curves
ggforest(cox_model, data = labelled_scort)
cox.zph(cox_model)

labelled_87211 <- read.csv("Supervised_Sets/87211.csv", row.names = 1)
clusters <- read.csv("12. gene list comp/87211_reduced_clusters.csv",row.names = 1)
labelled_87211$binary_canonical <- clusters=="Canonical"
clusters
cox_model_87211 <- coxph(Surv(OS, OS.Status==1) ~ binary_canonical, data = labelled_87211)
summary(cox_model_87211)
exp(cbind(coef(cox_model_87211),confint(cox_model_87211)))
## reproduce adams plot
surv_plot <- ggsurvplot(survfit(Surv(OS, OS.Status=="DECEASED") ~ cluster, data = labelled_scort), data = labelled_scort, pval = TRUE,xlim=c(0,120), ylim = c(0.5,1), xlab = "Time (months)", conf.int = TRUE,break.time.by = 12, risk.table = TRUE, ylab = "Overall Survival Probability", title = "Kaplan-Meier Curve for Overall Survival Different Clusters", legend.title = "Clusters", legend.labs = c("Atypical", "Canonical", "NonCanonical"))
surv_plot
head(labelled_scort$OS.Status)


### Do as pairs/atypical + noncanonical vs canonical




####### b) chi squared/log reg between cluster membership and TRG - assess performance with R2 measure/pvalues
clusters <- read.csv("12. gene list comp/Final_scort_Cluster_List.csv", row.names = 1)
labelled_scort$binary_canonical <- as.factor(clusters=="Canonical")
labelled_scort$Binary_response <- factor(labelled_scort$Binary_response, levels = c("poor response", "response"), labels =c(0,1))
trg_outcome_scort_model <- glm(Binary_response ~ binary_canonical, family = binomial, data = labelled_scort)

exp(cbind(coef(trg_outcome_scort_model),confint(trg_outcome_scort_model)))
tab_model(trg_outcome_scort_model)

multinom_model <- multinom(as.factor(TRG) ~ NFkB.Cluster, data = labelled_scort)
summary(multinom_model)
tab_model(multinom_model)
exp(coef(multinom_model))

contingency_table <- table(as.factor(labelled_scort$TRG), as.factor(labelled_scort$NFkB.Cluster))
chisq.test(contingency_table) # v strong evidence of relationship between clusters and therapy response



##################### Not sure what my conclusions are
 








