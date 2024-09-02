pheno_87211 <- read.csv("Processed_Data_Sets/pheno87211.csv", row.names = 1)
minipheno <- pheno_87211[,c("depth.of.invasion.after.rct.ch1","survival.time..month..ch1","disease.free.time..month..ch1","death.due.to.tumor.ch1","cancer.recurrance.after.surgery.ch1", "status")]
names(minipheno) <- c("TRG","OS","DFS", "OS.Status", "DFS.Status", "dummy")
table(!is.na(pheno_87211$depth.of.invasion.after.rct.ch1))
cc <- pheno_87211$depth.of.invasion.after.rct.ch1[complete.cases(minipheno)]
missing <- pheno_87211$depth.of.invasion.after.rct.ch1[!complete.cases(minipheno)]
table(minipheno$TRG[!complete.cases(minipheno)])
table(minipheno$OS.Status[!complete.cases(minipheno)])
table(minipheno$DFS.Status[!complete.cases(minipheno)])
table(is.na(minipheno$TRG))
table(is.na(minipheno$DFS))
library(naniar)
library(mice)
table(is.na(pheno_87211$status))
naniar::vis_miss(minipheno)

mcar_testsurv <- mcar_test(minipheno[,c("OS","DFS","OS.Status","DFS.Status")])
mcar_testtrg <- mcar_test(minipheno[,c(1,5)])
print(mcar_testsurv)
print(mcar_testtrg)

