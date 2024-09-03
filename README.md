# Repository for code used for MSc Dissertation
## Stage 1: Exploring different clustering in SCORT, GSE109057 and GSE145037
1a. Alternative - Alternate version of Chambers clusters with my own k-means on NF-kB for SCORT.
1a-1c. Cluster derivation for Each SCORT, GSE109057 and GSE145037 with 4 NF-kB sub-units. Run DEA.
1d-1e. Imposing SCORT cluster structure on GEO sets, running DEA.
2a-2c. Cluster derviation using k-means for whole dataset in each of the three datasets.
3a-3c. Cluster derivation using Hierarchical (PVClust) on the whole dataset for the three datasets.

## Stage 2: Compiling gene lists
4. Comparing lists from 1a-1c, 2a-2c and 3a-3c to see crossover between Differentially Expressed Genes between datasets.
5. Comparing genelists for clusters predicted from SCORT structure, such they are a more similar lists.

## Stage 3: Assessing Predictive Capacity of the Clusters and Genes (Mostly not reported)
6. Initial prediction of Survival and Therapy-Response from the clusters we derived (preliminary version)
7. Measure ability to prediuct cluster membership just from the genes in our signature.
8. Survival analysis for up and down-regulated categories for each genes in SCORT.
8a. 8 but applied to GSE87211. (Was added later)
9. Cluster membership for predicting Survival and Therapy-Response (reproducing Mr Chambers for SCORT, with our clusters, and novel work for GSE87211)
10. Attempt at a Neural Net to see overally predictive capacity of the genes we identified - ultimately failed, v poor accuracy.

## Stage 4: GSE87211 work
11. All analyses condcted in Stages 1 and 2 but for GSE87211.
12. Compiling gene lists - using clusters from assimilated methods in stage 1.

## Stage 5: Final stage of compiling for report
13. Compile final cluster versions using predictions from the SCORT structure.
14. Derive gene lists from the final clusters in 13. Then compile to final 20 genes.
15. Gene list derivation using original Adam's clusters for the sake of comparison

## Additional Folders
CB_Normalized_Data: Final datasets following harmonization
Processed_Data_sets: Earlier drafts of data sets, including phenotype data.
Supervised_Sets: Labelled datasets where scaled sets are combined with phenotype data.
Cluster_Survival_Plots_Scort_and_GSE87211: create the survival plots based on cluster membership.
Intial_Explorations: Initial work looking at the data, some of the pre-processing
Useful_functions: Functions used repeatedly throughout the project we wanted to store.
Some Plots for report: A few of the plots we want to include in the report
PrimeView.db: Annotation package for gse109057, gse145037

0.Cross_Batch_Normalization.R : Harmonization Script
Missiningess_test.R : Scrtipt to test if MCAR (if CCA is okay)

Notable Emissions:
- Script converting CEL to CSV for SCORT and GEO sets.








