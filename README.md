# Scene
__Inferring subtype-specific ceRNA modules in breast cancer__
## Background
The heterogeneity of breast cancer poses a fundamental challenge to clinical management, manifesting both in molecular subtype diversity and functionally distinct gene modules. Long non-coding RNAs (lncRNAs) acting as competing endogenous RNAs (ceRNAs) or microRNA (miRNA) sponges are emerging biomrakers of breast cancer. Accurate identification of subtype-specific ceRNA modules could contribute to precision medicine in breast cancer.<br><br>In this work, we propose a novel framework Scene (Subtype-specific ceRNA modules) to infer lncRNA-related breast cancer subtype-specific ceRNA modules from heterogeneous data, including gene expression data and priori information of miRNA targets.Firstly, the input data of Scene includes gene expression data and priori information of miRNA targets. For each breast cancer subtype, Scene identifies the ceRNA modules specific to it. In total, 22 subtype-specific ceRNA modules are inferred. To understand the biological significance of the identified subtype-specific ceRNA modules, Scene performs six types of downstream analysis, including heterogeneity analysis, distribution analysis, enrichment analysis, survival analysis, multi-class classification analysis, and immune infiltration analysis.<br><br>A schematic illustration of Scene is as follows.<br><br>![](https://github.com/YangHL24/Scene/blob/main/Scene.png)<br><br>For five breast cancer subtypes, most of ceRNA modules tend to be unique. Across 22 breast cancer-specific ceRNA modules, 20 ceRNA modules are significantly enriched in various biological processes or pathways, indicating that these modules play distinct biological functions in different subtypes. Survival analysis further indicates that all identified ceRNA modules serve as potential prognostic biomarkers capable of discriminating between high- and low-risk breast cancer groups. Moreover, classification analysis shows that all inferred ceRNA modules function as potential diagnostic biomarkers for distinguishing breast cancer subtypes. Finally, immune infiltration analysis reveals that all identified ceRNA modules show significant correlation with one or more immune cell types, suggesting their potential involvement in immune regulation within the tumor microenvironment.Collectively, our results demonstrate that Scene effectively elucidates both the specificity and commonality of ceRNA modules in breast cancer, advancing our understanding of its molecular pathogenesis.
## Description of file in R
* __Scence_framework.R__: Code for inferring and analyzing breast cancer subtype-specific ceRNA modules using the Scene framework.
## The usage of Scene
Paste all files into a single folder (set the folder as the directory of R environment). The users can simply run the scripts as follows.
```
source("R/Scence_case_study.R")
```
## Quick example to use Scene
To infer cancer subtype-specific ceRNA modules, users should prepare matching miRNA (optional), mRNA, ncRNA expression data as well as miRNA-target interactions. Paste the dataset and the Scene_case_study.R file into a single folder (set the folder to the directory of the R environment), users can use the following script to infer breast cancer-specific ceRNA modules. For convenience, we have prepared bulk breast cancer data with labeled sample subtypes for users.
```
# load bulk breast cancer dataset
load("lncRNA_se.RData")
load("mRNA_se.RData")
load("miRNA_se.RData")
load("miRTar_se.RData")

# load packages
library(reactome.db)
library(SummarizedExperiment)
library(igraph)
library(miRSM)
library(foreach)
library(doParallel)

# Inference of the ceRNA modules for all samples
set.seed(1234)
modulegenes_WGCNA_BRCA <- module_WGCNA(lncRNA_se,mRNA_se)
miRSM_WGCNA_SDC <- miRSM(miRNA_se, lncRNA_se, mRNA_se, miRTar_se,
                         modulegenes_WGCNA_BRCA,
                         num_shared_miRNAs = 3, pvalue.cutoff = 0.05,
                         method = "SDC", MC.cutoff = 0.8,
                         SMC.cutoff = 0.1, RV_method = "RV")

# Inference of breast cancer subtype-specific ceRNA modules
num.cores <- 3
cl <- makeCluster(num.cores)
registerDoParallel(cl)
nsamples <- 5
modulegenes_all_WGCNA<-modulegenes_WGCNA_BRCA
results_modulegenes_exceptk_WGCNA <- foreach(i = seq(nsamples), .packages = "miRSM") %dopar% {
  module_WGCNA(lncRNA_se[-index[[i]]$index, ], mRNA_se[-index[[i]]$index, ])
}
miRSM_SDC_exceptk_WGCNA<- foreach(i = seq(nsamples), .packages = "miRSM") %dopar% { miRSM(miRNA_se[-index[[i]]$index, ],
                                                                                          lncRNA_se[-index[[i]]$index,],
                                                                                          mRNA_se[-index[[i]]$index, ],
                                                                                          miRTar_se, results_modulegenes_exceptk_WGCNA[[i]],
                                                                                          method = "SDC",
                                                                                          SMC.cutoff = 0.1)}

stopCluster(cl)
stopImplicitCluster()

miRSM_SDC_all_WGCNA <- miRSM_WGCNA_SDC
Modulegenes_all_SDC_WGCNA <- miRSM_SDC_all_WGCNA[[2]]
Modulegenes_exceptk_SDC_WGCNA <- lapply(seq(nsamples), function(i) miRSM_SDC_exceptk_WGCNA[[i]][[2]])
Modules_SS_SDC_WGCNA <- miRSM_SS(Modulegenes_all_SDC_WGCNA, Modulegenes_exceptk_SDC_WGCNA)
Modules_SS_SDC_WGCNA
```
