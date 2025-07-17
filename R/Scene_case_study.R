## Input data
load("lncRNA_se.RData")
load("mRNA_se.RData")
load("miRNA_se.RData")
load("miRTar_se.RData")
load("BRCA_subtype_index.RData")
load('survival_result.RData')
load("lncREXP_final.RData")
load("mREXP_final.RData")
load('subtype_col')

# Install and load the required packages
library(org.Hs.eg.db)
library(reactome.db)
library(SummarizedExperiment)
library(GSEABase)
# Details of the miRSM R package can be found at http://bioconductor.org/packages/miRSM/ and https://github.com/zhangjunpeng411/miRSM
library(miRSM)
library(igraph)
# Details of the miRspongeR R package can be found at http://bioconductor.org/packages/miRspongeR/ and https://github.com/zhangjunpeng411/miRspongeR.
library(miRspongeR)
library(ggplot2)

############################## Inference of the ceRNA modules for all samples ###############################
set.seed(1234)
# lncRNA_se: SummarizedExperiment object of lncRNA: rows are sample, columns are lncRNA
# mRNA_se: SummarizedExperiment object of mRNA: rows are sample, columns are mRNA
# miRNA_se: SummarizedExperiment object of miRNA: rows are sample, columns are miRNA
# miRTar_se: SummarizedExperiment object of putative miRNA-target interactions
# Output: modulegenes_WGCNA_BRCA is a list of gene co-expression modules
modulegenes_WGCNA_BRCA <- module_WGCNA(lncRNA_se,mRNA_se)
# Output:miRSM_WGCNA_SDC is a list of ceRNA modules
miRSM_WGCNA_SDC <- miRSM(miRNA_se, lncRNA_se, mRNA_se, miRTar_se,
                         modulegenes_WGCNA_BRCA,
                         num_shared_miRNAs = 3, pvalue.cutoff = 0.05,
                         method = "SDC", MC.cutoff = 0.8,
                         SMC.cutoff = 0.1, RV_method = "RV")

################################ Inference of subtype-specific ceRNA modules for breast cancer ################################
library(foreach)
library(doParallel)
## Set the number of cores for parallel computing
num.cores <- 3
cl <- makeCluster(num.cores)
registerDoParallel(cl)
# nsamples: Number of breast cancer subtypes
nsamples <- 5
modulegenes_all_WGCNA<-modulegenes_WGCNA_BRCA

## Infer the gene co-expression module after knocking out the sample set
# lncRNA_se: SummarizedExperiment object of lncRNA: rows are sample, columns are lncRNA
# mRNA_se: SummarizedExperiment object of mRNA: rows are sample, columns are mRNA
# miRNA_se: SummarizedExperiment object of miRNA: rows are sample, columns are miRNA
# miRTar_se: SummarizedExperiment object of putative miRNA-target interactions
# index: A list of sample indices for each breast cancer subtype in the gene expression matrix
# Output: results_modulegenes_exceptk_WGCNA results_modulegenes_exceptk_WGCNA is a list of gene co-expression modules identified after removing subtype k
results_modulegenes_exceptk_WGCNA <- foreach(i = seq(nsamples), .packages = "miRSM") %dopar% {
  module_WGCNA(lncRNA_se[-index[[i]]$index, ], mRNA_se[-index[[i]]$index, ])
}
## Infer the ceRNA module after knocking out the sample set
# Output: miRSM_SDC_exceptk_WGCNA is a list of ceRNA modules identified after removing subtype k
miRSM_SDC_exceptk_WGCNA<- foreach(i = seq(nsamples), .packages = "miRSM") %dopar% { miRSM(miRNA_se[-index[[i]]$index, ],
                                                                                          lncRNA_se[-index[[i]]$index,], mRNA_se[-index[[i]]$index, ],
                                                                                          miRTar_se, results_modulegenes_exceptk_WGCNA[[i]],
                                                                                          method = "SDC",
                                                                                          SMC.cutoff = 0.1)}

stopCluster(cl)
stopImplicitCluster()

miRSM_SDC_all_WGCNA <- miRSM_WGCNA_SDC
Modulegenes_all_SDC_WGCNA <- miRSM_SDC_all_WGCNA[[2]]
Modulegenes_exceptk_SDC_WGCNA <- lapply(seq(nsamples), function(i) miRSM_SDC_exceptk_WGCNA[[i]][[2]])
# Output: Modules_SS_SDC_WGCNA is the list of ceRNA modules of subtype k
Modules_SS_SDC_WGCNA <- miRSM_SS(Modulegenes_all_SDC_WGCNA, Modulegenes_exceptk_SDC_WGCNA)
Modules_SS_SDC_WGCNA


####Visualization of the number of breast cancer-specific ceRNA modules
#Modules_SS_SDC_WGCNA: List of ceRNA modules identified for each breast cancer subtype
##Check whether the input length matches
subtype_names = c("Normal", "LumA", "Her2+", "LumB", "Basal")
if (length(subtype_names) != length(Modules_SS_SDC_WGCNA)) {
  stop("subtype_names length must match modules_list length")
}

ceRNA_modules <- data.frame(
  Subtype = subtype_names,
  "Number of ceRNA modules" = sapply(Modules_SS_SDC_WGCNA, length),
  row.names = NULL,
  check.names = FALSE
)
except_ceRNA_Modulegenes_plot <- ggplot(data = ceRNA_modules, mapping = aes(x = Subtype, y =`Number of ceRNA modules` , fill = Subtype, group = factor(1))) +
  ggtitle("BRCA Subtype-specific CeRNA Modules") +
  scale_fill_manual(values = c("#80FFFF", "#BFFFFF", "#FFD5FF", "#FFBFFF", "#FF80FF"))+
  theme_bw() +
  geom_bar(stat="identity",width=0.4)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +  # 向下扩展y轴范围
  theme(text = element_text(family = "serif"),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 21),  # 居中并增大标题字体大小
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13))+  # 调整图例文本大小为12
  geom_text(aes(label = `Number of ceRNA modules`), vjust = -0.5, color = "black", size = 7,family = "serif")
except_ceRNA_Modulegenes_plot

###################################################### Downstream analysis ###############################################
#### Heterogeneity analysis
# Define the subtype name and the corresponding sample name
subtype_mapping <- list("Normal" = "Sample 1",
                        "LumA" = "Sample 2",
                        "Her2+" = "Sample 3",
                        "LumB" = "Sample 4",
                        "Basal" = "Sample 5")
# The initialization list stores all the results
specific_ceRNA_module_list <- list()
for (subtype in names(subtype_mapping))
{
  ceRNA_data <- Modules_SS_SDC_WGCNA[[subtype_mapping[[subtype]]]]
  # Process each module and combine ceRNA and mRNA
  merged_data <- sapply(ceRNA_data, function(x) paste(c(x$ceRNA, x$mRNA), collapse = ""))
  specific_ceRNA_module_list[[subtype]] <- merged_data
}
library(UpSetR)
except_ceRM_upset<-upset(fromList(specific_ceRNA_module_list),
                         sets=names(specific_ceRNA_module_list),
                         point.size=3.5,
                         line.size=1.5,
                         mainbar.y.label="Number of shared modules",
                         main.bar.color = 'black',
                         matrix.color="black",
                         sets.x.label="Set size",
                         sets.bar.color=c('red', 'orange','green','blue','purple'),
                         order.by = "freq",
                         text.scale=c(3.2,3,2,2.5,2.5,3)
)
except_ceRM_upset

### Calculate the similarity of subtypes based on the ceRNA module
# Normal_ceRNA: the list of ceRNA modules of Normal subtype
# LumA_ceRNA: the list of ceRNA modules of LumA subtype
# Her2_ceRNA: the list of ceRNA modules of Her2 subtype
# LumB_ceRNA: the list of ceRNA modules of LumB subtype
# Basal_ceRNA: the list of ceRNA modules of Basal subtype
Normal_ceRNA<-Modules_SS_SDC_WGCNA$`Sample 1`
LumA_ceRNA<-Modules_SS_SDC_WGCNA$`Sample 2`
Her2_ceRNA<-Modules_SS_SDC_WGCNA$`Sample 3`
LumB_ceRNA<-Modules_SS_SDC_WGCNA$`Sample 4`
Basal_ceRNA<-Modules_SS_SDC_WGCNA$`Sample 5`

## Calculate the similarity between Normal subtype and LumA, Her2, LumB, and Basal subtypes respectively
sim_Normal_LumA<-module_group_sim(Normal_ceRNA, LumA_ceRNA, sim.method = "Simpson")
sim_Normal_Her2<-module_group_sim(Normal_ceRNA, Her2_ceRNA, sim.method = "Simpson")
sim_Normal_LumB<-module_group_sim(Normal_ceRNA, LumB_ceRNA, sim.method = "Simpson")
sim_Normal_Basal<-module_group_sim(Normal_ceRNA, Basal_ceRNA, sim.method = "Simpson")

## Calculate the similarity between LumA subtype and Her2, LumB, and Basal subtypes respectively
sim_LumA_Her2<-module_group_sim(LumA_ceRNA, Her2_ceRNA,sim.method = "Simpson")
sim_LumA_LumB<-module_group_sim(LumA_ceRNA, LumB_ceRNA,sim.method = "Simpson")
sim_LumA_Basal<-module_group_sim(LumA_ceRNA, Basal_ceRNA,sim.method = "Simpson")

## Calculate the similarity between Her2 subtype and LumB and Basal subtypes respectively
sim_Her2_LumB<-module_group_sim(Her2_ceRNA, LumB_ceRNA,sim.method = "Simpson")
sim_Her2_Basal<-module_group_sim(Her2_ceRNA, Basal_ceRNA,sim.method = "Simpson")

## Calculate the similarity between LumB subtype and Basal subtype separately
sim_LumB_Basal<-module_group_sim(LumB_ceRNA, Basal_ceRNA,sim.method = "Simpson")

#### Enrichment analysis
## Functional enrichment analysis of breast cancer-specific ceRNA module
library(openxlsx)
miRSM_Normal_FEA <- module_FA(Normal_ceRNA, Analysis.type = "FEA")
miRSM_LumA_FEA <- module_FA(LumA_ceRNA, Analysis.type = "FEA")
miRSM_Her_FEA <- module_FA(Her2_ceRNA, Analysis.type = "FEA")
miRSM_LumB_FEA <- module_FA(LumB_ceRNA, Analysis.type = "FEA")
miRSM_Basal_FEA <- module_FA(Basal_ceRNA, Analysis.type = "FEA")
# Normal
Nor_FEA_BP<-miRSM_Normal_FEA[[1]]
for (i in 1:length(Nor_FEA_BP)) {
  Nor_FEA_BP_data <- data.frame(Nor_FEA_BP[[i]])

  file_name <- paste0("Normal_", i, "FEA_BP_result.xlsx")
  write.xlsx(Nor_FEA_BP_data, file = file_name, rowNames = FALSE)
}
Nor_FEA_KEGG<-miRSM_Normal_FEA[[2]]
for (i in 1:length(Nor_FEA_KEGG)) {
  Nor_FEA_KEGG_data <- data.frame(Nor_FEA_KEGG[[i]])

  file_name <- paste0("Normal_", i, "FEA_KEGG_result.xlsx")
  write.xlsx(Nor_FEA_KEGG_data, file = file_name, rowNames = FALSE)
}
Nor_FEA_Reactome<-miRSM_Normal_FEA[[3]]
for (i in 1:length(Nor_FEA_Reactome)) {
  Nor_FEA_Reactome_data <- data.frame(Nor_FEA_Reactome[[i]])

  file_name <- paste0("Normal_", i, "FEA_Reactome_result.xlsx")
  write.xlsx(Nor_FEA_Reactome_data, file = file_name, rowNames = FALSE)
}
# Her2
Her_FEA_BP<-miRSM_Her_FEA[[1]]
for (i in 1:length(Her_FEA_BP)) {
  Her_FEA_BP_data <- data.frame(Her_FEA_BP[[i]])

  file_name <- paste0("Her_", i, "FEA_BP_result.xlsx")
  write.xlsx(Her_FEA_BP_data, file = file_name, rowNames = FALSE)
}
Her_FEA_KEGG<-miRSM_Her_FEA[[2]]
for (i in 1:length(Her_FEA_KEGG)) {
  Her_FEA_KEGG_data <- data.frame(Her_FEA_KEGG[[i]])

  file_name <- paste0("Her_", i, "FEA_KEGG_result.xlsx")
  write.xlsx(Her_FEA_KEGG_data, file = file_name, rowNames = FALSE)
}
Her_FEA_Reactome<-miRSM_Her_FEA[[3]]
for (i in 1:length(Her_FEA_Reactome)) {
  Her_FEA_Reactome_data <- data.frame(Her_FEA_Reactome[[i]])

  file_name <- paste0("Her_", i, "FEA_Reactome_result.xlsx")
  write.xlsx(Her_FEA_Reactome_data, file = file_name, rowNames = FALSE)
}
# LumB
LB_FEA_BP<-miRSM_LumB_FEA[[1]]
for (i in 1:length(LB_FEA_BP)) {
  LumB_FEA_BP_data <- data.frame(LB_FEA_BP[[i]])

  file_name <- paste0("LumB_", i, "FEA_BP_result.xlsx")
  write.xlsx(LumB_FEA_BP_data, file = file_name, rowNames = FALSE)
}
LB_FEA_KEGG<-miRSM_LumB_FEA[[2]]
for (i in 1:length(LB_FEA_KEGG)) {
  LumB_FEA_KEGG_data <- data.frame(LB_FEA_KEGG[[i]])

  file_name <- paste0("LumB_", i, "FEA_KEGG_result.xlsx")
  write.xlsx(LumB_FEA_KEGG_data, file = file_name, rowNames = FALSE)
}
LB_FEA_Reactome<-miRSM_LumB_FEA[[3]]
for (i in 1:length(LB_FEA_Reactome)) {
  LumB_FEA_Reactome_data <- data.frame(LB_FEA_Reactome[[i]])

  file_name <- paste0("LumB_", i, "FEA_Reactome_result.xlsx")
  write.xlsx(LumB_FEA_Reactome_data, file = file_name, rowNames = FALSE)
}
#Basal
Ba_FEA_BP<-miRSM_Basal_FEA[[1]]
for (i in 1:length(Ba_FEA_BP)) {
  Basal_FEA_BP_data <- data.frame(Ba_FEA_BP[[i]])

  file_name <- paste0("Basal_", i, "FEA_BP_result.xlsx")
  write.xlsx(Basal_FEA_BP_data, file = file_name, rowNames = FALSE)
}
Ba_FEA_KEGG<-miRSM_Basal_FEA[[2]]
for (i in 1:length(Ba_FEA_KEGG)) {
  Basal_FEA_KEGG_data <- data.frame(Ba_FEA_KEGG[[i]])

  file_name <- paste0("Basal_", i, "FEA_KEGG_result.xlsx")
  write.xlsx(Basal_FEA_KEGG_data, file = file_name, rowNames = FALSE)
}
Ba_FEA_Reactome<-miRSM_Basal_FEA[[3]]
for (i in 1:length(Ba_FEA_Reactome)) {
  Basal_FEA_Reactome_data <- data.frame(Ba_FEA_Reactome[[i]])

  file_name <- paste0("Basal_", i, "FEA_Reactome_result.xlsx")
  write.xlsx(Basal_FEA_Reactome_data, file = file_name, rowNames = FALSE)
}

## Disease enrichment analysis of breast cancer-specific ceRNA module
miRSM_Normal_DEA <- module_FA(Normal_ceRNA, Analysis.type = "DEA")
miRSM_LumA_DEA <- module_FA(LumA_ceRNA, Analysis.type = "DEA")
miRSM_Her_DEA <- module_FA(Her2_ceRNA, Analysis.type = "DEA")
miRSM_LumB_DEA <- module_FA(LumB_ceRNA, Analysis.type = "DEA")
miRSM_Basal_DEA <- module_FA(Basal_ceRNA, Analysis.type = "DEA")
###Normal
Nor_DEA_HDO<-miRSM_Normal_DEA[[1]]
for (i in 1:length(Nor_DEA_HDO)) {
  Nor_DEA_HDO_data <- data.frame(Nor_DEA_HDO[[i]])

  file_name <- paste0("Normal_", i, "DEA_HDO_result.xlsx")
  write.xlsx(Nor_DEA_HDO_data, file = file_name, rowNames = FALSE)
}
Nor_DEA_DisGeNET<-miRSM_Normal_DEA[[2]]
for (i in 1:length(Nor_DEA_DisGeNET)) {
  Nor_DEA_DisGeNET_data <- data.frame(Nor_DEA_DisGeNET[[i]])

  file_name <- paste0("Normal_", i, "DEA_DisGeNET_result.xlsx")
  write.xlsx(Nor_DEA_DisGeNET_data, file = file_name, rowNames = FALSE)
}
Nor_DEA_NCG<-miRSM_Normal_DEA[[3]]
for (i in 1:length(Nor_DEA_NCG)) {
  Nor_DEA_NCG_data <- data.frame(Nor_DEA_NCG[[i]])

  file_name <- paste0("Normal_", i, "DEA_NCG_result.xlsx")
  write.xlsx(Nor_DEA_NCG_data, file = file_name, rowNames = FALSE)
}
# LumA
LA_DEA_HDO<-miRSM_LumA_DEA[[1]]
for (i in 1:length(LA_DEA_HDO)) {
  LA_DEA_HDO_data <- data.frame(LA_DEA_HDO[[i]])

  file_name <- paste0("LumA_", i, "DEA_HDO_result.xlsx")
  write.xlsx(LA_DEA_HDO_data, file = file_name, rowNames = FALSE)
}
LA_DEA_DisGeNET<-miRSM_LumA_DEA[[2]]
for (i in 1:length(LA_DEA_DisGeNET)) {
  LA_DEA_DisGeNET_data <- data.frame(LA_DEA_DisGeNET[[i]])

  file_name <- paste0("LumA_", i, "DEA_DisGeNET_result.xlsx")
  write.xlsx(LA_DEA_DisGeNET_data, file = file_name, rowNames = FALSE)
}
LA_DEA_NCG<-miRSM_LumA_DEA[[3]]
for (i in 1:length(LA_DEA_NCG)) {
  LA_DEA_NCG_data <- data.frame(LA_DEA_NCG[[i]])

  file_name <- paste0("LumA_", i, "DEA_NCG_result.xlsx")
  write.xlsx(LA_DEA_NCG_data, file = file_name, rowNames = FALSE)
}
# Her2
Her_DEA_HDO<-miRSM_Her_DEA[[1]]
for (i in 1:length(Her_DEA_HDO)) {
  Her_DEA_HDO_data <- data.frame(Her_DEA_HDO[[i]])

  file_name <- paste0("Her_", i, "DEA_HDO_result.xlsx")
  write.xlsx(Her_DEA_HDO_data, file = file_name, rowNames = FALSE)
}
Her_DEA_DisGeNET<-miRSM_Her_DEA[[2]]
for (i in 1:length(Her_DEA_DisGeNET)) {
  Her_DEA_DisGeNET_data <- data.frame(Her_DEA_DisGeNET[[i]])

  file_name <- paste0("Her_", i, "DEA_DisGeNET_result.xlsx")
  write.xlsx(Her_DEA_DisGeNET_data, file = file_name, rowNames = FALSE)
}
Her_DEA_NCG<-miRSM_Her_DEA[[3]]
for (i in 1:length(Her_DEA_NCG)) {
  Her_DEA_NCG_data <- data.frame(Her_DEA_NCG[[i]])

  file_name <- paste0("Her_", i, "DEA_NCG_result.xlsx")
  write.xlsx(Her_DEA_NCG_data, file = file_name, rowNames = FALSE)
}
# LumB
LB_DEA_HDO<-miRSM_LumB_DEA[[1]]
for (i in 1:length(LB_DEA_HDO)) {
  LumB_DEA_HDO_data <- data.frame(LB_DEA_HDO[[i]])

  file_name <- paste0("LumB_", i, "DEA_HDO_result.xlsx")
  write.xlsx(LumB_DEA_HDO_data, file = file_name, rowNames = FALSE)
}
LB_DEA_DisGeNET<-miRSM_LumB_DEA[[2]]
for (i in 1:length(LB_DEA_DisGeNET)) {
  LumB_DEA_DisGeNET_data <- data.frame(LB_DEA_DisGeNET[[i]])

  file_name <- paste0("LumB_", i, "DEA_DisGeNET_result.xlsx")
  write.xlsx(LumB_DEA_DisGeNET_data, file = file_name, rowNames = FALSE)
}
LB_DEA_NCG<-miRSM_LumB_DEA[[3]]
for (i in 1:length(LB_DEA_NCG)) {
  LumB_DEA_NCG_data <- data.frame(LB_DEA_NCG[[i]])

  file_name <- paste0("LumB_", i, "DEA_NCG_result.xlsx")
  write.xlsx(LumB_DEA_NCG_data, file = file_name, rowNames = FALSE)
}
# Basal
Ba_DEA_HDO<-miRSM_Basal_DEA[[1]]
for (i in 1:length(Ba_DEA_HDO)) {
  Basal_DEA_HDO_data <- data.frame(Ba_DEA_HDO[[i]])

  file_name <- paste0("Basal_", i, "DEA_HDO_result.xlsx")
  write.xlsx(Basal_DEA_HDO_data, file = file_name, rowNames = FALSE)
}
Ba_DEA_DisGeNET<-miRSM_Basal_DEA[[2]]
for (i in 1:length(Ba_DEA_DisGeNET)) {
  Basal_DEA_DisGeNET_data <- data.frame(Ba_DEA_DisGeNET[[i]])

  file_name <- paste0("Basal_", i, "DEA_DisGeNET_result.xlsx")
  write.xlsx(Basal_DEA_DisGeNET_data, file = file_name, rowNames = FALSE)
}
Ba_DEA_NCG<-miRSM_Basal_DEA[[3]]
for (i in 1:length(Ba_DEA_NCG)) {
  Basal_DEA_NCG_data <- data.frame(Ba_DEA_NCG[[i]])

  file_name <- paste0("Basal_", i, "DEA_NCG_result.xlsx")
  write.xlsx(Basal_DEA_NCG_data, file = file_name, rowNames = FALSE)
}

#### Distribution analysis
## Extract the shared mirnas of each breast cancer-specific ceRNA module
# Nor_share: Shared miRNAs of each ceRNA module of Normal subtype
# LumA_share: Shared miRNAs of each ceRNA module of LumA subtype
# Her_share: Shared miRNAs of each ceRNA module of Her2 subtype
# LumB_share: Shared miRNAs of each ceRNA module of LumB subtype
# Basal_share: Shared miRNAs of each ceRNA module of Basal subtype
Nor_share<-share_miRs(miRNA_se,miRTar_se,Normal_ceRNA)
LumA_share<-share_miRs(miRNA_se,miRTar_se,LumA_ceRNA)
Her_share<-share_miRs(miRNA_se,miRTar_se,Her2_ceRNA)
LumB_share<-share_miRs(miRNA_se,miRTar_se,LumB_ceRNA)
Basal_share<-share_miRs(miRNA_se,miRTar_se,Basal_ceRNA)

Normal_share_upset<-upset(fromList(Nor_share),
                          sets=names(Nor_share),
                          point.size=3.5,
                          line.size=1.5,
                          mainbar.y.label="Number of shared miRNAs",
                          main.bar.color = 'black',
                          matrix.color="black",
                          sets.x.label="Set size",
                          sets.bar.color=c('red','blue'),
                          order.by = "freq",
                          text.scale=c(3,3,2,2.5,2.5,3)
)
Normal_share_upset
LumA_share_upset<-upset(fromList(LumA_share),
                        sets=names(LumA_share),
                        point.size=3.5,
                        line.size=1.5,
                        mainbar.y.label="Number of shared miRNAs",
                        main.bar.color = 'black',
                        matrix.color="black",
                        sets.x.label="Set size",
                        sets.bar.color=c('pink','red', 'orange','yellow','green','cyan','blue','purple','magenta'),
                        order.by = "freq",
                        text.scale=c(3.2,3,2,2.5,2.5,3)
)
LumA_share_upset
Her2_share_upset<-upset(fromList(Her_share),
                        sets=names(Her_share),
                        point.size=3.5,
                        line.size=1.5,
                        mainbar.y.label="Number of shared miRNAs",
                        main.bar.color = 'black',
                        matrix.color="black",
                        sets.x.label="Set size",
                        sets.bar.color=c('red','green','blue'),
                        order.by = "freq",
                        text.scale=c(3.2,3,2,2.5,2.5,3)
)
Her2_share_upset
LumB_share_upset<-upset(fromList(LumB_share),
                        sets=names(LumB_share),
                        point.size=3.5,
                        line.size=1.5,
                        mainbar.y.label="Number of shared miRNAs",
                        main.bar.color = 'black',
                        matrix.color="black",
                        sets.x.label="Set size",
                        sets.bar.color=c('red', 'orange','green','blue'),
                        order.by = "freq",
                        text.scale=c(3.2,3,2,2.5,2.5,3)
)
LumB_share_upset
Basal_share_upset<-upset(fromList(Basal_share),
                         sets=names(Basal_share),
                         point.size=3.5,
                         line.size=1.5,
                         mainbar.y.label="Number of shared miRNAs",
                         main.bar.color = 'black',
                         matrix.color="black",
                         sets.x.label="Set size",
                         sets.bar.color=c('red', 'orange','green','blue'),
                         order.by = "freq",
                         text.scale=c(3.2,3,2,2.5,2.5,3)
                         )
Basal_share_upset

####Survival analysis
lnc<-assay(lncRNA_se)
mR<-assay(mRNA_se)
miRNA<-assay(miRNA_se)
exp<-cbind(miRNA,mR,lnc)
row<-BRCA_sur$sample
row_index<-match(row,rownames(exp))
expdata<-exp[row_index,]

Normal_list <- list()
for (i in 1:length(Normal_ceRNA)) {
  Normal_list[[i]] <- c(Normal_ceRNA[[i]]$ceRNA, Normal_ceRNA[[i]]$mRNA)
}

LumA_list <- list()
for (i in 1:length(LumA_ceRNA)) {
  LumA_list[[i]] <- c(LumA_ceRNA[[i]]$ceRNA, LumA_ceRNA[[i]]$mRNA)
}

Her2_list <- list()
for (i in 1:length(Her2_ceRNA)) {
  Her2_list[[i]] <- c(Her2_ceRNA[[i]]$ceRNA, Her2_ceRNA[[i]]$mRNA)
}

LumB_list <- list()
for (i in 1:length(LumB_ceRNA)) {
  LumB_list[[i]] <- c(LumB_ceRNA[[i]]$ceRNA, LumB_ceRNA[[i]]$mRNA)
}

Basal_list <- list()
for (i in 1:length(Basal_ceRNA)) {
  Basal_list[[i]] <- c(Basal_ceRNA[[i]]$ceRNA, Basal_ceRNA[[i]]$mRNA)
}
# Normal_list: List of ceRNA modules of Normal subtype
# LumA_list: List of ceRNA modules of LumA subtype
# Her_sur: List of ceRNA modules of Her2 subtype
# LumB_sur: List of ceRNA modules of LumB subtype
# Basal_sur: List of ceRNA modules of Basal subtype
# expdata： Gene expression matrix including mRNA and lncRNA
Nor_sur<-moduleSurvival(Normal_list , expdata ,BRCA_sur,devidePercentage=.5, plot = TRUE)
LumA_sur<-moduleSurvival(LumA_list , expdata ,BRCA_sur,devidePercentage=.5, plot = TRUE)
Her_sur<-moduleSurvival(Her2_list , expdata ,BRCA_sur,devidePercentage=.5, plot = TRUE)
LumB_sur<-moduleSurvival(LumB_list , expdata ,BRCA_sur,devidePercentage=.5, plot = TRUE)
Basal_sur<-moduleSurvival(Basal_list , expdata ,BRCA_sur,devidePercentage=.5, plot = TRUE)


# sur_senlin：Data frame of survival analysis results
tabletext <- cbind(
  c("Module", sur_senlin$Module),
  c("Chi-square",round(sur_senlin$`Chi-square`,2)),
  c("HR (95% CI)", paste0(round(sur_senlin$`HR (95% CI)`, 2), " (", round(sur_senlin$HRlow95, 2), " - ", round(sur_senlin$HRup95, 2), ")")),
  c("p-value", sur_senlin$`p-value`)
)
library(forestplot)
forest_plot <- forestplot(
  labeltext = tabletext,
  #title = paste0(gene_name, " – Overall Survival"),
  #title_gp = grid::gpar(fontsize = 14, fontface = "bold", just = "left"),
  mean = c(NA, sur_senlin$`HR (95% CI)`),
  lower = c(NA, sur_senlin$HRlow95),
  upper = c(NA, sur_senlin$HRup95),
  hrzl_lines = list(`1` = grid::gpar(lwd = 2, col = "black"),
                    `2` = grid::gpar(lwd = 2, col = "black"),
                    `24` = grid::gpar(lwd = 2, col = "black")),
  is.summary = rep(FALSE, nrow(sur_senlin) + 1),
  zero = 1,
  boxsize = 0.2,
  lineheight = unit(0.7, "cm"),
  xlog = FALSE,
  col = fpColors(box = "#4DBBD5", lines = "black", zero = "gray50"),
  lwd.ci = 2.5,
  ci.vertices = TRUE,
  ci.vertices.height = 0.02,
  clip = c(0, 30),
  xticks = seq(0, 30, by = 5),
  graph.pos = 4,
  graphwidth = unit(6, "cm"),
  txt_gp = forestplot::fpTxtGp(ticks = grid::gpar(cex = 0.9)),
)
forest_plot

#### Classification analysis
library(e1071)
library(plyr)
library(mldr)
library(utiml)

#Organize the data frame for subtype grouping of the samples
BRCA_Exp<-rbind(lncREXP_final,mREXP_final)
BRCA_Exp<-t(BRCA_Exp)

subtype_type <- cbind(rownames(subtype_col), subtype_col)
colnames(subtype_type)[2]="type"
colnames(subtype_type)[1]="sample"
rownames(subtype_type) <- NULL
##Sort out the list of modules
names(Normal_list)[1] <- "Normal-ceRM1"
names(Normal_list)[2] <- "Normal-ceRM2"

names(LumA_list)[1]<-"LumA-ceRM1"
names(LumA_list)[2]<-"LumA-ceRM2"
names(LumA_list)[3]<-"LumA-ceRM3"
names(LumA_list)[4]<-"LumA-ceRM4"
names(LumA_list)[5]<-"LumA-ceRM5"
names(LumA_list)[6]<-"LumA-ceRM6"
names(LumA_list)[7]<-"LumA-ceRM7"
names(LumA_list)[8]<-"LumA-ceRM8"
names(LumA_list)[9]<-"LumA-ceRM9"

names(Her2_list)[1]<-"Her2-ceRM1"
names(Her2_list)[2]<-"Her2-ceRM2"
names(Her2_list)[3]<-"Her2-ceRM3"

names(LumB_list)[1]<-"LumB-ceRM1"
names(LumB_list)[2]<-"LumB-ceRM2"
names(LumB_list)[3]<-"LumB-ceRM3"
names(LumB_list)[4]<-"LumB-ceRM4"

names(Basal_list)[1]<-"Basal-ceRM1"
names(Basal_list)[2]<-"Basal-ceRM2"
names(Basal_list)[3]<-"Basal-ceRM3"
names(Basal_list)[4]<-"Basal-ceRM4"

## Function for multi-class classification analysis of ceRNA modules
# ceRExp: Gene expression matrix of ceRNA
# mRExp: Gene expression matrix of mRNA
# BRCA_subtype: A data frame with annotations of breast cancer subtypes for the samples
# Modulelist: A collection of breast cancer subtype-specific ceRNA module lists
# Output: module_classify is the matrix of multi-class classification analysis results of each ceRNA module
module.classify <- function(ceRExp, mRExp, BRCA_subtype, Modulelist, method = "br", base.algorith = "SVM", cv.folds = 10,
                            cv.sampling = "stratified", cv.seed = 12345) {

  module_ceRExp <- lapply(seq_along(Modulelist), function(i) ceRExp[, which(colnames(ceRExp) %in% Modulelist[[i]])])
  module_mRExp <- lapply(seq_along(Modulelist), function(i) mRExp[, which(colnames(mRExp) %in% Modulelist[[i]])])
  Basal <- as.numeric(BRCA_subtype[, 2] == "Basal")
  Her2 <- as.numeric(BRCA_subtype[, 2] == "Her2")
  LumA <- as.numeric(BRCA_subtype[, 2] == "LumA")
  LumB <- as.numeric(BRCA_subtype[, 2] == "LumB")
  Normal <- as.numeric(BRCA_subtype[, 2] == "Normal")
  module_classify <- list()

  for (i in seq_along(Modulelist)){

    temp <- as.data.frame(cbind(module_ceRExp[[i]], module_mRExp[[i]], Basal, Her2, LumA, LumB, Normal))
    Indices <- ncol(temp)
    temp_mldr <- mldr_from_dataframe(temp, labelIndices = c(Indices-4, Indices-3, Indices-2, Indices-1, Indices), name = "TEMPMLDR")
    temp_res <- cv(temp_mldr, method = method, base.algorith = base.algorith, cv.folds = cv.folds,
                   cv.sampling = cv.sampling, cv.seed = cv.seed)
    module_classify[[i]] <- temp_res

  }

  return(module_classify)
}

lncR<-t(lncREXP_final)
mR<-t(mREXP_final)
Modulelist<-c(Normal_list,LumA_list,Her2_list,LumB_list,Basal_list)

moduel_Classification_baseline<-do.call(cbind,module.classify(lncR, mR, subtype_type, Modulelist, method = "baseline"))
moduel_Classification<-do.call(cbind,module.classify(lncR, mR, subtype_type, Modulelist, method = "br"))
colnames(moduel_Classification)=names(Modulelist)
colnames(moduel_Classification_baseline)=names(Modulelist)

#### Immune infiltration analysis
library(GSVA)
# A list containing the expression matrix of genes in each Normal ceRNA module
Normal_exp_list <- list()
for (miRSM_name in names(Normal_ceRNA)) {
  lnc <- Normal_ceRNA[[miRSM_name]]$ceRNA
  lnc_in <- match(lnc, rownames(lncREXP_final))
  lncRExp <- lncREXP_final[lnc_in,]
  que <- any(is.na(lncRExp))

  mR <- Normal_ceRNA[[miRSM_name]]$mRNA
  mR_in <- match(mR, rownames(mREXP_final))
  mRExp <- mREXP_final[mR_in,]

  exp <- rbind(lncRExp, mRExp)

  Normal_exp_list[[miRSM_name]] <- exp
}
# Expression matrix of genes in the Normal ceRNA module
Normal_exp <- Normal_exp_list[[1]]
for (i in 2:length(Normal_exp_list)) {
  Normal_exp <- rbind(Normal_exp, Normal_exp_list[[i]])
}
# Normal_es: Enrichment score matrix of Normal ceRNA modules
Normal_para<-gsvaParam(Normal_exp,LumA_list)
Normal_es<-gsva(Normal_para,verbose=TRUE)
Normal_es_data<-as.data.frame(Normal_es)
rownames(Normal_es_data)=names(Normal_list)

# A list containing the expression matrix of genes in each LumA ceRNA module
LumA_exp_list <- list()
for (miRSM_name in names(LumA_ceRNA)) {
  Lu_lnc1 <- LumA_ceRNA[[miRSM_name]]$ceRNA
  Lu_lnc1_in <- match(Lu_lnc1, rownames(lncREXP_final))
  Lu_lncRExp1 <- lncREXP_final[Lu_lnc1_in,]
  que <- any(is.na(Lu_lncRExp1))

  Lu_mR1 <- LumA_ceRNA[[miRSM_name]]$mRNA
  Lu_mR1_in <- match(Lu_mR1, rownames(mREXP_final))
  Lu_mRExp1 <- mREXP_final[Lu_mR1_in,]

  Lu_exp1 <- rbind(Lu_lncRExp1, Lu_mRExp1)

  LumA_exp_list[[miRSM_name]] <- Lu_exp1
}
# Expression matrix of genes in the LumA ceRNA module
LumA_exp <- LumA_exp_list[[1]]
for (i in 2:length(LumA_exp_list)) {
  LumA_exp <- rbind(LumA_exp, LumA_exp_list[[i]])
}
# LumA_es: Enrichment score matrix of LumA ceRNA modules
LumA_para<-gsvaParam(LumA_exp,LumA_list)
LumA_es<-gsva(LumA_para,verbose=TRUE)
LumA_es_data<-as.data.frame(LumA_es)
rownames(LumA_es_data)=names(LumA_list)

# A list containing the expression matrix of genes in each Her2 ceRNA module
Her_exp_list <- list()
for (miRSM_name in names(Her2_ceRNA)) {
  Her_lnc1 <- Her2_ceRNA[[miRSM_name]]$ceRNA
  Her_lnc1_in <- match(Her_lnc1, rownames(lncREXP_final))
  Her_lncRExp1 <- lncREXP_final[Her_lnc1_in,]

  Her_mR1 <- Her2_ceRNA[[miRSM_name]]$mRNA
  Her_mR1_in <- match(Her_mR1, rownames(mREXP_final))
  Her_mRExp1 <- mREXP_final[Her_mR1_in,]

  Her_exp <- rbind(Her_lncRExp1, Her_mRExp1)

  Her_exp_list[[miRSM_name]] <- Her_exp
}
# Expression matrix of genes in the Her2 ceRNA module
Her_exp <- Her_exp_list[[1]]
for (i in 2:length(Her_exp_list)) {
  Her_exp <- rbind(Her_exp, Her_exp_list[[i]])
}
# Her2_es: Enrichment score matrix of Her2 ceRNA modules
Her_para<-gsvaParam(Her_exp,Her2_list)
Her2_es<-gsva(Her_para,verbose=TRUE)
Her2_es_data<-as.data.frame(Her2_es)
rownames(Her2_es_data)=names(Her2_list)

# A list containing the expression matrix of genes in each LumB ceRNA module
LumB_exp_list <- list()
for (miRSM_name in names(LumB_ceRNA)) {
  lnc <- LumB_ceRNA[[miRSM_name]]$ceRNA
  lnc_in <- match(lnc, rownames(lncREXP_final))
  lncRExp <- lncREXP_final[lnc_in,]

  mR <- LumB_ceRNA[[miRSM_name]]$mRNA
  mR_in <- match(mR, rownames(mREXP_final))
  mRExp <- mREXP_final[mR_in,]

  exp <- rbind(lncRExp, mRExp)

  LumB_exp_list[[miRSM_name]] <- exp
}
# Expression matrix of genes in the LumB ceRNA module
LumB_exp <- LumB_exp_list[[1]]
for (i in 2:length(LumB_exp_list)) {
  LumB_exp <- rbind(LumB_exp, LumB_exp_list[[i]])
}
# LumB_es: Enrichment score matrix of LumB ceRNA modules
LumB_para<-gsvaParam(LumB_exp,LumB_list)
LumB_es<-gsva(LumB_para,verbose=TRUE)
LumB_es_data<-as.data.frame(LumB_es)
rownames(LumB_es_data)=names(LumB_list)

# A list containing the expression matrix of genes in each Basal ceRNA module
Basal_exp_list <- list()
for (miRSM_name in names(Basal_ceRNA)) {
  lnc <- Basal_ceRNA[[miRSM_name]]$ceRNA
  lnc_in <- match(lnc, rownames(lncREXP_final))
  lncRExp <- lncREXP_final[lnc_in,]

  mR <- Basal_ceRNA[[miRSM_name]]$mRNA
  mR_in <- match(mR, rownames(mREXP_final))
  mRExp <- mREXP_final[mR_in,]

  exp <- rbind(lncRExp, mRExp)

  Basal_exp_list[[miRSM_name]] <- exp
}
# Expression matrix of genes in the Basal ceRNA module
Basal_exp <- Basal_exp_list[[1]]
for (i in 2:length(Basal_exp_list)) {
  Basal_exp <- rbind(Basal_exp, Basal_exp_list[[i]])
}
# Basal_es: Enrichment score matrix of Basal ceRNA modules
Basal_para<-gsvaParam(Basal_exp,Basal_list)
Basal_es<-gsva(Basal_para,verbose=TRUE)
Basal_es_data<-as.data.frame(Basal_es)
rownames(Basal_es_data)=names(Basal_list)

# module_es: Enrichment score matrix of ceRNA modules
module_es<-rbind(Normal_es_data,LumA_es_data,Her2_es_data,LumB_es_data,Basal_es_data)

library(readxl)
imm_data <- read_excel("Immune cell markers.xlsx", col_names = TRUE)
imm_data<-as.data.frame(imm_data)
immune_geneset <- lapply(1:ncol(imm_data), function(col_index) {
  imm_data[, col_index]
})
for (i in 1:length(immune_geneset)) {
  immune_geneset[[i]] <- na.omit(immune_geneset[[i]])
}
names(immune_geneset)<-colnames(imm_data)

EXP<-rbind(mREXP_final,lncREXP_final)
# EXP: Gene expression matrix: including lncRNA and mRNA
# immune_geneset: List of 24 immune cell markers
# exp_geneSet: Enrichment score matrix of 24 immune cells
exp_Param<-gsvaParam(EXP,immune_geneset,assay = NA_character_,annotation = NA_character_,
                     minSize = 1,maxSize = Inf,kcdf = "Gaussian",tau = 1,maxDiff = TRUE,absRanking = FALSE)
exp_geneSet<- gsva(exp_Param, verbose = TRUE)
exp_geneSet_t<-t(exp_geneSet)
module_es_t<-t(module_es)
## Calculate the correlation of two enrichment score matrixs
module_geneset_cor <- cor(module_es_t, exp_geneSet_t, use = "everything", method = "pearson")
library(pheatmap)
pheatmap(module_geneset_cor,cluster_rows = F,cluster_cols = F,
         color=colorRampPalette(c("#63baf8","white","#ff7c9c"))(100),
         show_colnames = T,
         show_rownames = T,
         fontsize = 13
         )
