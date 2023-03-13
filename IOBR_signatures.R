if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
depens<-c('tibble', 'survival', 'survminer', 'sva', 'limma', "DESeq2","devtools",
          'limSolve', 'GSVA', 'e1071', 'preprocessCore', 'ggplot2', "biomaRt",
          'ggpubr', "devtools", "tidyHeatmap", "caret", "glmnet", "ppcor","timeROC","pracma")
for(i in 1:length(depens)){
  depen<-depens[i]
  if (!requireNamespace(depen, quietly = TRUE)) BiocManager::install(depen)
}
if (!requireNamespace("remotes", quietly = TRUE)) install("remotes")
if (!requireNamespace("EPIC", quietly = TRUE))
  remotes::install_github("GfellerLab/EPIC", build_vignettes=TRUE)
if (!requireNamespace("MCPcounter", quietly = TRUE))
  remotes::install_github("ebecht/MCPcounter",ref="master", subdir="Source")
if (!requireNamespace("estimate", quietly = TRUE)){
  rforge <- "http://r-forge.r-project.org"
  install.packages("estimate", repos=rforge, dependencies=TRUE)
}  
if (!requireNamespace("IOBR", quietly = TRUE))
  remotes::install_github("IOBR/IOBR",ref="master")

library(IOBR)
library(EPIC)
library(estimate) 
library(MCPcounter)
library(tidyverse)
library(tidyHeatmap)
library(maftools)
library(ggpubr)
library(ggplot2)
library(survival)
library(GSVA)
###
#load("phe_met.Rdata")
load("metExp500.Rdata")
load("met500_GTEx.expSet.Rdata")
colnames(sig_tme)
load("/home/data/gaoyuzhen/Projects/MET500Projects/Signature/sig.TME.Rdata")

eset_stad<-metExp
sig_tme<-calculate_sig_score(pdata           = NULL,
                             eset            = eset_stad,
                             signature       = signature_tme,
                             method          = "pca",
                             mini_gene_count = 5)

sig_tme<-calculate_sig_score(pdata           = NULL,
                             eset            = eset_stad,
                             signature       = signature_collection,
                             method          = "ssgsea",
                             mini_gene_count = 5)
sig_meta<-calculate_sig_score(pdata           = NULL,
                              eset            = eset_stad,
                              signature       = signature_metabolism,
                              method          = "pca",
                              mini_gene_count = 2)


sig_res<-calculate_sig_score(pdata           = NULL,
                             eset            = eset_stad,
                             signature       = signature_collection,
                             method          = "integration",
                             mini_gene_count = 2)
##pathway
sig_hallmark<-calculate_sig_score(pdata           = NULL,
                                  eset            = eset_stad,
                                  signature       = hallmark,
                                  method          = "ssgsea",
                                  mini_gene_count = 2)

save(sig_tme,file ="/home/data/gaoyuzhen/Projects/MET500Projects/Signature/sig.TME.Rdata")