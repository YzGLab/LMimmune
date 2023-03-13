#
library(data.table)
library(GEOquery)
library(dplyr)
library(Seurat)
library(ggplot2)
library(umap)
library(hypeR)
library(msigdbr)
library(tidyverse)
library(pheatmap)
library(Seurat)
library(clustree)
library(NMF)##
###
phe<-read.delim2("MET500_geneExpression_M.meta.plus.txt")
##
exp<-fread("MET500_geneExpression_M.mx.txt")
metExp<-data.frame(exp)
rownames(metExp)<-metExp[,1]
metExp<-apply(metExp,2,as.numeric)
colnames(metExp)<-colnames(exp)[-1]
metExp<-metExp[,-1]
rownames(metExp)<-data.frame(exp)[,1]
metExp<-data.frame(metExp)
###
metExp[1:5,1:5]
####                                                                                                                                                                            
#ID2ENSEMBL<-read.delim2("probeMap_gencode.v23.annotation.gene.probemap")
#iOrd <- is.element(rownames(metExp),ID2ENSEMBL$id)
#metExp<- metExp[iOrd,]
#iOrd <- match(rownames(metExp),ID2ENSEMBL$id)
#rownames(metExp) <- ID2ENSEMBL$gene[iOrd]

probe2symbol<-ID2ENSEMBL
metExp <- metExp %>% 
  rownames_to_column(var="id") %>% 
  inner_join(probe2symbol,by="id") %>% 
  dplyr::select(-id) %>% 
  dplyr::select(gene,everything()) %>% 
  #mutate(rowMean =rowMeans(.[grep("ES", names(.))])) %>% 
  filter(gene != "NA") %>% 
  #arrange(desc(rowMean)) %>% 
  distinct(gene,.keep_all = T) %>% 
  #dplyr::select(-rowMean) %>% 
  column_to_rownames(var = "gene")

metExp[1:5,1:5]
colnames(metExp)
metExp<-metExp[,c(1:868)]
save(metExp,file="metExp500.Rdata")#############

##met500 
load("metExp500.Rdata")
###
metphe<-read_excel("MET500_geneExpression_M.meta.plus.xlsx")

##
metExp2<-metExp[,colnames(metExp) %in% metphe$Sample_id]


newexp_with_genesymobol[["metExp188"]]<-metExp2
metphe$tumor<-ifelse(metphe$tumor=="Colon","Colorectal",metphe$tumor)
allphe[["metExp188"]]<-data.frame(metphe)
