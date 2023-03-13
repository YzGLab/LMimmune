##
library(EPIC)
library(immunedeconv)
library(dplyr)
library(ggplot2)
library(tidyr)
#library(tibble)
library(MCPcounter)
library(xCell)

#####
load("phe_met.Rdata")
load("metExp500.Rdata")
load("met500_GTEx.expSet.Rdata")
########################################
genes <- data.table::fread("MCPcounter-master/Signatures/genes.txt",data.table = F)
probesets <- data.table::fread("MCPcounter-master/Signatures/probesets.txt",data.table = F,header = F)
### 成功！！！
colnames(metExp)
exprMat<-metExp
results_MCPcounter<- MCPcounter.estimate(exprMat,
                                         featuresType= "HUGO_symbols",
                                         probesets=probesets,
                                         genes=genes)
results_MCPcounter<-data.frame(results_MCPcounter)
colnames(results_MCPcounter)
results_MCPcounter<-t(results_MCPcounter)
#xcell
res_xcell <- xCellAnalysis(exprMat)
res_xcell<-t(res_xcell)
##quantiseq
res_quantiseq = deconvolute(exprMat, "quantiseq", tumor = TRUE)
res_quantiseq<-t(res_quantiseq)
colnames(res_quantiseq)<-res_quantiseq[1,]
res_quantiseq<-res_quantiseq[-1,]

#############
library(estimate)
estimateScores<-list()
for ( i in c("met500")){
  if (!dir.exists("split_sampl_Single")){
    dir.create("./split_sampl_Single")
  }
  print(i)
  mydata<-metExp
  #sample_split_list[[i]]<-mydata
  write.table(mydata,paste0("./split_sampl_Single/", i, ".txt"),sep="\t",quote=F,col.names=T)
}
for ( i in  c("met500")){
  print(i)
  in.file<-paste0("./split_sampl_Single/", i, ".txt")
  out.file<-paste0("./split_sampl_Single/", i,"outfile.gct")
  out.file2<-paste0("./split_sampl_Single/", i,"out.file_score.gct")
  filterCommonGenes(input.f= in.file , output.f=out.file, id="GeneSymbol")
  estimateScore(out.file, out.file2)
  #estimateScore[1:2,1:4]
  estimateScore<-read.table(out.file2, skip = 2, header = TRUE, sep = "\t")
  estimateScore<-data.frame(t(estimateScore))
  colnames(estimateScore)<- c('StromalScore','ImmuneScore', 'ESTIMATEScore',"TumorPurity")
  estimateScore<-estimateScore[-c(1,2),]
  #estimateScores[[i]]<-estimateScore
  message(paste0(i,"calculation of Immunescore is done!"))
}
##
library(ImmuneSubtypeClassifier)
calls <- callEnsemble(exprMat, geneids='symbol')
immunefeaturs<-list()
immunefeaturs[["MCP"]]<-results_MCPcounter
immunefeaturs[["xcell"]]<-res_xcell
immunefeaturs[["quantiseq"]]<-res_quantiseq
immunefeaturs[["estimate"]]<-estimateScore
immunefeaturs[["ImmuneSubtype"]]<-calls[,-1]
immunefeaturs<-do.call(cbind,immunefeaturs)
colnames(immunefeaturs)

immunefeaturs<-cbind(immunefeaturs,IPSs)
colnames(immunefeaturs)
###
library(readxl)### HCC genesignatures
Hep_score <- read_excel("~/Projects/LiverCancer/1-s2.0-S0092867420316135-mmc3.xlsx")
Hep_score<-t(Hep_score)
Hep_score<-split(Hep_score,rownames(Hep_score))
for( i in names(Hep_score)){
  Hep_score[[i]]<-Hep_score[[i]][Hep_score[[i]]!=" "]
  Hep_score[[i]]<-Hep_score[[i]][!is.na(Hep_score[[i]])]
}
metExp<-expSet[[1]]

##
save(immunefeaturs,file="immunefeaturs.Rdata")

load("immunefeaturs.Rdata")



#install.packages("DESeq")
#BiocManager::install("DESeq2")
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("GSVA")
######
library("ExperimentHub")
library("easierData")
score_signature_genes <- suppressMessages(easierData::get_scores_signature_genes())
#####
#save(score_signature_genes,file="score_signature_genes.Rdata")


#install.packages("~/Projects/MET500Projects/Signature/IMvigor210CoreBiologies_1.0.0.tar.gz", repos = NULL, type = "source")
#eh <- ExperimentHub()
#query(eh, "easierData")
#########################
library("IMvigor210CoreBiologies")
data(human_gene_signatures)
##############
ind_genes <- human_gene_signatures
# variables
ind_genes<-ind_genes[-c(8,9)]
names(ind_genes)[7]<-"Pan F TBRS"
##
#################################
Combined_signatures <- read_excel("Signature/Combined_signatures.xlsx")
Combined_signatures<-split(Combined_signatures$GeneSymbol,Combined_signatures$Name)
#Combined_signatures<-mapply(c, Combined_signatures,ind_genes, SIMPLIFY=FALSE)
ind_genes<-appendList(Combined_signatures,ind_genes)

#ind_genes<-appendList(Pathways,ind_genes)

goi <- names(ind_genes)



pd<-data.frame(colnames(m2))
# calculate gene set scores
for (sig in goi) {
  pd[,sig] <- NA
  genes <- ind_genes[[sig]]
  genes <- genes[genes %in% rownames(m2)]
  tmp <- m2[genes, , drop=FALSE]
  pd[, sig] <- gsScore(tmp)
}
#install.packages("ropenblas")
#ropenblas::rcompiler()

######合并两个list 方法
source("Signature/appendList.R")
appendList <- function (x, val) 
{
  stopifnot(is.list(x), is.list(val))
  xnames <- names(x)
  for (v in names(val)) {
    x[[v]] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]])) 
      appendList(x[[v]], val[[v]])
    else c(x[[v]], val[[v]])
  }
  x
}

colnames(immunefeaturs)
#immunefeaturs<-immunefeaturs[,-c(107:212)]

immunefeaturs<-cbind(immunefeaturs,pd[,-1])
####

########################
load("immunefeaturs.Rdata")
immunefeaturs<-cbind(immunefeaturs,gsva_es)
save(immunefeaturs,file="immunefeaturs.Rdata")
#####
load("~/Projects/MET500Projects/mets-immunecluster-master/data/pathways_for_gsva_mets.Rdata")
#Antigen processing machinery
library(clusterProfiler)
library(GSVA)
metaboslism<-read.gmt("Immunemetabolism_geneset.txt") %>% data.frame()
metaboslism<-metaboslism[metaboslism$gene!="",]
metaboslism<-split(metaboslism$gene,metaboslism$term)

Pathways<-appendList(Pathways,metaboslism)
m2<-metExp
norm.expMat<-as.matrix(m2)
gsva_es <- gsva(norm.expMat,Pathways,method="ssgsea",abs.ranking=F,min.sz = 2,ssgsea.norm=TRUE,kcdf="Poisson",parallel.sz=70)###Array data,ssgsea.norm=TRUE
gsva_es<-t(gsva_es) %>% data.frame()
colnames(immunefeaturs)
immunefeaturs<-immunefeaturs[,-c(233:344)]
immunefeaturs<-cbind(immunefeaturs,gsva_es)
#####