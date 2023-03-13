library(GSVA)
library(limma)
library(clusterProfiler)

###################################
load("~/Projects/PanCancerdata/panGTEX_mRNA_exprSet.Rdata")
load("metExp500.Rdata")
gmt_immune<-read.gmt("c7.all.v7.5.1.symbols.gmt")
gmt_immune<-split(gmt_immune$gene,gmt_immune$term)
###############
table(mRNA_exprSet$type)
####################
#mRNA_exprSet$type<-ifelse(mRNA_exprSet$type=="Liver","liver","nonLiver")
rownames(mRNA_exprSet)<-mRNA_exprSet[,1]
phe_mRNA_exprSet<-mRNA_exprSet[,c(1,2)]
##
table(phe_mRNA_exprSet$type)
####################################
mRNA_exprSet<-mRNA_exprSet[-c(1,2)]
colnames(mRNA_exprSet)
mRNA_exprSet<-t(mRNA_exprSet)
###
###############################
my_data_frame<-data.frame(mRNA_exprSet)
chunk <- 1000
n <- ncol(my_data_frame)
r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
r<-data.frame(r,colnames(my_data_frame))
#d <- split(my_data_frame,r)
sample_split_list<-list()
table(r)
for ( i in unique(r$r)){
  if (!dir.exists("split_sampl_Single_nromal")){
    dir.create("./split_sampl_Single_nromal")
  }
  print(i)
  mydata<-my_data_frame[,colnames(my_data_frame) %in% r[r$r==i,]$colnames.my_data_frame.]
  sample_split_list[[i]]<-mydata
  write.table(mydata,paste0("./split_sampl_Single_nromal/", i, ".txt"),sep="\t",quote=F,col.names=T)
}
###
normal_gsva<-list()
for(i in c(1:8)){
  mydata<-sample_split_list[[i]]
  norm.expMat<-as.matrix(mydata)
  gsva_es <- gsva(norm.expMat,gmt_immune,method="gsva",abs.ranking=F,min.sz = 5,ssgsea.norm=TRUE,kcdf="Poisson",parallel.sz=70)###Array data,ssgsea.norm=TRUE
  #gsva_es <- gsva(norm.expMat,gmt.list,method="ssgsea",abs.ranking=F,kcdf="Gaussian",parallel.sz=10)#RNA-Seq 
  gsva_es[1:5,1:5]
  normal_gsva[[i]]<-gsva_es
}
save(normal_gsva,file="GTExnormal_gsva.Rdata")
##########
load("GTExnormal_gsva.Rdata")
normal_gsvas<-do.call(cbind,normal_gsva)
normal_gsvas[1:5,1:5]
#############################################################################
norm.expMat<-as.matrix(metExp)
gsva_es <- gsva(norm.expMat,gmt_immune,method="gsva",abs.ranking=F,min.sz = 5,ssgsea.norm=TRUE,kcdf="Poisson",parallel.sz=70)###Array data,ssgsea.norm=TRUE
#gsva_es <- gsva(norm.expMat,gmt.list,method="ssgsea",abs.ranking=F,kcdf="Gaussian",parallel.sz=10)#RNA-Seq 
gsva_es[1:5,1:5]
metExp_gsva<-gsva_es
save(metExp_gsva,file="metExp_gsva.Rdata")
##################
load("GTExnormal_gsva.Rdata")
normal_gsvas<-do.call(cbind,normal_gsva)
normal_gsvas[1:5,1:5]
########
gsva_list<-list()
gsva_list[["normalGTEx"]]<-normal_gsvas
gsva_list[["met500"]]<-metExp_gsva
###
phe_met<-phe_met[,-c(2,3)]
phe_list<-list()
phe_list[["normalGTEx"]]<-phe_mRNA_exprSet
phe_list[["met500"]]<-phe_met
save(gsva_list,phe_list,file="new_gsva2.Rdata")
###
myphe<-phe_list[[1]]
table(myphe$biopsy_tissue)
myphe<-myphe[!myphe$tissue %in% c("bone_marrow","brain","liver", "lung","skin"),]
phe_list[[1]]<-myphe

#myphe<-phe_list[[2]]
#table(myphe$type)


######################
############################### JITC singautures for clusters.
load("~/Projects/MET500Projects/mets-immunecluster-master/data/pathways_for_gsva_mets.Rdata")
load("~/Projects/PanCancerdata/panGTEX_mRNA_exprSet.Rdata")
load("metExp500.Rdata")
###
load("new_gsva.Rdata")
norm.expMat<-as.matrix(metExp)
gsva_es <- gsva(norm.expMat,Pathways,method="gsva",abs.ranking=F,min.sz = 5,ssgsea.norm=TRUE,kcdf="Poisson",parallel.sz=70)###Array data,ssgsea.norm=TRUE
gsva_list[["met500_pathways"]]<-gsva_es 
##
#for nomral
normal_gsva_pathways<-list()
for(i in c(1:8)){
  print(i)
  tmp1 <- read.table(paste0("./split_sampl_Single_nromal/", i, ".txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
  norm.expMat<-as.matrix(tmp1)
  gsva_es <- gsva(norm.expMat,Pathways,method="gsva",abs.ranking=F,min.sz = 5,ssgsea.norm=TRUE,kcdf="Poisson",parallel.sz=70)###Array data,ssgsea.norm=TRUE
  #gsva_es <- gsva(norm.expMat,gmt.list,method="ssgsea",abs.ranking=F,kcdf="Gaussian",parallel.sz=10)#RNA-Seq 
  gsva_es[1:5,1:5]
  normal_gsva_pathways[[i]]<-gsva_es
}
normal_gsva_pathways<-do.call(cbind,normal_gsva_pathways)
normal_gsva_pathways[1:5,1:5]
gsva_list[["normalpanGTEX_pathways"]]<-normal_gsva_pathways
save(gsva_list,phe_list,file="new_gsva.Rdata")
