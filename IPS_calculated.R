####################################################
##
##   This R-script can be used to calculate Immunophenoscore (IPS) and generate Immunophenogram from "EXPR.txt" and "IPS_genes.txt"
##   (C) ICBI, Medical University of Innsbruck, Biocenter, Division of Bioinformatics
##   Version 1.0 08.07.2016
##   Needs packages ggplot2,grid,gridExtra
##
####################################################

library(ggplot2)
library(grid)
library(gridExtra)

## calculate Immunophenoscore
ipsmap<- function (x) {
  if (x<=0) {
    ips<-0
  } else {
    if (x>=3) {
      ips<-10
    } else {
      ips<-round(x*10/3, digits=0)
    }
  }
  return(ips)
}

## Assign colors 
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)
mapcolors<-function (x) {
  za<-NULL
  if (x>=3) {
    za=1000
  } else {
    if (x<=-3) {
      za=1
    } else {
      za=round(166.5*x+500.5,digits=0)
    }
  }
  return(my_palette[za])
}
my_palette2 <- colorRampPalette(c("black", "white"))(n = 1000)
mapbw<-function (x) {
  za2<-NULL
  if (x>=2) {
    za2=1000
  } else {
    if (x<=-2) {
      za2=1
    } else {
      za2=round(249.75*x+500.5,digits=0)
    }
  }
  return(my_palette2[za2])
}





#########################################

#
## Read expression data from tab-delimited text file, with official human gene symbols (HGNC) in the first columns
## and expression values (i.e. log2(TPM+1)) for each sample in the other columns
#gene_expression<-read.table("EXPR.txt",row.names=1,header=TRUE, sep="\t", dec = ".",check.names=FALSE)
colnames(metExp)
  exprMat<-metExp

  ###########
  #exprMat<-data.frame(exprMat)
  gene_expression<-exprMat
  sample_names<-names(gene_expression)
  #
  ## Read IPS genes and corresponding weights from tab-delimited text file "IPS_genes.txt"
  # For different 
  IPSG<-read.table("/home/data/gaoyuzhen/Projects/MET500Projects/Immunophenogram-master/IPS_genes.txt",header=TRUE, sep="\t", dec = ".",check.names=FALSE)
  unique_ips_genes<-as.vector(unique(IPSG$NAME))
  #
  IPS<-NULL
  MHC<-NULL
  CP<-NULL
  EC<-NULL
  SC<-NULL
  AZ<-NULL
  
  # Gene names in expression file
  GVEC<-row.names(gene_expression)
  # Genes names in IPS genes file
  VEC<-as.vector(IPSG$GENE)
  # Match IPS genes with genes in expression file
  ind<-which(is.na(match(VEC,GVEC)))
  # List genes missing or differently named
  MISSING_GENES<-VEC[ind]
  dat<-IPSG[ind,]
  if (length(MISSING_GENES)>0) {
    cat("differently named or missing genes: ",MISSING_GENES,"\n")
  }
  for (x in 1:length(ind)) {
    print(IPSG[ind,])
  }
  
  for (i in 1:length(sample_names)) {	
    GE<-gene_expression[[i]]
    mGE<-mean(GE)
    sGE<-sd(GE)
    Z1<-(gene_expression[as.vector(IPSG$GENE),i]-mGE)/sGE
    W1<-IPSG$WEIGHT
    WEIGHT<-NULL
    MIG<-NULL
    k<-1
    for (gen in unique_ips_genes) {
      MIG[k]<- mean(Z1[which (as.vector(IPSG$NAME)==gen)],na.rm=TRUE)
      WEIGHT[k]<- mean(W1[which (as.vector(IPSG$NAME)==gen)])
      k<-k+1
    }
    WG<-MIG*WEIGHT
    MHC[i]<-mean(WG[1:10]);if(is.na(MHC[i])){MHC[i]=0}
    CP[i]<-mean(WG[11:20]);if(is.na(CP[i])){CP[i]=0}
    EC[i]<-mean(WG[21:24]);if(is.na(EC[i])){EC[i]=0}
    SC[i]<-mean(WG[25:26]);if(is.na(SC[i])){SC[i]=0}
    AZ[i]<-sum(MHC[i],CP[i],EC[i],SC[i])
    IPS[i]<-ipsmap(AZ[i])
  }
  
  DF<-data.frame(SAMPLE=sample_names,MHC=MHC,EC=EC,SC=SC,CP=CP,AZ=AZ,IPS=IPS)
  IPSs<-DF
  #write.table(DF,file="IPS.txt",row.names=FALSE, quote=FALSE,sep="\t")
  #barplot(DF[[7]],DF[[1]],col=c(rep("yellow",7),rep("green",4)))
}
#################################################################
