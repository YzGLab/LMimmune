####heatmap for subtype and its score. 
library(readxl)
load("metExp500.Rdata")
load("new_gsva2.Rdata")
load("met500_GTEx.expSet.Rdata")
load(file.path(Score.path,"phe_total.Rdata"))
var_cat<-read.csv(file.path(Score.path,"phe_total_cat.csv"))
## first step  different
tmp<-metExp
###
comparedsignatures<-read_excel("/home/data/gaoyuzhen/Projects/MET500Projects/Signature/Combined_signatures.xlsx")
##
table(comparedsignatures$Name)
tmp1<-tmp[rownames(tmp) %in% comparedsignatures$GeneSymbol,]
####
table(phe_total1$group_cluster)
rt<-phe_total1
colnames(tmp1)<-gsub("-",".",colnames(tmp1))
exprSet<-tmp1[,colnames(tmp1) %in% rt$Sample_id]
deglist<-list()
for (i in unique(phe_total1$group_cluster)){
  print(i)
  #rt<-phe_list[[i]]
  #exprSet<-gsva_list[[i]]
  rt$groups<-ifelse(rt$group_cluster==i,"yes",'no')
  # high risk vs others
  subt <- data.frame(condition =rt$groups,
                     row.names = colnames(exprSet))
  twoclasslimma(subtype  = subt, # subtype information (must contain a column named 'condition')
                featmat  = exprSet, # expression file (fill detect data scale automatically)
                treatVar = "yes", # name of treatment group
                ctrlVar  = "no", # name of control group
                prefix   = i, # prefix of the DE file
                overwt   = TRUE, # whether over write files that already existed
                sort.p   = TRUE, # if sorting the results by adjusted p value
                verbose  = FALSE, # if showing verbose result
                res.path = Score.path) # path for result
  tmp1 <- read.table(file.path(Score.path,paste0(i,"_limma_test_result.yes_vs_no.txt")),
                     sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
  diffSig <-  tmp1[which(tmp1$log2fc >2 & tmp1$padj < 0.05),]
  diffSig$gene<-rownames(diffSig)
  deglist[[i]]<-diffSig
}####

diffSigs<-do.call(rbind,deglist)
#diffSigs<-sapply(deglist,c)
unique(diffSigs$gene)

#write.csv(diffSig_liver,file=file.path(Score.path,"GTEX_diffSig_liver.csv"))
