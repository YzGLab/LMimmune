################################
library(GSVA)
library(limma)
library(clusterProfiler)
library(sva)
library(dplyr)
library(ggplot2)
library(clustree)
library(data.table)
library(gplots)
library(ComplexHeatmap)
library(ggpubr)
library(scales)
library(RColorBrewer)
library(ggvenn)
library(fpc)
library(data.table)
library(pheatmap)
library(ComplexHeatmap)
# step one 临床资料整理
#######################
##for met500
#
######挑选一些必要的组织
metasite<-c("adrenal", "bladder",'bone_marrow','brain','breast','cervix', 'colon','liver','lung',"pancreas",'prostate','skin','thyroid') 
normaltissue<-c("Adrenal Gland","Bladder","Bone Marrow","Brain","Breast","Uterus","Colon","Liver","Lung","Pancreas","Prostate","Skin","Thyroid")
###
metasite<-c('bone_marrow','brain','liver','lung','skin') 
normaltissue<-c("Bone Marrow","Brain","Liver","Lung","Skin")

#####

load("metExp500.Rdata")
phe_met<-read.delim2("MET500_geneExpression_M.meta.plus.txt")
phe_met$tissue<-ifelse(phe_met$tissue=="","other",phe_met$tissue)
phe_met<-phe_met[phe_met$tissue!="other",]
phe_met<-phe_met[!phe_met$tissue %in% metasite,]###########删除 liver 原位转移的内容
phe_met<-phe_met[!(phe_met$tissue=="liver" & phe_met$biopsy_tissue %in% "liver"),]###########删除 liver 原位转移的内容
phe_met<-phe_met[!(phe_met$tissue=="bone_marrow" & phe_met$biopsy_tissue=="bone_marrow"),]###########删除 liver 原位转移的内容
phe_met<-phe_met[!(phe_met$tissue=="brain" & phe_met$biopsy_tissue=="brain"),]###########删除 liver 原位转移的内容
phe_met<-phe_met[!(phe_met$tissue=="lung" & phe_met$biopsy_tissue=="lung"),]###########删除 liver 原位转移的内容
phe_met<-phe_met[!(phe_met$tissue=="skin" & phe_met$biopsy_tissue=="skin"),]###########删除 liver 原位转移的内容

table(phe_met$tissue)
table(phe_met$biopsy_tissue)
table(phe_met$biopsy_tissue,phe_met$tissue)
#################################################
phe_met<-phe_met[phe_met$biopsy_tissue %in% metasite,]
phe_met$biopsy_tissue_bio<-ifelse(phe_met$biopsy_tissue=="liver","liver","nonliver")
rownames(phe_met)<-phe_met$Sample_id
table(phe_met$tissue)
table(phe_met$biopsy_tissue)
####pie
pie(prop.table(table(phe_met$biopsy_tissue)))
pie(prop.table(table(phe_met$tissue)))
#######
load("~/Projects/PanCancerdata/panGTEX_mRNA_exprSet.Rdata")
table(mRNA_exprSet$type)
####################
#mRNA_exprSet$type<-ifelse(mRNA_exprSet$type=="Liver","liver","nonLiver")
rownames(mRNA_exprSet)<-mRNA_exprSet[,1]
phe_mRNA_exprSet<-mRNA_exprSet[,c(1,2)]
phe_mRNA_exprSet<-phe_mRNA_exprSet[phe_mRNA_exprSet$type %in% normaltissue,]
table(phe_mRNA_exprSet$type)####
phe_mRNA_exprSet$biopsy_tissue_bio<-ifelse(phe_mRNA_exprSet$type=="Liver","liver","nonliver")
rownames(phe_mRNA_exprSet)<-gsub("-",".",phe_mRNA_exprSet[,1])
###pie 图构造
pie(prop.table(table(phe_mRNA_exprSet$type)))
####
#####################
barplot(prop.table(table(phe_met$tissue)))
barplot(prop.table(table(phe_met$biopsy_tissue)))



##############如何展示 ，ggplot2 , top plot 2 top metastasis

results<-"Results"
Piepath    <-"Results/Piepath"
if (!file.exists(results)) { dir.create(results) }
if (!file.exists(Piepath)) { dir.create(Piepath) }

####################
tissuesite<-unique(phe_met$tissue)
subsite<-tissuesite[2]

for (subsite in tissuesite) {
  print(subsite)
  sub_phe_met<-phe_met[phe_met$tissue %in% subsite,]
  pdf(file.path(Piepath,paste0(subsite,".pdf")),width = 5,height = 5)
  pie(prop.table(table(sub_phe_met$biopsy_tissue)),col=mycol)
  legend("topleft",subsite,cex =1)
  dev.off()
}
######
ggplot(phe_met, aes(x=biopsy_tissue,fill=tissue)) +
  #geom_bar(stat="identity",position="fill") + 
  geom_bar(position="fill") +
  scale_y_continuous(labels = percent) +
  #theme_light()+
  ylab("%Metastasis")+
  #facet_wrap(~tissue,scales="free_x")+
  theme_classic()+
  scale_fill_manual(values = mycol)+
  theme(axis.text.x = element_text(angle =30, hjust = 1,size = 10), 
        axis.text.y = element_text(angle = 90, size = 12),
        axis.title.y = element_text(angle = 90, size = 15)) +
   geom_text(aes( y=..count../tapply(..count.., ..x.. ,sum)[..x..], 
               label=paste(..count..,paste("(",percent(..count../tapply(..count.., ..x.. ,sum)[..x..]),")",sep=""),sep=" ")),
          stat="count",color="white", position=position_fill(0.5), vjust=0.5,size=3) +
  theme(legend.position = "right")

ggsave(file.path(Piepath,"tissue_site_count_percent.pdf"),
       width = 14.38, height =6.85,
       bg = "transparent", dpi = 300, useDingbats=FALSE)

ggplot(phe_met) +
  geom_bar( aes(x=biopsy_tissue,fill=tissue)) + 
  #geom_bar(position="fill") +
  #scale_y_continuous(labels = percent) +
  #theme_light()+
  ylab("% sample")+
  #facet_wrap(~tissue,scales="free_x")+
  theme_classic()+
  scale_fill_manual(values = mycol)+
  theme(axis.text.x = element_text(angle =30, hjust = 1,size = 10), 
        axis.text.y = element_text(angle = 90, size = 12),
        axis.title.y = element_text(angle = 90, size = 15)) +
  theme(legend.position =  "right")
ggsave(file.path(Piepath,"tissue_site_count.pdf"),
       width = 14.38, height =6.85,
       bg = "transparent", dpi = 300, useDingbats=FALSE)

###
phe_GTEx<-phe_mRNA_exprSet
expSet<-list()
expSet[["metExp500"]]<-metExp
expSet[["panGTEX"]]<-data.frame(t(mRNA_exprSet[,-c(1,2)]))


save(expSet,phe_GTEx,phe_met,file="met500_GTEx.expSet.Rdata")
save.image(file = "liver_metastasis_2022_06_10.Rdata")
