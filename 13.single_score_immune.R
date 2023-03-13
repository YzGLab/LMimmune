#####immune predicted###
### for immune score 
load(file.path(DEG_cluster.path_gene ,"DEGenes_special_liver_metastasis.Rdata" ))
load("~/Projects/ImmunePro/immunecohorts_allRNASeqandPhe20220224.Rdata")
#
#save(immunecohorts,immunephe,file="~/Projects/ImmunePro/immunecohorts_allRNASeqandPhe20220224.Rdata")

names(immunecohorts) 
LMiSscoreimmunelist<-list()
for(set in names(immunecohorts)[-c(5,10,16)]){
  myexp<-immunecohorts[[set]]
  rt=as.matrix(myexp)
  dimnames=list(rownames(rt),colnames(rt))
  data=matrix(as.numeric(as.matrix(rt)),nrow=nrow(rt),dimnames=dimnames)
  data=avereps(data)
  myexp=data[rowMeans(data)>0,]
  #####Score
  geness<-DEG
  myexp<-data.frame(myexp)
  exprSet<-myexp[geness,]
  exprSet<-data.frame(exprSet)
  exprSet<-na.omit(exprSet)
  pca.results <- prcomp(t(exprSet), center = TRUE, scale. = FALSE)
  coef <- pca.results$rotation[rownames(exprSet), c("PC1", "PC2", "PC3")] %>% as.data.frame()
  coef$PC <- coef$PC1 + coef$PC2 + coef$PC3 
  #write.csv(coef %>% rownames_to_column("gene"), "PCA_coef.csv", row.names = F)##save the coeficient
  l <- lapply(names(exprSet), FUN = function(x) {
    coef$PC %*% exprSet[, x]
  })
  #score <- do.call(c, l)
  score<-sapply(l,cbind)
  score<-data.frame(score)
  rownames(score)<-colnames(exprSet)
  LMiSscore <- cbind(rownames(score),score)
  colnames(LMiSscore)<-c("geo_accession","LMiScore")
  ###
  myphe<-immunephe[[set]]
  myphe$geo_accession<-rownames(myphe)
  LMiSscore<-merge(LMiSscore,myphe,by="geo_accession")
  LMiSscoreimmunelist[[set]]<-LMiSscore
  message(set,"GSVA is done!")
}

##
library(pROC)
library(tidyverse)
library(data.table)
library(tidyr)




#######################################################
LMiSscore<-LMiSscoreimmunelist[[1]]
####
LMiSscore<-LMiSscore[LMiSscore$Best.Confirmed.Overall.Response!="NE",]
LMiSscore$LMiScore<-as.numeric(LMiSscore$LMiScore)
boxplot(LMiSscore$LMiScore~ LMiSscore$Cat_Response)
LMiSscore$Cat_Response<-as.factor(LMiSscore$Cat_Response)

c<- ggplot(LMiSscore,
           aes_string(x="Best.Confirmed.Overall.Response", y="LMiScore", 
                      fill ="Best.Confirmed.Overall.Response",
                      color ="Best.Confirmed.Overall.Response")) +
  geom_point(aes(shape=Best.Confirmed.Overall.Response), 
             position = position_jitter(width = .15,height=-0.7),
             size=1)+
  xlab("")+
  #scale_shape_manual(values=shapes) +
  scale_color_manual(values=mycol) +
  scale_fill_manual(values=mycol) +
  geom_boxplot(notch = F, alpha = 0.65, width=0.2,
               outlier.size = 0.65)+
  geom_half_violin(position = position_nudge(x = -.15),alpha = 0.65)+
  #scale_y_continuous(limits = quantile(rt[,gene], c(0.05, 0.95)))+
  #scale_fill_manual(values= c("#E31A1C","#2E9FDF")) +
  #ggtitle("CPSS score", "stratified by InternStudent") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle =30, hjust = 1,size = 10), 
        axis.text.y = element_text(colour="grey",angle = 90, size = 12),
        axis.title.y = element_text(angle = 90, size = 15)) +
  theme(legend.position = "right")
c+stat_compare_means(method = 'anova') #kruskal.test
c+stat_compare_means(method = 't.test',
                     #label.y = c(2.5,3,3.5),
                     #label = "p.signif",
                     symnum.args=list(),
                     comparisons = list(c("C4", "C3"),c("C4","C2"),c("C4", "C1")))
###

names(LMiSscoreimmunelist)
###
LMiSscore<-LMiSscoreimmunelist[["Liu_NatMedicine_Met_Melanoma_2019"]]
LMiSscore<-LMiSscoreimmunelist[["Braun_RCC_Nat_medicine2020"]]
cox.sig<-c()
##################
library(survminer)
library(survival)
colnames(LMiSscore)
#LMiSscore<-LMiSscore[!is.na(LMiSscore$Met.Disease.Status),]
LMiSscore[,"LMiScore"]<-as.numeric(as.character(LMiSscore[,"LMiScore"]))
LMiSscore<-LMiSscore[LMiSscore$met.liver.cumul=="Y",]
sur.cut<-surv_cutpoint(LMiSscore,time = "OS", event = "OS.Event", variables = "LMiScore")
LMiSscore$group=ifelse(LMiSscore[,"LMiScore"]>(sur.cut$cutpoint[,1]),1,0)
LMiSscore$group=ifelse(LMiSscore[,"LMiScore"]>median(LMiSscore[,"LMiScore"]),1,0)
table(LMiSscore$group)
Gcox1<-coxph(Surv(OS,OS.Event)~group,data=LMiSscore)
GSum<-summary(Gcox1)
HR<-round(GSum$coefficients[,2],3)
Pvalue<-round(GSum$coefficients[,5],4)
CI<-paste0(round(GSum$conf.int[,3:4],3),collapse='-')
coeff<-round(GSum$coefficients[1],3)
se<-GSum$coefficients[1,3]
low<-round(GSum$conf.int[,3],3)
up<-round(GSum$conf.int[,4],3)
cox.p<-data.frame('Characteristics'= set,
                  'Hazard Ratio'= HR,
                  'CI95'=CI,
                  "coeff"=coeff,
                  "se"=se,
                  "low"=low,
                  "up"=up,
                  'P-value'=Pvalue,
                  "set"=set,
                  "numberofpatients"=nrow(rt))
cox.sig=rbind(cox.sig,cox.p)
message(set,"cox regression is Done!")
#LMiSscore$liver<-ifelse(LMiSscore$Met.Disease.Status=="Liver",1,0)
kmfit<- survfit(Surv(OS,OS.Event)~group,data=LMiSscore)
p.val<-Pvalue
HR <- paste("Hazard Ratio= ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low,2), round(up,2), sep = " - "), sep = "")
print(paste0(set," p= ",round(p.val,2)," ",HR))
ggsurvplot(kmfit,#surv.median.line = "hv", # Add medians survival
           # Change legends: title & labels
           #main = "Survival curve",
           #legend.title =set,
           #legend.labs = c("Low","High"),
           #xlim = c(0, 200),
           xlab="Time(months)",
           ylab="Overall Survival",
           size = 0.5,
           #fun="cumhaz",
           #fun='event',
           # Add p-value and tervals
           #pval = TRUE,
           #test.for.trend = TRUE,###group more than 2groups
           #break.time.by = 30,
           #conf.int = TRUE,
           #group.by=,
           # Add risk table
           #risk.table = TRUE,
           tables.height = 0.185,
           tables.theme = theme_cleantable(),
           palette = c("#1F78B4", "#E31A1C"),
           #ggtheme = theme_bw(), # Change ggplot2 theme
           #font.title="OS",
           font.main =15,
           font.x =  15,
           font.y = 15,
           font.tickslab =25,
           #在左下???标出pvalue、HR???95% CI
           #???小的p value标为p < 0.001
           pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                      paste("p = ",round(p.val,3), sep = "")), HR, CI, sep = "\n"))


dev.off()
###




























clin<-read.delim2("/home/data/gaoyuzhen/Projects/MET500Projects/mets-immunecluster-master/data/liu_clinical.txt")
####
gmt_immune<-read.gmt("c7.all.v7.5.1.symbols.gmt")
gmt_immune<-split(gmt_immune$gene,gmt_immune$term)

load(file.path(DEG.path,"DEGs_pathway.Rdata"))
###
names(LMiSscoreimmunelist)
DEGs_pathway_list<-gmt_immune[DEGs_pathway]
#myexp<-read.delim2("/home/data/gaoyuzhen/Projects/MET500Projects/mets-immunecluster-master/data/liu_expression_logTPM.txt")
#row.names<-rownames(myexp)
#myexp<-apply(myexp,2,as.numeric,na.rm = TRUE)
#rownames(myexp)<-row.names
LMiSscore<-LMiSscoreimmunelist[["Liu_NatMedicine_Met_Melanoma_2019"]]
myexp<-immunecohorts[["Liu_NatMedicine_Met_Melanoma_2019"]]
norm.expMat<-as.matrix(myexp)
gsva_es <- gsva(norm.expMat,DEGs_pathway_list,method="gsva",abs.ranking=F,min.sz = 5,ssgsea.norm=TRUE,kcdf="Poisson",parallel.sz=10)###Array data,ssgsea.norm=TRUE
##cluster for it 
dist.e=dist(t(gsva_es),method='euclidean')
hclust(dist.e, method = "complete", members = NULL)
tree <- hclust(dist.e, method = "ward.D2")
plot(tree,labels = FALSE, hang = -3, main = "MetaSubtype")
#plot(tree, hang = -5)
###############################
group<-cutree(tree, k =3)
group[group==3] = "High"
group[group==2] = "Low"
group[group==1] = "Medium"
group<-factor(group,levels = c("Low","Medium","High"))
#tree<-cutree(tree,h=0.7)
###########################################
group_cluster<-group
group_cluster<-data.frame(group_cluster)
group_cluster <-cbind(geo_accession=rownames(group_cluster),group_cluster)

LMiSscore<-merge(group_cluster,LMiSscore,by="geo_accession")

kmfit<- survfit(Surv(OS,OS.Event)~group_cluster,data=LMiSscore)
p.val<-Pvalue
HR <- paste("Hazard Ratio= ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low,2), round(up,2), sep = " - "), sep = "")
print(paste0(set," p= ",round(p.val,2)," ",HR))
ggsurvplot(kmfit,#surv.median.line = "hv", # Add medians survival
           # Change legends: title & labels
           #main = "Survival curve",
           #legend.title =set,
           #legend.labs = c("Low","High"),
           #xlim = c(0, 200),
           xlab="Time(months)",
           ylab="Overall Survival",
           size = 0.5,
           #fun="cumhaz",
           #fun='event',
           # Add p-value and tervals
           pval = TRUE,
           #test.for.trend = TRUE,###group more than 2groups
           #break.time.by = 30,
           #conf.int = TRUE,
           #group.by=,
           # Add risk table
           #risk.table = TRUE,
           tables.height = 0.185,
           tables.theme = theme_cleantable(),
           #palette = c("#1F78B4", "#E31A1C"),
           #ggtheme = theme_bw(), # Change ggplot2 theme
           #font.title="OS",
           font.main =15,
           font.x =  15,
           font.y = 15,
           font.tickslab =15
           #在左下???标出pvalue、HR???95% CI
           #???小的p value标为p < 0.001
           )
