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
  geness<-rownames(diffsig)
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

for (set in names(LMiSscoreimmunelist)[-2]){
  print(set)
  myscore<-LMiSscoreimmunelist[[set]]
  myscore$LMiScore<-as.numeric(myscore$LMiScore)
  pdf(paste0(set,"responseofimmunotherapy.pdf"))
  rocobj1=roc(myscore$Cat_Response,myscore$LMiScore,
              plot=TRUE,
              ci=TRUE,
              print.thres=FALSE, print.auc=TRUE,
              xlim=c(1,0),ylim=c(0,1),
              #smooth=TRUE,
              col='#E41A1C')
  tables<-data.frame(table(myscore$Cat_Response))
  prot<-prop.table(table(myscore$Cat_Response)) %>% melt()
  prot$value<-round(prot$value,2)
  means<-aggregate(LMiScore~Cat_Response,data=myscore,mean)
  #sd<-aggregate(LMiScore~Cat_Response,data=myscore,sd)
  prot<-cbind(prot,means)
  prot$LMiScore<-round(prot$LMiScore,2)
  prot<-cbind(prot,tables$Freq)
  prot$responserate<-paste0(prot$`tables$Freq`,"(",prot$value,")")
  prot<-prot[,c(1,4,6)]
  #colnames(prot)<-paste0(set,"_",colnames(prot))
  title(set)
  legend("bottomright",paste(prot$Var1,prot$responserate,prot$LMiScore))
  dev.off()
  
}

names(LMiSscoreimmunelist)
LMiSscore<-LMiSscoreimmunelist[[1]]
####
LMiSscore$Best.Confirmed.Overall.Response
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
c+stat_compare_means(method = 'kruskal.test') 
c+stat_compare_means(method = 't.test',
                     #label.y = c(2.5,3,3.5),
                     #label = "p.signif",
                     symnum.args=list(),
                     comparisons = list(c("C4", "C3"),c("C4","C2"),c("C4", "C1"),c("C4", "C5")))
###


###

cox.sig<-c()
##################
library(survminer)
library(survival)
colnames(LMiSscore)
LMiSscore[,"LMiScore"]<-as.numeric(as.character(LMiSscore[,"LMiScore"]))
sur.cut<-surv_cutpoint(LMiSscore,time = "OS", event = "OS.Event", variables = "LMiScore")
LMiSscore$group=ifelse(LMiSscore[,"LMiScore"]>(sur.cut$cutpoint[,1]),1,0)##bug 


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
