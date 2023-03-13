###
library(readxl)
PANCAN_tpm <- readRDS("~/Projects/PanCancerdata/PANCAN_tpm.RDS")
tcga<-read_excel("/home/data/gaoyuzhen/Projects/MET500Projects/PANCANCER_information.xlsx",1)
tcga<-data.frame(tcga)
tcga<-tcga[tcga$A8_New_Event_Tissue=="Liver",]
tcga$Sample<-paste0(tcga$Sample,"-01A")
colnames(PANCAN_tpm)[1]
tcga_exp<-PANCAN_tpm[,colnames(PANCAN_tpm) %in% tcga$Sample]
tcga<-tcga[tcga$Sample %in% colnames(tcga_exp),]
identical(colnames(tcga_exp),tcga$Sample)
#

tmp1=as.matrix(tcga_exp)
dimnames=list(rownames(tmp1),colnames(tmp1))
data=matrix(as.numeric(as.matrix(tmp1)),nrow=nrow(tmp1),dimnames=dimnames)
tmp1=avereps(tmp1)
tmp1=tmp1[rowMeans(tmp1)>0,]
#####Score
geness<-DEG
tmp1<-data.frame(tmp1)
exprSet<-tmp1[geness,]
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
colnames(LMiSscore)<-c("Sample","LMiScore")
###
tcga$Sample<-gsub("-",".",tcga$Sample)
LMiSscore<-merge(LMiSscore,tcga,by="Sample")
message("TCGA"," GSVA is done!")
save(LMiSscore,file=file.path(TCGAresponse.path,'TCGA_liver_LMiSscore.Rdata'))
#####

##
load(file.path(TCGAresponse.path,'TCGA_liver_LMiSscore.Rdata'))
LMiSscore$LMiScore<-scale(LMiSscore$LMiScore)
write.csv(LMiSscore,file.path(TCGAresponse.path,'TCGA_liver_LMiSscore.csv'))

write.csv(LMiSscore,file = 'TCGA_liver_LMiSscore.csv')
TCGAresponse.path <-"Results/TCGAresponse.path"
if (!file.exists(TCGAresponse.path )) { dir.create(TCGAresponse.path) }


cox.sig<-c()
##################
library(survminer)
library(survival)
colnames(LMiSscore)
LMiSscore$OS.time<-as.numeric(LMiSscore$OS.time)
LMiSscore$OS<-ifelse(LMiSscore$OS=="Dead",1,0)
#LMiSscore<-LMiSscore[!is.na(LMiSscore$Met.Disease.Status),]
LMiSscore[,"LMiScore"]<-as.numeric(as.character(LMiSscore[,"LMiScore"]))
#LMiSscore<-LMiSscore[LMiSscore$Livermeta=="Yes",]
sur.cut<-surv_cutpoint(LMiSscore,time = "OS.time", event = "OS", variables = "LMiScore")
LMiSscore$group=ifelse(LMiSscore[,"LMiScore"]>(sur.cut$cutpoint[,1]),1,0)
#LMiSscore$group=ifelse(LMiSscore[,"LMiScore"]>median(LMiSscore[,"LMiScore"]),1,0)
table(LMiSscore$group)
Gcox1<-coxph(Surv(OS.time,OS)~group,data=LMiSscore)
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
#LMiSscore$liver<-ifelse(LMiSscore$Met.Disease.Status=="Liver",1,0)
colnames(LMiSscore)
kmfit<- survfit(Surv(OS.time/30,OS)~group,data=LMiSscore)
p.val<-Pvalue
HR <- paste("Hazard Ratio= ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low,2), round(up,2), sep = " - "), sep = "")
print(paste0(" p= ",round(p.val,4)," ",HR))





pdf(file.path(TCGAresponse.path,"repsonse_OS_TCGA_liver.pdf"),width = 6.49, height = 6.58,onefile = F )
ggsurvplot(kmfit,#surv.median.line = "hv", # Add medians survival
           # Change legends: title & labels
           #main = "Survival curve",
           #legend.title =set,
           legend.labs =c("Low-LCMS","High-LCMS"),
           #xlim = c(0, 200),
           xlab="Time(months)",
           ylab="Overall Survival",
           size = 0.5,
           #fun="cumhaz",
           #fun='event',
           # Add p-value and tervals
           #pval = TRUE,
           #test.for.trend = TRUE,###group more than 2groups
           break.time.by =20,
           #conf.int = TRUE,
           #group.by=,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.185,
           tables.theme = theme_cleantable(),
           palette = c("#1F78B4", "#E31A1C"),
           #ggtheme = theme_bw(), # Change ggplot2 theme
           #font.title="OS",
           font.main =15,
           font.x =  15,
           font.y = 15,
           font.tickslab =15,
           #在左下???标出pvalue、HR???95% CI
           #???小的p value标为p < 0.001
           pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                      paste("p = ",round(p.val,4), sep = "")), HR, CI, sep = "\n"))


dev.off()

