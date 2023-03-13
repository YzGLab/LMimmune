source("Rcode_meta/twoclasslimma.R")
library(gghalves)
library(ggpubr)
######claculate the score for the subtype .
load(file.path(DEG_cluster.path_gene ,"DEGenes_special_liver_metastasis.Rdata" ))
##
library(readxl)
load("metExp500.Rdata")
load("new_gsva2.Rdata")
load("met500_GTEx.expSet.Rdata")
## first step  different
mycol<-c( '#D95F02','#7570B3','#E7298A', '#66A61E','#E6AB02','#A6761D','#666666',
          '#377EB8','#4DAF4A', '#984EA3','#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999',
          '#8DD3C7', '#FFFFB3','#BEBADA', '#FB8072','#80B1D3','#FDB462','#B3DE69', '#FCCDE5',
          '#D9D9D9','#BC80BD','#CCEBC5', '#FFED6F','#66C2A5', '#FC8D62','#8DA0CB','#E78AC3',
          '#A6D854','#FFD92F', '#E5C494', '#B3B3B3','#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C',
          '#1B9E77','#7FC97F','#BEAED4','#FDC086','#FFFF99','#386CB0',
          '#FB9A99', '#E31A1C','#FDBF6F', '#FF7F00', '#CAB2D6', '#6A3D9A', '#FFFF99', '#B15928')
################################


########计算评分 并计算差异
Score.path <-"Results/ScoreFeatures"
if (!file.exists(Score.path )) { dir.create(Score.path) }

################
myexp<-metExp
rt=as.matrix(myexp)
dimnames=list(rownames(rt),colnames(rt))
data=matrix(as.numeric(as.matrix(rt)),nrow=nrow(rt),dimnames=dimnames)
data=avereps(data)
myexp=data[rowMeans(data)>0,]
#####Score
geness<-DEG
myexp<-data.frame(data)
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
colnames(LMiSscore)<-c("Sample_id","LMiScore")
###############
# save(coef,file="LMiScore_coef.Rdata")
write.csv(coef,file=file.path(Score.path ,"DEGenes_coef.csv"))####
#################################################

##simple test 
load("immunefeaturs.Rdata")
colnames(immunefeaturs)
load(file.path(Cluster.path,"phe_met_cluster_info.Rdata"))
###########################################
phe_met$group_cluster<-paste0("C",phe_met$group_cluster)
#immunefeaturs<-cbind(Sample_id=rownames(immunefeaturs),immunefeaturs)
colnames(phe_met)
phe_met$Sample_id<-gsub("-",".",phe_met$Sample_id)
phe_total<-merge(phe_met,immunefeaturs,by="Sample_id")
##az  asZXaszwaq
#phe_met$group_cluster<-factor(phe_met$group_cluster,levels=c("1liver_very_high","2liver_high","3liver_median","4liver_low"))
#################################
table(phe_met$group_cluster)

##计算score 特征 ggplot2
phe_total1<-merge(LMiSscore,phe_total,by="Sample_id")
phe_total1$group_cluster<-factor(phe_total1$group_cluster,levels=paste0("C",a))

phe_total1$LMiScore<-as.numeric(phe_total1$LMiScore)
boxplot(phe_total1$LMiScore~ phe_total1$group_cluster)
#,ylim=c(-500,1000)
save(phe_total1,file=file.path(Score.path,"phe_total.Rdata"))
write.csv(phe_total1,file=file.path(Score.path,"phe_total.csv"))




##############################  
shape_level <- nlevels(factor(phe_total1$tissue))
if (shape_level < 15){
  shapes = (0:shape_level) %% 15
} else{
  shapes = c(0:14,c((15:shape_level) %% 110 + 18))
}
##
colnames(phe_total1)
phe_total1$LMiScore<-log(phe_total1$LMiScore+2000)
c<- ggplot(phe_total1,
           aes_string(x="group_cluster", y="LMiScore", 
                      fill ="group_cluster",
                      color ="group_cluster")) +
  geom_point(aes(shape=biopsy_tissue), 
             position = position_jitter(width = .15,height=-0.7),
             size=1)+
  xlab("")+
  #scale_shape_manual(values=shapes) +
  scale_color_manual(values=mycol) +
  scale_fill_manual(values=mycol) +
  geom_boxplot(notch = F, alpha = 0.65, width=0.2,
               outlier.size = 0.65)+
  geom_half_violin(position = position_nudge(x = -.15),alpha = 0.65)+
  #scale_y_continuous(limits = quantile(phe_total1[,"LMiScore"], c(0.10, 0.99)))+
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
                     comparisons = list(c("C4", "C3"),c("C4","C2"),c("C4","C5"),c("C4","C1"),c("C4","C6")))

ggsave(file.path(Score.path,"meta_immune_score.pdf"),
       width = 7, height =5,
       bg = "transparent", dpi = 300, useDingbats=FALSE)



###tissue site
colnames(phe_total1)
b<-data.frame(table(phe_total1$biopsy_tissue,phe_total1$tissue))
b<-b[b$Freq>21,]
vars<-unique(b$Var2)[1:4]
subphe_total1<-phe_total1[phe_total1$biopsy_tissue=="liver" & phe_total1$tissue %in% vars ,]
c<- ggplot(subphe_total1,
           aes_string(x="tissue", y="LMiScore", 
                      fill ="tissue",
                      color ="tissue")) +
  geom_point(aes(shape=tissue), 
             position = position_jitter(width = .15,height=-0.7),
             size=1)+
  xlab("")+
  #scale_shape_manual(values=shapes) +
  scale_color_manual(values=mycol) +
  scale_fill_manual(values=mycol) +
  geom_boxplot(notch = F, alpha = 0.65, width=0.2,
               outlier.size = 0.65)+
  geom_half_violin(position = position_nudge(x = -.15),alpha = 0.65)+
  #scale_y_continuous(limits = quantile(phe_total1[,"LMiScore"], c(0.10, 0.99)))+
  #scale_fill_manual(values= c("#E31A1C","#2E9FDF")) +
  #ggtitle("CPSS score", "stratified by InternStudent") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle =30, hjust = 1,size = 10), 
        axis.text.y = element_text(colour="grey",angle = 90, size = 12),
        axis.title.y = element_text(angle = 90, size = 15)) +
  theme(legend.position = "right")
c+stat_compare_means(method = 'kruskal.test') 


ggsave(file.path(Score.path,"meta_site_immune_score.pdf"),
       width =5, height =5,
       bg = "transparent", dpi = 300, useDingbats=FALSE)
####
c<- ggplot(phe_total1,
           aes_string(x="biopsy_tissue", y="LMiScore", 
                      fill ="biopsy_tissue", 
                      color ="biopsy_tissue")) +
  geom_jitter(aes_string(shape="tissue"), 
              position = position_jitter(width = .15,height=-0.7),
              size=1)+
  #ylim(0,4000)+
  xlab("")+
  scale_shape_manual(values=shapes) +
  scale_color_manual(values=mycol) +
  scale_fill_manual(values=mycol) +
  geom_boxplot(notch = F, alpha = 0.65, width=0.2,
               outlier.size = 0.65)+
  geom_half_violin(position = position_nudge(x = -.15),alpha = 0.65)+
  #scale_y_continuous(limits = quantile(rt[,gene], c(0.05, 0.95)))+
  #ggtitle("CPSS score", "stratified by InternStudent") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle =30, hjust = 1,size = 10), 
        axis.text.y = element_text(colour="grey",angle = 90, size = 12),
        axis.title.y = element_text(angle = 90, size = 15)) +
  theme(legend.position = "right")
c+stat_compare_means(method = 'kruskal.test')
c+stat_compare_means(method = 'wilcox.test', 
                     #label.y = c(2.5,3,3.5),
                     #label = "p.signif",
                     #symnum.args=list(),
                     comparisons = list(c("liver", "bone_marrow"),
                                        c("liver", "brain"),
                                        c("liver", "lung"),
                                        c("liver", "skin")))
ggsave(file.path(Score.path,"Liver_meta_site_immune_score.pdf"),
       width =7, height =5,
       bg = "transparent", dpi = 300, useDingbats=FALSE)

################ gsva 方法

norm.expMat<-as.matrix(myexp)
gsva_es <- gsva(norm.expMat,list(DEG),method="ssgsea",abs.ranking=F,min.sz = 5,ssgsea.norm=TRUE,kcdf="Poisson",parallel.sz=70)###Array data,ssgsea.norm=TRUE
#gsva_es <- gsva(norm.expMat,gmt.list,method="ssgsea",abs.ranking=F,kcdf="Poisson",parallel.sz=10)#RNA-Seq 
gsva_es<-t(gsva_es)
gsva_es<-cbind(Sample_id=gsub("-",".",rownames(gsva_es)),gsva_es)##
colnames(gsva_es)<-c("Sample_id","LMiScore")
phe_total2<-merge(gsva_es,phe_total,by="Sample_id")
#phe_total2$group_cluster<-factor(phe_total2$group_cluster,levels=c("C3","C1","C4","C2"))
phe_total2$group_cluster<-factor(phe_total2$group_cluster,levels=paste0("C",a))
phe_total2$LMiScore<-scale(as.numeric(phe_total2$LMiScore))
boxplot(phe_total2$LMiScore~ phe_total2$group_cluster)
phe_total2$LMiScore<-log(phe_total1$LMiScore)
c<- ggplot(phe_total2,
           aes_string(x="group_cluster", y="LMiScore", 
                      fill ="group_cluster",
                      color ="group_cluster")) +
  geom_point(aes_string(shape="group_cluster"), 
             position = position_jitter(width = .15,height=-0.7),
             size=1)+
  xlab("")+
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
                     comparisons = list(c("C3", "C1"),c("C3","C4"),c("C3", "C2")))

#save(phe_total,phe_total1,phe_total2,file="phe_total.Rdata")
#######################
colnames(phe_total1)
###
table(phe_total1$biopsy_tissue)
subphe<-phe_total1[phe_total1$biopsy_tissue=="liver",]
#####
c<- ggplot(subphe,
           aes_string(x="group_cluster", y="LMiScore", 
                      fill ="group_cluster",
                      color ="group_cluster")) +
  geom_point(aes(shape=biopsy_tissue), 
             position = position_jitter(width = .15,height=-0.7),
             size=1)+
  xlab("")+
  #scale_shape_manual(values=shapes) +
  scale_color_manual(values=mycol) +
  scale_fill_manual(values=mycol) +
  geom_boxplot(notch = F, alpha = 0.65, width=0.2,
               outlier.size = 0.65)+
  geom_half_violin(position = position_nudge(x = -.15),alpha = 0.65)+
  #scale_y_continuous(limits = quantile(phe_total1[,"LMiScore"], c(0.10, 0.99)))+
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
                     comparisons = list(c("C4", "C3"),c("C4","C2"),c("C4","C5"),c("C4","C1"),c("C4","C6")))


ggsave(file.path(Score.path,"meta_immune_score.pdf"),
       width = 7, height =5,
       bg = "transparent", dpi = 300, useDingbats=FALSE)

###corrrelationship 


library(readxl)
library(ggcorrplot)
library(ggthemes)
library(Hmisc)
display.brewer.all()
library(RColorBrewer)
library(corrplot)#
library("ggsci")
library("scales")
library("ggplot2")
library("gridExtra")
library(data.table)
library(GSVA)
library(ggplot2)
library(ggcor)
###
Degexp<-metExp[DEG,gsub("-",".",colnames(metExp)) %in% phe_met$Sample_id]
corrplot::corrplot(cor(t(Degexp)))
rt<-data.frame(t(Degexp))
head(rt) 
correlation<-rcorr(as.matrix(rt))
correlation_r<-correlation$r
correlation_p_values<-correlation$P
#write.csv(correlation_r,file = "E_correlation.csv")
#write.csv(correlation_p_values,file = "E_correlation_p_values.csv")
ggcorrplot(correlation_r, method = "circle",
           outline.col = "white", ggtheme = theme_bw(),
           type="upper",
           colors = c("#6D9EC1", "white", "#E46726"),
           lab = TRUE, lab_size = 2,) 

ggcorrplot(correlation_r, method = "circle",
           outline.col = "white", ggtheme = theme_bw(),
           colors = c("#6D9EC1", "white", "#E46726"),
           lab = TRUE, lab_size = 2,
           p.mat = correlation_p_values, insig = "blank") 
?ggcorrplot
ggcorrplot(correlation_r, method = "c", type = "upper")
ggcorrplot(correlation_r, add = TRUE, type = "lower", method = "number", diag = FALSE, tl.pos = "n", cl.pos = "n")
###############


pdf(file=file.path(DEG_cluster.path_gene,"genes_cor.pdf"),width = 8,height = 8)
corrplot(as.matrix(correlation_r), 
         method="number", 
         type="upper", 
         #col=colorRampPalette(c("yellow2","goldenrod","green"))(8), 
         col=brewer.pal(n=9, name="YlGnBu"), 
         tl.col="black", 
         #tl.srt=45, 
         #diag=FALSE,
         tl.pos='l',
         pch.cex = 1, 
         pch = 2,
         number.cex = 0.5,
         number.font = 1
         #p.mat = as.matrix(correlation_p_values)
         #sig.level = 0.05 
         #insig = "label_sig"
)
corrplot(as.matrix(correlation_r), 
         method="pie", 
         type="lower", 
         #method = "number",
         #col=colorRampPalette(c("yellow2","goldenrod","green"))(8), 
         col=brewer.pal(n=9, name="YlGnBu"),
         tl.col="black", 
         tl.srt=45, 
         tl.pos='n',
         cl.pos = 'n',
         #diag=FALSE,
         p.mat = as.matrix(correlation_p_values), 
         sig.level = 0.05, 
         add = TRUE, 
         insig = "blank")
dev.off()
##

immPath.score<-cbind(phe_total1[,c(1,2)],rt)
immPath.score
immCorLMiScore <- NULL
for (i in colnames(immPath.score)[-1]) {
  cr <- cor.test(as.numeric(immPath.score[,i]),
                 as.numeric(immPath.score[,"LMiScore"]),
                 method = "pearson")
  immCorLMiScore <- rbind.data.frame(immCorLMiScore,
                                     data.frame(gene = "LMiScore",
                                                path = i,
                                                r = cr$estimate,
                                                p = cr$p.value,
                                                stringsAsFactors = F),
                                     stringsAsFactors = F)
}
immCorLMiScore$sign <- ifelse(immCorLMiScore$r > 0,"pos","neg")
immCorLMiScore$absR <- abs(immCorLMiScore$r)
immCorLMiScore$rSeg <- as.character(cut(immCorLMiScore$absR,c(0,0.25,0.5,0.75,1),labels = c("0.25","0.50","0.75","1.00"),include.lowest = T))
immCorLMiScore$pSeg <- as.character(cut(immCorLMiScore$p,c(0,0.001,0.01,0.05,1),labels = c("<0.001","<0.01","<0.05","ns"),include.lowest = T))
immCorLMiScore[nrow(immCorLMiScore),"pSeg"] <- "Not Applicable"

immCorLMiScore$rSeg <- factor(immCorLMiScore$rSeg, levels = c("0.25","0.50","0.75","1.00"))
immCorLMiScore$pSeg <- factor(immCorLMiScore$pSeg, levels = c("<0.001","<0.01","<0.05","Not Applicable","ns"))
immCorLMiScore$sign <- factor(immCorLMiScore$sign, levels = c("pos","neg"))

p1 <- quickcor(immPath.score[,-1], 
               type = "lower",
               show.diag = TRUE) + 
  geom_colour() +
  add_link(df = immCorLMiScore, 
           mapping = aes(colour = pSeg, size = rSeg, linetype = sign),
           spec.key = "gene",
           env.key = "path",
           diag.label = FALSE) +
  scale_size_manual(values = c(0.5, 1, 1.5, 2)) +
  scale_color_manual(values = mycol[6:12]) +
  scale_fill_gradient2(low = "#9483E1",mid = "white",high = "#E31A1C",midpoint=0) +
  remove_axis("x")
p1

ggsave(file=file.path(DEG_cluster.path_gene,"genes_cor_score.pdf"),width = 8,height = 8)

