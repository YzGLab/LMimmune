##GSE131418 大规模测序 结肠癌转移 （原位组织） 
GEOliverscore.path <-"Results/GEOliverscore.path"
if (!file.exists(GEOliverscore.path )) { dir.create(GEOliverscore.path) }
##
liver_exp<-read.delim2("/home/data/gaoyuzhen/Projects/MET500Projects/mets-immunecluster-master/data/ex_gse51244.1.txt")
clin_liver<-read.delim2("/home/data/gaoyuzhen/Projects/MET500Projects/mets-immunecluster-master/data/clin_gse51244.1.txt")

#sub1exp<-liver_exp["SUB1",] %>% t() %>% data.frame() %>% cbind(GEO_ID=rownames(.),.)
colnames(clin_liver)
#clin_liver<-merge(sub1exp,clin_liver,by="GEO_ID")
#clin_liver$SUB1<-as.numeric(clin_liver$SUB1)
 # boxplot(clin_liver$SUB1~clin_liver$TISSUE)

  liver_exp2<-apply(liver_exp,2,as.numeric)
  rownames(liver_exp2)<-rownames(liver_exp)
  tmp1=as.matrix(liver_exp2)
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
  colnames(LMiSscore)<-c("GEO_ID","LMiScore")
  ###
  LMiSscore<-merge(LMiSscore,clin_liver,by="GEO_ID")

  save(LMiSscore,file=file.path(GEOliverscore.path,'GEOliverscore.Rdata'))
####################
  colnames(LMiSscore)
  c<- ggplot(LMiSscore,
             aes_string(x="TISSUE", y="LMiScore", 
                        fill ="TISSUE",
                        color ="TISSUE")) +
    geom_point(aes(shape=TISSUE), 
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