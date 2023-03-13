
###features of lmscore 
LiverscoreFeatures.path <-"Results/LiverscoreFeatures"
ComparisonResults_score_bio <-"Results/LiverscoreFeatures/ComparisonResults_score_bio"
ComparisonResults_score<-"Results/LiverscoreFeatures/ComparisonResults_score"
if (!file.exists(LiverscoreFeatures.path )) { dir.create(LiverscoreFeatures.path) }
if (!file.exists(ComparisonResults_score_bio )) { dir.create(ComparisonResults_score_bio ) }
if (!file.exists(ComparisonResults_score)) { dir.create(ComparisonResults_score) }

 ##subtype  features
load("immunefeaturs.Rdata")
load(file.path(Cluster.path,"phe_met.Rdata"))
load("metExp500.Rdata")
##########################
load(file.path(Score.path,"phe_total.Rdata"))

colnames(phe_total1)
boxplot(phe_total1$AZ~ phe_total1$group_cluster)
phe_total1$estimate.ImmuneScore<-as.numeric(phe_total1$estimate.ImmuneScore)
phe_total_score<-phe_total1
phe_total_score<-data.frame(phe_total_score)

for(gene in comparVariables){
  print(gene)
  phe_total_score[,gene]<-as.numeric(phe_total_score[,gene])
c<- ggplot(phe_total_score,
           aes_string(x="group_cluster", y=gene, 
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
  #scale_y_continuous(limits = quantile(rt[,gene], c(0.05, 0.95)))+
  #scale_fill_manual(values= c("#E31A1C","#2E9FDF")) +
  #ggtitle("CPSS score", "stratified by InternStudent") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle =30, hjust = 1,size = 10), 
        axis.text.y = element_text(colour="grey",angle = 90, size = 12),
        axis.title.y = element_text(angle = 90, size = 15)) +
  theme(legend.position = "right")
#c+stat_compare_means(method = 'kruskal.test') 
c+stat_compare_means(method = 'wilcox.test',
                     symnum.args=list(),
                     comparisons=list(c("C4", "C3"),c("C4","C2"),c("C4","C5"),c("C4","C1"),c("C4","C6")))

ggsave(file.path(LiverscoreFeatures.path,paste0("ComparisonResults_score/",gene,".pdf")),height=5.5,width=6.24)
ggsave(file.path(LiverscoreFeatures.path,paste0("ComparisonResults_score/",gene,".png")),height=5.5,width=6.24)
}
##"'
gene<-comparVariables[1]
phe_total_score$group_cluster_bio<-ifelse(phe_total_score$group_cluster=="C4","1","0")
for(gene in comparVariables){
  print(gene)
  phe_total_score[,gene]<-as.numeric(phe_total_score[,gene])
  c<- ggplot(phe_total_score,
             aes_string(x="group_cluster_bio", y=gene, 
                        fill ="group_cluster_bio",
                        color ="group_cluster_bio")) +
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
    #scale_y_continuous(limits = quantile(rt[,gene], c(0.05, 0.95)))+
    #scale_fill_manual(values= c("#E31A1C","#2E9FDF")) +
    #ggtitle("CPSS score", "stratified by InternStudent") + 
    theme_classic() +
    theme(axis.text.x = element_text(angle =30, hjust = 1,size = 10), 
          axis.text.y = element_text(colour="grey",angle = 90, size = 12),
          axis.title.y = element_text(angle = 90, size = 15)) +
    theme(legend.position = "right")
  c+stat_compare_means(method = 't.test') 
  #c+stat_compare_means(method = 't.test',
                 #      symnum.args=list(),
                     #  comparisons=list(c("C4", "C3"),c("C4","C2"),c("C4","C5"),c("C4","C1"),c("C4","C6")))
  
  ggsave(file.path(LiverscoreFeatures.path,paste0("ComparisonResults_score_bio/",gene,".pdf")),height=5.5,width=6.24)
  ggsave(file.path(LiverscoreFeatures.path,paste0("ComparisonResults_score_bio/",gene,".png")),height=5.5,width=6.24)
}