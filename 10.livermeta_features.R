

###features of lmscore 
LiverFeatures.path <-"Results/LiverFeatures"
ComparisonResults_bio <-"Results/LiverFeatures/ComparisonResults_bio"
ComparisonResults<-"Results/LiverFeatures/ComparisonResults"
if (!file.exists(LiverFeatures.path )) { dir.create(LiverFeatures.path) }
if (!file.exists(ComparisonResults_bio )) { dir.create(ComparisonResults_bio) }
if (!file.exists(ComparisonResults )) { dir.create(ComparisonResults) }
### phe analysis lmiscore.
########################################
library(gghalves)
library(ggpubr)
library(RColorBrewer)
##
################
tables<-data.frame(table(phe_total1$biopsy_tissue))
tables<-tables[tables$Freq>5,] 
phe_total3<-phe_total1[phe_total1$biopsy_tissue %in% tables$Var1,]

chisq.test(table(phe_total3$ImmuneSubtype.BestCall,phe_total3$group_cluster))
phe_total3<-data.frame(phe_total3)
phe_total3$tissue<-ifelse(phe_total3$tissue=="","other",phe_total3$tissue)
colnames(phe_total3)
comparVariables<-colnames(phe_total3)[15:220]
#
shape_level <- nlevels(factor(phe_total3$tissue))
if (shape_level < 15){
  shapes = (0:shape_level) %% 15
} else{
  shapes = c(0:14,c((15:shape_level) %% 110 + 18))
}
# gene<-comparVariables[1]
for(gene in comparVariables){
  print(gene)
  phe_total3[,gene]<-as.numeric(phe_total3[,gene])
  ################
  c<- ggplot(phe_total3,
             aes_string(x="biopsy_tissue", y=gene, 
                        fill ="biopsy_tissue", 
                        color ="biopsy_tissue")) +
    geom_point(aes(shape=tissue), 
               position = position_jitter(width = .15,height=-0.7),
               size=1)+
    xlab("")+
    scale_shape_manual(values=shapes) +
    scale_color_manual(values=c(mycol)) +
    scale_fill_manual(values=c(mycol)) +
    geom_boxplot(notch = F, alpha = 0.65, width=0.2, outlier.size = 0.65)+
    geom_half_violin(position = position_nudge(x = -.15),alpha = 0.65)+
    #ggtitle("CPSS score", "stratified by InternStudent") + 
    theme_light() +
    theme(axis.text.x = element_text(angle =30, hjust = 1,size = 10), 
          axis.text.y = element_text(colour="grey",angle = 90, size = 12),
          axis.title.y = element_text(angle = 90, size = 15)) +
    theme(legend.position = "right")
  #c+stat_compare_means(method = 'kruskal.test')
  c+stat_compare_means(method = 't.test', 
                       #label.y = c(2.5,3,3.5),
                       #label = "p.signif",
                       #symnum.args=list(),
                       comparisons = list(c("liver", "bone_marrow"),
                                          c("liver", "brain"),
                                          c("liver", "lung"),
                                          c("liver", "skin")))
  ggsave(file.path(LiverFeatures.path ,paste0("ComparisonResults/",gene,".pdf")),height=6.5,width=6.5)
  ggsave(file.path(LiverFeatures.path ,paste0("ComparisonResults/",gene,".png")),height=6.5,width=6.5)
}
###

for(gene in comparVariables){
  print(gene)
  phe_total3[,gene]<-as.numeric(phe_total3[,gene])
  ################
  c<- ggplot(phe_total3,
             aes_string(x="biopsy_tissue_bio", y=gene, 
                        fill ="biopsy_tissue_bio", 
                        color ="biopsy_tissue_bio")) +
    geom_jitter(aes_string(shape="tissue"), 
                position = position_jitter(width = .15,height=-0.7),
                size=1)+
    #ylim(0,4000)+
    xlab("")+
    scale_shape_manual(values=shapes) +
    scale_color_manual(values=c(mycol)) +
    scale_fill_manual(values=c(mycol)) +
    geom_boxplot(notch = F, alpha = 0.65, width=0.2, outlier.size = 0.65)+
    geom_half_violin(position = position_nudge(x = -.15),alpha = 0.65)+
    #scale_y_continuous(limits = quantile(rt[,gene], c(0.05, 0.95)))+
    #scale_fill_manual(values= c("#E31A1C","#2E9FDF")) +
    #ggtitle("CPSS score", "stratified by InternStudent") + 
    theme_classic() +
    theme(axis.text.x = element_text(angle =30, hjust = 1,size = 10), 
          axis.text.y = element_text(colour="grey",angle = 90, size = 12),
          axis.title.y = element_text(angle = 90, size = 15)) +
    #theme(legend.position = "top")
    theme(legend.position = "right")
  c+stat_compare_means(method = 't.test')
  ggsave(file.path(LiverFeatures.path,paste0("ComparisonResults_bio/",gene,".pdf")),height=5.5,width=6.24)
  ggsave(file.path(LiverFeatures.path,paste0("ComparisonResults_bio/",gene,".png")),height=5.5,width=6.24)
}
dev.off()

#######
####比较差异
library(plyr)
library(survival)
library(ComplexHeatmap)
load(file.path(Score.path,"phe_total.Rdata"))
colnames(phe_total1)
rt<-phe_total1
#rt<-phe_total1[phe_total1$group_cluster %in% c("C4","C3"),]
########
k.testfuc<-function(x) {
  rt[,x]<-as.numeric(rt[,x])
  if(mean(as.numeric(rt[,x]))!=0 & !is.na(mean(as.numeric(rt[,x]))) ){
    k.test<-kruskal.test(as.numeric(rt[,x])~rt[,"group_cluster"])
    data.frame(k.test$p.value,k.test$statistic,x)
  }
}
k.test<-kruskal.test(as.numeric(rt[,"LMiScore"])~rt[,"group_cluster"])
##
heatmap_var<- colnames(rt)[c(15:382)]
##
heatmap_vars<-lapply(heatmap_var,k.testfuc) %>% ldply(.,cbind)###################
#heatmap_vars<-lapply(heatmap_vars,cbind)
#write.csv(heatmap_vars,file=file.path(Score.path,"phe_total_cat3.csv"))
heatmap_vars<-heatmap_vars[heatmap_vars$k.test.p.value<0.001,]

##complexheatmap for subtype
#save(diffSig_liver,file=file.path(DEG.path_gene,"GTEX_diffSig_liver.Rdata"))
rt<-phe_total1
#rt$LMiScore<-scale(rt$LMiScore)
rt$LMiScore<-scale(rt$LMiScore)
rt<-rt[order(rt$LMiScore,decreasing = F),]
#rt<-rt[order(rt$group_cluster,decreasing = T),]


###
var_cat<-read.csv(file.path(Score.path,"phe_total_cat2.csv"))
var_cat<-var_cat[!var_cat$Class %in% c("quantiseq","ImmuneSubtype","Pathway","gene","IPS","Metabolism","ImmnePathway","xcell.Epithelial","xcell.HSC"),]
var_cat<-var_cat[var_cat$Var %in% heatmap_vars$x,]
heatmap_vars<-heatmap_vars[heatmap_vars$x %in% var_cat$Var,]
##
rownames(var_cat)<-var_cat$Var
var_cat<-var_cat[heatmap_vars$x,]

unique(var_cat$Class)
##
cats<-c("CancerPathway","xcell.Lymphoid","xcell.Stroma", "ImmuneResponse", "xcell.Myeloid","MCP" , "EstimateScore","xcell.estimate")
unique(rt$cohort)
#show_col(mycol)
col = list(LMiScore=colorRamp2(c(-2,0,8), c("navyblue", "white", "red")),
  group_cluster= c('C4'="#E7298A",'C3'="#7570B3",'C2' ="#66A61E","C5"="#66A61E","C6"="#FF7F00","C1"="#A65628"),
  biopsy_tissue= c("bone_marrow" = "#7570B3", "brain" = "#66A61E", "liver" = "#E7298A", "lung" = "#66A61E", "skin" ="#A65628"),
  biopsy_tissue_bio=c("nonliver" = "White", "liver" = "black"),
  cohort= c("OV"="#D95F02","BRCA"= "#7570B3","ESCA"="#E7298A", "SARC"="#66A61E", "SECR"="#E6AB02",
            "KDNY"="#A6761D", "CHOL"="#666666","PAAD"="#377EB8", "PRAD"="#4DAF4A", "STAD"="#FFFFB3",
            "COLO"= "#BEBADA", "ACC"="#FB8072","HNSC"="#80B1D3", "THYM"="#FDB462", "TGCT"="#B3DE69", "BLCA"= "#FCCDE5" )
)
############
right_annotation = rowAnnotation(foo2 = 18:1, bar2 = anno_barplot(runif(18)))
top_annotation = HeatmapAnnotation(foo1 = nrow(rt):1, bar1 = anno_points(rt$LMiScore))
ha = rowAnnotation(foo = anno_horizon(rt[,c(25,50)], negative_from_top = TRUE))
ha4 = rowAnnotation(foo = anno_horizon(apply(rt[,var_cat$Var],2,as.numeric),
                                      gp = gpar(pos_fill = "orange", neg_fill = "darkgreen")),annotation_width =unit(1, "cm"))

type1<-rt[,c(3,8,13,14)]
ha1= HeatmapAnnotation(df=data.frame(type1),
                       col=col,
                       #height = unit(30, "mm"),
                       gap=unit(0.4, "mm"),
                       #bar1 = anno_points(rt$LMiScore),
                       #foo = anno_horizon( negative_from_top = TRUE),
                       #bar2 = anno_barplot(rt$LMiScore),
                       #gp = gpar(fontsize=8),
                       #annotation_width =unit(4, "cm"),
                       annotation_name_side ="left"
                       #annotation_height = unit(3, "cm")
                       #width = 1
)
type2<-rt[,c(33,35,36,37,41,43,44,29,21,24)]
ha2= HeatmapAnnotation(df=data.frame(type2),col=col,
                      #gp = gpar(fontsize=8),
                      gap = unit(0, "points"),
                      annotation_name_side ="left")


#cats<-unique(var_cat$Class) %>% data.frame()
#cats<-cats %>%
 # dplyr::arrange(-dplyr::row_number())## change the arrange

matrix<-c(cats)
ha33 = rowAnnotation(foo = anno_empty(border = TRUE, 
                                      width = max_text_width(unlist(matrix)) + unit(0, "mm")))
###establish the heatmap
library(circlize)
heatcluster<-t(scale(apply(rt[,var_cat$Var],2,as.numeric)))####gao
heatcluster[heatcluster>1]=1
heatcluster[heatcluster< -1]=-1
p1<-Heatmap(heatcluster,
            name = "Z-score",
            cluster_columns = F,
            cluster_rows =T,
            #jitter = TRUE,
            #column_order=c("1high","2median","3low"),
            row_split = var_cat$Class,
            column_split = rt$group_cluster,
            #row_labels = splits,
            gap = unit(1, "mm"),
            #column_split =2,
            #row_split =4,
            #cluster_rows = color_branches(row_dend, k = 2),
            #cluster_columns = color_branches(col_dend, k = 2),
            col = colorRamp2(c(-1,0,1), c("navyblue", "white", "red")),
            show_heatmap_legend = TRUE,
            row_names_gp = gpar(fontsize =7),
            #column_names_gp = gpar(fontsize = 8),
            row_names_side = "left",
            #row_names_rot = 30,
            #row_names_max_width = unit(3, "cm"),
            show_column_names = FALSE,
            show_row_names = TRUE,
            width = unit(120, "mm"),
            #height = unit(140, "mm"),
            top_annotation = ha1,
            #bottom_annotation=ha,
            #right_annotation = ha4,
            left_annotation=ha33,
            row_dend_side = "right",
            show_row_dend=FALSE,
            show_column_dend = FALSE,
            #name = "ht2",
            border = TRUE,
            na_col = "grey",
            column_title_gp = gpar(fill = "WHITE", col = "BLACK", border = "WHITE"),
            row_title_gp = gpar(fill = "WHITE", col = "BLACK", border = "WHITE"),
            #row_title = "TME cells and HiF genes",
            #column_title = "Combinded TME with HIF "
)

pdf(file.path(LiverFeatures.path ,"Immune_liverfeatures_heatmap.pdf"),width=10.68,height = 8.16)

draw(p1, newpage = TRUE, 
     #column_title = "Representative Genes and cells", 
     column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
     heatmap_legend_side = "right")


for(i in 1:8) {
  decorate_annotation("foo", slice = i, {
    grid.rect(x =0, width = unit(1, "mm"), gp = gpar(fill = i, col = NA), just = "left")
    grid.text(paste(matrix[[i]], collapse = "\n"), x = unit(1, "mm"), just = "left")
  })
}
##
dev.off()
#######################################
##

