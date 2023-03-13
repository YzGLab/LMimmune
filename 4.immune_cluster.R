library(GSVA)
library(limma)
library(clusterProfiler)
library(sva)
library(dplyr)
library(ggplot2)
library(clustree)
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
library(GEOquery)
library(ggplot2)
library(umap)
library(hypeR)
library(msigdbr)
library(tidyverse)
library(Seurat)
library(clustree)
library(survival)
library(survminer)
source("Rcode_meta/twoclasslimma.R")

#############################
#load("new_gsva.Rdata")
load("new_gsva2.Rdata")
load("met500_GTEx.expSet.Rdata")
#load("Meta_liver_met500_20220610.Rdata")

Cluster.path<-"Results/cluster_DEG"
if (!file.exists(Cluster.path)) { dir.create(Cluster.path) }

mycol<-c( '#D95F02','#7570B3','#E7298A', '#66A61E','#E6AB02','#A6761D','#666666',
           '#377EB8','#4DAF4A', '#984EA3','#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999',
           '#8DD3C7', '#FFFFB3','#BEBADA', '#FB8072','#80B1D3','#FDB462','#B3DE69', '#FCCDE5',
           '#D9D9D9','#BC80BD','#CCEBC5', '#FFED6F','#66C2A5', '#FC8D62','#8DA0CB','#E78AC3',
           '#A6D854','#FFD92F', '#E5C494', '#B3B3B3','#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C',
           '#1B9E77','#7FC97F','#BEAED4','#FDC086','#FFFF99','#386CB0',
           '#FB9A99', '#E31A1C','#FDBF6F', '#FF7F00', '#CAB2D6', '#6A3D9A', '#FFFF99', '#B15928')
show_col()
################################
names(gsva_list)
#save(phe_list,gsva_list,file="new_gsva.Rdata")
#########
#######计算差异！
deglist<-list()
for (i in names(gsva_list)[c(1,2)]){
  print(i)
  rt<-phe_list[[i]]
  exprSet<-gsva_list[[i]]
  exprSet<-exprSet[,colnames(exprSet) %in% rownames(rt)]
  rt$groups<-ifelse(rt$biopsy_tissue_bio=="liver","yes",'no')
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
                res.path = Cluster.path ) # path for result
  tmp1 <- read.table(file.path(Cluster.path,paste0(i,"_limma_test_result.yes_vs_no.txt")),
                     sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
  diffSig <- tmp1[which(tmp1$t > 6 & tmp1$padj < 0.001),]
  deglist[[i]]<-diffSig
}####
##

DEGs_pathway1<-deglist[[2]][!rownames(deglist[[2]]) %in% rownames(deglist[[1]]),]#
##
DEGs_pathway1

write.csv(DEGs_pathway1,file=file.path(Cluster.path,"DEGs_pathway_sig.csv"))

colnames(DEGs_pathway1)
ortdf<-DEGs_pathway1
ortdf<-ortdf[order(ortdf$t,decreasing = F),]
ortdf$ID<-rownames(ortdf)
ortdf$score<-ortdf$t
library(ggplot2)
library(data.table)
library(ggpubr)
library(ggsci)

make_plot <- function(data,x,y){
  ggplot(data) + 
    geom_col(aes(x={{x}},y={{y}}, fill = {{x}} > 0),
             size = 0.25, color = "white")+
    geom_point(aes(x={{x}},y={{y}},
                   color=ifelse({{x}} > 0,"#BA7A70","#829BAB")),size=4.1)+
    geom_text(aes(x = ifelse({{x}} > 0, -.005, .005),y = {{y}}, 
                  label = gene_2,
                  color=ifelse({{x}} > 0,"#BA7A70","#829BAB"),
                  hjust = ifelse({{x}} > 0, 1, 0)),size = 3.8)+
    geom_vline(xintercept=0,size=1,color="grey40")+
    scale_y_discrete(expand = c(.025,.025))+
    scale_fill_manual(values = c("TRUE" = "#BA7A70","FALSE" = "#829BAB"))+
    scale_color_manual(values = c("#829BAB","#BA7A70"))+
    coord_cartesian(clip = "off") +  
    theme_minimal() + 
    theme(panel.grid = element_blank(),
          plot.background = element_rect(fill="Aliceblue",color="Aliceblue"),
          axis.text.y =  element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          axis.text.x = element_text(face = "bold", size =rel(1), color = "black"))
}
make_plot(p,cor,gene_2)

ortdf %>% 
  mutate(ID=fct_relevel(ID,ortdf$ID)) %>%
  ggplot(aes(ID, score,fill=mycol[1:29])) + 
  geom_bar(stat = 'identity',width = 0.4) + 
  scale_colour_manual(values = mycol) +
  scale_fill_manual(values = mycol) +
  coord_cartesian(clip = "off") +
  coord_flip() + 
  #scale_fill_manual(values = c('palegreen3', 'snow3', 'dodgerblue4'), guide = FALSE) + 
  #画2条虚线
  #geom_hline(yintercept = c(-cutoff,cutoff), 
           #  color="white",
         #    linetype = 2, #画虚线
           #  size = 0.3) + #线的粗细
  #写label
  geom_text(data = ortdf,
           aes(x=ID, y= 0, label= paste0(" ", ID)),#bar跟坐标轴间留出间隙
           size = 3, #字的大小
            hjust = "inward" ) +  #字的对齐方式
 # geom_text(data = subset(df, score > 0),
          #  aes(x=ID, y= -0.1, label=ID, color = group),
          #  size = 3, hjust = "outward") +  

  theme_minimal() + 
  theme(panel.grid = element_blank(),
        plot.background = element_rect(fill="Aliceblue",color="Aliceblue"),
        axis.text.y =  element_blank(),
        #axis.title = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(face = "bold", size =rel(1), color = "black"))+
  xlab("") + ylab("t value of limma, liver metastasis \n versus liver metastasis")+
  theme(legend.position = "none")  #去除网格线

ggsave(file.path(Cluster.path,"ggvenn_pathway_barplot.pdf"),width = 5.4,height = 8.5)


colnames(DEGs_pathway1)

library(EnhancedVolcano)
res<-DEGs_pathway1
D1<-EnhancedVolcano(res,
                    lab = rownames(res),
                    x = 't',
                    y = 'padj',
                    #pCutoff =0.05,
                    FCcutoff = 6,
                    #xlim = c(-5.5, 5.5),
                    #ylim = c(0, -log10(10e-12)),
                    pointSize = 1.5,
                    labSize = 2.5,
                    #title = 'Differential expression(normal liver vs other tissue)',
                    # subtitle = '',
                    # caption = 'FC cutoff, 1.333; p-value cutoff, 10e-4',
                    legendPosition = "bottom",
                    legendLabSize = 14,
                    col = mycol,
                    colAlpha = 0.9,
                    drawConnectors = TRUE,
                    #hline = c(10e-8),
                    widthConnectors = 0.5)
D1
ggsave(file.path(Cluster.path,"Liver_special_pathway.pdf"),
       width = 7.5, height =6.73,
       bg = "transparent", dpi = 300, useDingbats=FALSE)
####

##
targetList<-list()
for (set in names(deglist)){
  diffSig<-deglist[[set]]
  targetList[[set]]<-rownames(diffSig)
}
dev.new()
##save it
ggvenn::ggvenn(targetList,fill_color = mycol)
ggsave(file.path(Cluster.path,"ggvenn_pathway.pdf"))
venn(targetList)

#names(targetList)
DEGs_pathway<-targetList[[2]][!targetList[[2]] %in% targetList[[1]]]#
save(DEGs_pathway,file=file.path(Cluster.path,"DEGs_pathway.Rdata"))
write.csv(DEGs_pathway,file=file.path(Cluster.path,"DEGs_pathway.csv"))
###
gmt_immune<-read.gmt("Signature/c7.all.v7.5.1.symbols.gmt")
gmt_immune<-split(gmt_immune$gene,gmt_immune$term)

DEGs_pathwaylist<-gmt_immune[DEGs_pathway]
DEGs_pathwaylists<-do.call(rbind,DEGs_pathwaylist) %>% data.frame()
write.csv(DEGs_pathwaylists,file=file.path(Cluster.path,"DEGs_pathwaylists.csv"))

################################################################

set.seed(12345)
Cluster.path<-"Results/cluster_DEG"
if (!file.exists(Cluster.path)) { dir.create(Cluster.path) }


set<-"met500"
rt<-phe_list[[set]]
#rt<-rt[rt$biopsy_tissue=="liver",]
gsva_es<-gsva_list[[set]]
gsva_es<-gsva_es[,colnames(gsva_es) %in% rownames(rt)]
gsva_es<-gsva_es[rownames(gsva_es) %in% DEGs_pathway,]
write.csv(gsva_es,file =file.path(Cluster.path,"DEGs_pathwaylistsGSVAscoreMeta_stastasis.csv"))
#clusterK = 3
#M <- cor(t(gsva_es), method = "pearson")
#用corrplot画图，便于标出各cluster的黑色方框
#?corrplot
#corrplot::corrplot(M, 
      #   method = "pie", #用颜色展示相关系数，还可以改为"circle" (default), "square", "ellipse", "number", "pie", "shade"
      #   order = "hclust", 
      #   tl.cex = 0.4,
       #  hclust.method = "ward.D2", 
       #  addrect = clusterK, #画黑色方框
      #   tl.pos = "l", #"n" means don't add textlabel
       #  col = rev(brewer.pal(n = 11, name = "RdBu"))) # 运行?brewer.pal查看更多配色方案
#gsva_es<-gsva_list[[3]]
#colnames(gsva_es)
#gsva_es<-data.frame(gsva_es)
#gsva_es<-gsva_es[,colnames(gsva_es) %in% rownames(rt)]
#cluster for Liver metastasis
### to cluster in method 
##cluster for it 
dist.e=dist(t(gsva_es),method='euclidean')
hclust(dist.e, method = "complete", members = NULL)
tree <- hclust(dist.e, method = "ward.D2")
plot(tree,labels = FALSE, hang = -3, main = "MetaSubtype")
#plot(tree, hang = -5)
# save figures
pdf(file.path(Cluster.path,paste0(set,"hcluster_fig.pdf")))
plot(tree,labels = FALSE, hang = -3, main = "MetaSubtype")
dev.off()
###############################
#cutree 
group<-cutree(tree, k =6)
#tree<-cutree(tree,h=0.7)
####
#group<-factor(group,levels = c("Low","Medium","High"))
group_cluster<-group
group_cluster<-data.frame(group_cluster)
group_cluster <-cbind(Sample_id=rownames(group_cluster),group_cluster)
##合并临床数据
phe_met<-phe_list[[set]]
phe_met<-merge(group_cluster,phe_met,by="Sample_id")
#phe_met$group_cluster<-factor(phe_met$group_cluster,levels=c("3","2","1","4"))
prop.table(table(phe_met$group_cluster,phe_met$biopsy_tissue),1)
table(phe_met$group_cluster)
#phe_met$group_cluster[phe_met$group_cluster %in% c(4,6)]<-4

########################
phe_met$group_cluster<-factor(phe_met$group_cluster,levels=unique(phe_met$group_cluster))
a<-as.numeric(order(prop.table(table(phe_met$group_cluster,phe_met$biopsy_tissue),1)[,3]))
phe_met$group_cluster<-factor(phe_met$group_cluster,levels=a)
##################
write.csv(phe_met,file=file.path(Cluster.path,"phe_met_cluster_info.csv"))
save(phe_met,file=file.path(Cluster.path,"phe_met_cluster_info.Rdata"))

##########################
pro<-prop.table(table(phe_met$group_cluster,phe_met$biopsy_tissue),1) %>% melt()
colnames(pro)<-c("immuneclass","Tissue","precent")

pro$immuneclass<-factor(pro$immuneclass,levels=a)
write.csv(pro,file=file.path(Cluster.path,"phe_met_cluster_info_prob.csv"))
###
#pro$immuneclass<-ifelse(pro$immuneclass=="3","high",ifelse(pro$immuneclass=="2","low","2Median"))
p1 <- ggplot(pro, aes(x=immuneclass,y=precent,fill=Tissue)) +
  geom_bar(stat="identity") + 
  #geom_bar(position="fill") +
  scale_y_continuous(labels = percent) +
  ylab("% sample")+
  #theme_light()+
  facet_wrap(~Tissue,scales="free_x",nrow = 5)+
  theme_classic()+
  scale_fill_manual(values = mycol)+
  theme(axis.text.x = element_text(angle =30, hjust = 1,size = 10), 
        axis.text.y = element_text(angle = 90,hjust = 0.2, size = 10),
        axis.title.y = element_text(angle = 90, size = 15)) +
  theme(legend.position = "none")
p1
####
ggsave(file.path(Cluster.path,paste0(set,"immunesuppression_score_Liver_immuneclass.pdf")),
       width = 6, height =10.5,
       bg = "transparent", dpi = 300, useDingbats=FALSE)
#####
p2 <- ggplot(pro, aes(x=Tissue,y=precent,fill=Tissue)) +
  geom_bar(stat="identity") + 
  #geom_bar(position="fill") +
  scale_y_continuous(labels = percent) +
  #theme_light()+
  ylab("% sample")+
  facet_wrap(~immuneclass,scales="free_x",nrow = 1)+
  #theme_classic()+
  scale_fill_manual(values = c(mycol))+
  theme_light()+
  theme(axis.text.x = element_text(angle =30, hjust = 1,size = 10), 
        axis.text.y = element_text(angle = 90,hjust = 0.2, size = 10),
        axis.title.y = element_text(angle = 90, size = 15)) +
  theme(legend.position = "none")
p2 

ggsave(file.path(Cluster.path,paste0(set,"_Liver_specail_high.pdf")),
       width = 11, height =4.73,
       bg = "transparent", dpi = 300, useDingbats=FALSE)
#######经过鉴定，
#phe_met$group_cluster<-factor(phe_met$group_cluster,levels=c("3","2","1","4"))
save(phe_met,file=file.path(Cluster.path,"phe_met.Rdata"))


#######contruction of the heatmap data
colnames(gsva_es)
rt<-t(gsva_es)
rt<-cbind(Sample_id=rownames(rt),rt)
phe_met<-merge(phe_met,rt,by="Sample_id")
phe_met$group_cluster<-as.factor(phe_met$group_cluster)
######
mycol<-c( '#D95F02','#7570B3','#E7298A', '#66A61E','#E6AB02')
ggplot(phe_met, aes(x=tissue,fill=group_cluster)) +
  #geom_bar(stat="identity",position="fill") + 
  geom_bar(position="fill") +
  scale_y_continuous(labels = percent) +
  #theme_light()+
  ylab("% sample")+
  #facet_wrap(~tissue,scales="free_x")+
  theme_classic()+
  scale_fill_manual(values = c('#D95F02','#7570B3', '#66A61E','#E6AB02','#A6761D','#E7298A'))+
  theme(axis.text.x = element_text(angle =30, hjust = 1,size = 10), 
        axis.text.y = element_text(angle = 90, size = 12),
        axis.title.y = element_text(angle = 90, size = 15)) +
  geom_text(aes( y=..count../tapply(..count.., ..x.. ,sum)[..x..], 
                 label=paste(..count..,paste("(",percent(..count../tapply(..count.., ..x.. ,sum)[..x..]),")",sep=""),sep=" ")),
            stat="count",color="white", position=position_fill(0.5), vjust=0.5,size=3) +
  theme(legend.position = "right")



ggsave(file.path(Cluster.path,"group_cluster_tissue_site_count_percent.pdf"),
       width = 14.38, height =6.85,
       bg = "transparent", dpi = 300, useDingbats=FALSE)

ggplot(phe_met) +
  geom_bar( aes(x=tissue,fill=group_cluster)) + 
  #geom_bar(position="fill") +
  #scale_y_continuous(labels = percent) +
  #theme_light()+
  ylab("% sample")+
  #facet_wrap(~tissue,scales="free_x")+
  theme_classic()+
  scale_fill_manual(values = c('#D95F02','#7570B3', '#66A61E','#E6AB02','#A6761D','#E7298A'))+
  theme(axis.text.x = element_text(angle =30, hjust = 1,size = 10), 
        axis.text.y = element_text(angle = 90, size = 12),
        axis.title.y = element_text(angle = 90, size = 15)) +
  theme(legend.position =  "right")
ggsave(file.path(Cluster.path,"group_cluster_tissue_site_count.pdf"),
       width = 14.38, height =6.85,
       bg = "transparent", dpi = 300, useDingbats=FALSE)


##

#"#FF7F00","#A6CEE3","#E31A1C","#33A02C","#999999"
col = list(
  group_cluster= c('4'="#E31A1C",
                   '6'=mycol[6],
                   '5'=mycol[3],
                   '2'=mycol[2],
                   '1'=mycol[4],
                   "3"="grey"),
  biopsy_tissue_bio= c("nonliver"="white","liver"="black")
  #tissue= c("Bone"="#FF7F00","Brain"="#A6CEE3", "Liver"="#E31A1C", "Lung"="#33A02C"),
  #tumor= c("Breast"="#A6CEE3", "Colorectal"="#E31A1C", "Kidney"=mycol[4], "Lung"="#33A02C" ,"Prostate"=mycol[5],"Skin"=mycol[9] )
)
############

type<-phe_met[,c(2,6,12,13)]
ha1= HeatmapAnnotation(df=data.frame(type),
                       col=col,
                       #height = unit(30, "mm"),
                       gap=unit(0.4, "mm"),
                       #gp = gpar(fontsize=8),
                       #annotation_width =unit(4, "cm"),
                       annotation_name_side ="left"
                       #annotation_height = unit(3, "cm")
                       #width = 1
)
colnames(phe_met)
plotdata<-phe_met[,c(14:42)]
plotdata<-apply(plotdata,2,as.numeric)
rownames(plotdata)<-phe_met[,1]
heatcluster<-t(scale(plotdata))##
heatcluster[heatcluster >2]  <- 2 # 
heatcluster[heatcluster < -2] <- -2 # 
#display.brewer.all()
p3<-Heatmap(heatcluster,
            name = "Z-score",
            cluster_columns = F,
            cluster_rows =T,
            #split = 4,
            #row_km = 2,
            #jitter = TRUE,
            #column_order=c("1","2","3","4"),
            #row_split = splits,
            column_split = phe_met$group_cluster,
            #row_labels = splits,
            gap = unit(1, "mm"),
            show_heatmap_legend = TRUE,
            row_names_gp = gpar(fontsize =9),
            #column_names_gp = gpar(fontsize = 8),
            col=colorRampPalette(c('blue','white','#E7298A'))(100), 
            #col=brewer.pal(n=9, name="YlGnBu"), 
            row_names_side = "left",
            #row_names_rot = 30,
            #row_names_max_width = unit(3, "cm"),
            show_column_names = FALSE,
            show_row_names = FALSE,
            #width = unit(120, "mm"),
            #height = unit(140, "mm"),
            top_annotation = ha1,
            #bottom_annotation=ha,
            #right_annotation = ha33,
            #left_annotation=ha33,
            row_dend_side = "left",
            show_row_dend=TRUE,
            show_column_dend = FALSE,
            #name = "ht2",
            border = TRUE,
            na_col = "grey",
            column_title_gp = gpar(fill = "WHITE", col = "BLACK", border = "WHITE"),
            row_title_gp = gpar(fill = "WHITE", col = "BLACK", border = "WHITE"),
            row_title = "Liver metastasis special Immune signatures",
            #column_title = "Combinded TME with HIF "
)
p3

pdf(file.path(Cluster.path,paste0(set,"Immunosuppression_heatmap.pdf")),width=10.68,height = 6.16)
p3
dev.off()

save.image(file="Meta_liver_met500_20220610.Rdata")
###

table(phe_met$group_cluster,phe_met$cohort)
table(phe_met$group_cluster,phe_met$tissue)

table(phe_met$group_cluster,phe_met$biopsy_tissue_bio)
table1<-cbind(table(phe_met$group_cluster,phe_met$biopsy_tissue),
              table(phe_met$group_cluster,phe_met$biopsy_tissue_bio),
              table(phe_met$group_cluster,phe_met$tissue)) %>% t()
table2<-cbind(prop.table(table(phe_met$group_cluster,phe_met$biopsy_tissue),2),
              prop.table(table(phe_met$group_cluster,phe_met$biopsy_tissue_bio),2),
              prop.table(table(phe_met$group_cluster,phe_met$tissue),2)) %>% t()
tables<-cbind(table1,table2)
write.csv(tables,file=file.path(Cluster.path,"tables_phe_met_cluster_info.csv"))