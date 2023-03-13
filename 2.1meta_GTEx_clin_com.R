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
library(survminer)
library(survival)
library(readxl)
library(tibble)
library(ggrepel)
library(forcats)
library(reshape2)
library(cowplot)
library(data.table)
library(tidyr)

##################################
mycol<-c('#D95F02','#7570B3','#E7298A', '#66A61E','#E6AB02','#A6761D','#666666','#377EB8',
          '#4DAF4A', '#984EA3','#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999','#8DD3C7',
           '#FFFFB3','#BEBADA', '#FB8072','#80B1D3','#FDB462','#B3DE69', '#FCCDE5',
          '#D9D9D9','#BC80BD','#CCEBC5', '#FFED6F','#66C2A5', '#FC8D62','#8DA0CB','#E78AC3',
          '#A6D854','#FFD92F', '#E5C494', '#B3B3B3','#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C',
          '#1B9E77','#7FC97F','#BEAED4','#FDC086','#FFFF99','#386CB0',
          '#FB9A99', '#E31A1C','#FDBF6F', '#FF7F00', '#CAB2D6', '#6A3D9A', '#FFFF99', '#B15928')
show_col(mycol)
##

results<-"Results"
Piepath    <-"Results/Nature_Piepath"
if (!file.exists(results)) { dir.create(results) }
if (!file.exists(Piepath)) { dir.create(Piepath) }

###
library(ggplot2)
piecol<-c('#66A61E','#E6AB02','#E31A1C','#A6761D','#666666')
###
load("new_gsva2.Rdata")
###

names(phe_list)
phe<-phe_list[[1]]
######
prop.table(table(phe$biopsy_tissue))
pro<-cbind(prop.table(table(phe$biopsy_tissue)) %>% melt(),data.frame(table(phe$biopsy_tissue))[2])
#colnames(pro)
colnames(pro)<-c('tissue',"perc","Freq")
pro$lable<-paste0(pro$tissue,",",pro$Freq,"(",round(pro$perc*100),"%",")")
ggplot(pro, aes(x = "", y = perc, fill = tissue)) +
  geom_col(color = "black") +
  geom_text_repel(aes(label =   lable), 
             #color = mycol[1:5],
             position = position_stack(vjust = 0.5),
             show.legend = FALSE) +
  scale_fill_manual(values =piecol)+
  guides(fill = guide_legend(title = "Metastasis site")) +
  #scale_fill_viridis_d("D") +
  coord_polar(theta = "y") + 
  theme_void()
ggsave(file.path(Piepath,"1.total_metastasis site.pdf"),
       width = 6, height =6,
       bg = "transparent", dpi = 300, useDingbats=FALSE)
###
dev.off()



for(site in unique(phe$biopsy_tissue) ) {
  print(site)
  subphe<-phe[phe$biopsy_tissue==site,]
  subpro<-cbind(prop.table(table(subphe$tissue)) %>% melt(),data.frame(table(subphe$tissue))[2])
  #subpro<-prop.table(table(subphe$tissue)) %>% melt()
  colnames(subpro)<-c('tissue',"perc","Freq")
  subpro$lable<-paste0(subpro$tissue,",",subpro$Freq,"(",round(subpro$perc*100),"%",")")
  ggplot(subpro, aes(x = "", y = perc, fill = tissue)) +
    geom_col(color = "black") +
    geom_text_repel(aes(label =   lable), 
               #color = mycol[1:5],
               position = position_stack(vjust = 0.5),
               show.legend = FALSE) +
    scale_fill_manual(values =mycol)+
    guides(fill = guide_legend(title = "Metastasis site")) +
    #scale_fill_viridis_d("D") +
    coord_polar(theta = "y") + 
    theme_void()
  ggsave(file.path(Piepath,paste0(site,"metastasis site.pdf")),
         width = 6, height =6,
         bg = "transparent", dpi = 300, useDingbats=FALSE)
}

######
phe<-phe_list[[2]]
#aggregate(x =phe$type, by = list(phe$type), FUN = "count")
#phe %>% count(gender, wt = runs)
colnames(phe)
pro<-cbind(data.frame(prop.table(table(phe$type))),data.frame(table(phe$type))[2])
colnames(pro)<-c('tissue',"perc","Freq")
pro$lable<-paste0(pro$tissue,",",pro$Freq,"(",round(pro$perc*100),"%",")")

ggplot(pro, aes(x = "", y = perc, fill =tissue)) +
  geom_col(color = "black") +
  geom_text_repel(aes(label =   lable), 
                  #color = mycol[1:5],
                  position = position_stack(vjust = 0.5),
                  show.legend = FALSE) +
  scale_fill_manual(values =piecol)+
  guides(fill = guide_legend(title = "Normal site")) +
  #scale_fill_viridis_d("D") +
  coord_polar(theta = "y") + 
  theme_void()
ggsave(file.path(Piepath,"2.total_Normal_site.pdf"),
       width = 6, height =6,
       bg = "transparent", dpi = 300, useDingbats=FALSE)

#####################################clinical data outcome...


