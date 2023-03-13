#####
DEG_cluster.path_gene <-"Results/DEG_cluster_gene"
if (!file.exists(DEG_cluster.path_gene )) { dir.create(DEG_cluster.path_gene ) }

#Figure 4 
library(limma)
source("Rcode_meta/twoclasslimma.R")
######claculate the score for the subtype .
load("immunefeaturs.Rdata")
load(file.path(Cluster.path,"phe_met.Rdata"))
load("metExp500.Rdata")
load(file.path(DEG.path_gene,"GTEX_diffSig_liver.Rdata"))
##
################
mycol<-c( '#D95F02','#7570B3','#E7298A', '#66A61E','#E6AB02','#A6761D','#666666',
          '#377EB8','#4DAF4A', '#984EA3','#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999',
          '#8DD3C7', '#FFFFB3','#BEBADA', '#FB8072','#80B1D3','#FDB462','#B3DE69', '#FCCDE5',
          '#D9D9D9','#BC80BD','#CCEBC5', '#FFED6F','#66C2A5', '#FC8D62','#8DA0CB','#E78AC3',
          '#A6D854','#FFD92F', '#E5C494', '#B3B3B3','#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C',
          '#1B9E77','#7FC97F','#BEAED4','#FDC086','#FFFF99','#386CB0',
          '#FB9A99', '#E31A1C','#FDBF6F', '#FF7F00', '#CAB2D6', '#6A3D9A', '#FFFF99', '#B15928')
################################
###########################################
phe_met$group_cluster<-paste0("C",phe_met$group_cluster)
#immunefeaturs<-cbind(Sample_id=rownames(immunefeaturs),immunefeaturs)
phe_met$Sample_id<-gsub("-",".",phe_met$Sample_id)
phe_total<-merge(phe_met,immunefeaturs,by="Sample_id")
##az  asZXaszwaq
#################################
table(phe_met$group_cluster)
#phe_met    $group_cluster<-ifelse(phe_met$group_cluster %in% c("C1","C6","C4"),"C4",phe_met$group_cluster)

###############################################################################
 ##直接算的 C4 
  rt<-phe_met[phe_met$group_cluster %in% c("C4","C3"),]
  #rt<-rt[rt$group_cluster!=i,]
  exprSet<-metExp
  exprSet<-exprSet[,colnames(exprSet) %in% rt$run.id]
  rt$groups<-ifelse(rt$group_cluster=="C4","yes",'no')
  #####
  subt <- data.frame(condition =rt$groups,
                     row.names = colnames(exprSet))
  twoclasslimma(subtype  = subt, # subtype information (must contain a column named 'condition')
                featmat  = exprSet, # expression file (fill detect data scale automatically)
                treatVar = "yes", # name of treatment group
                ctrlVar  = "no", # name of control group
                prefix   = "Top", # prefix of the DE file
                overwt   = TRUE, # whether over write files that already existed
                sort.p   = TRUE, # if sorting the results by adjusted p value
                verbose  = FALSE, # if showing verbose result
                res.path = DEG_cluster.path_gene) # path for result
  tmp1 <- read.table(file.path(DEG_cluster.path_gene ,"Top_limma_test_result.yes_vs_no.txt"),
                     sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
  diffSig_cluster <- tmp1[which(tmp1$log2fc>2 & tmp1$padj < 0.05),]
  
  

  
 #####  ?EnhancedVolcano
  library(EnhancedVolcano)
  res<-tmp1
  colnames(res)
  d2<-EnhancedVolcano(res,
                      lab = rownames(res),
                      x = 'log2fc',
                      y = 'padj',
                      pCutoff =0.05,
                      FCcutoff = 2,
                      xlim = c(-4.5, 4.5),
                      ylim = c(0, -log10(10e-12)),
                      pointSize = 1.5,
                      labSize = 2.5,
                      title = 'Differential expression(C4 vs C3 tissue)',
                      # subtitle = '',
                      caption = 'FC cutoff, 1.333; p-value cutoff, 10e-4',
                      legendPosition = "bottom",
                      legendLabSize = 14,
                      col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                      colAlpha = 0.9,
                      drawConnectors = TRUE,
                      hline = c(10e-8),
                      widthConnectors = 0.5)
  d2
  met<-"met500"
  ggsave(file.path(DEG_cluster.path_gene,"_Liver_met_gene.pdf"),
         width = 7.5, height =6.73,
         bg = "transparent", dpi = 300, useDingbats=FALSE)
  ##

  ###########
  show_col(mycol)
  mycol
  ##venn plot for ex
  targetList<-list()
  targetList[["Liver"]]<-rownames(diffSig_liver)
  targetList[["Meta-Liver"]]<-rownames(diffSig_cluster)
  ggvenn::ggvenn(targetList,fill_color = c('#8DA0CB','#E78AC3'))
  #venn(targetList),  "#377EB8",    "#E31A1C"
  ggsave(file.path(DEG_cluster.path_gene ,"gene_ggvenn_pathway.pdf"))
  DEG<-rownames(diffSig_cluster)[!rownames(diffSig_cluster) %in% rownames(diffSig_liver)]
  
  diffsig_lm_gene<-diffSig_cluster[rownames(diffSig_cluster) %in% DEG,]
  
  ######
  save(DEG,diffsig_lm_gene,file =file.path(DEG_cluster.path_gene ,"DEGenes_special_liver_metastasis.Rdata" ))
  write.csv(diffsig_lm_gene,file=file.path(DEG_cluster.path_gene ,"DEGenes_special.csv"))####

##generated the score for LIMISCORe
  

  
###
  load(file.path(DEG_cluster.path_gene ,"DEGenes_special_liver_metastasis.Rdata" ))
  library(ggnewscale)
  library("clusterProfiler")
  options(connectionObserver = NULL)
  library("org.Hs.eg.db")
  library("enrichplot")
  library("ggplot2")
  
  genes=DEG
  entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
  gene<-entrezIDs <- as.character(entrezIDs)
  
  
  #---Go????
  Go<- enrichGO(gene = gene,
                keyType = "ENTREZID",
                OrgDb = org.Hs.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.1,
                ont="all",
                readable =T)
  write.table(Go,file.path(DEG_cluster.path_gene,file="Cluster_GO-up.txt"),sep="\t",quote=F,row.names = F)
  
  pdf(file=file.path(DEG_cluster.path_gene,"Cluster_Go-up.pdf"),width = 6,height = 8)
  dotplot(Go,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free') + scale_color_continuous(low='purple', high='green')+ aes(shape=GeneRatio > 0.1)
  dev.off()
  
  pdf(file=file.path(DEG_cluster.path_gene,"Cluster_Go-net.pdf"),width = 9,height = 8)
  Go2<- pairwise_termsim(Go)
  emapplot(Go2, showCategory = 10,color = "p.adjust")+ scale_color_continuous(low='purple', high='green')
  dev.off()
  
  pdf(file=file.path(DEG_cluster.path_gene,"Cluster_Go-circ.pdf"),width = 8,height = 5)
  cnetplot(Go2, showCategory = 10,categorySize="count",circular = FALSE, colorEdge = TRUE,cex_gene = 0.4,cex_label_gene =0.5)
  dev.off()
  
  pdf(file=file.path(DEG_cluster.path_gene,"Cluster_Go-circ.pdf"),width = 8,height = 5)
  setReadable(Go2, 'org.Hs.eg.db', 'ENTREZID') 
  test<-c(grep("gly",Go$Description,value=TRUE),grep("glu",Go$Description,value=TRUE))
  #test<-c("ATP generation from ADP","ATP metabolic process","glycolytic process","generation of precursor metabolites and energy")
  #cnetplot(Go2, showCategory = test, categorySize="count",circular = TRUE, colorEdge = TRUE)
  cnetplot(Go2, showCategory = test, categorySize="count",circular = FALSE, colorEdge = TRUE,cex_gene = 0.4,cex_label_gene =0.5)
  dev.off()
  ####
  
  
  
  df %>% 
    ggplot() + 
    geom_col(aes(x = rating_diff, y = company_location, fill = rating_diff > 0),
             size = 0.25, color = "white")+
    geom_point(aes(x = rating_diff,y = company_location,
                   color=ifelse(rating_diff > 0,"#E11B4D","#8456BA")),size=5)+
    geom_text(aes(x = ifelse(rating_diff > 0, -.005, .005),y = company_location, 
                  label = company_location,
                  color=ifelse(rating_diff > 0,"#E11B4D","#8456BA"),
                  hjust = ifelse(rating_diff > 0, 1, 0)),size = 4)+
    geom_vline(xintercept=0,size=1,color="grey40")+
    scale_x_continuous(expand = expansion(add = c(0,.2)),
                       breaks = seq(-.4,.2, by = .2)) + 
    scale_y_discrete(expand = c(.025,.025))+
    scale_fill_manual(values = c("TRUE" = "#E11B4D","FALSE" = "#8456BA"))+
    scale_color_manual(values = c("#8456BA","#E11B4D"))+
    coord_cartesian(clip = "off") +  
    theme_minimal() + 
    theme(panel.grid = element_blank(),
          plot.background = element_rect(fill="Aliceblue",color="Aliceblue"),
          axis.text.y =  element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          axis.text.x = element_text(face = "bold", size =rel(1), color = "black"))
  
  
  ###
  

  ###
  diffSig_clusters<-list()
  for ( i in unique(phe_met$group_cluster)[-3]){
    print(i)
    rt<-phe_met[phe_met$group_cluster %in% c(i,"C3"),]
    #rt<-rt[rt$group_cluster!=i,]
    exprSet<-metExp
    exprSet<-exprSet[,colnames(exprSet) %in% rt$run.id]
    rt$groups<-ifelse(rt$group_cluster==i,"yes",'no')
    #####
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
                  res.path = DEG_cluster.path_gene) # path for result
    tmp1 <- read.table(file.path(DEG_cluster.path_gene ,paste0(i,"_limma_test_result.yes_vs_no.txt")),
                       sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    diffSig_cluster <- tmp1[which(tmp1$log2fc>2 & tmp1$padj < 0.05),]
    diffSig_cluster<-diffSig_cluster[!rownames(diffSig_cluster) %in% rownames(diffSig_liver),]
    diffSig_clusters[[i]]<-diffSig_cluster
  }
  
  targetList<-list()
  for (set in names(diffSig_clusters)){
    diffSig<-diffSig_clusters[[set]]
    targetList[[set]]<-rownames(diffSig)
  }
  ##save it
  
  ggvenn::ggvenn(targetList,fill_color = mycol)
  #venn(targetList)
  DEG<-targetList[["C4"]][!targetList[["C4"]] %in% sapply(targetList[c(1:2,3,5,6)],c)]
  DEG<-targetList[["C4"]]
  #DEG<-rownames(DEG)[!rownames(DEG) %in% rownames(diffSig_liver)]
 