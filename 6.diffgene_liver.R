##
source("Rcode_meta/twoclasslimma.R")
library(limma)

load("met500_GTEx.expSet.Rdata")

DEG.path_gene <-"Results/gene-liver_DEG"
if (!file.exists(DEG.path_gene )) { dir.create(DEG.path_gene  ) }

###
names(expSet)


  rt<-phe_GTEx
  exprSet<-expSet[['panGTEX']]
  exprSet<-exprSet[,colnames(exprSet) %in% rownames(rt)]
  rt$groups<-ifelse(rt$biopsy_tissue_bio=="liver","yes",'no')
  # high risk vs others
  subt <- data.frame(condition =rt$groups,
                     row.names = colnames(exprSet))
  twoclasslimma(subtype  = subt, # subtype information (must contain a column named 'condition')
                featmat  = exprSet, # expression file (fill detect data scale automatically)
                treatVar = "yes", # name of treatment group
                ctrlVar  = "no", # name of control group
                prefix   = "GTEx", # prefix of the DE file
                overwt   = TRUE, # whether over write files that already existed
                sort.p   = TRUE, # if sorting the results by adjusted p value
                verbose  = FALSE, # if showing verbose result
                res.path = DEG.path_gene ) # path for result
  tmp1 <- read.table(file.path(DEG.path_gene,"GTEx_limma_test_result.yes_vs_no.txt"),
                     sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
  diffSig_liver <- tmp1[which(tmp1$log2fc >2 & tmp1$padj < 0.05),]
  
  write.csv(diffSig_liver,file=file.path(DEG.path_gene,"GTEX_diffSig_liver.csv"))
  save(diffSig_liver,file=file.path(DEG.path_gene,"GTEX_diffSig_liver.Rdata"))
  ##
library(EnhancedVolcano)
  ?EnhancedVolcano
  res<-tmp1
  colnames(res)
  d1<-EnhancedVolcano(res,
                  lab = rownames(res),
                  x = 'log2fc',
                  y = 'padj',
                  pCutoff =0.05,
                  FCcutoff = 2,
                  xlim = c(-5.5, 5.5),
                  ylim = c(0, -log10(10e-12)),
                  pointSize = 1.5,
                  labSize = 2.5,
                  title = 'Differential expression(normal liver vs other tissue)',
                 # subtitle = '',
                 # caption = 'FC cutoff, 1.333; p-value cutoff, 10e-4',
                  legendPosition = "bottom",
                  legendLabSize = 14,
                  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                  colAlpha = 0.9,
                  drawConnectors = TRUE,
                  hline = c(10e-8),
                  widthConnectors = 0.5)
  d1
  ggsave(file.path(DEG_cluster.path_gene,"_Liver_normal_gene.pdf"),
         width = 7.5, height =6.73,
         bg = "transparent", dpi = 300, useDingbats=FALSE)
  
  
###one method to calculate the top genes. 
  source("Rcode_meta/twoclasslimma.R")
  ######claculate the score for the subtype .
  load("immunefeaturs.Rdata")
  load("phe_met.Rdata")
  load("metExp500.Rdata")
  ##
  
  
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
  colnames(phe_met)
  phe_met$Sample_id<-gsub("-",".",phe_met$Sample_id)
  phe_total<-merge(phe_met,immunefeaturs,by="Sample_id")
  ##az  asZXaszwaq
  #phe_met$group_cluster<-factor(phe_met$group_cluster,levels=c("1liver_very_high","2liver_high","3liver_median","4liver_low"))
  #################################
  table(phe_met$group_cluster)
  
  ####
  fig.path_gene <-"gene_DEG"
  if (!file.exists(fig.path_gene )) { dir.create(fig.path_gene ) }
  
  deglist<-list()
  for(i in unique(phe_met$group_cluster)){
    print(i)
    rt<-phe_met
    #rt<-rt[rt$group_cluster!=i,]
    exprSet<-metExp
    exprSet<-exprSet[,colnames(exprSet) %in% rt$run.id]
    rt$groups<-ifelse(rt$group_cluster==i,"yes",'no')
    group_list<-rt[,'groups']
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
                  res.path = fig.path_gene ) # path for result
    tmp1 <- read.table(file.path(fig.path_gene,paste0(i,"_limma_test_result.yes_vs_no.txt")),
                       sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    diffSig <- tmp1[which(tmp1$log2fc >2 & tmp1$padj < 0.05),]
    deglist[[i]]<-diffSig
  }  
  ##
  
  
  
  #####
  ##
  #save(deglist,file="liver_special_deglist.Rdata")
  ##
  #load("liver_special_deglist.Rdata")
  
  targetList<-list()
  for (set in names(deglist)){
    diffSig<-deglist[[set]]
    targetList[[set]]<-rownames(diffSig)
  }
  library(ggvenn)
  ggvenn(targetList,fill_color = mycol)##只能四组
  #library(VennDiagram)
  venn(targetList)
  
  #
  names(targetList)
  #DEGs_genes=Reduce(intersect,targetList[2:5])
  DEG<-targetList[[4]]
  DEG<-DEG[!DEG %in% rownames(diffSig_liver)]
  #DEG<-targetList[[4]][!targetList[[4]] %in% do.call(c,targetList[c(1:3)])]
  #DEG<-targetList[[4]][!targetList[[4]] %in% do.call(c,targetList[c(1:3,5)])]
  #DEG<-targetList[[5]][!targetList[[5]] %in% do.call(c,targetList[c(1:4)])]
  save(DEG,file="LMiScore.Rdata")
  ##save it
  
  ##
  DEGs_phe<-deglist[[4]]
  DEGss<-DEGs_phe[rownames(DEGs_phe) %in% DEG,]
  ###
