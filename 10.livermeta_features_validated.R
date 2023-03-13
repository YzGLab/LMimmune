#
set<-"validated"
#validate the features fro liver
 val_exp<-read.delim2("/home/data/gaoyuzhen/Projects/MET500Projects/mets-immunecluster-master/data/mets_expression_combat.txt")
 val_clin<-read.delim2("/home/data/gaoyuzhen/Projects/MET500Projects/mets-immunecluster-master/data/mets_infoClin_cluster.txt",sep = " ") 
###
 colnames(val_clin)
 ###
 m2<-val_exp
 m3<-apply(m2,2,as.numeric)
 rownames(m3)<-rownames(m2)
 m2<-m3;m3<-c()
 ##
load(file.path(Cluster.path,"DEGs_pathway.Rdata"))
     ###
     gmt_immune<-read.gmt("Signature/c7.all.v7.5.1.symbols.gmt")
     gmt_immune<-split(gmt_immune$gene,gmt_immune$term)
     
     DEGs_pathwaylist<-gmt_immune[DEGs_pathway]
 
     library(clusterProfiler)
 
      library(GSVA)
 norm.expMat<-as.matrix(m2)
 gsva_es <- gsva(norm.expMat,DEGs_pathwaylist,method="ssgsea",abs.ranking=F,min.sz = 2,ssgsea.norm=TRUE,kcdf="Poisson",parallel.sz=70)###Array data,ssgsea.norm=TRUE
# gsva_es<-t(gsva_es) %>% data.frame()

 ### to cluster in method 
 ##cluster for it 
 dist.e=dist(t(gsva_es),method='euclidean')
 hclust(dist.e, method = "complete", members = NULL)
 tree <- hclust(dist.e, method = "ward.D2")
 plot(tree,labels = FALSE, hang = -3, main = "MetaSubtype")
 #plot(tree, hang = -5)
 # save figures
 pdf(file.path(Cluster.path,paste0(set,"hcluster_fig_vali_JITC.pdf")))
 plot(tree,labels = FALSE, hang = -3, main = "MetaSubtype")
 dev.off()
################################################
 #cutree 
 group<-cutree(tree, k =6)
 #tree<-cutree(tree,h=0.7)
 ####
 #group<-factor(group,levels = c("Low","Medium","High"))
 group_cluster<-group
 group_cluster<-data.frame(group_cluster)
 group_cluster <-cbind(GEO_ID=rownames(group_cluster),group_cluster)
 ##合并临床数据
 colnames(val_clin)
 val_clin<-merge(group_cluster,val_clin,by="GEO_ID")
 #val_clin$group_cluster<-factor(val_clin$group_cluster,levels=c("3","2","1","4"))
 prop.table(table(val_clin$group_cluster,val_clin$MET_SITE),1)
 table(val_clin$group_cluster)
 
 ########################
 val_clin$group_cluster<-factor(val_clin$group_cluster,levels=unique(val_clin$group_cluster))
 a<-as.numeric(order(prop.table(table(val_clin$group_cluster,val_clin$MET_SITE),1)[,3]))
 val_clin$group_cluster<-factor(val_clin$group_cluster,levels=a)
 ##################
 write.csv(val_clin,file=file.path(Cluster.path,"val_clin_cluster_info.csv"))
 save(val_clin,file=file.path(Cluster.path,"val_clin_cluster_info.Rdata"))
 
 ##########################
 pro<-prop.table(table(val_clin$group_cluster,val_clin$MET_SITE),1) %>% melt()
 colnames(pro)<-c("immuneclass","Tissue","precent")
 
 pro$immuneclass<-factor(pro$immuneclass,levels=a)
 write.csv(pro,file=file.path(Cluster.path,"val_clin_cluster_info_prob.csv"))
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
 
 
 
 ##
 library(EPIC)
 library(immunedeconv)
 library(dplyr)
 library(ggplot2)
 library(tidyr)
 #library(tibble)
 library(MCPcounter)
 library(xCell)

 ########################################
 genes <- data.table::fread("TMEpackages/MCPcounter-master/Signatures/genes.txt",data.table = F)
 probesets <- data.table::fread("TMEpackages/MCPcounter-master/Signatures/probesets.txt",data.table = F,header = F)
 ### 成功！！！
 colnames(val_exp)
 exprMat<-val_exp
 results_MCPcounter<- MCPcounter.estimate(exprMat,
                                          featuresType= "HUGO_symbols",
                                          probesets=probesets,
                                          genes=genes)
 results_MCPcounter<-data.frame(results_MCPcounter)
 colnames(results_MCPcounter)
 results_MCPcounter<-t(results_MCPcounter)
 #xcell
 res_xcell <- xCellAnalysis(exprMat)
 res_xcell<-t(res_xcell)
 ##quantiseq
 exprMat<-data.frame(exprMat)
 res_quantiseq = deconvolute(exprMat, "quantiseq", tumor = TRUE)
 res_quantiseq<-t(res_quantiseq)
 colnames(res_quantiseq)<-res_quantiseq[1,]
 res_quantiseq<-res_quantiseq[-1,]
 
 #############
 library(estimate)
 estimateScores<-list()
 for ( i in c("met500_2")){
   if (!dir.exists("split_sampl_Single")){
     dir.create("./split_sampl_Single")
   }
   print(i)
   mydata<-val_exp
   #sample_split_list[[i]]<-mydata
   write.table(mydata,paste0("./split_sampl_Single/", i, ".txt"),sep="\t",quote=F,col.names=T)
 }
 for ( i in  c("met500_2")){
   print(i)
   in.file<-paste0("./split_sampl_Single/", i, ".txt")
   out.file<-paste0("./split_sampl_Single/", i,"outfile.gct")
   out.file2<-paste0("./split_sampl_Single/", i,"out.file_score.gct")
   filterCommonGenes(input.f= in.file , output.f=out.file, id="GeneSymbol")
   estimateScore(out.file, out.file2)
   #estimateScore[1:2,1:4]
   estimateScore<-read.table(out.file2, skip = 2, header = TRUE, sep = "\t")
   estimateScore<-data.frame(t(estimateScore))
   colnames(estimateScore)<- c('StromalScore','ImmuneScore', 'ESTIMATEScore',"TumorPurity")
   estimateScore<-estimateScore[-c(1,2),]
   #estimateScores[[i]]<-estimateScore
   message(paste0(i,"calculation of Immunescore is done!"))
 }
 ##
 library(ImmuneSubtypeClassifier)
 calls <- callEnsemble(exprMat, geneids='symbol')
 immunefeaturs<-list()
 immunefeaturs[["MCP"]]<-results_MCPcounter
 immunefeaturs[["xcell"]]<-res_xcell
 immunefeaturs[["quantiseq"]]<-res_quantiseq
 immunefeaturs[["estimate"]]<-estimateScore
 immunefeaturs[["ImmuneSubtype"]]<-calls[,-1]
 immunefeaturs<-do.call(cbind,immunefeaturs)
 colnames(immunefeaturs)
 
 immunefeaturs<-cbind(immunefeaturs,IPSs)
 colnames(immunefeaturs)
 ###
 library(readxl)### HCC genesignatures
 Hep_score <- read_excel("~/Projects/LiverCancer/1-s2.0-S0092867420316135-mmc3.xlsx")
 Hep_score<-t(Hep_score)
 Hep_score<-split(Hep_score,rownames(Hep_score))
 for( i in names(Hep_score)){
   Hep_score[[i]]<-Hep_score[[i]][Hep_score[[i]]!=" "]
   Hep_score[[i]]<-Hep_score[[i]][!is.na(Hep_score[[i]])]
 }

 library("ExperimentHub")
 library("easierData")
 score_signature_genes <- suppressMessages(easierData::get_scores_signature_genes())
 #####
 #save(score_signature_genes,file="score_signature_genes.Rdata")
 
 
 #install.packages("~/Projects/MET500Projects/Signature/IMvigor210CoreBiologies_1.0.0.tar.gz", repos = NULL, type = "source")
 #eh <- ExperimentHub()
 #query(eh, "easierData")
 #########################
 library("IMvigor210CoreBiologies")
 data(human_gene_signatures)
 ##############
 ind_genes <- human_gene_signatures
 # variables
 ind_genes<-ind_genes[-c(8,9)]
 names(ind_genes)[7]<-"Pan F TBRS"
 ##
 #################################
 Combined_signatures <- read_excel("Signature/Combined_signatures.xlsx")
 Combined_signatures<-split(Combined_signatures$GeneSymbol,Combined_signatures$Name)
 #Combined_signatures<-mapply(c, Combined_signatures,ind_genes, SIMPLIFY=FALSE)
 ind_genes<-appendList(Combined_signatures,ind_genes)
 
 #ind_genes<-appendList(Pathways,ind_genes)
 
 goi <- names(ind_genes)[!names(ind_genes) %in% c("IgG_19272155","IL12_score_21050467", "IL13_score_21050467" ,
                                                  "IL2_score_21050467","IL4_score_21050467","IL8_21978456","PD1_data","PD1_PDL1_score","PDL1_data")]
 
 

 pd<-data.frame(colnames(m2))
 # calculate gene set scores
 for (sig in goi) {
   pd[,sig] <- NA
   genes <- ind_genes[[sig]]
   genes <- genes[genes %in% rownames(m2)]
   tmp <- m2[genes, , drop=FALSE]
   pd[, sig] <- gsScore(tmp)
 }
 #install.packages("ropenblas")
 #ropenblas::rcompiler()
 colnames(pd)
 ######合并两个list 方法
 source("Signature/appendList.R")
 appendList <- function (x, val) 
 {
   stopifnot(is.list(x), is.list(val))
   xnames <- names(x)
   for (v in names(val)) {
     x[[v]] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]])) 
       appendList(x[[v]], val[[v]])
     else c(x[[v]], val[[v]])
   }
   x
 }

 immunefeaturs<-cbind(immunefeaturs,pd[,-1])
 ########################
 load("~/Projects/MET500Projects/mets-immunecluster-master/data/pathways_for_gsva_mets.Rdata")
 #Antigen processing machinery
 library(clusterProfiler)
 library(GSVA)
 norm.expMat<-as.matrix(m2)
 gsva_es <- gsva(norm.expMat,Pathways,method="ssgsea",abs.ranking=F,min.sz = 2,ssgsea.norm=TRUE,kcdf="Poisson",parallel.sz=70)###Array data,ssgsea.norm=TRUE
 gsva_es<-t(gsva_es) %>% data.frame()
 colnames(immunefeaturs)
 immunefeaturs<-cbind(immunefeaturs,gsva_es)
 #####
 save(immunefeaturs,file="val_immunefeaturs.Rdata")
 
 colnames(immunefeaturs)
 immunefeaturs<-cbind(GEO_ID=rownames(immunefeaturs),immunefeaturs)
  colnames( val_clin)
  colnames(val_clin)
  val_clin$MET_SITE_bio<-ifelse(val_clin$MET_SITE=="liver","liver","nonliver")
 val_clin<-merge( val_clin,immunefeaturs,by="GEO_ID")
####
 library(plyr)
 library(survival)
 library(ComplexHeatmap)
 rt<-val_clin
 #rt<-phe_total1[phe_total1$group_cluster %in% c("C4","C3"),]
 ########
 k.testfuc<-function(x) {
   rt[,x]<-as.numeric(rt[,x])
   if(mean(as.numeric(rt[,x]))!=0 & !is.na(mean(as.numeric(rt[,x]))) ){
     k.test<-kruskal.test(as.numeric(rt[,x])~rt[,"group_cluster"])
     data.frame(k.test$p.value,k.test$statistic,x)
   }
 }
 #k.test<-kruskal.test(as.numeric(rt[,"LMiScore"])~rt[,"group_cluster"])
 ##
 heatmap_var<- colnames(rt)[c(13:209)]
 ##
 heatmap_vars<-lapply(heatmap_var,k.testfuc) %>% ldply(.,cbind)###################
 #heatmap_vars<-lapply(heatmap_vars,cbind)
 #write.csv(heatmap_vars,file=file.path(Score.path,"phe_total_cat3.csv"))
 heatmap_vars<-heatmap_vars[heatmap_vars$k.test.p.value<0.001,]
 ##
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


 #show_col(mycol)
 col = list(
   group_cluster= c('1'="#E31A1C",'2' ="#33A02C","6"=mycol[5],"3"=mycol[8],"4"=mycol[9],'5'="grey"),
   #tissue= c("Bone"="#FF7F00","Brain"="#A6CEE3","Liver"="#E31A1C", "Lung"="#33A02C" ),
   MET_SITE_bio=c("liver"="black","nonliver"="white"),
   tumor= c("Breast"="#A6CEE3","Colorectal"="#E31A1C","Kidney"=mycol[4], "Lung"="#33A02C" ,"Prostate"=mycol[5],"Skin"=mycol[9] )
 )
 ############

 ha4 = rowAnnotation(foo = anno_horizon(apply(rt[,var_cat$Var],2,as.numeric),
                                        gp = gpar(pos_fill = "orange", neg_fill = "darkgreen")),annotation_width =unit(1, "cm"))
 
 colnames(rt)
 type1<-rt[,c(2,3,4,6,9,210)]
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
 
 pdf(file.path(LiverFeatures.path ,"Immune_liverfeatures_heatmap_val.pdf"),width=10.68,height = 8.16)
 
 draw(p1, newpage = TRUE, 
      #column_title = "Representative Genes and cells", 
      column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
      heatmap_legend_side = "right")
 
 
 for(i in 1:7) {
   decorate_annotation("foo", slice = i, {
     grid.rect(x =0, width = unit(1, "mm"), gp = gpar(fill = i, col = NA), just = "left")
     grid.text(paste(matrix[[i]], collapse = "\n"), x = unit(1, "mm"), just = "left")
   })
 }
 ##
 dev.off()
 