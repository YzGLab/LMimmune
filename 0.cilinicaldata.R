####Figure2 clinical plot for liver metastasis

#####
results<-"Results"
Clinicloutcome   <-"Results/Clinicloutcome"
if (!file.exists(results)) { dir.create(results) }
if (!file.exists(Clinicloutcome)) { dir.create(Clinicloutcome) }



##
library(survminer)
library(survival)
library(readxl)
library(ggplot2)
library(tibble)
library(scales)
library(ggrepel)
library(forcats)
library(reshape2)
library(cowplot)
library(data.table)
library(tidyr)

CellData <- read_excel("/home/data/gaoyuzhen/Projects/MET500Projects/Results/ClinicalDataForMeta/CellData_1.xlsx")
###
colnames(CellData)
#######
rt<-CellData
rt<-data.frame(rt)
rt<-rt[rt$sample_type=="Metastasis",]
#rt<-rt[,c("os_status","os_days","MetaLiver")]
#rt<-rt[!rt$os_days<30,]
  rt$os_days<-as.numeric(as.character(rt$os_days))
  table(rt$os_status)
  rt$os_status<-ifelse(rt$os_status=="alive",0,1)
  rt$os_status<-as.numeric(as.character(rt$os_status))
  #rt$os_status<-as.factor(rt$os_status)
  rt<-rt[which(rt$os_days!='NA' & rt$os_status!="NA"),]
  kmfit<- survfit(Surv(os_days,os_status)~MetaLiver,data=rt)
pdf(file.path(Clinicloutcome,"Total_OS.pdf"),width = 6.49, height = 6.58,onefile = F )
  ggsurvplot(kmfit,#surv.median.line = "hv", # Add medians survival
             # Change legends: title & labels
             main = "Survival curve",
             #legend.title = paste0(set,"-",sig),
             legend.labs = c("OnlyLiver","IncludeLiver","ExculdeLiver"),
             xlim = c(0, 2000),
             xlab="Time(days)",
             ylab="Overall Survival",
             size = 0.1,
             #fun="cumhaz",
             #fun='event',
             # Add p-value and tervals
             pval = TRUE,
             #test.for.trend = TRUE,###group more than 2groups
             break.time.by = 200,
             #conf.int = TRUE,
             #group.by=,
             # Add risk table
             risk.table = TRUE,
             tables.height = 0.185,
             tables.theme = theme_cleantable(),
             palette = c('#E6AB02',"#E31A1C", "#1F78B4"),
             #ggtheme = theme_bw(), # Change ggplot2 theme
             #font.title="OS",
             font.main =15,
             font.x =  15,
             font.y = 15,
             font.tickslab =15
             #在左下???标出pvalue、HR???95% CI
             #???小的p value标为p < 0.001
             #pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
             # paste("p = ",round(p.val,3), sep = "")), HR, CI, sep = "\n")
  )
dev.off()
  
table(sub_rt$MetaLiver)
  cox.sig<-c()
  cancernames<-unique(rt$cancer_type)
  #cancernames<-c("breast cancer","colorectal cancer","prostate cancer","renal cell carcinoma","non-small cell lung cancer","melanoma")
  for(set in cancernames){
  sub_rt<-rt[rt$cancer_type==set,]
  sub_rt$met_count.Liver<-as.numeric(sub_rt$met_count.Liver)
  sub_rt$met_count.Liver<-ifelse(sub_rt$met_count.Liver==0,0,1)
      Gcox1<-coxph(Surv(os_days,os_status)~met_count.Liver,data=sub_rt)
      GSum<-summary(Gcox1)
      HR<-round(GSum$coefficients[,2],3)
      Pvalue<-round(GSum$coefficients[,5],3)
      CI<-paste0(round(GSum$conf.int[,3:4],3),collapse='-')
      coeff<-round(GSum$coefficients[1],3)
      se<-GSum$coefficients[1,3]
      low<-round(GSum$conf.int[,3],3)
      up<-round(GSum$conf.int[,4],3)
      cox.p<-data.frame('dataset'= set,
                        'Hazard Ratio'= HR,
                        'CI95'=CI,
                        "coeff"=coeff,
                        "se"=se,
                        "low"=low,
                        "up"=up,
                        'P-value'=Pvalue,
                        "nofliver"=data.frame(table(sub_rt$MetaLiver))[2,2],
                        "prec"=prop.table(table(sub_rt$MetaLiver))[2],
                        "numberofpatients"=nrow(sub_rt))
      cox.sig=rbind(cox.sig,cox.p)
      message(set,"cox regression is Done!")
      ##plot the km curves
      kmfit<- survfit(Surv(os_days,os_status)~met_count.Liver,data=sub_rt)
      p.val<-round(Pvalue,3)
      HR <- paste("Hazard Ratio= ", round(HR,2), sep = "")
      CI <- paste("95% CI: ", paste(round(low,2), round(up,2), sep = " - "), sep = "")
      #if (!dir.exists("metastasis_COHORT")){
     #  dir.create("./metastasis_COHORT")
     # }
      #if(p.val<0.05){
        pdf(file.path(Clinicloutcome,paste0(set,"_",HR,".pdf")),width = 6.49, height = 6.58,onefile = F )
        print(ggsurvplot(kmfit,#surv.median.line = "hv", # Add medians survival
                         # Change legends: title & labels
                         main = "Survival curve",
                         #legend.title = paste0(set,"-",sig),
                         legend.labs = c("NoLiverMeta","LiverMeta"),
                         xlim = c(0, 2000),
                         xlab="Time(days)",
                         ylab="Overall Survival",
                         size = 1,
                         #fun="cumhaz",
                         #fun='event',
                         # Add p-value and tervals
                         #pval = TRUE,
                         #test.for.trend = TRUE,###group more than 2groups
                         break.time.by = 200,
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
                                                    paste("p = ",round(p.val,3), sep = "")), HR, CI, sep = "\n")) )
        invisible(dev.off())
     # }
    }
####
  write.csv(cox.sig,file =file.path(Clinicloutcome,"LiverMeta_cox.sig.csv" ))
  
  cox.sig
  cox.sig<-cox.sig[order(cox.sig$Hazard.Ratio,decreasing = T),]
  data<-cox.sig
  colnames(data)
  data$dataset<-paste0(data$dataset,ifelse(data$P.value<0.01,"**",ifelse(data$P.value<0.05,"*","")))
  data$direction<-ifelse(data$Hazard.Ratio>1,1,0)
  data$direction<-as.factor(data$direction)
  ##
  data$id<-as.numeric(as.factor(data$dataset))
  label_data <- data
  number_of_bar <- nrow(label_data)
  angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust <- ifelse( angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)
  head(label_data)


  p <- data %>% 
    mutate(dataset==fct_relevel(dataset,data$dataset)) %>%
    ggplot(aes(x=dataset, y=Hazard.Ratio, fill=direction)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
    geom_bar(stat="identity", alpha=0.5) +
    scale_fill_manual(values = mycol[c(9,3)])+
    #ylim(-100,120) +
    theme_minimal() +
    #scale_fill_brewer(palette = "Set1") +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm") 
    ) +
    coord_polar()+ 
    geom_text(data=label_data, aes(x=dataset, y=Hazard.Ratio+0.5, 
                                 label=dataset, hjust=hjust), 
  color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) 
 p
ggsave(file.path(Clinicloutcome,"prognsisof Livermetastasis.pdf"),width = 7,height = 7)
  
  
  ##plot for cox.sig
  colnames(cox.sig)
  
  
  library(metafor)

  coxmeta<-c()

  res.rma_each<- metafor::rma(yi = cox.sig[,4], sei = cox.sig[,5],method="FE")
  results<-round(c(res.rma_each$pval,exp(res.rma_each$b),exp(res.rma_each$ci.lb),exp(res.rma_each$ci.ub)),3)

  coxmeta<-rbind(cox.sig[,c(7,1,5,6)],results)

  colnames(coxmeta)<-c("pvalue","HR","95CI_L","95CI_U")

  head(coxmeta)
  library(DT)
  datatable(coxmeta, filter = "top", 
            options = list(pageLength = 5))
  datatable(coxmeta, 
            extensions = 'Buttons',
            options = list(dom = 'Bfrtip', 
                           buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))

  #####forest
  library("forestplot")
  library("magrittr")
  library("checkmate")
  library(grid)
  ##########
  rt<-cox.sig
  rownames(rt)<-rt$dataset
  colnames(rt)
  rt$HR<-round(rt$Hazard.Ratio,2)
  rt$up<-round(rt$up,2)
  rt$low<-round(rt$low,2)
  rt$pvalue<-round(rt$P.value,2)
 # rt$dataset<- fct_reorder(rt$dataset, rt$HR)
  rt<-rt[order(rt$nofliver,decreasing = T),]
  #colnames(rt)
  #rt<-rt[order(rt$Hazard.Ratio,decreasing = T),]
  ######
  tabletext <- cbind(c("\nSub_Cancer",NA,rownames(rt),NA),
                     c("Hazard Ratio\n(95% CI)",NA, 
                       paste0(format(rt$HR,nsmall=2),
                              " (",format(rt$low,nsmall = 2),"-",format(rt$up,nsmall = 2),")",sep=""),NA),
                     c("p-value",NA,rt$pvalue,NA),
                     c("\nprec",NA,round(rt$prec,2),NA),
                     c("N",NA,rt$numberofpatients,NA))
  pdf(file.path(Clinicloutcome,"forestplot_OS.pdf"),width = 8,height = 6,onefile=FALSE)
  forestplot(labeltext=tabletext, #
             mean=c(NA,NA,rt$HR,NA),#HR
             lower=c(NA,NA,rt$low,NA), #95%
             upper=c(NA,NA,rt$up,NA),#95%
             #title="Hazard Ratio",
             graph.pos=2,#�֦�
             graphwidth = unit(.3,"npc"),#
             #fn.ci_norm="fpDrawDiamondCI",#
             col=fpColors(box="steelblue", lines="black", zero = "black"),#
             #col=fpColors(box="black", lines="black", zero = "black"),#
             #boxsize=c(NA,NA,NA,rt$numberofpatients,NA)/200,#
             boxsize=0.3,
             #lwd.ci=2,ci.vertices.height = 0.1,ci.vertices=TRUE,#
             zero=1,#zero
             #xlog = TRUE,
             lwd.zero=1,#zero
             #grid = structure(c(rt[1,]$Hazard.Ratio), gp = gpar(col = "black", lty=2,lwd=2)),#()
             #xticks = c(0.5,0.75, 1,1.25,1.5),#
             #clip = c(0.1,2.5), 
             #lwd.xaxis=2,#X
             #xlab="     <-Favour Combination  Therapy       Favour Target  Therapy->",#X
             hrzl_lines=list("3" = gpar(lwd=2, col="black"),#
                             #"4" = gpar(lwd=60,lineend="butt", columns=c(1:4), col="#99999922"),#
                             "31" = gpar(lwd=2, col="black")),#"nrow(rt)+5
             txt_gp=fpTxtGp(label=gpar(cex=1.25),#
                            ticks=gpar(cex=1.25),
                            xlab=gpar(cex = 1.25),
                            title=gpar(cex = 1.25)),
             #is.summary = c(T,rep(F,27)),#
             #lineheight = unit(.75,"cm"),#
             align=c("l","c","c"),#
             #cex=10,
             colgap = unit(0.1,"cm"),#
             #mar=unit(rep(0.25, times = 4), "cm"),#
             new_page = T#
  )
    dev.off()
################################################clinical 

   mycol<-c( '#D95F02','#7570B3','#E7298A', '#66A61E','#E6AB02','#A6761D','#666666',
            '#377EB8','#4DAF4A', '#984EA3','#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999',
            '#8DD3C7', '#FFFFB3','#BEBADA', '#FB8072','#80B1D3','#FDB462','#B3DE69', '#FCCDE5',
            '#D9D9D9','#BC80BD','#CCEBC5', '#FFED6F','#66C2A5', '#FC8D62','#8DA0CB','#E78AC3',
            '#A6D854','#FFD92F', '#E5C494', '#B3B3B3','#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C',
            '#1B9E77','#7FC97F','#BEAED4','#FDC086','#FFFF99','#386CB0',
            '#FB9A99', '#E31A1C','#FDBF6F', '#FF7F00', '#CAB2D6', '#6A3D9A', '#FFFF99', '#B15928')
  ##
  
  ##figure2 clinical of metalive
rt<-CellData
rt<-data.frame(rt)
#rt<-rt[rt$sample_type=="Metastasis",]
 pie( table(rt$sample_type))
 pie( table(rt$cancer_type))
 rt<-data.frame(rt)

 #x=1:length(rs),
 colnames(rt)
 #rt$MetaLiver[rt$sample_type!="Metastasis"]="0"
 num<-data.table(table(rt$cancer_type,rt$MetaLiver))
 num<-num %>%
   pivot_wider(
     names_from =c("V2"),
     values_from = N
   )
 
 rownames(num)<-paste0(num$V1,"(",num$`1`+num$`2`+num$`3`,")")
 rownames(prop)<-rownames(num)
 
 
 prop<-as.matrix(prop.table(table(rt$cancer_type,rt$MetaLiver),1))
 prop<-data.frame(prop)
 
 prop<-prop %>%
   pivot_wider(
     names_from =c("Var2"),
     values_from = Freq
   )
 prop<-data.frame(prop)

  colnames(prop)<-c("cancer_type","onlyLivermeta_precent","Livermeta_precent","non_Livermeta_precent")

 prop$liver_precent<-prop$onlyLivermeta_precent+prop$Livermeta_precent
 prop$liver_precent<-as.numeric(prop$liver_precent)

  g3<-ggplot(prop,aes(x= reorder(rownames(prop),-liver_precent),liver_precent,fill=cancer_type))+
   geom_bar(stat ="identity")+
   scale_color_manual(values=c(mycol)) +
   scale_fill_manual(values=c(mycol)) +
   geom_text_repel(aes(label = cancer_type), 
                   #color = mycol[1:5],
                   position = position_stack(vjust = 0.5),
                   show.legend = FALSE)+
   theme_classic() +xlab("")+
   #coord_polar(theta = "y") + 
   theme(axis.text.x = element_text(angle =70, hjust = 0.8,size = 10), 
         axis.text.y = element_text(colour="grey",angle = 30, size = 8),
         axis.title.y = element_text(angle = 90, size = 15))+
   theme(legend.position = "none")+coord_flip()
 g3

 #
 g3<-ggplot(prop,aes(x= reorder(rownames(prop),-liver_precent),liver_precent,fill=rownames(prop)))+
   geom_bar(stat ="identity")+
   geom_bar(stat="identity", alpha=0.5) +
   scale_fill_manual(values = mycol[c(9,3)])+
   #ylim(-100,120) +
   theme_bw() +
   #scale_fill_brewer(palette = "Set1") +
   coord_polar()+ 
   scale_color_manual(values=c(mycol)) +
   scale_fill_manual(values=c(mycol)) +
   #geom_text_repel(aes(label = cancer_type), 
                   #color = mycol[1:5],
                  # position = position_stack(vjust = 0.5),
                   #show.legend = FALSE)+
   theme_classic() +xlab("")+
   #coord_polar(theta = "y") + 
   theme(axis.text.x = element_text(angle =30, hjust = 0.7,size = 10), 
         axis.text.y = element_text(colour="grey",angle = 30, size = 8),
         axis.title.y = element_text(angle = 90, size = 15))+
   theme(legend.position = "none")
 g3

 ggsave(file.path(Clinicloutcome,"onlyLivermeta_precent_full.pdf"),width = 10,height =10)
 
 
 prop$liver_precent<-prop$onlyLivermeta_precent
 g1<-ggplot(prop,aes(x= reorder(rownames(prop),-liver_precent),liver_precent,fill=cancer_type))+
   geom_bar(stat ="identity")+
   scale_color_manual(values=c(mycol)) +
   scale_fill_manual(values=c(mycol)) +
   geom_text_repel(aes(label = cancer_type), 
                   #color = mycol[1:5],
                   position = position_stack(vjust = 0.5),
                   show.legend = FALSE)+
   theme_classic() +xlab("")+
   theme(axis.text.x = element_text(angle =70, hjust = 0.8,size = 10), 
         axis.text.y = element_text(colour="grey",angle = 30, size = 8),
         axis.title.y = element_text(angle = 90, size = 15))+
   theme(legend.position =  "none")+coord_flip()
 
 g1

 prop$liver_precent<-prop$Livermeta_precent
 g2<-ggplot(prop,aes(x= reorder(rownames(prop),-liver_precent),liver_precent,fill=cancer_type))+
   geom_bar(stat ="identity")+
   scale_color_manual(values=c(mycol)) +
   scale_fill_manual(values=c(mycol)) +
   theme_classic() +
   geom_text_repel(aes(label = cancer_type), 
                   #color = mycol[1:5],
                   position = position_stack(vjust = 0.5),
                   show.legend = FALSE)+
  xlab("")+
   theme(axis.text.x = element_text(angle =70, hjust = 0.8,size = 10), 
         axis.text.y = element_text(colour="grey",angle =30, size = 8),
         axis.title.y = element_text(angle = 90, size = 15))+
   theme(legend.position = "none")+coord_flip()
 g2
# ggsave("onlyLivermeta_precent.pdf",width = 13.83,height = 7.48)
 
 plot_grid(g1, g2, g3,
           labels = c("A", "B","C"), #
           rel_heights = c(1,1,1), # 3??ͼ?ı???
           #label_x=0,
           #label_y=1,
           align = 'v',ncol = 1, axis="lr", scale = c(1,1,1), greedy = F)
 ggsave(file.path(Clinicloutcome,"onlyLivermeta_precent.pdf"),width = 6,height = 18)

 
 # "tmb" # "msi_score"                 "msi_type" 
  #######################
 colnames(rt)
 table(rt$met_count.Live)
 site<-"met_count.Liver"
 metasites<-c("met_count.Liver","met_count.Bone","met_count.Skin","met_count.CNS.Brain","met_count.Lung" )
 metasites<-colnames(rt)[c(38:57)]
 for(site in metasites ) {
   print(site)
   subphe<-rt[rt[,site]!=0,]
   subpro<-cbind(prop.table(table(subphe$cancer_type)) %>% melt(),data.frame(table(subphe$cancer_type))[2])
   #subpro<-prop.table(table(subphe$tissue)) %>% melt()
   colnames(subpro)<-c('tissue',"perc","Freq")
   subpro$tissue <- fct_reorder(subpro$tissue, subpro$Freq)
   subpro$lable<-paste0(subpro$tissue,",",subpro$Freq,"(",round(subpro$perc*100,1),"%",")")
   ggplot(subpro, aes(x = "", y = perc, fill = tissue)) +
     geom_col(color = "black") +
     geom_label_repel(aes(label =   lable), max.overlaps=50,
                     position = position_stack(vjust = 0.5),
                     show.legend = FALSE) +
     scale_fill_manual(values =mycol)+
     guides(fill = guide_legend(title = "Metastasis site")) +
     #scale_fill_viridis_d("D") +
     coord_polar(theta = "y") + 
     theme_void()+
     theme(legend.position = "none")
   
   ggsave(file.path(Clinicloutcome,paste0(site,"orignial_cellDATA_site.pdf")),
          width =10, height =10,
          bg = "transparent", dpi = 300, useDingbats=FALSE)
 }
 

 
###TMB
 CellData<- read_excel("/home/data/gaoyuzhen/Projects/MET500Projects/Results/ClinicalDataForMeta/CellData_1.xlsx")
 rt<-CellData
 rt<-data.frame(rt)
 rt<-rt[rt$sample_type=="Metastasis",]
 rt<-rt[rt$MetaLiver!="3",]
###
 colnames(rt)
 rt$tmb<-as.numeric(rt$tmb)
 rt$MetaLiver<-as.factor(rt$MetaLiver)
 table(rt$MetaLiver)
 boxplot(rt$tmb~rt$MetaLiver)
# ggplot(rt, aes(x = MetaLiver, y =tmb, fill =MetaLiver)) +
 library(gghalves)
 library(ggpubr)
 
   c<- ggplot(rt,
              aes(x= MetaLiver, y=msi_score, 
                         fill =MetaLiver,
                         color =MetaLiver)) +
   geom_point(aes(shape=MetaLiver), 
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
