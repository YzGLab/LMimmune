##drug resisitence

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
phe_met$group_cluster<-paste0("C",phe_met$group_cluster)
#immunefeaturs<-cbind(Sample_id=rownames(immunefeaturs),immunefeaturs)
phe_met$Sample_id<-gsub("-",".",phe_met$Sample_id)
phe_total<-merge(phe_met,immunefeaturs,by="Sample_id")
##az  asZXaszwaq
#################################
table(phe_met$group_cluster)
#phe_met    $group_cluster<-ifelse(phe_met$group_cluster %in% c("C1","C6","C4"),"C4",phe_met$group_cluster)
phe_drug<-phe_met[phe_met$group_cluster %in% c("C4","C3"),]

#
phe_drug$group<-ifelse(phe_drug$group_cluster=="C4","1","2")
head(phe_drug)
rownames(phe_drug)<-phe_drug$Sample_id

##
#("Erlotinib", "Rapamycin", "Sunitinib","PHA-665752", "MG-132", "Paclitaxel", "Cyclopamine", "AZ628", "Sorafenib", 
#"VX-680", "Imatinib", "TAE684", "Crizotinib", "Saracatinib", "S-Trityl-L-cysteine", "Z-LLNle-CHO", "Dasatinib",
#GNF-2, CGP-60474, CGP-082996, A-770041, WH-4-023, WZ-1-84, BI-2536, BMS-536924, BMS-509744, CMK, 
#Pyrimethamine, JW-7-52-1, A-443654, GW843682X, MS-275, Parthenolide, KIN001-135, TGX221, Bortezomib, 
#XMD8-85, Roscovitine, Salubrinal, Lapatinib, GSK269962A, Doxorubicin, Etoposide, Gemcitabine,
#Mitomycin C, Vinorelbine, NSC-87877, Bicalutamide, QS11, CP466722, Midostaurin, CHIR-99021, AP-24534, 
#AZD6482, JNK-9L, PF-562271, HG-6-64-1, JQ1, JQ12, DMOG, FTI-277, OSU-03012, Shikonin, AKT inhibitor VIII, 
#Embelin, FH535, PAC-1, IPA-3, GSK-650394, BAY 61-3606, 5-Fluorouracil, Thapsigargin, Obatoclax Mesylate, 
#BMS-754807, Lisitinib, Bexarotene, Bleomycin, LFM-A13, GW-2580, AUY922, Phenformin, Bryostatin 1, Pazopanib, 
#LAQ824, Epothilone B, GSK1904529A, BMS345541, Tipifarnib, BMS-708163, Ruxolitinib, AS601245, Ispinesib Mesylate,
#TL-2-105, AT-7519, TAK-715, BX-912, ZSTK474, AS605240, Genentech Cpd 10, GSK1070916, KIN001-102, LY317615, GSK429286A,
#FMK, QL-XII-47, CAL-101, UNC0638, XL-184, WZ3105, XMD14-99, AC220, CP724714, JW-7-24-1, NPK76-II-72-1, STF-62247,
#NG-25, TL-1-85, VX-11e, FR-180204, Tubastatin A, Zibotentan, YM155, NSC-207895, VNLG/124, AR-42, CUDC-101, Belinostat, 
#I-BET-762, CAY10603, Linifanib , BIX02189, CH5424802, EKB-569, GSK2126458, KIN001-236, KIN001-244, KIN001-055, KIN001-260,
##KIN001-266, Masitinib, MP470, MPS-1-IN-1, BHG712, OSI-930, OSI-027, CX-5461, PHA-793887, PI-103, PIK-93, 
#SB52334, TPCA-1, TG101348, Foretinib, Y-39983, YM201636, Tivozanib, GSK690693, SNX-2112, QL-XI-92, XMD13-2, QL-X-138, XMD15-27)

drugsnames<-read.delim2("/home/data/gaoyuzhen/Projects/MET500Projects/Drugs.txt",header=F)
drugsnames<-drugsnames$V1
library(MOVICS)
phe_drug_exp<-data.frame(metExp)
pseudo.moic.total <- list("clust.res" = phe_drug,
                          "mo.method" = "test")
pseudo.moic.total$clust.res$samID <- rownames(pseudo.moic.total$clust.res)
pseudo.moic.total$clust.res$clust <- pseudo.moic.total$clust.res$group
drug.prad<- compDrugsen(moic.res    = pseudo.moic.total,
                        norm.expr   = phe_drug_exp[,pseudo.moic.total$clust.res$samID]+1, # double guarantee sample order
                        #drugs       = c("Paclitaxel","Cisplatin","5-Fluorouracil","Sorafenib", "Sunitinib"), # a vector of names of drug in GDSC
                        #drugs=c("Erlotinib","Rapamycin","Sunitinib", "PHA-665752","MG-132","Paclitaxel","Cyclopamine", "AZ628","Sorafenib" ,"VX-680","Imatinib", "TAE684"),
                        drugs=c(drugsnames[-1]),
                        #tissueType  = "digestive_system", # choose specific tissue type to construct ridge regression model
                        clust.col   = c("#0088A2","lightgrey"),
                        fig.path="Results",
                        test.method = "parametric", # statistical testing method
                        seed = 123456,
                        prefix      = "BOXVIOLIN OF ESTIMATED IC50") 
