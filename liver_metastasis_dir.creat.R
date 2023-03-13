##lib and path set
results<-"Results"
Clinicloutcome   <-"Results/Clinicloutcome"
if (!file.exists(results)) { dir.create(results) }
if (!file.exists(Clinicloutcome)) { dir.create(Clinicloutcome) }

results<-"Results"
Piepath    <-"Results/Nature_Piepath"
if (!file.exists(results)) { dir.create(results) }
if (!file.exists(Piepath)) { dir.create(Piepath) }


Cluster.path<-"Results/cluster_DEG"
if (!file.exists(Cluster.path)) { dir.create(Cluster.path) }

DEG.path_gene <-"Results/gene-liver_DEG"
if (!file.exists(DEG.path_gene )) { dir.create(DEG.path_gene  ) }

DEG_cluster.path_gene <-"Results/DEG_cluster_gene"
if (!file.exists(DEG_cluster.path_gene )) { dir.create(DEG_cluster.path_gene ) }
################
Score.path <-"Results/ScoreFeatures"
if (!file.exists(Score.path )) { dir.create(Score.path) }

###features of lmscore 
LiverFeatures.path <-"Results/LiverFeatures"
ComparisonResults_bio <-"Results/LiverFeatures/ComparisonResults_bio"
ComparisonResults<-"Results/LiverFeatures/ComparisonResults"
if (!file.exists(LiverFeatures.path )) { dir.create(LiverFeatures.path) }
if (!file.exists(ComparisonResults_bio )) { dir.create(ComparisonResults_bio) }
if (!file.exists(ComparisonResults )) { dir.create(ComparisonResults) }

###features of lmscore 
LiverscoreFeatures.path <-"Results/LiverscoreFeatures"
ComparisonResults_score_bio <-"Results/LiverscoreFeatures/ComparisonResults_score_bio"
ComparisonResults_score<-"Results/LiverscoreFeatures/ComparisonResults_score"
if (!file.exists(LiverscoreFeatures.path )) { dir.create(LiverscoreFeatures.path) }
if (!file.exists(ComparisonResults_score_bio )) { dir.create(ComparisonResults_score_bio ) }
if (!file.exists(ComparisonResults_score)) { dir.create(ComparisonResults_score) }
######
TCGAresponse.path <-"Results/TCGAresponse.path"
if (!file.exists(TCGAresponse.path )) { dir.create(TCGAresponse.path) }
##GSE131418 大规模测序 结肠癌转移 （原位组织） 
GEOliverscore.path <-"Results/GEOliverscore.path"
if (!file.exists(GEOliverscore.path )) { dir.create(GEOliverscore.path) }
###SingleScoreFeatures
SingleScore.path <-"Results/SingleScoreFeatures"
if (!file.exists(SingleScore.path )) { dir.create(SingleScore.path) }