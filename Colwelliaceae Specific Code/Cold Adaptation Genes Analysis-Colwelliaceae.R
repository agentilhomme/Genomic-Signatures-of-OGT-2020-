source("Cold Adaptation Genes Calc.R")

#---------------------------------------------------------load files---------------------------------------------------------------------
CladeID = read.table("./Data/Sub-Clade ID.tsv",fill=T,head=T,sep="\t",colClasses="character")

#---------------------------------------------------------format files ------------------------------------------------------------------
colnames(CladeID)[1] = "strn"
#set up whole dataframe to analyze by cluster
CladeID$patric.id = genome_ID$genome.genome_id[match(CladeID$strn,genome_ID$genome.genome_name)]

#--------------------------------------------------Analyze cold and warm clusters of hornerea ----------------------------------------------
# isolate the 3 different hornerea clades from the 
CH_cold_id = subset.data.frame(CladeID, ClusterName == "Clade A2");rownames(CH_cold_id)=CH_cold_id$patric.id
CH_cold = GeneP1[,unique(CH_cold_id$strn)]
CH_cold = t(CH_cold>0)+0

CH_warm_1_id = subset.data.frame(CladeID, ClusterName == "Clade A3");rownames(CH_warm_1_id) = CH_warm_1_id$patric.id
CH_warm_1 = GeneP1[,unique(CH_warm_1_id$strn),drop = F]
CH_warm_1 = t(CH_warm_1>0)+0

CH_warm_id = subset.data.frame(CladeID, ClusterName == "Clade A1"); rownames(CH_warm_id) = CH_warm_id$patric.id
CH_warm = GeneP1[,unique(CH_warm_id$strn),drop = F]
CH_warm = t(CH_warm>0)+0

# remove protein families that are not present in any strain in the cluster
CH_cold_PFG = CH_cold[,colSums(CH_cold)>0] 
CH_cold_PFG = t(CH_cold_PFG)

CH_warm_PFG = CH_warm[,colSums(CH_warm)>0]
CH_warm_PFG = t(CH_warm_PFG)

CH_warm_1_PFG = CH_warm_1[,colSums(CH_warm_1)>0, drop = F]
CH_warm_1_PFG = t(CH_warm_1_PFG)

# keep protein present in all strains
CH_cold_commonPFG = CH_cold_PFG[rowSums(CH_cold_PFG) == ncol(CH_cold_PFG),] 
CH_cold_commonPFG = rownames(CH_cold_commonPFG)

CH_warm_commonPFG = CH_warm_PFG[rowSums(CH_warm_PFG) == ncol(CH_warm_PFG),]
CH_warm_commonPFG = rownames(CH_warm_commonPFG)

CH_warm_1_commonPFG = rownames(CH_warm_1_PFG)

# set diff of the the different unique global protein families 
diffcw1 = setdiff(CH_cold_commonPFG,CH_warm_1_commonPFG) # present in the cold but not in the warm 1
diffcw = setdiff(CH_cold_commonPFG,CH_warm_commonPFG) # present in the cold but not in the warm cluster
diffw1c = setdiff(CH_warm_1_commonPFG,CH_cold_commonPFG) # present in the warm 1 but not in the cold
diffwc = setdiff(CH_warm_commonPFG,CH_cold_commonPFG) # present in the warm cluster but not in the cold 

#------------------------------------------------create cluster different protein datasets----------------------------------------------------
colnames(blast_output_feat)[5] = "genome.genome_id"

# get blast information of  protein families present in the cold cluster that are cold adaptive
ccc.1=blast_output_feat[blast_output_feat$feature.pgfam_id %in% diffcw1,]
ccc.all=blast_output_feat[blast_output_feat$feature.pgfam_id %in% diffcw,]

# get first (best) blast hit only of the protein present in each cold-adaptive protein family
ccc.1=ccc.1[!duplicated(ccc.1$feature.patric_id),]
ccc.all=ccc.all[!duplicated(ccc.all$feature.patric_id),]

# get all features (PATRIC) for the cold cluster genomes
ddd.1=feat[feat$genome.genome_id %in% rownames(CH_cold_id),]
ddd.all=feat[feat$genome.genome_id %in% rownames(CH_cold_id),]

# get only protein families (PATRIC) that match protein families that are present in the cold but absent in the warm
fff.1=ddd.1[ddd.1$feature.pgfam_id %in% diffcw1,]
fff.all=ddd.all[ddd.all$feature.pgfam_id %in% diffcw,]

# merge blast hits and features so that there are no discrepancies between the blast outputs and the PATRIC data 
ggg.1=merge.data.frame(ccc.1,fff.1,by="feature.patric_id",all.y=T)
ggg.all=merge.data.frame(ccc.all,fff.all,by="feature.patric_id",all.y=T)

# merge blast hits, features, and patric global family names because some protein families still don't have descriptors
hhh.1=merge.data.frame(ggg.1,cfg,by.x="feature.pgfam_id.y",by.y="feature.pgfam_id",all=F)
hhh.all=merge.data.frame(ggg.all,cfg,by.x="feature.pgfam_id.y",by.y="feature.pgfam_id",all=F)

# create customary data frame of  protein IDs ( present in Clade A2, Absent Clade A3)
diffcw1_data = hhh.1[,c(1,2,5,9,12)]
diffcw1_data = filter(diffcw1_data,e.val < 1e-20 | is.na(e.val))
diffcw1_data$e.val = NULL
diff1PFG = diffcw1_data[!duplicated(diffcw1_data$feature.pgfam_id.y),] # get the number of unique protein families 
diff1PFG = diff1PFG[,-c(2,3)]
colnames(diff1PFG) = c("Global Protein Family ID","GPF Description")

# data frame with unique protein IDs and genome IDs
diffcw1_data=diffcw1_data[!duplicated(diffcw1_data$feature.patric_id),]
colnames(diffcw1_data)= c("Global Protein Family (GPF) ID","Protein ID","Genome ID", "GPF Description")
diffcw1_data$`Protein ID` = NULL; diffcw1_data = unique(diffcw1_data)

# create customary data frame of  protein IDs ( present in Clade A2, Absent Clade A1)
diffcw_data = hhh.all[,c(1,2,5,9,12)]
diffcw_data = filter(diffcw_data,e.val < 1e-20 | is.na(e.val)); diffcw_data$e.val = NULL
diffPFG = diffcw_data[!duplicated(diffcw_data$feature.pgfam_id.y),] # get the number of unique protein families 
diffPFG = diffPFG[,-c(2,3)]
colnames(diffPFG) = c("Global Protein Family ID","GPF Description")

# data frame with unique protein IDs and genome IDs
diffcw_data=diffcw_data[!duplicated(diffcw_data$feature.patric_id),]
colnames(diffcw_data)= c("Global Protein Family (GPF) ID","Protein ID","Genome ID", "GPF Description")
diffcw_data$`Protein ID` = NULL; diffcw_data = unique(diffcw_data)

# merge the dataframes to print
diffcw_data$`Clade A1` = "Absent"
diffcw1_data$`Clade A3` = "Absent"
diff_name = full_join(diffcw_data,diffcw1_data,by = c("Global Protein Family (GPF) ID","Genome ID","GPF Description"))

diff1PFG$`Clade A1` = "Absent"
diffPFG$`Clade A3` = "Absent"
diff_PFG = full_join(diff1PFG,diffPFG,by = c("Global Protein Family ID","GPF Description"))

#-----------------------------------------------------------------------export dataframes----------------------------------------------------------------------
write.xlsx(diff_name, file = "./data/Differential PFG in Colwellia A Sub-Clusters.xlsx", sheetName = "Cold Present, Warm Cluster Absent, w/ Genome",row.names = F)
write.xlsx(diff_PFG, file = "./data/Differential PFG in Colwellia A Sub-Clusters.xlsx", sheetName = "Unique PFG warm absent",append = T,row.names = F)
