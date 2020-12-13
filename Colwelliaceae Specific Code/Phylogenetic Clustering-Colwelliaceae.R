# Seperate family tree into main genera and into sub-genera clusters based on OGT calculations
library(ape) # read.tree()
library(phytools) # midpoint.root()
library(treeio) # rename_taxa()
library(DECIPHER) # ReadDendrogram()
library(ggtree)
library(ggpubr) # stat_cor

#-----------------------------------------------------------------------load files------------------------------------------------------------------------
# load PATRIC phylogenetic files
treefile="./Data/Colwelliaceae Rooted Tree.nwk"
genomes = read.table("./Data/Colwelliaceae Rooted PATRIC Data.tsv",fill=T,head=T,sep="\t",colClasses="character")

# load growth rates, OGT, and OGT parameters of Colwelliaceae strains (Literature and novel)
GR = read.table("./Data/growthrates.tsv",fill=T,head=T,sep="\t",colClasses="character")
OGT = read.table("./Data/Strain Ratk OGT Final.tsv",fill=T,head=T,sep="\t",colClasses="character")
biospec = read.csv("./References/Literature Colwellia Growthrates.csv")

#nodes of interest to split the tree
main_nodes = c(169,145,142,116) # order of shallow to deep branched clade
temp_nodes = c(116,142,145,154,169,172,173,179)

#-----------------------------------------------------------------------set up files----------------------------------------------------------------------- 
#read tree
tre.phylo = read.tree(treefile)

# remove any genomes that had poor assembly
tre.phylo = drop.tip(tre.phylo,"1380381.3")

# subset genomes that are only used in the phylogenetic tree 
common_ids = intersect(tre.phylo$tip.label,genomes$genome.genome_id)
tre.phylo = keep.tip(tre.phylo, common_ids) # keeping the tips that are in the dataframe 

#set up tree and rename tips with strain names 
tre.phylo = midpoint.root(tre.phylo)
tre.phylo = ladderize(tre.phylo)
tre.phylo = rename_taxa(tre.phylo,genomes,genome.genome_id,genome.genome_name) # rename the ends. 

# convert to dendrogram
dend = tre.phylo
dend$node.label = NULL
dend <- ReadDendrogram(textConnection(write.tree(dend)))

#set up column names and format to numeric 
colnames(OGT)[1] = "strain"
OGT$OGT = as.numeric(OGT$OGT)
GR = GR[,c(1,3,4)]
GR$mumax = as.numeric(GR$mumax)

#----------------------------------------------------------------------set up tree splits and plot--------------------------------------------------------
#set up clusters for the main groups
tre.phylo.main = groupClade(tre.phylo,main_nodes)
tree.DF.main = as.data.frame(as_tibble(tre.phylo.main))
tree.DF.main = tree.DF.main[1:length(tre.phylo.main$tip.label),] # get just the groups associated with the tips

#set up clusters for the temperature phenotype groups
tre.phylo.temp = groupClade(tre.phylo,temp_nodes)
tree.DF.temp = as.data.frame(as_tibble(tre.phylo.temp))
tree.DF.temp = tree.DF.temp[1:length(tre.phylo.temp$tip.label),]

#--------------------------------------------------------------------set up ID files for each cluster type ------------------------------------------------
# set up the names of main clusters based on the group number
Main_ClusterID = data.frame(group = 0:4,ClusterName = c("Other Colwelliaceae","Clade A","Clade B","Other Colwelliaceae","Clade C"))
tree.DF.main = merge.data.frame(tree.DF.main,Main_ClusterID, by = "group", all = T)
tree.DF.main = tree.DF.main[,5:6]

# set up the names of temperature clusters based on the group number
Temp_ClusterID = data.frame(group = c(0:8),ClusterName = c("Other Colwelliaceae","Clade C","Other Colwelliaceae","Clade B2","Clade B1","Other Colwelliaceae","Clade A3","Clade A2","Clade A1"))
tree.DF.temp = merge.data.frame(tree.DF.temp,Temp_ClusterID, by = "group", all = T)
tree.DF.temp = tree.DF.temp[,5:6]

#set up color scheme
Temp.Clust.colors = data.frame(Clade = c("Clade A1","Clade A2","Clade A3","Clade B1","Clade B2","Clade C","Other Colwelliaceae"),
                               clustCol = c("seagreen3","olivedrab4","chartreuse3","darkorange2","darkgoldenrod3",
                                            "purple","black"))
col <- as.character(Temp.Clust.colors$clustCol)
names(col) = as.character(Temp.Clust.colors$Clade)

# --------------------------------------------- write files -----------------------------------------------------------------------------
write.table(tree.DF.main,"./Data/Genera Clade ID.tsv",sep = "\t",quote = F,row.names = F)
write.table(tree.DF.temp,"./Data/Sub-Clade ID.tsv",sep = "\t",quote = F,row.names = F)

#-------------------------------------Calculate average pairwise differences by cluster type --------------------------------------------------------
#phylo pairwise distances by genera
DistM.main = as.matrix(cophenetic.phylo(tre.phylo.main))
DistM.main = as.data.frame(DistM.main)

Main_Clade_Avg = data.frame()
for(i in 1:length(unique(tree.DF.main$ClusterName))) {
  
  #get clade ID
  clade_ID = sort(unique(tree.DF.main$ClusterName))[i]
  print(as.character(clade_ID))
  
  #get strain names within the clade
  names = tree.DF.main[tree.DF.main$ClusterName == clade_ID,][["label"]]
  
  #subset pairwise distance matrix by names and upper tri for calculations
  data = DistM.main[names,names]
  x = data*upper.tri(data)
  x[x==0] = NA
  
  #find average
  avg = mean(as.vector(t(x)),na.rm = T)
  
  #create dataframe
  avg_df = data.frame(clade = paste(clade_ID),Dist_avg = avg,strain_num = length(names))
  Main_Clade_Avg = rbind(Main_Clade_Avg,avg_df)
}
#print
Main_Clade_Avg

#phylo pairwise distances by genera sub-clades
DistM.temp = as.matrix(cophenetic.phylo(tre.phylo.temp))
DistM.temp = as.data.frame(DistM.temp)
Temp_Clade_Avg = data.frame()
for(i in 1:length(unique(tree.DF.temp$ClusterName))) {
  
  #get clade ID
  clade_ID = sort(unique(tree.DF.temp$ClusterName))[i]
  print(as.character(clade_ID))
  
  #get strain names within the clade
  names = tree.DF.temp[tree.DF.temp$ClusterName == clade_ID,][["label"]]
  
  #subset pairwise distance matrix by names and upper trie for calculations
  data = DistM.temp[names,names]
  x = data*upper.tri(data)
  x[x==0] = NA
  
  #find average
  avg = mean(as.vector(t(x)),na.rm = T)
  
  # create dataframe
  avg_df = data.frame(clade = paste(clade_ID),Dist_avg = avg,strain_num = length(names))
  Temp_Clade_Avg = rbind(Temp_Clade_Avg,avg_df)
}
#print
Temp_Clade_Avg

#--------------------------------Calculate average and standard devation of growth rate and OGT by Main Clade and Sub-Clade---------------------------
# merge growth rate data, OGT data, temperature, and main cluster ID
GR_OGT_Clade = merge.data.frame(GR,OGT,by = "strain", all = F)
GR_OGT_Clade$`Temp ClusterName` = tree.DF.temp$ClusterName[match(GR_OGT_Clade$strain,tree.DF.temp$label)]
GR_OGT_Clade$`Main ClusterName` = tree.DF.main$ClusterName[match(GR_OGT_Clade$strain,tree.DF.main$label)]
GR_OGT_Clade$Type = NULL
GR_OGT_Clade = GR_OGT_Clade[!is.na(GR_OGT_Clade$`Temp ClusterName`),]

# calculate by Temperature Clade the OGT
OGT_TC = GR_OGT_Clade %>% dplyr::group_by(`Temp ClusterName`) %>% dplyr::summarize(`Average OGT` = paste0(signif(mean(OGT),2)," ± ",signif(sd(OGT),2)),
                                                                  `Min OGT` = signif(min(OGT)),`Max OGT` = signif(max(OGT)), N = length(unique(strain)))
# calculate by Main Clade the OGT
OGT_MC = GR_OGT_Clade %>% dplyr::group_by(`Main ClusterName`) %>% dplyr::summarize(`Average OGT` = paste0(signif(mean(OGT),2)," ± ",signif(sd(OGT),2)),
                                                                            `Min OGT` = signif(min(OGT)),`Max OGT` = signif(max(OGT)), N = length(unique(strain)))
#calculate by Temperature Clade the growth rate
GR_TC = GR_OGT_Clade %>% dplyr::group_by(`Temp ClusterName`) %>% dplyr::summarize(`Average GR` = paste0(signif(mean(mumax),2)," ± ",signif(sd(mumax),2)),
                                                                            `Min GR` = signif(min(mumax)),`Max GR` = signif(max(mumax)), N = length(unique(strain)))
#calculate by Main Clade the growth rate
GR_MC = GR_OGT_Clade %>% dplyr::group_by(`Main ClusterName`) %>% dplyr::summarize(`Average GR` = paste0(signif(mean(mumax),2)," ± ",signif(sd(mumax),2)),
                                                                           `Min GR` = signif(min(mumax)),`Max GR` = signif(max(mumax)), N = length(unique(strain)))
# ANOVA of growthrate between temperature tested
GR_Temp_Tested_aov = aov(mumax~temp,data = GR_OGT_Clade)
summary(GR_Temp_Tested_aov)

#---------------------------------------------------------------- Plot Phylogenetic Tree -------------------------------------------------------------
# merge temp ID with OGT and clade colors and make the strains rownames
TC_df = merge(OGT, tree.DF.temp, by.x = "strain", by.y = "label", all.y= TRUE)
TC_df$`SubColor` = Temp.Clust.colors$clustCol[match(TC_df$ClusterName,Temp.Clust.colors$Clade)]
rownames(TC_df) = TC_df$strain
TC_df = as.matrix(TC_df)

#set the strain names in the order of the phylogenetic tree
TC_df = TC_df[tre.phylo$tip.label,]
TC_df = as.data.frame(TC_df)
TC_df$SubColor = as.character(TC_df$SubColor)

#set Sub-Clade colors as a factor list
SubClade_color = setNames(TC_df$SubColor,TC_df$ClusterName)

# create a dataframe for colors based on bootstrap values 
nodeColor = data.frame(node = tre.phylo$node.label,color = NA)
nodeColor$node = as.character(nodeColor$node)
nodeColor$nodenum = as.numeric(nodeColor$node)
for (i in 1:nrow(nodeColor)){
  if((is.na(nodeColor$nodenum[i]) & nodeColor$node[i] == "Root") == T){
    nodeColor$color[i] = "white"
  }else if(nodeColor$node[i] == "100"){
    nodeColor$color[i] = "dark grey"
  } else if(is.na(nodeColor$nodenum[i]) & nodeColor$node[i] == ""){
    nodeColor$color[i] = "brown"
  }else if(nodeColor$nodenum[i] < 100 && nodeColor$nodenum[i] > 90){
    nodeColor$color[i] = "dark grey"
  }else {
    nodeColor$color[i] = "brown"
  }
}

#plain tree with strains as tiplabels
p = ggtree(tre.phylo, size = 2) + 
  theme_tree2() +
  geom_nodepoint(color = nodeColor$color,size = 4)

# add an X where the main clades are, and their main cluster names
p1 = p + geom_nodepoint(aes(subset = node == main_nodes[4], x = x - branch.length * 0.5),shape = 4,size = 7) +
  geom_text2(aes(subset = node == main_nodes[4],x = x - branch.length * 0.5, label = "Clade \nC", fontface = "italic"),vjust = -0.5,size = 8, color = "purple") + 
  geom_nodepoint(aes(subset = node == main_nodes[2], x = x - branch.length * 0.5),shape = 4,size = 7) + 
  geom_text2(aes(subset = node == main_nodes[2],x = x - branch.length * 0.5, label = "Clade \nB",fontface = "italic"),vjust = -0.5,size = 8, hjust = 0.3,color = "orange") +
  geom_nodepoint(aes(subset = node == main_nodes[1], x = x - branch.length * 0.5),shape = 4,size = 7) + 
  geom_text2(aes(subset = node == main_nodes[1],x = x - branch.length * 0.5, label = "Clade \nA",fontface = "italic"),vjust = -0.5,size = 8,color = "darkgreen")

#have colored strain names according to sub-clade and all types of OGT source 
p2 = p1 %<+% TC_df + 
  geom_tiplab(aes(color = ClusterName),cex = 3.5,offset = .05) + 
  geom_tippoint(aes(subset = !is.na(Type),shape = Type,x = x + 0.03),size = 2.5) + 
  scale_color_manual(values = SubClade_color) + 
  theme(legend.position = "none") + 
  xlim(0,3.2) 

# add bar to designate the subclades 
p3 = p2 +
  geom_strip("Colwellia BRX8.6","Colwellia hornerae strain IC036",barsize = 2,label = "Clade A1",offset = 0.8,fontsize = 4,hjust = -0.1) + 
  geom_strip("Colwellia MB02u-12","Colwellia MB02u-6",barsize = 2,label = "Clade B1", offset = 0.8,fontsize = 4,hjust = -0.1) + 
  geom_strip("Colwellia sp. C1TZA3","Colwellia mytili strain KCTC 52417",barsize = 2, label = "Clade B2 (BRX10.3 & Other Clade B) ",offset = 0.8,fontsize = 3,hjust = -0.1)+
  geom_strip("Colwellia Bg11-28","Colwellia echini strain A3T",barsize = 2, label = "Clade C",offset = 0.8,fontsize = 4,hjust = -0.1) + 
  geom_strip("Colwellia MB3u-64","Colwellia MB02u-14",barsize = 2,label = "Clade A2",offset = 0.8,fontsize = 4,hjust = -0.1) + 
  geom_strip("Colwellia Bg11-12",'Colwellia Bg11-12',barsize = 2, label = "Clade A3",offset = 0.8,fontsize = 4,hjust = -0.1)

#save the figure 
tiff("./Plots/Phylogenetic Tree.tiff",width = 12,height = 17, units = "in", res = 500)
p3
dev.off()

#------------------------------------------ Plot growth rate vs. OGT  by temperature tested colored by Sub-Clade -------------------------------------
# add the colors to the dataframe and filter it by OGTs under 25 degrees C and set up facet names
GR_OGT_Clade$clustCol = Temp.Clust.colors$clustCol[match(GR_OGT_Clade$`Temp ClusterName`,Temp.Clust.colors$Clade)]
GR_OGT_Clade$temp = as.numeric(GR_OGT_Clade$temp)
Temp.labs = c("(A) Growth at -1°C","(B) Growth at 4°C","(C) Growth at 11°C","(D) Growth at 17°C")
names(Temp.labs) = c("-1","4","11","17")

#plot and save the figure
tiff("./Plots/Avg GR vs. OGT by Temp & Clade.tiff",width = 15,height = 10, units = "in", res = 1200)
ggplot(GR_OGT_Clade,aes(x = OGT,y = mumax)) +
  geom_point(aes(color = `Temp ClusterName`)) +
  stat_smooth(method = "lm", se = F,aes(group = 1))+
  stat_cor(method = "pearson",aes(label = paste(..rr.label..,..r.label..,..p.label.., sep = "~`,`~")),label.x = 6)+
  labs(x = "Optimal Growth Temperature (°C)", y = expression("Growth Rate" ~ (d^{-1})))+
  scale_color_manual("Taxonomic Clades",values = col) + 
  facet_grid(.~temp,scales = "free",labeller = labeller(temp = Temp.labs)) +
  theme_classic(base_line_size = 1, base_size = 16, base_rect_size = 1) + 
  theme(panel.grid.major = element_line(color = "lightgrey", size = 0.1), panel.border = element_rect(colour = "black", size = 1, fill = NA), legend.position = "bottom")
dev.off()

#---------------------------------------------------- Box plot of OGT by Temperature Sub-Clade ---------------------------------------------------------
#plot and save the figure 
tiff("./Plots/OGT Boxplot by Temp Cluster.tiff",width = 11,height = 8, units = "in", res = 1200)
ggplot(GR_OGT_Clade,aes(x = `Temp ClusterName`,y = OGT)) + 
  geom_boxplot(lwd = 1) +
  labs(x = "Taxonomic Clades", y = "Optimal Growth Temperature (°C)") +
  theme_classic(base_line_size = 1, base_size = 16) + 
  theme(panel.grid.major = element_line(color = "lightgrey", size = 0.1))
dev.off()


