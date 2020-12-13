library(xlsx)
library(ape)
library(phytools)
library(phylogram)
library(dplyr)

#---------------------------------------------------------load files---------------------------------------------------------------------
# load PATRIC phylogenetic and protein family files
treefile="./Data/Colwelliaceae Rooted Tree.nwk"
genomes = read.table("./Data/Colwelliaceae UnRooted PATRIC Data.tsv",fill=T,head=T,sep="\t",colClasses="character")
feat = read.table("./Data/Colwelliaceae PF Features.tsv",fill = T,head = T,sep = "\t",colClasses = "character")
cfg = read.table("./Data/Colwelliaceae_families_global.tsv",fill=T,head=F,sep="\t",colClasses="character",quote="")

#load presence/absence dataframe calculated from Protein Family Indices.R
Gene_Prescence = read.csv("./Data/Gene Prescence.csv",fill = T,header = T, check.names = F)

#load OGT dataframe calculated from Strain OGT.R
ratk = read.table("./Data/Strain Ratk OGT Final.tsv",fill=T,head=T,sep="\t",colClasses="character")

#blast output of identified cold adaptation genes against strains of interst
blast_output = read.table("./Data/blast_CA_Colwelliaceae.txt",fill = T, header = T, sep = '\t',quote="")

#---------------------------------------------------------format files ------------------------------------------------------------------
rownames(genomes) = genomes$genome.genome_id
colnames(feat)[1] = "genome.genome_id"
rownames(Gene_Prescence) = Gene_Prescence[,1]; Gene_Prescence[,1] = NULL
colnames(ratk) = c("strn","topt","Type")

rownames(cfg)=cfg$V1
cfg$V1=NULL
colnames(cfg) = c("feature.pgfam_id","type","description")

#set up tree
tre.phylo = read.tree(treefile)
tre.phylo = drop.tip(tre.phylo,"1380381.3") #drop long tip (poor genome assembly)
common_ids = intersect(tre.phylo$tip.label,rownames(genomes))
tre.phylo = keep.tip(tre.phylo, common_ids) # keeping the tips that are in the dataframe 
tre.phylo$node.label = NULL
tre.phylo = midpoint.root(tre.phylo)
tre.phylo = ladderize(tre.phylo)

dend <- read.dendrogram(textConnection(write.tree(tre.phylo))) # convert phylo to dend
dend = as.dendrogram(dend)

#----------------------------------------create dataframe of cold adaptation proteins take from the family of interest-------------------------------------------------------
# set up blast information dataframe 
blast_output$sid = gsub(":","|",blast_output$sid)
colnames(blast_output)[4] = "feature.patric_id"
blast_output_feat = merge.data.frame(blast_output,feat,by = "feature.patric_id", all = F)
blast_output_feat = blast_output_feat[,-c(2,4,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21)]

# add the patric ids to the OGT strains
genome_ID = genomes[,c(2,3)]
genome_ID = subset(genome_ID,genome.genome_id %in% tre.phylo$tip.label)
ratk$id.patric = genome_ID$genome.genome_id[match(ratk$strn,genome_ID$genome.genome_name)]
ratk$topt = as.numeric(ratk$topt)
ratk = na.omit(ratk)
rownames(ratk) = ratk$id.patric
blast_match = as.character(unique(blast_output_feat$feature.pgfam_id))

#set up the Gene Presence dataframe
Gene_Prescence = Gene_Prescence[,labels(dend)] # reorders the matrix according to the dendrogram tip labels
GeneP1 = as.matrix(Gene_Prescence)
GeneP1[is.infinite(GeneP1)] = NA
GeneP1[is.nan(GeneP1)] = NA
GeneP1 = na.omit(GeneP1)

# subset the data from just the uniqhits
GeneP0 = GeneP1 # find the pure global protein family differences between the OGT clusters( not just the U.S)
GeneP1 = GeneP1[rownames(GeneP1) %in% blast_match,] 

#set up the strain names instead of the Patric ID numbers 
colnames(GeneP1) = genomes[colnames(GeneP1),"genome.genome_name"]
GeneP2 = GeneP1[rowSums(GeneP1>0)>0,] #  removes protein families that are only present in one/ 5 genome 

