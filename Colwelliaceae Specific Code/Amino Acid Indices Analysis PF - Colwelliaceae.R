# create pangenome vs core genome plot, volcano plot, and pairwise differences of protein family and OGT
source("./Plain Code/f-pangenome.R")
library(DECIPHER) # ReadDendrogram()
library(ape) # read.tree()
library(reshape2) # melt()
library(RColorBrewer) # brewer.ap()
library(phytools) #midpoint.root()
library(dendextend)  # set()
library(colorspace) # rainbow_hcl()
library(phylogram) # read.dendrogram()
library(rlang) # is_empty()
library(dplyr) # filter()
library(ggplot2) # ggplot()
library(ggridges) # geom_density_ridges()
theme_set(theme_ridges())

# ------------------------------------------------------------load files -------------------------------------------------------------
# load PATRIC phylogenetic and protein family files
treefile="./Data/Colwelliaceae Rooted Tree.nwk"
genomes = read.table("./Data/Colwelliaceae Rooted PATRIC Data.tsv",fill=T,head=T,sep="\t",colClasses="character")
colfeat = read.table("./Data/Colwelliaceae PF Features.tsv",fill = T,head = T,sep = "\t",colClasses = "character") # patric protein features
plfdat = read.table("./Data/Colwelliaceae_families_local.tsv",sep="\t",row.names = 1, stringsAsFactors=F, quote="",comment="")[,3,drop=F]
pgfdat = read.table("./Data/Colwelliaceae_families_global.tsv",sep="\t",row.names=1, stringsAsFactors=F, quote="",comment="")[,3,drop=F]

# load Protein family amino acid indice calculations
colind = read.table("./Data/Protein Indices Data.tsv",row.names = 1, sep="\t",stringsAsFactors = F) # Protein Amino Acid Index
colsave = readRDS("./Data/Colwelliaceae PF AA Calc.RDS")

# load OGT by strain and  Sub-Clade IDs
ratk = read.table("./Data/Strain Ratk OGT Final.tsv",sep="\t",head=T,stringsAsFactors=FALSE) # OGT table
Temp.Clust.ID = read.table("./Data/Sub-Clade ID.tsv",fill=T,head=T,sep="\t",colClasses="character")

# set up list of literature based Colwellia growth rates 
ex_col = paste(c("Colwellia psychrerythraea ACAM 605","Colwellia hornerae strain ACAM 607","Colwellia piezophila ATCC BAA-637","Colwellia demingiae strain ACAM 459",
                 "Colwellia psychotropica"),collapse = "|")
taxaname = "Colwellia (PGF)"
#---------------------------------------------------------- Format Files -----------------------------------------------------------------
#read whole-genome tree from PATRIC
tre.dend = ReadDendrogram(treefile)
tre.phylo = read.tree(treefile)

#set up patric metatable
rownames(genomes) = genomes$genome.genome_id

#set up protein features
rownames(colfeat) = colfeat$feature.patric_id
colfeat = colfeat[,-2]
colnames(colfeat)[1] = "genome.genome_id"

# limit datasets to common genomes and merge protein features and amino acid indices into dataframe
common_ids = intersect(rownames(colfeat), rownames(colind))
colfeat = colfeat[common_ids, ]
colind = colind[common_ids, ]
colout = data.frame(colfeat, colind)

# add the patric ids to the OGT strains
rownames(ratk) = ratk$strn
genome_ID = genomes[,c(2,3)]
genome_ID = subset(genome_ID,genome.genome_id %in% tre.phylo$tip.label)
ratk$id.patric = genome_ID$genome.genome_id[match(ratk$strn,genome_ID$genome.genome_name)]

#-------------------------------------------------------set up phylogeny -----------------------------------------------------------------
#use the global protein family
mat = as.matrix(colsave[[1]])
mat[is.infinite(mat)] = NA
mat[is.nan(mat)] = NA
mat[mat==0] = NA

# keep only genomes that are in both the tree and the amino acid indices calculations
mattmp = mat[,intersect(tre.phylo$tip.label, colnames(mat))]
dend = keep.tip(tre.phylo, intersect(tre.phylo$tip.label,colnames(mattmp)))

# clean up and reorder tree
dend$node.label = NULL # Need to remove node labels
dend = midpoint.root(dend)
dend = ladderize(dend)
dend = drop.tip(dend,"1380381.3") #drop long tip (poor genome assembly)

# convert from 'phylo' to 'dendrogram'
dend1 = ReadDendrogram(textConnection(write.tree(dend)))
#dend1 = set(dend1, "branches_lwd", 5)

#--------------------------------------------------------------set up color scheme----------------------------------------------------------
#set up OGT colors
cols = rev(brewer.pal(11, "RdBu"))
pal = colorRampPalette(c("blue", "red")) # Define colour pallete
pal = colorRampPalette(cols) # Use the following line with RColorBrewer

#Rank variable for colour assignment
rowcols = ratk$OGT[match(labels(dend1),ratk$id.patric)]
rc.order = findInterval(rowcols, sort(rowcols))
rowcols = pal(length(sort(rowcols)))[rc.order]

#set up Sub-Clade Colors
Temp.Clust.colors = c("Clade A1" = "seagreen3","Clade A2" = "olivedrab4","Clade A3" = "chartreuse3", 
                      "Clade B1" = "darkorange2", "Clade B2" = "darkgoldenrod3", 
                      "Clade C" = "purple")

#------------------------------------------------Amino Acid Index Distribution by protein in each genome---------------------------------------
# set up amino acid indices by protein dataframe
PFAA = colind

# subset the dataframe by each amino acid indices and remove the bad data to see the bulk of the data for each
Arg.Lys_PF = PFAA[,c(1,7)]; rownames(Arg.Lys_PF) = NULL; Arg.Lys_PF$variable = "Arginine-Lysine Ratio";colnames(Arg.Lys_PF) = c("value","name","variable")
Arg.Lys_PF = subset(Arg.Lys_PF, value > 0)
Arg.Lys_PF = subset(Arg.Lys_PF, value < 0.75)

Acidic_PF = PFAA[,c(2,7)]; rownames(Acidic_PF) = NULL; Acidic_PF$variable = "Acidic Residue Proportion";colnames(Acidic_PF) = c("value","name","variable")
Acidic_PF = subset(Acidic_PF, value > 0.05 )
Acidic_PF = subset(Acidic_PF, value < 0.2 )

GRAVY_PF = PFAA[,c(3,7)]; rownames(GRAVY_PF) = NULL; GRAVY_PF$variable = "GRAVY";colnames(GRAVY_PF) = c("value","name","variable")
GRAVY_PF = subset(GRAVY_PF,value > -1 )
GRAVY_PF = subset(GRAVY_PF,value < 0.5)

Proline_PF = PFAA[,c(4,7)];rownames(Proline_PF) = NULL; Proline_PF$variable = "Proline Residue Proportion";colnames(Proline_PF) = c("value","name","variable")
Proline_PF = subset(Proline_PF,value > 0 )
Proline_PF = subset(Proline_PF,value < 0.08)

Aroma_PF = PFAA[,c(5,7)];rownames(Aroma_PF) = NULL; Aroma_PF$variable = "Aromaticity Index";colnames(Aroma_PF) = c("value","name","variable")
Aroma_PF = subset(Aroma_PF, value > 0 )
Aroma_PF = subset(Aroma_PF, value < 0.2)

Alipha_PF = PFAA[,c(6,7)];rownames(Alipha_PF) = NULL; Alipha_PF$variable = "Aliphatic Index";colnames(Alipha_PF) = c("value","name","variable")
Alipha_PF = subset(Alipha_PF, value > 50 )
Alipha_PF = subset(Alipha_PF, value < 150)

PFAA_melt = do.call("rbind", list(Arg.Lys_PF, Acidic_PF, GRAVY_PF,Proline_PF,Aroma_PF,Alipha_PF))

tiff("./Plots/Genome AA Distribution.tiff",width = 11,height = 8, units = "in", res = 500)
ggplot(PFAA_melt,aes(x = value, y = name)) + geom_density_ridges() + 
  facet_wrap(.~variable,nrow = 1, ncol = 6,scales = "free_x",labeller = label_wrap_gen(width = 20)) +
  labs(y = "Strain")+
  theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 15, face = "bold", hjust = 0.5),strip.text = element_text(size = 10, face = "bold"),
                           axis.text.y = element_text(size = 4.5),axis.text.x = element_text(size = 8)) 
dev.off()

#---------------------------------------------------------------Create volcano plot ------------------------------------------------------------
# set up local and global protein family names
famdat = rbind(plfdat,pgfdat)
colnames(famdat) = "description"

# set up phylogenetic tree for the volcano plot
dend.clust=1:length(labels(dend1))
names(dend.clust) = labels(dend1)
dend2 = set(dend1, "branches_lwd", 3)

pdf(file="./Plots/Volcano.pdf", width=10,height=6)  

for(j in 1:length(colsave)) {
  print(names(colsave)[j])
  
  #isolate dataframe from list
  mat = as.matrix(colsave[[j]])
  mat[is.infinite(mat)] = NA
  mat[is.nan(mat)] = NA
  
  # subset strains of interest with calculated or identified OGT
  common_ids = intersect(colnames(mat),ratk$id.patric)
  ratmat = mat[,common_ids]
  ogt = ratk[match(colnames(ratmat),ratk$id.patric),"OGT"]
  
  ratmat = ratmat[rowSums(!is.na(ratmat)) > 0,] 
  ratmat = ratmat[rowSums(!is.na(ratmat)) > 1,] 
  ratmat = ratmat[rowSums(!is.na(ratmat)) > 2,] 
  ratmat = ratmat[rowSums(!is.na(ratmat)) > 3,]
  
  # calculate the correlation values ( r & p.value) by protein family
  cormat = apply(ratmat,1,function(x) (cor.test(x,ogt)[c("estimate","p.value")]))
  corout = data.frame(matrix(unlist(cormat),ncol=2,byrow=T))
  rownames(corout) = names(cormat)
  colnames(corout) = c("r","p.value")
  corout = corout[which(abs(corout$r)<1),] #remove perfect matches (artifacts)
  
  # set up bon ferroni signficance threshold
  pmin = 0.05/nrow(corout)
  
  # plot volcano plot
  plot(corout$r,log10(corout$p.value),pch=19,col=rgb(0,0,1,0.3), main=names(colsave)[j])
  lines(x=c(-1,1),y=c(log10(pmin),log10(pmin)))
  
  # enrichment values of protein family quantiies that have positive or negative correlation to OGT
  print(sum(corout$p.value<pmin & corout$r > 0)/sum(corout$p.value<pmin & corout$r < 0)) # enrichment value of protein families that are significantly correlated
  print(sum(corout$r > 0)/sum(corout$r < 0)) # enrichment value of all protein families
  
  ## make correlation plots of only significantly correlated protein families
  corout = corout[order(corout$r,corout$p.value,decreasing=c(TRUE,FALSE)),]
  getgenes = rownames(corout)[corout$p.value < pmin]
  
  # stipulation if there is only one significant protein family
  if(length(getgenes)<1) next
  
  #get the protein family name
  genenames = famdat[getgenes,"description"]
  
  # set up color scheme for which strain do the protein family belong to
  colpick = dend.clust[colnames(ratmat)]
  names(colpick) = colnames(ratmat)
  cols = rainbow_hcl(length(dend.clust),l=50,c=100)[colpick]
  
  # loop to plot correlation of significant protein families
  for(i in 1:length(getgenes)) {
    
    k = getgenes[i]
    cols[is.na(cols)] = "#888888"
    
    par(mfcol=c(1,2),mar = c(4,4,4,4))
    plot(ogt,ratmat[k,],
         xlab = "Optimal Growth Temperature", ylab=names(colsave)[j],
         main=(paste0(names(colsave)[j],"\n",k,"\n",genenames[i],"\n","r = ",round(corout[k,"r"],3),", p = 10^",round(log10(corout[k,"p.value"]),3))),
         cex.main=0.7, pch=19, col=cols)
    abline(lm(ratmat[k,]~ogt))
    
    plot(dend2,horiz=T,xaxt="n",yaxt="n",leaflab="none")
    dend.cols = rainbow_hcl(length(dend.clust),l=50,c=100)
    dend.cols[!(dend.clust %in% colpick[!is.na(ratmat[k,])])] = NA
    dend.cols = cbind(dend.cols, rowcols)
    colored_bars(dend.cols,horiz=T,rowLabels = NA, add=T)
  }
}
dev.off()

# ----------------------------------------Principle Component Analysis of Protein Families by OGT and Sub-Clade----------------------------------------
Temp.Clust.ID$Temp_colors = Temp.Clust.colors[match(Temp.Clust.ID$ClusterName, names(Temp.Clust.colors))]
Temp.Clust.ID[is.na(Temp.Clust.ID)] = "black"

# plot principle component analysis for all amino acid indices ( global and local protein families) by mean and length
for(i in 1:length(colsave)) {
  
  #get dataframe
  cs = colsave[[i]]
  cs = as.matrix(cs)
  cs[is.nan(cs)] = NA
  
  #limit the dataframe to the core protein families
  cs = cs[rowSums(!is.na(cs)) > 0.8*ncol(cs),]
  
  #if na remains convert it to the average amino acid index
  cs[is.na(cs)] = mean(cs,na.rm=T)
  
  #calculate pca
  pc = prcomp(t(cs),center=T)
  
  #set up names and color schemes
  csn = genomes$genome.genome_name[match(colnames(cs),genomes$genome.genome_id)]
  css = as.factor(Temp.Clust.ID$ClusterName[match(csn,Temp.Clust.ID$label)]) # as.factor
  plot(pc$x[,1],pc$x[,2],pch=19,col=css,main=names(colsave)[i])
}

# Example: plot pca according to phylogeny ( Aliphatic Index, global protein)
cs = colsave[[13]] 
cs = as.matrix(cs)
cs[is.nan(cs)] = NA
cs = cs[rowSums(!is.na(cs)) > 0.8*ncol(cs),]
cs[is.na(cs)] = mean(cs,na.rm=T)
pc = prcomp(t(cs),center=T)
csn = genomes$genome.genome_name[match(colnames(cs),genomes$genome.genome_id)]
csc = as.factor(Temp.Clust.ID$Temp_colors[match(csn,Temp.Clust.ID$label)]) # get the cluster color
cst = as.factor(Temp.Clust.ID$ClusterName[match(csn,Temp.Clust.ID$label)]) # get the clustertype
cso = as.numeric(ratk$OGT[match(csn,ratk$strn)])

# set up color scheme for Sub_Genus Clusters and OGT 
rowcols.temp = cso
rc.order.temp = findInterval(rowcols.temp, sort(rowcols.temp))
rowcols.temp = pal(length(sort(rowcols.temp)))[rc.order.temp]

data = as.data.frame(pc$x)
data$ClusterName = cst
data$OGT = cso
data$OGT_color = rowcols.temp

tiff("./Plots/PCA by Cluster.tiff",width = 8,height = 5.5, units = "in", res = 1200)
ggplot(data,aes(PC1,PC2,color = ClusterName)) + 
  geom_point(size = 3) + 
  scale_color_manual(name = "Clade",labels = c("Clade A1","Clade A2","Clade A3","Clade B1","Colwellia BRX10.3 & Other Clade B (Clade B2)","Clade C","Other Colwelliaceae"),
                     values = c("seagreen3","olivedrab4","chartreuse3","darkorange2","darkgoldenrod3","purple","black")) +
  guides(fill = guide_legend(title = "Clades")) +
  ggtitle("(B)  PCA of Aliphatic Index by Protein Family Colored by Sub-Clades") + 
  theme_classic(base_line_size = 1, base_size = 11) + 
  theme(panel.grid.major = element_line(color = "lightgrey", size = 0.1),legend.position = "bottom", plot.title = element_text(hjust = 0.5))
dev.off()

tiff("./Plots/PCA by OGT.tiff",width = 8,height = 5.5, units = "in", res = 1200)
ggplot(data,aes(PC1,PC2, colour = OGT)) + 
  geom_point(color = as.character(data$OGT_color), size = 3) + 
  ggtitle("(A)  PCA of Aliphatic Index by Protein Family Colored by OGT") + 
  theme_classic(base_line_size = 1, base_size = 11) + 
  theme(panel.grid.major = element_line(color = "lightgrey", size = 0.1), plot.title = element_text(hjust = 0.5))
dev.off()


#------------------------------------------------------------------ pangenome vs. core genome plot -----------------------------------------------------
# use global protein families dataframe
mat = as.matrix(colsave[[1]])
mat[is.infinite(mat)] = NA
mat[is.nan(mat)] = NA
mat[mat==0] = NA

#subset data
matsmall = mat[,colnames(mattmp)]
matcols = colnames(matsmall)

#calculate gene family frequency spectrum
mat = !is.na(as.matrix(colsave[[1]])) # use global families instead of local families
mat = mat + 0
mat = mat[,matcols]

Gk <- f.getspectrum(mat)

genomesize <- median(colSums(mat>0)) # mean genome size measured in gene families
ng <- dim(mat)[2] #number of genomes

# Calculate 100 permutations each of the pangenome and core genome
perm.pangenome <- f.pangenome(mat,100)
perm.core <- f.core(mat,100)

# Calculate the exact mean pan and core genome curves
# from the gene frequency spectrum G(k)
mean.pangenome <- f.meanpancore(Gk)$pan
mean.core <- f.meanpancore(Gk)$core
pancore <- c(mean.pangenome,mean.core)

# Calculate the RMS value for the permutations
rms.perm <- mean(f.rms(c(mean.pangenome,mean.core),rbind(perm.pangenome,perm.core)))

tiff(file="./Plots/pancore.tiff",width = 12,height = 8, units = "in", res = 500)
par(mar =c(5,5,1,1))
# Prepare a new plot window
plot(1:ng,xlim=c(1,ng),ylim=c(0.9*min(mean.core), 1.1*max(mean.pangenome)),log="",
     xlab="Genomes added", ylab="PATRIC Global Protein Family",pch='',cex.lab = 1.5,cex.axis = 1.3)

# Plot polygons outlining permutations
polygon(c(1:ng, rev(1:ng)), c(apply(perm.core,1,min), rev(apply(perm.core,1,max))), col="gray88",border=NA)
polygon(c(1:ng, rev(1:ng)), c(apply(perm.pangenome,1,min), rev(apply(perm.pangenome,1,max))), col="gray88",border=NA)

# Add the mean pan and core genome curves to the plot
points(1:ng,mean.pangenome,type='l')
points(1:ng,mean.core,type='l')

dev.off()

# -------------------------------------Compare the percent shared & unshared to the pairwise OGT differences------------------------------
#GLOBAL and LOCAL FAMILIES differential gene abundance
fams=c("global")

protfam=fams[1]
colfams = as.matrix(colsave[[1]]>0)

coldiff = matrix(nrow=ncol(colfams),ncol=ncol(colfams))
colnames(coldiff) = colnames(colfams)
rownames(coldiff) = colnames(colfams)

colsim = coldiff
coluniq = coldiff

for(j in 1:nrow(coldiff)) {
  for(k in 1:nrow(coldiff)) {
    dfs = rownames(colfams)[ setdiff( which(colfams[,j]), which(colfams[,k])) ]
    sim = rownames(colfams)[ intersect( which(colfams[,j]), which(colfams[,k])) ]
    # unq = rownames(colfams)[ intersect( which(colfams[,j]), which(colfams[,k])) & !]
    if(!is_empty(dfs)) coldiff[j,k] = length(dfs)
    if(!is_empty(sim)) colsim[j,k] = length(sim)
  }
}

# sum matrices
coldiffs = coldiff + t(coldiff)
keepcol = intersect(labels(dend1),colnames(coldiffs))

coldiffs = coldiffs[keepcol,keepcol]
colsim = colsim[keepcol, keepcol]

coldiffs = coldiffs * lower.tri(coldiffs)
colsim = colsim * lower.tri(colsim)

coldiffs = coldiffs + t(colsim)

genomesizes = colSums(colfams,na.rm=T)[rownames(coldiffs)]
sizemat = outer(genomesizes,genomesizes,'+')
rownames(sizemat) = genomes$genome.genome_name[match(rownames(sizemat),genomes$genome.genome_id)]
colnames(sizemat) = rownames(sizemat)

rownames(coldiffs) = genomes$genome.genome_name[match(rownames(coldiffs),genomes$genome.genome_id)]
colnames(coldiffs) = rownames(coldiffs)

OGT_diff = outer(ratk$OGT,ratk$OGT,`-`)
rownames(OGT_diff) = rownames(ratk)
colnames(OGT_diff) = rownames(ratk)

# ### distances
# dist.tre = cophenetic.phylo(tre.phylo)
# dist.plot = dist.tre[rownames(coldiffs),rownames(coldiffs)]
# # 
common_ids = intersect(rownames(OGT_diff),rownames(coldiffs))
coldiffs = coldiffs[rownames(coldiffs) %in% common_ids,colnames(coldiffs) %in% common_ids]
OGT_diff = OGT_diff[rownames(OGT_diff) %in% common_ids,colnames(OGT_diff) %in% common_ids]
sizemat = sizemat[rownames(sizemat) %in% common_ids,colnames(sizemat) %in% common_ids]

coldiffs = coldiffs[rownames(OGT_diff),colnames(OGT_diff)]
sizemat = sizemat[rownames(OGT_diff),colnames(OGT_diff)]
OGT_diff = abs(OGT_diff)

# shared fams
x = OGT_diff*upper.tri(OGT_diff)
x[x==0] = NA
y = coldiffs*upper.tri(coldiffs)
y[y==0] = NA
z = 100*(coldiffs*upper.tri(coldiffs))/(sizemat*upper.tri(sizemat))

# non-shared fams
x1 = OGT_diff*lower.tri(OGT_diff)
x1[x1==0] = NA
y1 = coldiffs*lower.tri(coldiffs)
y1[y==0] = NA
z1 = 100*(coldiffs*lower.tri(coldiffs))/(sizemat*lower.tri(sizemat))

tiff("./Plots/UnShared PFG Percent vs.OGT Differences.tiff",width = 8,height = 11, units = "in", res = 500)
plot(x1,z1,pch=19,col=rgb(0,0,1,0.3), xlab="Pairwise OGT differences", ylab="Percent Unshared Protein Families")
dev.off()

tiff("./Plots/Shared PFG Percent vs.OGT Differences.tiff",width = 8,height = 11, units = "in", res = 500)
plot(x,z,pch=19,col=rgb(0,0,1,0.3), xlab="Pairwise OGT Differences", ylab="Percent Shared Protein Families")
dev.off()

