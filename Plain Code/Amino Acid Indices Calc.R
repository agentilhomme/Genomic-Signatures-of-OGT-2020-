library(phylogram) # read.dendrogram()
library(ape) # read.tree()
library(Rcpi) # readFASTA()
library(protr) # protcheck()
library(stringr) # str_remove()
library(seqinr) # c2s()
library(DECIPHER) # ReadDendrogram()
library(reshape2) # melt()

#---------------------------------------------------------load files---------------------------------------------------------------------
# load PATRIC phylogenetic files and metadata
treefile="./Data/Colwelliaceae Rooted Tree.nwk"
genomes = read.table("./Data/Colwelliaceae Rooted PATRIC Data.tsv",fill=T,head=T,sep="\t",colClasses="character")
colfeat = read.table("./Data/Colwelliaceae PF Features.tsv",fill = T,head = T,sep = "\t",colClasses = "character") # patric protein features

# create vector list of genome faa files from PATRIC
orgfiles = list.files("./Genomes", pattern = "*.faa") 

# load amino acid information 
aminos = read.table("./References/Amino Acid Info.txt", sep = "\t", header = T, stringsAsFactors = F)

#-------------------------------------------------------set up files----------------------------------------------------------------------- 
#read phylogenetic tree from PATRIC in dendrogram and phylogenetic format
tre.dend = ReadDendrogram(treefile)
tre.phylo = read.tree(treefile)
tre.phylo = drop.tip(tre.phylo,"1380381.3") #drop long tip (poor genome assembly)

rownames(genomes) = genomes$genome.genome_id
rownames(aminos) = aminos$abb1
aminos =  aminos[order(aminos$amino.acid),] # set the amino acid alphabetically

# get scientific names for PATRIC ID of the phylogenetic tree
tree_ids = intersect(tre.phylo$tip.label,rownames(genomes))
tre.phylo = keep.tip(tre.phylo,tree_ids)

tree_files = paste(tree_ids,".faa",sep = "")
orgfiles = intersect(tree_files,orgfiles)
genomes$genome.genome_id.1 = paste(genomes$genome.genome_id,".faa", sep = "")

#set up protein features
rownames(colfeat) = colfeat$feature.patric_id
colfeat = colfeat[, -2]
colnames(colfeat)[1] = "genome.genome_id"

#-----------------------------------------------------indice functions----------------------------------------------------------------------
arg_lys_ratio = function(x) {
  x["R",]/ (x["R",] + x["K",] )
}
aliphatic_index = function(x) {
  #normalize the protein to 1
  mole_fraction = sweep(x, MARGIN = 1, STATS = aminos$mw, FUN = `*`)
  total_mole_fraction = colSums(mole_fraction)
  mole_percent = (mole_fraction/total_mole_fraction)*100	
  mole_percent["A",] + 2.9*mole_percent["V",] + 3.9*(mole_percent["I",] + mole_percent["L",])
}
aromaticity = function(x) {
  colSums(sweep(x, MARGIN = 1, STATS = aminos$aromatic, FUN = `*`))	
}
acidic_residue = function(x) {
  x["D",] + x["E",]
}
gravy = function(x) {
  colSums(sweep(x, MARGIN = 1, STATS = aminos$hydropathicity, FUN = `*`))	
}
proline_residue = function(x) {
  x["P",]
}

#---------------------------------------------------Calculate the amino acid indices for the entire genome------------------------------------------
output = data.frame()
for( file in orgfiles) {
  
  #read fasta file
  fasta_in = readFASTA(paste0("./Genomes/",file))# set the directory and read the fasta file 
  file_name = paste(file)
  file_name = str_remove(file_name,".faa")
  fasta_good = fasta_in[(sapply(fasta_in, protcheck))] # check if all the protein sequences are valid
  fasta_good = c2s(paste(fasta_good,sep = "",collapse = ""))
  aac = data.frame(extractAAC(fasta_good)) # gives amino acid composition in each amino acid sequence in genome, gives amino acid frequencies
  
  #calculate amino acid indice and input them into dataframe
  output[file_name,1] = arg_lys_ratio(aac)
  output[file_name,2] = aliphatic_index(aac)
  output[file_name,3] = aromaticity(aac)
  output[file_name,4] = acidic_residue(aac)
  output[file_name,5] = gravy(aac)
  output[file_name,6] = proline_residue(aac)
}
# add column titles
colnames(output) <- c("Arg-Lys Ratio","Aliphatic Index","Aromaticity Index","Acidic Residue Proportion","GRAVY","Proline Residue Proportion")

#---------------------------------------------Calculate the amino acid indices for each protein in each genome---------------------------------------
IndicesData <- data.frame()
for(file in orgfiles){
  
  # read fasta file
  fasta_in = readFASTA(paste0("./Genomes/",file))
  fasta_good = fasta_in[(sapply(fasta_in, protcheck))] # check if all the protein sequences are valid
  aac = sapply(fasta_good, extractAAC) # gives amino acid composition in each amino acid sequence in genome, gives amino acid frequencies
  
  #get scientific name
  print(file)
  name = genomes$genome.genome_name[match(file,genomes$genome.genome_id.1)]
  print(name)
  
  # calculate amino acid ratios by protein family
  ArgLysRatio <- as.data.frame(arg_lys_ratio(aac))
  AcidicResRatio <- as.data.frame(acidic_residue(aac)) 
  GravyRatio <- as.data.frame(gravy(aac)) 
  ProlResRatio <- as.data.frame(proline_residue(aac)) 
  AromaRatio <- as.data.frame(aromaticity(aac)) 
  AliphaRatio <- as.data.frame(aliphatic_index(aac))
  
  # bind the data together into one dataframe
  IndicesTotal <- cbind(ArgLysRatio,AcidicResRatio,GravyRatio,ProlResRatio,AromaRatio,AliphaRatio)
  colnames(IndicesTotal) = c("Arg-Lys Ratio","Acidic Residue Proportion","GRAVY","Proline Residue Proportion","Aromaticity Index","Aliphatic Index")
  IndicesTotal$name = name 
  IndicesData <- rbind(IndicesData,IndicesTotal)
}

# ---------------------Calculate the average and length of the amino acid indices by local and global protein family within each genome---------------
colind = IndicesData

# limit datasets to common genomes
common_ids = intersect(rownames(colfeat), rownames(colind))
colfeat = colfeat[common_ids, ]
colind = colind[common_ids, ]

#merge protein features and protein amino acid indicees
colout = data.frame(colfeat, colind)

# initialize list
colsave = list()

# for each index
for (var in 4:9) {
  #for each statistical summary
  for (fun in c("mean","length")) {
    if(var > 4 & fun=="length") next
    #for each type of protein family
    for (fam in 2:3) {
      
      #subset data
      coltmp = colout[, c(1, fam, var)]
      coltmp = coltmp[which(coltmp[,2] != ""), ]
      savename = paste0(colnames(colout)[fam], ".", colnames(colout)[var], ".", fun)
      print(savename)
      
      # summarize global families
      if (fam == 2) {
        
        tmp.mat = reshape2::dcast(coltmp,formula = feature.pgfam_id ~ genome.genome_id,
                                  fun.aggregate = eval(parse(text=fun)))
        
        rownames(tmp.mat) = tmp.mat[,1]
        colsave[[savename]] = tmp.mat[,-1]
        
      }
      
      #summarize local families
      if (fam == 3) {
        tmp.mat = reshape2::dcast(coltmp,formula = feature.plfam_id ~ genome.genome_id,
                                  fun.aggregate = eval(parse(text=fun)))
        
        rownames(tmp.mat) = tmp.mat[,1]
        colsave[[savename]] = tmp.mat[,-1]
      }
    }
  }
}

# ------------------------------------------------------------------write files -----------------------------------------------------------------------
write.table(output,file="./Data/WG Indice Data.tsv",quote=F, sep = "\t",row.names = T) # amino acid calculations by genome
write.table(IndicesData,file = "./Data/Protein Indices Data.tsv",row.names = T, sep = '\t') # amino acid calculations by protein
write.csv(colsave[["feature.pgfam_id.Arg.Lys.Ratio.length"]],"./Data/Gene Prescence.csv",row.names = T) # presence or absence of global protein families
saveRDS(colsave,file="./Data/Colwelliaceae PF AA Calc.RDS") # save the protein family amino acid calculations in a list format

