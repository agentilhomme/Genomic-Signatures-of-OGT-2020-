library(readxl) # read XLS files
source("./Plain Code/Strain OGT Ratkowsky Calc.R")

# ----------------------------------------------------load files -----------------------------------------------------------------------
lit_col = read_excel("./References/Literature Colwellia OGT.xlsx", sheet = "OGT Data") # other OGT data from the literature

#------------------------------------------------- Format strain names-------------------------------------------------------------------
# format dataframe
Final_OGT = as.data.frame(Final_OGT)
Final_OGT2 = Final_OGT
Final_OGT2$strn = as.character(Final_OGT2$strn)

# set strain names of literature based growth rate Colwellia to actually strain names (PATRIC format)
Final_OGT2$strn = gsub("Colwellia psychrerythraea","Colwellia psychrerythraea ACAM 605",Final_OGT2$strn) 
Final_OGT2$strn = gsub("Colwellia hornerae","Colwellia hornerae strain ACAM 607",Final_OGT2$strn)
Final_OGT2$strn = gsub("Colwellia piezophila","Colwellia piezophila ATCC BAA-637",Final_OGT2$strn) 
Final_OGT2$strn = gsub("Colwellia demingiae","Colwellia demingiae strain ACAM 459",Final_OGT2$strn) 

#-------------------------- Merge dataframe with OGTs of other Colwelliaceae found in the literature ------------------------------------
Final_OGT2 = rbind(Final_OGT2,lit_col)

#----------------------------------------write table of final Colwelliaceae OGT----------------------------------------------------------
write.table(Final_OGT2,"./Data/Strain Ratk OGT Final.tsv",sep = "\t",quote = F,row.names = F)
