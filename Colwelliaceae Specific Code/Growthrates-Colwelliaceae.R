### Calculate bacterial growth rates from plate reader data and filters them according to
### the absolute median deviation

library(growthrates) # calculates growth rates
library(reshape2) # melt/cast
library(readxl) # read XLS files
library(ggforce)
library(stats)

#-------------------------------------load files ---------------------------------------------------
source("./Plain Code/Growthrates Calc.R")
#--------------------------------- format strain names ---------------------------------------------
grd = as.data.frame(grd)
grd$strain = as.character(grd$strain)
grd$`mumax/day` = as.numeric(grd$`mumax/day`)

for(i in 1:nrow(grd)){
  strain = as.character(grd[i,1])
  grd[i,1] = paste0("Colwellia ",strain)# if colwellia is not in the string
}

# -----------------------remove outliers using median absolute deviation------------------------------
# robust outlier detection the goal of robust statistics is to find a fit that is close
# to the fit we would have found without the outliers.
# We can then identify the outliers by their large deviation from that robust fit.
temp_strn_ls = split(grd,list(grd$strain,grd$temp))
mad_out_detect <- function(x, thres, na.rm = TRUE) {
  return(ifelse(abs(x - median(x, na.rm = na.rm)) >= thres * mad(x, na.rm = na.rm),TRUE,FALSE))
} # false is a non outlier and true is an outlier

out_test_df = data.frame()
for (i in 1:length(temp_strn_ls)){
  data = temp_strn_ls[[i]]
  temp = as.character(paste(unique(data$temp)))
  strain = as.character(paste(unique(data$strain)))
  rep = data$replicate
  r2 = data$r2
  data = data$`mumax/day`
  
  mad_out_7 = mad_out_detect(data,7)
  
  out_test = data.frame(strain = strain,replicate = rep,temp = temp,mumax = data,r2 = r2,mad_out_7 = mad_out_7)
  out_test_df = rbind(out_test_df,out_test)
}
grd = subset.data.frame(out_test_df,mad_out_7 == F)
grd$mad_out_7 = NULL

#--- calculate average and standard deviation of the growth rates by strain and temperature tested ---------
# calculate growth rate average by temperature
GR_temp = grd %>% group_by(temp) %>% dplyr::summarize(growthrate = paste0(signif(mean(mumax),2)," ± ",signif(sd(mumax),2)))

#create supplementary table for growth rate average by temperature and strain 
GR_strain = grd %>% group_by(temp,strain) %>% dplyr::summarize(growthrate = paste0(signif(mean(mumax),2)," ± ",signif(sd(mumax),2)))
GR_strain = as.data.frame(GR_strain)
GR_strain = recast(GR_strain,strain~temp, id.var = c("strain","temp"))

# ----------------------------------------- write files----------------------------------------
write.table(grd,"./Data/growthrates.tsv",sep = "\t",quote = F,row.names = F)
write.csv(GR_temp,"./Plots/growthrates_temp.csv",row.names = T)
write.csv(GR_strain,"./Plots/growthrates_strain.csv",row.names = T)

