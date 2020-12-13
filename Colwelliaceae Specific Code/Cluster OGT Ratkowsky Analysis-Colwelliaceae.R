library(scales)
library(dplyr)
library(reshape2)
# ----------------------------------------------------read in files---------------------------------------------------------
ratk_TC = read.table("./Data/Clade Ratk Params.tsv",fill=T,head=T,sep="\t",colClasses="character")
ratk_TC_med = read.table("./Data/Clade Ratk Med Params.tsv",fill=T,head=T,sep="\t",colClasses="character")
ratk_TC_sd = read.table("./Data/Clade Ratk SD Params.tsv",fill=T,head=T,sep="\t",colClasses="character")
ratk = read.table("./Data/Strain Ratk OGT Outputs Clean.tsv",sep = "\t",header = T,colClasses = "character")
Temp.Clust.ID = read.table("./Data/Sub-Clade ID.tsv",fill=T,head=T,sep="\t",colClasses="character")
GR = read.table("./Data/growthrates.tsv",fill=T,head=T,sep="\t",colClasses="character")

#-------------------------------------set up color schemes ------------------------------------------------------------
Temp.Clust.colors = c("Clade A1" = "seagreen3","Clade A2" = "olivedrab4","Clade A3" = "chartreuse3", 
                      "Clade B1" = "darkorange2", "Clade B2" = "darkgoldenrod3", 
                      "Clade C" = "purple")
#----------------------------------------------------set up functions ------------------------------------------------------------------------
# load functions
rat83 = function(x, b, T1, cc, T2) {
  rT = b * (x - T1) * (1 - exp(cc * (x - T2)))
  return(rT)
}

rat82a = function(x, b, T0) {
  r = b*(x - T0)
}

#------------------------------------------------------set up  files -------------------------------------------------------------------------
colnames(Temp.Clust.ID)[1] = "strain"

#set up ratk output file
cols.num_topt = colnames(ratk)[2:7]
ratk[cols.num_topt] = sapply(ratk[cols.num_topt],as.numeric)

# set up growth rate file
cols.num_GR = colnames(GR)[2:5]
GR[cols.num_GR] = sapply(GR[cols.num_GR],as.numeric)

# set up the ratkowsky fit by temperature cluster 
init.plot = as.data.frame(ratk_TC_med)
rownames(init.plot) = init.plot[,1]
cols.num = colnames(init.plot)[2:6]
init.plot[cols.num] = sapply(init.plot[cols.num],as.numeric)
ratk_TC_sd[cols.num] = sapply(ratk_TC_sd[cols.num],as.numeric)

# add color column to dataframe
init.plot$col = Temp.Clust.colors[match(rownames(init.plot),names(Temp.Clust.colors))]
ratk_TC$col = Temp.Clust.colors[match(ratk_TC$Clade,names(Temp.Clust.colors))]
grd = merge.data.frame(Temp.Clust.ID,GR,by = "strain", all = F)

#set up growth rates data to plot
colnames(grd) = c("strain_original","strain","rep","temp","rate","r2") #'strain' is actually 'clade'
grd$strain_original = as.factor(grd$strain_original)
grd$temp = as.numeric(grd$temp) + 273
grd$rep = as.factor(grd$rep)
grd$rate = as.numeric(grd$rate)
grd = na.omit(grd)
grd$sqperday = sqrt(grd$rate) 

#----------------------------------------plot Clade ratkowsky fits by Individual Clade-------------------------------------------------------
tiff("./Plots/Clade Ratk Individual.tiff",width = 11,height = 8, units = "in", res = 500)
op = par(mfcol=c(3,2), mar = c(3,3,2,1),oma = c(2,3,0,0))
for(i in 1:length(unique(ratk_TC$Clade))) {

  #set up clade name and temperature range
  clade = sort(unique(ratk_TC$Clade))[i] 
  print(as.character(clade))
  xs = seq(268,298,0.1)
  
  # subset the dataframes of growth rates and ratkowksy parameters by Clade
  sgrd = grd[grd$strain == clade,c("temp","sqperday")]
  sgrd = rbind(sgrd,aggregate(sgrd,by = list(sgrd$temp), FUN=mean)[,c("temp","sqperday")])
  Clade_df = ratk_TC[ratk_TC$Clade == clade,c("b","Tmin","c","Tmax","topt","col")]
  
  # set up clade parameters to numeric and melt the dataframe
  cols.num = colnames(Clade_df)[1:5]
  Clade_df[cols.num] = sapply(Clade_df[cols.num],as.numeric)
  Clade_df_melt = melt(Clade_df,id.vars = "col")
  
  # identify the median, mean,sd, each parameter for each clade and make a dataframe of each
  Clade_df_med = Clade_df_melt %>% group_by(col,variable) %>% dplyr::summarize(med = median(value))
  Clade_df_med = recast(Clade_df_med,col~variable)
  Clade_df_mean = Clade_df_melt %>% group_by(col,variable) %>% dplyr::summarize(med = mean(value))
  Clade_df_mean = recast(Clade_df_mean,col~variable)
  Clade_df_sd = Clade_df_melt %>% group_by(col,variable) %>% dplyr::summarize(med = sd(value))
  Clade_df_sd = recast(Clade_df_sd,col~variable)
 
  #Calculate the median fit of the clade 
  rt.fit_med = rat83(xs, Clade_df_med$b, Clade_df_med$Tmin, Clade_df_med$c, Clade_df_med$Tmax)
  rt.fit_med[rt.fit_med < 0] = NA
  
  # plot original growth rates, individual ratkowsky fits, and median ratkowski fits
  plot(sgrd$temp - 273, sgrd$sqperday^2, pch=19, col=Clade_df_med$col,cex =0.5, # remove the ^2 to look at the shape, not for the paper
       xlim=c(-1,17), 
       ylim=c(0,3),
       cex.axis=1.3)
  for ( j in 1:nrow(Clade_df)){
    rt.fit = rat83(xs, Clade_df$b[j], Clade_df$Tmin[j], Clade_df$c[j], Clade_df$Tmax[j])
    rt.fit[rt.fit < 0] = NA
    lines(xs-273, rt.fit^2, type="l",lwd=0.5, col = alpha(Clade_df_med$col, 0.03))
  }
  lines(xs-273, rt.fit_med^2, type="l",lwd=3, col = Clade_df_med$col)
  
  # Find the OGT from the maximum growth rate of the median line
  topt_med = round(xs[which.max(rt.fit_med)],2)
  GR_med = max(rt.fit_med,na.rm = T)^2
  
  # plot titles,mean and sd of parameters 
  title(paste0(clade,"   OGT = ",round(topt_med-273,2),"°C ± ",round(Clade_df_sd$topt,1)),cex.main=1.2)
  b_text = paste0("b =  ",round(Clade_df_med$b,2)," ± ",round(Clade_df_sd$b,2))
  c_text = paste0("c =  ",round(Clade_df_med$c,2)," ± ",round(Clade_df_sd$c,2))
  Tmin_text = paste0("Tmin =  ",round(Clade_df_med$Tmin,2),"K ± ",round(Clade_df_sd$Tmin,2))
  Tmax_text = paste0("Tmax =  ",round(Clade_df_med$Tmax,2),"K ± ",round(Clade_df_sd$Tmax,2))
  text(x = 12, y = 2.2,paste(b_text,c_text,Tmin_text,Tmax_text,sep = "\n") ,cex=1,adj = c(0,0))
}
title(xlab = "Temperature (°C)",ylab = expression("Average Growth Rate" ~ (d^{-1})),outer = T, line = 1, cex.lab = 1.8)
par(op)
dev.off()

#----------------------------------------Plot Clade ratkowsky fits for Clades of Interest (Clade A2 & Clade A3)-------------------------------------------------------
tiff("./Plots/Clade Ratk A2 & A3.tiff",width = 11,height = 8, units = "in", res = 500)
par(mfrow=c(1,1),mar = c(4,5,2,1))
#set up empty plot
plot(NA,NA,
     xlim=c(-1,17),
     ylim=c(0,3),
     xlab = "Temperature (°C)", ylab = expression("Average Growth Rate" ~ (d^{-1})),
     cex.axis=1.3, cex.lab=1.8)
for(i in 2:3) {
  
  #set up clade name and temperature range
  clade = sort(unique(ratk_TC$Clade))[i] #i
  print(as.character(clade))
  xs = seq(268,298,0.1)
  
  # subset the dataframes of growth rates and ratkowksy parameters by Clade
  sgrd = grd[grd$strain == clade,c("temp","sqperday")]
  sgrd = rbind(sgrd,aggregate(sgrd,by = list(sgrd$temp), FUN=mean)[,c("temp","sqperday")])
  Clade_df = ratk_TC[ratk_TC$Clade == clade,c("b","Tmin","c","Tmax","topt","col")]
  
  # set up clade parameters to numeric and melt the dataframe
  cols.num = colnames(Clade_df)[1:5]
  Clade_df[cols.num] = sapply(Clade_df[cols.num],as.numeric)
  Clade_df_melt = melt(Clade_df,id.vars = "col")
  
  # identify the median,each clade and make a dataframe of each
  Clade_df_med = Clade_df_melt %>% group_by(col,variable) %>% dplyr::summarize(med = median(value))
  Clade_df_med = recast(Clade_df_med,col~variable)
  
  #Calculate the median fit of the clade 
  rt.fit_med = rat83(xs, Clade_df_med$b, Clade_df_med$Tmin, Clade_df_med$c, Clade_df_med$Tmax)
  rt.fit_med[rt.fit_med < 0] = NA
  
  if ( i == 2){
    for ( j in 1:nrow(Clade_df)){
      rt.fit = rat83(xs, Clade_df$b[j], Clade_df$Tmin[j], Clade_df$c[j], Clade_df$Tmax[j])
      rt.fit[rt.fit < 0] = NA
      lines(xs-273, rt.fit^2, type="l",lwd=0.5, col = alpha(Clade_df_med$col, 0.3))
    }
  } else {
    for ( j in 1:nrow(Clade_df)){
      rt.fit = rat83(xs, Clade_df$b[j], Clade_df$Tmin[j], Clade_df$c[j], Clade_df$Tmax[j])
      rt.fit[rt.fit < 0] = NA
      lines(xs-273, rt.fit^2, type="l",lwd=0.5, col = alpha(Clade_df_med$col, 0.03))
    }
    
  }
  lines(xs-273, rt.fit_med^2, type="l",lwd=5) #,col = Clade_df_med$col)
  
}
legend("topleft",legend = c("Clade A2", "Clade A3"), bty = "n",
       lwd = 4, col =Temp.Clust.colors[2:3],cex = 1)
dev.off()

#----------------------------------------Plot Overlapped Clade ratkowsky fits -------------------------------------------------------
tiff("./Plots/Clade Ratk All.tiff",width = 11,height = 8, units = "in", res = 500)
plot(NA,NA,
     xlim=c(-1,17),
     ylim=c(0,3),
     xlab = "Temperature (°C)", ylab = expression("Average Growth Rate" ~ (d^{-1})),
     cex.axis=1.3, cex.lab=1.8)

for(i in 1:nrow(init.plot)) {
  #i = 1
  strn = rownames(init.plot)[i]
  print(as.character(strn))

  xs = seq(0,1000,0.1)
  rt.fit = rat83(xs, init.plot$b[i], init.plot$Tmin[i], init.plot$c[i], init.plot$Tmax[i])
  rt.fit[rt.fit < 0] = NA
  
  xs = xs - 273
  lines(xs, rt.fit^2, type="l", col=init.plot$col[i],lwd=6, lty = i)
}
legend("topright",legend = rownames(init.plot), bty = "n",
       lwd = 4, col = init.plot$col,lty = 1:6 ,cex = 1)

dev.off()
