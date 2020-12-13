### Calculate bacterial growth rates from plate reader data and filters them according to
### the absolute median deviation

library(growthrates) # calculates growth rates
library(reshape2) # melt/cast
library(readxl) # read XLS files
library(ggforce)
library(stats)
# -----------------------------read in files---------------------------------------------
dat.in = read_excel("./Data/OD_Raw_Data.xlsx", sheet = "Data")

#------------------------------ set up file ---------------------------------------------
dat.in$`Start Time` = rep(dat.in[[1,3]],nrow(dat.in))
dat.in$`TD` = as.numeric(difftime(dat.in$`Time Taken`,dat.in$`Start Time`,units = "hours"))
dat.in$`Time Taken` = NULL
dat.in$`Start Time` = NULL

# cut out bad data
dat.in = dat.in[!(dat.in$Temperature == 4 & dat.in$`Plate Number`== 3 & dat.in$TD > 40),]
dat.in = dat.in[!(dat.in$Temperature == 11 & dat.in$`Plate Number` == 8 & dat.in$TD > 40),]
dat.in = dat.in[!(dat.in$Temperature == 17 & dat.in$TD > 40),]
dat.in$`Plate Number` = NULL

# --------------------------------convert OD 600 to cells/ml ------------------------------
#A.Huston Colwellia equation: log10(cells/ml) = 9.9 -4/(1+(OD600/0.14)^0.8)
OD600_to_cellconc = function(value){
  cellperml = 10^(9.9 -4/(1+(value/0.14)^0.8))
}
test_dat.in = data.frame(apply( dat.in[,grep("OD",colnames(dat.in))],2,OD600_to_cellconc))
dat.in[,4:11] = test_dat.in
colnames(dat.in)[4:11] = 1:8

#--------------------------------calculate maximum growth rates ---------------------------
# melt data
dat.all = reshape2::melt(dat.in,id = c("TD","Temperature","Time Point","Strain"))
colnames(dat.all) = c("time","temp","timepoint","strain","replicate","value")
dat.all = dat.all[,c("strain","replicate","temp","time","value")]

# fit splines to all data
many_spline_fits = all_splines(value ~ time | replicate + temp + strain,data = dat.all, spar = 0.5)
many_spline_res = results(many_spline_fits)
grd = many_spline_res[,c("strain","replicate","temp","mumax","r2")]
rownames(grd) = NULL

grd = as.data.frame(grd)
grd$strain = as.character(grd$strain)

#output data
grd = grd[,c("strain","replicate","temp","mumax","r2")]

rownames(grd) = NULL

#convert growthrates from cells/hour to cells/day
grd$mumax = grd$mumax*24
colnames(grd)[4] = "mumax/day"

#--- calculate average and standard deviation of the growth rates by strain and temperature tested ---------
# calculate growth rate average by temperature
GR_temp = grd %>% group_by(temp) %>% dplyr::summarize(growthrate = paste0(signif(mean(`mumax/day`),2)," ± ",signif(sd(`mumax/day`),2)))


#create supplementary table for growth rate average by temperature and strain 
GR_strain = grd %>% group_by(temp,strain) %>% dplyr::summarize(growthrate = paste0(signif(mean(`mumax/day`),2)," ± ",signif(sd(`mumax/day`),2)))
GR_strain = as.data.frame(GR_strain)
GR_strain = recast(GR_strain,strain~temp, id.var = c("strain","temp"))

# ----------------------------------------- write files-----------------------------------------------------
write.table(grd,"./Data/growthrates.tsv",sep = "\t",quote = F,row.names = F)
write.csv(GR_temp,"./Plots/growthrates_temp.csv",row.names = T)
write.csv(GR_strain,"./Plots/growthrates_strain.csv",row.names = T)
# -----------------------------------------plot data--------------------------------------------------------
pdf(file="./Plots/xyplot_strain_temps_cells.pdf", width=24, height=24)
xyplot(value ~ time|strain+as.factor(temp), data = dat.all, groups = replicate, pch = 16, cex = 0.5)
dev.off()


pdf(file="./Plots/many_spline_fits_cells.pdf",width=16,height=8)
par(mfrow = c(4, 8))
par(mar = c(2.5, 4, 2, 1))

plot(many_spline_fits,xlab = "Hours",ylab = "Cells/ml")
xyplot(mumax ~ temp|strain, data = many_spline_res, layout = c(7, 4),scales = "free")

dev.off()



