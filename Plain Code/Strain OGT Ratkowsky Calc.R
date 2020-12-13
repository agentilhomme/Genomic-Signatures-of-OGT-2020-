# use Ratkowksy 1983 model to predict optimal growth temperatures

library(dplyr)
library(minpack.lm)

# ----------------------------------------------------read in files---------------------------------------------------------
genus = "Colwellia"
GR = read.table("./Data/growthrates.tsv",fill=T,head=T,sep="\t",colClasses="character")
biospec = read.csv("./References/Literature Colwellia Growthrates.csv")

#-----------------------------------------------------set up files ---------------------------------------------------------
colnames(GR) = c("strain","replicate","temp","mumax","r2") # by day 

# load possible proxy format of Colwellia species growth rates 
biospec$rate.per.hour = 60*biospec$rate.per.minute
genome.list = as.character(unique(biospec$binomial.name[grepl(genus,biospec$binomial.name)]))
biogrd = biospec[biospec$binomial.name %in% genome.list,]
biogrd = data.frame("strain"=biogrd$binomial.name,
                    "replicate"="OD1",
                    "temp"=biogrd$T.C,
                    "mumax"=biogrd$rate.per.hour,
                    "r2"=0.99)
#------------------------------------create dataframes from weighted and Monte Carlo approach-------------------------------
#create dataframe for plain weighted boostrap
grd = rbind(GR,biogrd)
colnames(grd) = c("strain","rep","temp","rate","r2")
grd$temp = as.numeric(as.character(grd$temp)) + 273 #calculate Kelvin temperature
grd$rate = as.numeric(as.character(grd$rate))
grd$r2 = as.numeric(as.character(grd$r2))
grd$sqperday = sqrt(grd$rate) #Ratkowsky uses square root of rate to 1/day
grd = na.omit(grd)
grd$weights = -1*log10(1-grd$r2) #calculate weights based on R^2 values

# data frame to create monte-carlo distribution for just my strains
GR$replicate = NULL
colnames(GR)[1] = "strn"
GR$strn = as.factor(GR$strn);GR$temp = as.factor(GR$temp);GR$mumax = as.numeric(GR$mumax);GR$r2 = as.numeric(GR$r2)
GR_average = GR %>% dplyr::group_by(strn,temp) %>%  dplyr::summarize(GR_average = mean(mumax),GR_SD = sd(mumax))  
GR_average = as.data.frame(GR_average)

#---------------------------------------------------set up Ratkowsky Function------------------------------------------------
# rT = (cc * (T - T1)(1 - exp( k * (T - T2))))^2
# Ratkowsky et al. 1983 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC217594
# rT is the growth rate
# T the temperature in Kelvin
# T1 and T2 the minimum and maximum temperatures at which rate of growth is zero
# sqrt(cc) * k1 is the slope of the regression 
# k is a constant

rat83 = function(x, cc, T1, k, T2) {
  rT = cc * (x - T1) * (1 - exp(k * (x - T2)))
  return(rT)
}

#----------------------------------------------------------Calculate OGT-----------------------------------------------------
init.avg = list("b.est"=0.01, "Tmin.est"=182, "c.est"=0.06, "Tmax.est"=293) # initial parameters
jmax = 1000 # number of bootstraps
init.save.topt = matrix(NA,nrow=jmax*length(unique(grd$strain)), ncol=9) # set up empty dataframe 
for(i in 1:length(unique(grd$strain))) {
  
  for(j in 1:jmax) {
    
    #initialize
    rt = NULL
    topt = NA 
    
    strn = sort(unique(grd$strain))[i] #i
    print(as.character(strn))
    
    # -------------------------------------------------------------Find Parameters with weights------------------------------------------------------------------------
    sgrd = grd[grd$strain == strn,c("temp","sqperday","weights")]
    sgrd = rbind(sgrd,aggregate(sgrd,by = list(sgrd$temp), FUN=mean)[,c("temp","sqperday","weights")]) 
    
    #estimate parameters before fitting for non-linear weighted boostrap using directions in Ratkowsky 1983
    mintemp3 = sort(unique(sgrd$temp))[1]
    midtemp3 = sort(unique(sgrd$temp))[3]
    maxtemp3 = max(sgrd$temp)
    
    Tmin.est3 = (sgrd$sqperday[sgrd$temp == midtemp3] * sgrd$temp[sgrd$temp == mintemp3] -
                   sgrd$sqperday[sgrd$temp == mintemp3] * sgrd$temp[sgrd$temp == midtemp3])/
      (sgrd$sqperday[sgrd$temp == midtemp3] - sgrd$sqperday[sgrd$temp == mintemp3])
    
    b.est3 = (sgrd$sqperday[sgrd$temp == mintemp3] - sgrd$sqperday[sgrd$temp == midtemp3])/(mintemp3 - midtemp3)
    
    lhs1_3 = log(1 - sgrd$sqperday[sgrd$temp == mintemp3]/(b.est3 * (maxtemp3 - Tmin.est3)))
    lhs2_3 = log(1 - sgrd$sqperday[sgrd$temp == mintemp3]/(b.est3 * (midtemp3 - Tmin.est3)))
    
    Tmax.est3 = (lhs1_3*midtemp3 - lhs2_3*maxtemp3)/(lhs1_3-lhs2_3)
    
    c.est3 = (lhs1_3 - lhs2_3)/(maxtemp3 - midtemp3)
    
    Tmin.est3 = median(Tmin.est3[Tmin.est3 > 0 & Tmin.est3 < mintemp3 & !is.nan(Tmin.est3)])
    if(is.na(Tmin.est3)) Tmin.est3 = 250
    b.est3 = median(b.est3[b.est3 > 0 & !is.nan(b.est3)])
    Tmax.est3 = median(Tmax.est3[Tmax.est3>0 & !is.nan(Tmax.est3)])
    c.est3 = median(c.est3[c.est3>0 & !is.nan(c.est3)])
    
    #--------------------------------------------------------------Find Parameters with Monte Carlo-------------------------------------------------------------------------- 
    #subset the data dand create distribution around average and SD for monte carlo analysis
    sgrd2 = rbind(data.frame("temp"=-1,"sqperday"=sqrt(rnorm(1, GR_average[GR_average$temp==-1 & GR_average$strn == strn,][["GR_average"]], GR_average[GR_average$temp==-1 & GR_average$strn == strn,][["GR_SD"]]))),
                  data.frame("temp"= 4,"sqperday"=sqrt(rnorm(1, GR_average[GR_average$temp== 4 & GR_average$strn == strn,][["GR_average"]], GR_average[GR_average$temp== 4 & GR_average$strn == strn,][["GR_SD"]]))),
                  data.frame("temp"=11,"sqperday"=sqrt(rnorm(1, GR_average[GR_average$temp==11 & GR_average$strn == strn,][["GR_average"]], GR_average[GR_average$temp==11 & GR_average$strn == strn,][["GR_SD"]]))),
                  data.frame("temp"=17,"sqperday"=sqrt(rnorm(1, GR_average[GR_average$temp==17 & GR_average$strn == strn,][["GR_average"]], GR_average[GR_average$temp==17 & GR_average$strn == strn,][["GR_SD"]]))))
    sgrd2$temp = sgrd2$temp + 273
    
    #estimate parameters before fitting for monte carlo distribution boostrap using directions in Ratkowsky 1983
    mintemp2 = sort(unique(sgrd2$temp))[1]
    midtemp2 = sort(unique(sgrd2$temp))[3]
    maxtemp2 = max(sgrd2$temp)
    
    Tmin.est2 = (sgrd2$sqperday[sgrd2$temp == midtemp2] * sgrd2$temp[sgrd2$temp == mintemp2] - 
                   sgrd2$sqperday[sgrd2$temp == mintemp2] * sgrd2$temp[sgrd2$temp == midtemp2])/
      (sgrd2$sqperday[sgrd2$temp == midtemp2] - sgrd2$sqperday[sgrd2$temp == mintemp2])
    
    b.est2 = (sgrd2$sqperday[sgrd2$temp == mintemp2] - sgrd2$sqperday[sgrd2$temp == midtemp2])/(mintemp2 - midtemp2)
    
    lhs1_2 = log(1 - sgrd2$sqperday[sgrd2$temp == mintemp2]/(b.est2 * (maxtemp2 - Tmin.est2)))
    lhs2_2 = log(1 - sgrd2$sqperday[sgrd2$temp == mintemp2]/(b.est2 * (midtemp2 - Tmin.est2)))
    
    Tmax.est2 = (lhs1_2*midtemp2 - lhs2_2*maxtemp2)/(lhs1_2-lhs2_2)
    
    c.est2 = (lhs1_2 - lhs2_2)/(maxtemp2 - midtemp2)
    
    Tmin.est2 = median(Tmin.est2[Tmin.est2 > 0 & Tmin.est2 < mintemp2 & !is.nan(Tmin.est2)])
    if(is.na(Tmin.est2)) Tmin.est2 = 250
    b.est2 = median(b.est2[b.est2 > 0 & !is.nan(b.est2)])
    Tmax.est2 = median(Tmax.est2[Tmax.est2>0 & !is.nan(Tmax.est2)])
    c.est2 = median(c.est2[c.est2>0 & !is.nan(c.est2)])
    
    #---------------------------------------------------------fit the models for the weighted and monte carlo --------------------------------------------
    # do fitting using nonlinear weighted least squares for plain boostrap
    try({
      rt_weight = nlsLM(sqperday ~ b * (temp - Tmin) * (1 - exp(cc * (temp - Tmax))),
                        data = sgrd,
                        weights = sgrd$weights,
                        start = list(b = b.est3, Tmin = Tmin.est3, cc = c.est3, Tmax = Tmax.est3),
                        algorithm = "port", trace = FALSE, control = list(maxiter=1000, tol=1e-6, minFactor=1e-5, printEval = TRUE, warnOnly = TRUE))
    },silent = T)
    
    try({
      if(is.null(rt_weight)) {
        rt_weight = nlsLM(sqperday ~ b * (temp - Tmin) * (1 - exp(cc * (temp - Tmax))),
                          data = sgrd,
                          weights = sgrd$weights,
                          start = list(b = init.avg$b.est, Tmin = init.avg$Tmin.est, cc = init.avg$c.est, Tmax = init.avg$Tmax.est),
                          algorithm = "port", trace = FALSE, control = list(maxiter=1000, tol=1e-6, minFactor=1e-5, printEval = TRUE, warnOnly = TRUE))
        
      }
    }, silent=TRUE)
    
    
    # Non-weighted nonlinear model based on ratkosky model for monte carlo distribution 
    try({
      rt_MC = nlsLM(sqperday ~ b * (temp - Tmin) * (1 - exp(cc * (temp - Tmax))),
                    data = sgrd2,
                    start = list(b = b.est2, Tmin = Tmin.est2, cc = c.est2, Tmax = Tmax.est2),
                    algorithm = "port", trace = FALSE, control = list(maxiter=1000, tol=1e-6, minFactor=1e-5, printEval = TRUE, warnOnly = TRUE))
    },silent = T)
    
    try({
      # if that fails, try using previous good estimates as starting values  
      if(is.null(rt_MC)) {
        rt_MC = nlsLM(sqperday ~ b * (temp - Tmin) * (1 - exp(cc * (temp - Tmax))),
                      data = sgrd2,
                      start = list(b = init.avg$b.est, Tmin = init.avg$Tmin.est, cc = init.avg$c.est, Tmax = init.avg$Tmax.est),
                      algorithm = "port", trace = FALSE, control = list(maxiter=1000, tol=1e-6, minFactor=1e-5, printEval = TRUE, warnOnly = TRUE))
        
      }
    }, silent=TRUE)
    
    #--------------------------------------------------------------FInd and save the OGTs-----------------------------------------------------------------------
    xs = seq(0,1000,0.1)
    xs = xs - 273
    if(!is.null(rt_weight)) {
      rtco_weight = as.data.frame(t(coef(rt_weight)))
      rt.fit_weight = rat83(xs, rtco_weight$b, rtco_weight$Tmin, rtco_weight$cc, rtco_weight$Tmax)
      rt.fit_weight[rt.fit_weight < 0] = NA
      Weight_OGT = xs[which.max(rt.fit_weight)]-273
    } else {
      Weight_OGT = NA
    }
    
    if(!is.null(rt_MC)) {
      rtco.MC = as.data.frame(t(coef(rt_MC)))
      rt.fit.MC = rat83(xs, rtco.MC$b, rtco.MC$Tmin, rtco.MC$cc, rtco.MC$Tmax)
      rt.fit.MC_fit = rat83(c(272,277,284,290), rtco.MC$b, rtco.MC$Tmin, rtco.MC$cc, rtco.MC$Tmax)
      rt.fit.MC[rt.fit.MC < 0] = NA
      rt.fit.MC_fit[rt.fit.MC_fit < 0] = NA
      
      # where m is the fitted and o is the observed square root of growth( from the 4 data points of the original MC distribution)
      RMSE = function(m,o){ 
        sqrt(mean((m-o)^2))
      }
      resid1 = RMSE(rt.fit.MC_fit,sgrd2$sqperday)
      
      MC_NL_OGT = xs[which.max(rt.fit.MC)]-273
    } else {
      resid1 = NA
      MC_NL_OGT = NA
    }
    
    init.save.topt[(i-1)*jmax+j,] = c(as.character(strn),j,Weight_OGT,rtco.MC$b,rtco.MC$Tmin,rtco.MC$cc,rtco.MC$Tmax,MC_NL_OGT,resid1)
  }
}
colnames(init.save.topt) = c("strn","j","Weight_OGT","MC_b","MC_Tmin","MC_cc","MC_Tmax","MC_OGT","MC.RMSE")
init.save.topt = as.data.frame(init.save.topt)
#----------------------------------------------- format final OGT table --------------------------------------------------------------
#format data table
OGT = init.save.topt
OGT$j = NULL
cols.num_OGT = colnames(OGT)[2:8]
OGT[cols.num_OGT] = sapply(OGT[cols.num_OGT],as.character)
OGT[cols.num_OGT] = sapply(OGT[cols.num_OGT],as.numeric)

#subset out literature strains and take the weighted OGT of strain as their final OGT
genome.list = paste(genome.list,collapse = "|")
ex_topt = subset.data.frame(OGT,grepl(genome.list,strn))
ex_topt = ex_topt %>% group_by(strn) %>% dplyr::summarize(OGT = median(Weight_OGT))
ex_topt = as.data.frame(ex_topt)
ex_topt$Type = "lit_ratk"

OGT = OGT %>% filter(!grepl(genome.list,strn))
OGT = as.data.frame(OGT)

# filter monte carlo outputs by RSME and take the median of the OGT by strain
OGT$Weight_OGT = NULL
OGT_params = subset(OGT,OGT$`MC.RMSE` < exp(-30))
Final_OGT = OGT_params %>% group_by(strn) %>% dplyr::summarize(OGT = median(MC_OGT, OGT_sd = sd(MC_OGT)))

Final_OGT$Type = "ratk"
Final_OGT = rbind(Final_OGT,ex_topt)

#------------------------------------------------- write table ------------------------------------------------------------------------
write.table(init.save.topt,"./Data/Strain Ratk OGT Outputs Raw.tsv",sep = "\t",quote = F,row.names = F) # writes OGT file based on growth rates file
write.table(Final_OGT,"./Data/Strain Ratk OGT Final.tsv",sep = "\t",quote = F,row.names = F)
write.table(OGT_params,"./Data/Strain Ratk OGT Outputs Clean.tsv",sep = "\t",quote = F,row.names = F)

