library(dplyr)
library(minpack.lm)
library(reshape2)

# ----------------------------------------------------read in files---------------------------------------------------------
GR = read.table("./Data/growthrates.tsv",fill=T,head=T,sep="\t",colClasses="character")
Temp.Clust.ID = read.table("./Data/Sub-Clade ID.tsv",fill=T,head=T,sep="\t",colClasses="character")

#---------------------------------------------------set up Ratkowsky files -------------------------------------------------
colnames(Temp.Clust.ID)[1] = "strain"
grd = merge.data.frame(Temp.Clust.ID,GR,by = "strain", all = F)

# read in data
colnames(grd) = c("strain_original","strain","rep","temp","rate","r2") #'strain' is actually 'clade'
grd$strain_original = as.factor(grd$strain_original)
grd$temp = as.numeric(grd$temp) + 273
grd$rep = as.factor(grd$rep)
grd$rate = as.numeric(grd$rate)
grd = na.omit(grd)
grd$sqperday = sqrt(grd$rate) 

clade_avg = grd %>% group_by(strain,temp) %>% dplyr::summarize(GR_average = mean(rate),GR_SD = sd(rate))
clade_avg = as.data.frame(clade_avg)
clade_avg$temp = as.factor(clade_avg$temp)

#---------------------------------------------------set up Ratkowsky Function------------------------------------------------
# rT = (cc * (T - T1)(1 - exp( k * (T - T2))))^2
# Ratkowsky et al. 1983 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC217594
# rT is the growth rate
# T the temperature in Kelvin
# T1 and T2 the minimum and maximum temperatures at which rate of growth is zero
# sqrt(cc) * k1 is the slope of the regression 
# k is a constant

rat83 = function(x, b, T1, cc, T2) {
  rT = b * (x - T1) * (1 - exp(cc * (x - T2)))
  return(rT)
}

rat82a = function(x, b, T0) {
  r = b*(x - T0)
}

#-----------------------------------plot the ratkowsky fits of each cluster from a monte carlo distribution-------------------
jmax = 5000
init.save = matrix(NA,nrow=jmax*length(unique(grd$strain)), ncol=8)
for(i in 1:length(unique(grd$strain))) {
  for(j in 1:jmax) {
    
    #initialize
    rt = NULL
    topt = NA
    resid1 = NA
    
    strn = sort(unique(grd$strain))[i]
    print(as.character(strn))
    
    #---------------------------------------------------Find Parameters with Monte Carlo----------------------------------------- 
    #subset data
    sgrd = rbind(data.frame("temp"=-1,"sqperday"=sqrt(rnorm(1, clade_avg[clade_avg$temp== 272 & clade_avg$strain == strn,][["GR_average"]], clade_avg[clade_avg$temp== 272 & clade_avg$strain == strn,][["GR_SD"]]))),
                 data.frame("temp"= 4,"sqperday"=sqrt(rnorm(1, clade_avg[clade_avg$temp== 277 & clade_avg$strain == strn,][["GR_average"]], clade_avg[clade_avg$temp== 277 & clade_avg$strain == strn,][["GR_SD"]]))),
                 data.frame("temp"=11,"sqperday"=sqrt(rnorm(1, clade_avg[clade_avg$temp==284 & clade_avg$strain == strn,][["GR_average"]], clade_avg[clade_avg$temp==284 & clade_avg$strain == strn,][["GR_SD"]]))),
                 data.frame("temp"=17,"sqperday"=sqrt(rnorm(1, clade_avg[clade_avg$temp==290 & clade_avg$strain == strn,][["GR_average"]], clade_avg[clade_avg$temp==290 & clade_avg$strain == strn,][["GR_SD"]]))))
    sgrd$temp = sgrd$temp + 273
    
    
    mintemp = sort(unique(sgrd$temp))[1]
    midtemp = sort(unique(sgrd$temp))[3]
    maxtemp = max(grd$temp)
    
    #estimate parameters before fitting, using directions in Ratkowsky 1983
    Tmin.est = (sgrd$sqperday[sgrd$temp == midtemp] * sgrd$temp[sgrd$temp == mintemp] - 
                  sgrd$sqperday[sgrd$temp == mintemp] * sgrd$temp[sgrd$temp == midtemp])/
      (sgrd$sqperday[sgrd$temp == midtemp] - sgrd$sqperday[sgrd$temp == mintemp])
    
    b.est = (sgrd$sqperday[sgrd$temp == mintemp] - sgrd$sqperday[sgrd$temp == midtemp])/
      (mintemp - midtemp)
    
    lhs1 = log(1 - sgrd$sqperday[sgrd$temp == mintemp]/(b.est * (maxtemp - Tmin.est)))
    lhs2 = log(1 - sgrd$sqperday[sgrd$temp == mintemp]/(b.est * (midtemp - Tmin.est)))
    
    Tmax.est = (lhs1*midtemp - lhs2*maxtemp)/(lhs1-lhs2)
    c.est = (lhs1 - lhs2)/(maxtemp - midtemp)
    Tmin.est = median(Tmin.est[Tmin.est > 0 & Tmin.est < mintemp & !is.nan(Tmin.est)])
    if(is.na(Tmin.est)) Tmin.est = 250
    b.est = median(b.est[b.est > 0 & !is.nan(b.est)])
    Tmax.est = median(Tmax.est[Tmax.est>0 & !is.nan(Tmax.est)])
    c.est = median(c.est[c.est>0 & !is.nan(c.est)])
    
    #---------------------------------------------------------fit the models for the weighted and monte carlo ------------------------
    try({
      rt = nlsLM(sqperday ~ b * (temp - Tmin) * (1 - exp(cc * (temp - Tmax))),
                 data = sgrd,
                 start = list(b = b.est, Tmin = Tmin.est, cc = c.est, Tmax = Tmax.est),
                 algorithm = "port", trace = FALSE, control = list(maxiter=1000, tol=1e-6, minFactor=1e-5, printEval = TRUE, warnOnly = TRUE));
    }, silent=T)
    
    try({
      # if that fails, try using previous good estimates as starting values  
      if(is.null(rt)) {
        rt = nlsLM(sqperday ~ b * (temp - Tmin) * (1 - exp(cc * (temp - Tmax))),
                   data = sgrd,
                   start = list(b = init.avg$b.est, Tmin = init.avg$Tmin.est, cc = init.avg$c.est, Tmax = init.avg$Tmax.est),
                   algorithm = "port", trace = FALSE, control = list(maxiter=1000, tol=1e-6, minFactor=1e-5, printEval = TRUE, warnOnly = TRUE));
      }
    }, silent=T)
    
    #--------------------------------------------------------------FInd and save the OGTs------------------------------------------------
    xs = seq(0,1000,0.1)
    xs = xs - 273
    if(!is.null(rt)) {
      # finish plotting
      rtco = as.data.frame(t(coef(rt)))
      rt.fit = rat83(xs, rtco$b, rtco$Tmin, rtco$cc, rtco$Tmax)
      rt.fit_fit = rat83(c(272,277,284,290), rtco$b, rtco$Tmin, rtco$cc, rtco$Tmax)
      rt.fit[rt.fit < 0] = NA
      rt.fit_fit[rt.fit_fit < 0] = NA
      
      RMSE = function(m,o){ # where m is the fitted and o is the observed square root of growth( from the 4 data points of the original MC distribution)
        sqrt(mean((m-o)^2))
      }
      resid1 = RMSE(rt.fit_fit,sgrd$sqperday)
      
      topt = round(xs[which.max(rt.fit)],2)-273
    } else {
      topt = sgrd$temp[which.max(sgrd$sqperday)]-273
      resid1 = NA
    } 
    init.save[(i-1)*jmax+j,] = c(as.character(strn),j,rtco$b, rtco$Tmin, rtco$cc, rtco$Tmax, topt,resid1)
  }
}

#-------------------------------------------------------------- format final OGT table -----------------------------------------------------
# set up raw ratkowsky parameters
colnames(init.save) = c("Clade","j","b", "Tmin", "c", "Tmax", "topt","resid")
init.save = as.data.frame(init.save)
init.save$resid = as.numeric(as.character(init.save$resid))

# filter sets of parameters with a certain RSME value
clade_ratk_filt = subset(init.save,init.save$resid < exp(-30))

# melt filtered ratkowsky table
ratk_TC_melt = melt(clade_ratk_filt,id.vars = c("Clade","j","resid")) #ratk_TC is clade_ratk_filt
ratk_TC_melt$value = as.numeric(ratk_TC_melt$value)

# take the median and standard deviation
init.plot_RSME_med = ratk_TC_melt %>% group_by(Clade,variable) %>% dplyr::summarize(med = median(value))
init.plot_RSME_sd = ratk_TC_melt %>% group_by(Clade,variable) %>% dplyr::summarize(sd = sd(value))
init.plot_RSME_med_re = recast(init.plot_RSME_med,Clade~variable)
init.plot_RSME_sd_re = recast(init.plot_RSME_sd,Clade~variable)

#--------------------------------------------------------------Write Files--------------------------------------------------------------------------------------
write.table(x=clade_ratk_filt,file="./Data/Clade Ratk Params.tsv",sep="\t",quote=F)
write.table(x=init.save,file="./Data/Clade Ratk Param Raw.tsv",sep="\t",quote=F)
write.table(x=init.plot_RSME_med_re,file="./Data/Clade Ratk Med Params.tsv",sep="\t",quote=F)
write.table(x=init.plot_RSME_sd_re,file="./Data/Clade Ratk SD Params.tsv",sep="\t",quote=F)
