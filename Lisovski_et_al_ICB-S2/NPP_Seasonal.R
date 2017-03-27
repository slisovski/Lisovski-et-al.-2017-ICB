#############################################################################
### Project: Defining the Degree of Seasonality (Lisovski et al 2017 ICB)  ##
### Author: Simeon Lisovski                                                ##
### Date: March 2017                                                       ##
#############################################################################
### Seasonal Components: Montly Net Primary Productivity (1x1 degree grid) ##
#############################################################################

library(doParallel)
library(biwavelet)
library(raster)
library(RNetCDF)
library(maptools)
  data(wrld_simpl)

load("Revision_Analysis/Output/NPP/npp_monthly_array.RData")
# plot(raster(nppA_MNT[,,3]))

month.ind <- seq(as.POSIXct("2006-01-15"), as.POSIXct("2015-12-15"), by = "month")

########################  
#### cosine function ###
########################  

leastS.cos <- function(params, sd) {
  fit  <- params[1]*cos(pi*((1: length(y))/(length(y)/((length(y)/12)*2))) + (pi+params[2]))
  if(is.matrix(Mx)) {
    -sum(apply(Mx, 2, function(x) dnorm(x = x, mean = fit, sd = sd, log = T)), na.rm = T)
  } else {
    -sum(dnorm(x= Mx[!is.na(Mx)], mean=fit[!is.na(Mx)], sd=sd, log=TRUE))
  }
}

leastS.cos.5 <- function(params, sd) {
  fit  <- params[1]*cos(pi*((1: length(y))/(length(y)/((length(y)/(12/2))*2))) + (pi+params[2]))
  if(is.matrix(Mx)) {
    -sum(apply(Mx, 2, function(x) dnorm(x = x, mean = fit, sd = sd, log = T)), na.rm = T)
  } else {
    -sum(dnorm(x= Mx[!is.na(Mx)], mean=fit[!is.na(Mx)], sd=sd, log=TRUE))
  }
} 

## parallel init ########
cl0 <- detectCores()    #
cl <- makeCluster(cl0-1)#
registerDoParallel(cl)  #
# stopImplicitCluster() #
#########################

npp.s  <- list(Amp = array(dim = c(dim(nppA_MNT)[1], dim(nppA_MNT)[2], 10)),
               Pre = array(dim = c(dim(nppA_MNT)[1], dim(nppA_MNT)[2], 10)),
               wt  = array(dim = c(dim(nppA_MNT)[1], dim(nppA_MNT)[2], 2)))
# load("Revision_Analysis/Output/NPP/npp_monthly_seasonality.RData")

for(r in 1:dim(npp.s$Amp)[1]) {
  for(c in 1:dim(npp.s$Amp)[2]) {
    
    cat(paste(rep("\b\b\b\b\b\b", 6), collapse = ""))
    cat(r, " - ", c, "\r")
    
    y <- nppA_MNT[r, c, ]
    # plot(y, type = "o")
    
    if(all(is.na(npp.s$Amp[r,c,])) & (!all(is.na(y)) & abs(diff(quantile(y, probs = c(0.9, 0.1), na.rm = T)))>0.15)) { 
      
      ### Wavelet analysis to determine significant periodicity
      
      wt  <- wt(cbind(1:length(y),y))
      power  <- log2(wt$power.corr)
      
      time   <- wt$t 
      period <- wt$period/12
      
      tmp.pow <- apply(power[,-c(1:12, (ncol(power)-11):ncol(power))], 1, median, na.rm = T) 
      tmp.sig <- apply(wt$signif[,-c(1:12, (ncol(power)-11):ncol(power))], 1, median, na.rm = T)
      
      # plot(period, tmp.pow, type = "o")
      # points(period, tmp.pow, pch = 16, col = ifelse(tmp.sig>1, "red", "grey90"))
      
      ind <- c(1, which((!is.na(tmp.pow[-length(tmp.pow)]) &  is.na(tmp.pow[-1])) |
                          ( is.na(tmp.pow[-length(tmp.pow)]) & !is.na(tmp.pow[-1]))), length(tmp.pow))
      
      tmp02 <- split(data.frame(period, tmp.pow, tmp.sig), f = cut(1:length(tmp.pow), breaks = unique(ind)))
      
      pck.periods <- function(p) {
        if(sum(!is.na(p[,2]))>1) {
          ind01 <- c(ifelse((p[1, 2]>p[2, 2]), 1, NA), which((p[,2]>c(NA, p[-nrow(p),2]) & c(p[-1,2], NA)<p[,2])),
                     ifelse(p[nrow(p),2]>p[(nrow(p)-1),2], nrow(p), NA))
          p[ind01[!is.na(ind01)],]
        } else p[which.max(p[,2]),]
      }
      
      tmp03 <- do.call(rbind, lapply(tmp02, pck.periods))[,1:3]
        wt.out <- subset(tmp03[order(tmp03$tmp.pow, decreasing = T),], tmp.sig>1)
      
      if(!is.null(wt.out) & (any(wt.out[,1]>0.85 & wt.out[,1]<1.15) | any(wt.out[,1]>0.45 & wt.out[,1]<0.65))) {
    
        dat <- data.frame(t = 1:dim(nppA_MNT)[3], 
                        years = as.numeric(format(month.ind, "%Y")), 
                        y = y,
                        y.norm = scale(y))
      

        ## normalize x
        Mx <- y - mean(y, na.rm = T)
        
        fit0 <- optim(fn = leastS.cos, par = c(a = 50, b = 0), sd = 0.001)
        if(fit0$par[1]==0){
          Mx <- suppressWarnings(predict(loess(1:nrow(dat) ~ dat$y, span = 0.02)))
          fit0 <- optim(fn = leastS.cos, par = c(a = 50, b = 0), sd = 0.001)
          Mx <- y - mean(y, na.rm = T)
        }	
        curve.1 <- fit0$par[1]*cos(pi*((1:length(y))/(length(y)/((length(y)/12)*2))) + 
                                     (pi+fit0$par[2])) +  mean(y, na.rm=T)
        
        # plot(month.ind, y, type = "o")
        #   lines(month.ind, curve.1, lwd = 2, col = "orange")
        
        spl  <- which(diff(curve.1[-length(curve.1)])<0 & diff(curve.1[-1])>0) +1
        tmp3 <- split(data.frame(dat), f = cut(1:nrow(dat), breaks = c(0, spl, nrow(dat))))
        
        ind2 <- sapply(tmp3, function(x) nrow(x))
        if(ind2[1]<7) tmp3 <- tmp3[-1]
        ind3 <- sapply(tmp3, function(x) nrow(x))
        if(ind3[length(ind3)]<7) {tmp3 <- tmp3[-length(ind3)]; ind3 <- as.vector(ind3[-length(ind3)])}
        
        tmp4 <- do.call('rbind', tmp3)
          rownames(tmp4) <- 1:nrow(tmp4)
        yrs <- as.vector(apply(cbind(1:length(ind3), ind3), 1, function(x) rep(x[1], x[2])))
        if(is.list(yrs)) tmp4$yrs <- unlist(yrs) else tmp4$yrs <- yrs
        
        ### Predictability
        ##  normalize x
        ind.years <- suppressWarnings(cbind(1:(length(tmp3)-4), 4:(length(tmp3)-1), as.numeric(ind3)[-c(1:4)]))
        
        fc <- foreach(loop=1:nrow(ind.years), .combine = rbind, .packages = "forecast") %dopar% {
          ts <- ts(tmp4$y.norm[tmp4$yrs%in%c(ind.years[loop,1]:ind.years[loop,2])], frequency = ind.years[loop,3])
          ets <- stlf(ts, method = "arima", h = ind.years[loop,3])
          # plot(ets)
          cbind(c(ets$mean), as.matrix(ets$lower), as.matrix(ets$upper))}
        
        outTempR <- data.frame(years = tmp4$yrs[min(which((tmp4$yrs==5))):nrow(tmp4)], 
                               y     = tmp4$y[min(which((tmp4$yrs==5))):nrow(tmp4)], 
                               y.norm= tmp4$y.norm[min(which((tmp4$yrs==5))):nrow(tmp4)], as.matrix(fc))
        
        outTempR$t <- 1:nrow(outTempR)
        
        fcy <- foreach(loop=c(min(outTempR$years):max(outTempR$years)), .combine = rbind) %dopar% {
          tmp <- outTempR[outTempR$years==loop,]
          RSS <- tmp$y.norm - mean(tmp$y.norm, na.rm = T)
          SSE <- sum((tmp$y.norm-tmp[,4])^2, na.rm = T)
          
          cbind(yrs  = loop,
                amp  = as.numeric(abs(diff(quantile(tmp$y, probs = c(0.975, 0.025), na.rm = T)))),
                pred = (sum(RSS^2, na.rm = T)-SSE)/sum(RSS^2, na.rm = T))
        }
        
        npp.s$Amp[r,c,fcy[,1]] <- fcy[,2]
        npp.s$Pre[r,c,fcy[,1]] <- ifelse(fcy[,3]<0, 0, fcy[,3])
        npp.s$wt[r,c,1] <- wt.out[1,1]
        npp.s$wt[r,c,2] <- wt.out[1,2]
        
      }
    }
  }
  
  if(floor(r/5)==r/5) {
    png("NPP.png", width = 400, height = 700)
    opar <- par(mfrow = c(3,1), mar = c(2,2,1,1))
    plot(raster(apply(npp.s$Amp, 1:2, mean, na.rm = T), xmn = -180, xmx = 180, ymn = -90, ymx = 90))
    plot(wrld_simpl, add = T)
    plot(raster(apply(npp.s$Pre, 1:2, mean, na.rm = T), xmn = -180, xmx = 180, ymn = -90, ymx = 90))
    plot(wrld_simpl, add = T)
    plot(raster(npp.s$wt[,,2], 1:2, xmn = -180, xmx = 180, ymn = -90, ymx = 90))
    plot(wrld_simpl, add = T)
    par(opar)
    dev.off()
    
    save(npp.s, file = "Revision_Analysis/Output/NPP/npp_monthly_seasonality.RData")
  }
} # end r


amp <- raster(apply(npp.s$Amp[,,6:10], 1:2, mean, na.rm = T), xmn = -180, xmx = 180, ymn = -90, ymx = 90)
pre <- raster(apply(npp.s$Pre[,,6:10], 1:2, mean, na.rm = T), xmn = -180, xmx = 180, ymn = -90, ymx = 90)
pow <- raster(npp.s$wt[,,2], 1:2, xmn = -180, xmx = 180, ymn = -90, ymx = 90)
per <- raster(npp.s$wt[,,1], 1:2, xmn = -180, xmx = 180, ymn = -90, ymx = 90)

plot(per, breaks = c(0, 0.75, 1.5), col = c("goldenrod4", "darkgreen"))
map(interior = F, add = T)


plot(pow, breaks = c(-3, seq(0, 9, length = 50)), col = rev(terrain.colors(45)))
  map(interior = F, add = T)

  pre[] <- ifelse(amp[]<=1, NA, pre[])
  
plot(pre)
  map(interior = F, add = T)
amp[] <- ifelse(amp[]<=1, NA, amp[])
plot(amp)
map(interior = F, add = T)



pre[] <- ifelse(pow[]<=0, NA, pre[])









