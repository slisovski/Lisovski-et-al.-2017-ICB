## RNCEP II: Temperature (seasonality)

library(raster)
library(zoo)
library(biwavelet)
library(forecast)
library(RNetCDF)
library(doParallel)
library(maptools)
  data(wrld_simpl)

  
load("Results/Temperature/tempWeek_array.RData")
range(tempWeek-273.15)
  # plot(raster(tempWeek[,,4]))



########################  
#### cosine function ###
########################    
leastS.cos <- function(params, sd) {
  fit  <- params[1]*cos(pi*((1: length(x))/(length(x)/((length(x)/52)*2))) + (pi+params[2]))
  if(is.matrix(Mx)) {
    -sum(apply(Mx, 2, function(x) dnorm(x = x, mean = fit, sd = sd, log = T)), na.rm = T)
  } else {
    -sum(dnorm(x = Mx[!is.na(Mx)], mean=fit[!is.na(Mx)], sd=sd, log=TRUE))
  }
}    


## parallel init ########
cl0 <- detectCores()    #
cl <- makeCluster(cl0-1)#
registerDoParallel(cl)  #
# stopImplicitCluster() #
#########################

temp.s  <- list(Amp = array(dim = c(dim(tempWeek)[1], dim(tempWeek)[2], 35)),
                Pre = array(dim = c(dim(tempWeek)[1], dim(tempWeek)[2], 35)))

# load("Results/Temperatue/temp_s.RData")


for(r in 1:dim(tempWeek)[1]) {
  for(c in 1:dim(tempWeek)[2]) {
    
    cat(paste(rep("\b\b\b\b\b\b", 8), collapse = ""))
    cat(r, " - ", c, "\r")
    
    if(all(is.na(temp.s[[1]][r,c,]))) {  
      y <- tempWeek[r, c, ] - 273.15
      # plot(y, pch = 16, cex = 0.5, type = 'o')
      
      
      dat <- data.frame(t = 1:dim(tempWeek)[3], 
                        years = as.numeric(substring(dimnames(tempWeek)[3][[1]],1,4)), 
                        y = y)

      # plot(NA, xlim = c(0,  1040), ylim = c(min(dat$y), 26),
      #      xlab = "", ylab = "Temperature", xaxt = "n", las = 1, bty = "n")
      # with(dat[dat$years%in%c(1993:2012),], lines(1:1040, y, type = "l", 
      #                                             col = "orange", lwd = 2))
      # axis(1, at = which(!duplicated(dat[dat$years%in%c(1993:2012),"years"])),
      #      labels = 1993:2012)
      # 
      ### Years
      x  <- dat$y 
      Mx <- x - mean(x, na.rm = T)
      fit0 <- optim(fn = leastS.cos, par = c(a = 50, b = 0), sd = 0.001)
      
      curve.1 <- fit0$par[1]*cos(pi*((1:length(x))/(length(x)/((length(x)/52)*2))) + 
                                   (pi+fit0$par[2])) +  mean(x, na.rm=T)
      
      spl  <- which(diff(curve.1[-length(curve.1)])<0 & diff(curve.1[-1])>0) +1
      tmp3 <- split(data.frame(dat, x), f = cut(1:nrow(dat), breaks = c(0, spl, nrow(dat))))
      
      ind2 <- sapply(tmp3, function(x) nrow(x))
      if(ind2[1]<48) tmp3 <- tmp3[-1]
      ind3 <- sapply(tmp3, function(x) nrow(x))
      if(ind3[length(ind3)]<48) {tmp3 <- tmp3[-length(ind3)]; ind3 <- as.vector(ind3[-length(ind3)])}
      
      tmp4 <- do.call('rbind', tmp3)
        rownames(tmp4) <- 1:nrow(tmp4)
      yrs <- as.vector(apply(cbind(1:length(ind3), ind3), 1, function(x) rep(x[1], x[2])))
      if(is.list(yrs)) tmp4$yrs <- unlist(yrs) else tmp4$yrs <- yrs
      
      

     ### (Overall) Predictability
     ## normalize x
     ind.years <- suppressWarnings(cbind(1:(length(tmp3)-4), 4:(length(tmp3)-1), as.numeric(ind3)[-c(1:4)]))
  
      fc <- foreach(loop=1:nrow(ind.years), .combine = rbind, .packages = "forecast") %dopar% {
          ts <- ts(tmp4$x[((ind.years[loop,1]-1)*52):(((ind.years[loop,2])*52))], frequency = 52)
          ets <- stlf(ts, method = "arima", h = ind.years[loop,3])
          # plot(ets)
          cbind(c(ets$mean), as.matrix(ets$lower), as.matrix(ets$upper))}
  
      outTempR <- data.frame(years = tmp4$years[((4*52)+1):nrow(tmp4)][1:nrow(fc)], 
                             y     = tmp4$y[((4*52)+1):nrow(tmp4)][1:nrow(fc)], as.matrix(fc))
  
      outTempR$t <- 1:nrow(outTempR)
      outTempR0 <- subset(outTempR, years%in%c(1992:2015))
      outTempR1 <- subset(outTempR, years%in%c(1995:2015))
      
      # 
      # png("pred03.png", res = 100, width = 1400, height = 600)
      # par(bg = "transparent", mar = c(0,0,0,0), bty = "n")
      # plot(NA, xlim = c(1, max(outTempR1$t)), 
      #      ylim = c(range(outTempR[,-c(1:2, ncol(outTempR))], na.rm = T)),
      #      xaxt = "n", yaxt = "n", xlab = "", ylab = "")
      # # lines(1:nrow(outTempR), outTempR[,3], col = "darkblue", lwd = 3)
      # polygon(c(outTempR1$t, rev(outTempR1$t)),
      #         c(outTempR1[,5], rev(outTempR1[,7])), col = rgb(.4,.6,.92, alpha = 0.5), border  = NA)
      # polygon(c(outTempR1$t, rev(outTempR1$t)),
      #         c(outTempR1[,4], rev(outTempR1[,6])), col = rgb(.4,.6,.92, alpha = 0.8), border  = NA)      
      # lines(outTempR0$t, outTempR0$y, col = "black", lwd =2)
      # dev.off()
      
      
      outTempR$ind <- tmp4$yrs[tmp4$yrs>4]
  
      fcy <- foreach(loop=5:max(outTempR$ind), .combine = rbind) %dopar% {
        tmp <- outTempR[outTempR$ind==loop,]
        RSS <- tmp$y - mean(tmp$y, na.rm = T)
        SSE <- sum((tmp$y-tmp[,3])^2, na.rm = T)
  
        cbind(yrs  = as.numeric(names(table(tmp$years)))[order(table(tmp$years), decreasing = T)][1],
              amp  = as.numeric(abs(diff(quantile(tmp$y, probs = c(0.975, 0.025), na.rm = T)))),
              pred = (sum(RSS^2, na.rm = T)-SSE)/sum(RSS^2, na.rm = T))
      }
  
      temp.s$Amp[r,c,match(fcy[,1], c(1982:2016))] <- fcy[,2]
      temp.s$Pre[r,c,match(fcy[,1], c(1982:2016))] <- fcy[,3]
    
    }
    
    if(floor(r/5)==r/5) {
      png("Temp.png", width = 400, height = 700)
      opar <- par(mfrow = c(2,1), mar = c(2,2,1,1))
      plot(raster(apply(temp.s$Amp, 1:2, mean, na.rm = T), xmn = -180, xmx = 180, ymn = -90, ymx = 90))
      plot(wrld_simpl, add = T)
      plot(raster(apply(temp.s$Pre, 1:2, mean, na.rm = T), xmn = -180, xmx = 180, ymn = -90, ymx = 90))
      plot(wrld_simpl, add = T)
      par(opar)
      dev.off()
      
      save(temp.s, file = "Results/Temperature/temp_s.RData")
    }
}
}

