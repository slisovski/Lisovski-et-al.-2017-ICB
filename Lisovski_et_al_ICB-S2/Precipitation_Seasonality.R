## GPCP: Precipitation (seasonality)

library(raster)
library(zoo)
library(forecast)
library(RNetCDF)
library(doParallel)
library(maptools)
  data(wrld_simpl)


load("Results/Precipitation/precWeek_array.RData")
plot(raster(precW[,,1], xmn = -180, xmx = 180, ymn = -90, ymx = 90))
plot(wrld_simpl, add = T)
# hist(precW, xlim = c(0, 100), breaks = seq(min(precW, na.rm = T), max(precW, na.rm = T)+10, by = 10))


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

prec.s  <- list(Amp = array(dim = c(dim(precW)[1], dim(precW)[2], 35)),
                Pre = array(dim = c(dim(precW)[1], dim(precW)[2], 35)))

# load("Results/Precipitation/prec_s.RData")


for(r in 125:dim(precW)[1]) {
  for(c in 1:dim(precW)[2]) {
    
    cat(paste(rep("\b\b\b\b\b\b", 9), collapse = ""))
    cat(r, " - ", c, "\r")
    
    if(all(is.na(prec.s[[1]][r,c,]))) {  
      y <- precW[r,c,]
      # plot(y, pch = 16, cex = 0.5, type = 'o')
      y <- na.approx(y, rule = 2)
      
      dat <- data.frame(t = 1:dim(precW)[3], 
                        years = as.numeric(substring(dimnames(precW)[3][[1]],1,4)),
                        y = y)
      
      # plot(NA, xlim = c(0,  988), ylim = c(min(dat$y), 26),
      #      xlab = "", ylab = "Precipitation", xaxt = "n", las = 1, bty = "n")
      # with(dat[dat$years%in%c(1997:2016),], lines(1:988, y, type = "h", 
      #                                             col = "darkblue", lwd = 2))
      # axis(1, at = which(!duplicated(dat[dat$years%in%c(1997:2016),"years"])),
      #      labels = 1997:2015)
      
      
      
      ### Power ~ Years
      x  <- dat$y 
      Mx <- x - mean(x, na.rm = T)
      
      fit0 <- optim(fn = leastS.cos, par = c(a = 50, b = 0), sd = 0.001)
      
      curve.1 <- fit0$par[1]*cos(pi*((1:length(x))/(length(x)/((length(x)/52)*2))) + 
                                   (pi+fit0$par[2])) +  mean(x, na.rm=T)
      
      spl <- which(diff(curve.1[-length(curve.1)])<0 & diff(curve.1[-1])>0) +1
      tmp3 <- split(data.frame(dat), f = cut(1:nrow(dat), breaks = c(0, spl, nrow(dat))))
      
      ind2 <- sapply(tmp3, function(x) nrow(x))
      if(ind2[1]<48) tmp3 <- tmp3[-1]
      ind3 <- sapply(tmp3, function(x) nrow(x))
      if(ind3[length(ind3)]<48) {tmp3 <- tmp3[-length(ind3)]; ind3 <- as.vector(ind3[-length(ind3)])}
      
      tmp4 <- do.call('rbind', tmp3)
        rownames(tmp4) <- 1:nrow(tmp4)
      yrs <- as.vector(apply(cbind(1:length(ind3), ind3), 1, function(x) rep(x[1], x[2])))
      if(is.list(yrs)) tmp4$yrs <- unlist(yrs) else tmp4$yrs <- yrs
      
    ### (Overall) Predictability
    ind.years <- suppressWarnings(cbind(1:(length(tmp3)-4), 4:(length(tmp3)-1), as.numeric(ind3)[-c(1:4)]))
       
    fc <- foreach(loop=1:nrow(ind.years), .combine = rbind, .packages = "forecast") %dopar% {
       ts <- ts(tmp4$y[((ind.years[loop,1]-1)*52):(((ind.years[loop,2])*52))], frequency = 52)
       ets <- stlf(ts, method = "arima", h = ind.years[loop,3])
       cbind(c(ets$mean), as.matrix(ets$lower), as.matrix(ets$upper))
     }
       
    outPrecR <- data.frame(years = tmp4$years[((4*52)+1):nrow(tmp4)][1:nrow(fc)],
                           y     = tmp4$y[((4*52)+1):nrow(tmp4)][1:nrow(fc)], as.matrix(fc))
     
    # plot(NA, xlim = c(1, nrow(outPrecR)), ylim = c(range(outPrecR[,-c(1:2)], na.rm = T)))
    # lines(1:nrow(outPrecR), outPrecR[,3], col = "darkblue", lwd = 3)
    # lines(1:nrow(outPrecR), outPrecR$y, col = "violetred1")
     
    tmp3 <- split(data.frame(outPrecR), f = tmp4$yrs[-which(tmp4$yrs<5)])
       
    fcy <- foreach(loop=1:length(tmp3), .combine = rbind) %dopar% {
       tmp <- tmp3[[loop]]
       RSS <- tmp$y - mean(tmp$y, na.rm = T)
       SSE <- sum((tmp$y-tmp[,3])^2, na.rm = T)
         
       cbind(yrs  = as.numeric(names(table(tmp$years)))[order(table(tmp$years), decreasing = T)][1],
             amp  = as.numeric(abs(diff(quantile(tmp$y, probs = c(0.975, 0.025), na.rm = T)))),
             pred = (sum(RSS^2, na.rm = T)-SSE)/sum(RSS^2, na.rm = T))
      }
      
      prec.s$Amp[r,c,match(fcy[,1], c(1982:2016))] <- fcy[,2]
      prec.s$Pre[r,c,match(fcy[,1], c(1982:2016))] <- fcy[,3]
       
  }
  if(floor(r/5)==r/5) {
    png("Prec.png", width = 400, height = 700)
    opar <- par(mfrow = c(2,1), mar = c(2,2,1,1))
    plot(raster(apply(prec.s$Amp, 1:2, mean, na.rm = T), xmn = -180, xmx = 180, ymn = -90, ymx = 90),
         legend = F)
    plot(wrld_simpl, add = T)
    plot(raster(apply(prec.s$Pre, 1:2, mean, na.rm = T), xmn = -180, xmx = 180, ymn = -90, ymx = 90),
         breaks = seq(0, 0.8, by = .1), col = rev(heat.colors(9)), legend = F)
    plot(wrld_simpl, add = T)
    par(opar)
    dev.off()
    
    save(prec.s, file = "Results/Precipitation/prec_s.RData")
  }
  }
}

