################################
### NDVI into 1x1 Grid Array ###
################################

library(raster)
library(zoo)
library(biwavelet)
library(forecast)
library(RNetCDF)
library(maptools)
  data(wrld_simpl)

# flsN <- list.files("~/Desktop/VHP_SM_SMN")
#   yearsN <- as.numeric(substring(flsN, 17, 20))
#   weekN  <- as.numeric(substring(flsN, 21, 23))
# 
# index <- data.frame(Date = seq(as.POSIXct(paste0(yearsN[1], "-01-01")),
#                                as.POSIXct(paste0(yearsN[length(flsN)], "-12-31")), by = "day"))
#   
# ndvi.tmp <- data.frame(Year = as.numeric(format(index$Date, "%Y")))
# index$Week = unlist(lapply(split(ndvi.tmp, f = ndvi.tmp$Year), function(x) rep(1:53, each = 7)[1:nrow(x)]))   
# 
# index$NDVI <- merge(data.frame(Year = as.numeric(format(index$Date, "%Y")), Week = index$Week),
#                     data.frame(Ind = 1:length(flsN), Year = yearsN, Week = weekN), by = c("Year", "Week"), all.x = T)$Ind
#   
# index <- aggregate(index, by = list(paste(format(index$Date, "%Y"), index$Week)), max)[,2:4]
# index <- index[order(index$Date),]
# index <- index[index$Week<53,]
#   
# index$yday <- as.numeric(format(index$Date, "%j"))
# index <- index[,c(1,2,4,3)]
#   
# index <- subset(index, format(Date, "%Y")%in%c(1982:2015))
# 
# vhA <- array(dim = c(180, 360, nrow(index)))
# 
    r.out <- raster(extent(c(-180, 180, -90, 90)), res = c(1,1))  
      proj4string(r.out) <- proj4string(wrld_simpl)
#   
# for(i in 1:nrow(index)) {
#   
#   cat("\b\b\b\b\b\b\b\b\b\b\b\b")
#   cat(i, " of ", nrow(index), "\r")
#   
#   if(!is.na(index$NDVI[i])){
#   tmp <- raster(paste0("~/Desktop/VHP_SM_SMN/", flsN[index$NDVI[i]]))
#     tmp[tmp[]<0] <- NA
#     tmp1 <- raster:::resample(tmp, r.out, method = "ngb")
# 
#   vhA[,,i] <-  as.matrix(tmp1)
#   }
# }
# save(vhA, file = "vh_array.RData")
load("vh_array.RData")

#######################  
######## ocean mask ###
#######################
  r0 <- rasterize(wrld_simpl, r.out)
  t1 <- r.out
    t1[] <- vhA[,,2]  
    t1[] <- ifelse(r0[]>0, t1[], NA)
  plot(t1)

  ind1 <- apply(vhA, 1:2, function(x) !all(is.na(x)))
  ind1[is.na(as.matrix(r0))] <- NA
  r1 <- raster(ind1)
    extent(r1) <- extent(r.out)

    
#######################  
######## snow mask  ###
#######################    
# fid <- open.nc("~/Downloads/nhsce_v01r01_19661004_20160801.nc")
# #print.nc(fid)
#     
# dat<-read.nc(fid)
#     t<-dat$time
#     z<-dat$snow_cover_extent
#     ylat<-dat$longitude
#     xlon<-dat$latitude
#     close.nc(fid)
# 
#     start <- as.POSIXct("1966-10-03", "GMT")
#     date <- start + (t*24*60*60)
#     
# ### Projection and raster extent
# lon <- raster(ylat)
# lat <- raster(xlon)
# snow <- brick(z)
#     
# ## Projection
# projSnow <- "+proj=stere +lat_0=90 +lon_0=10"
# xy <- project(cbind(values(lon), values(lat)), projSnow)
#     
# # Rasterize
# r <- raster(xmn = min(xy[,1]), xmx = max(xy[,1]), ymn = min(xy[,2]), ymx = max(xy[,2]), nrow=88, ncol=88)
# 
# snow.r <- r.out  
#   
# for(s in which(format(date, "%Y")=="2005")[1]:max(which(format(date, "%Y")=="2015"))) {
#   tmp1 <- rasterize(xy, y = r, field = values(snow[[s]]))    
#     proj4string(tmp1) <- projSnow    
#   tmp2 <- projectRaster(tmp1, r.out)
#   tmp3 <- raster:::resample(tmp2, r.out)
#   tmp3[] <- ifelse(tmp3[]>0, 1, NA)
#   
#   snow.r[] <- apply(cbind(snow.r[],tmp3[]), 1, sum, na.rm = T)
# }
# snow.r[] <- snow.r[]/11
# save("snow.r", file = "snowR.RData")
load("snowR.RData")

########################  
#### cosine function ###
########################    
leastS.cos <- function(params, sd) {
  fit  <- params[1]*cos(pi*((1: length(x1))/(length(x1)/((length(x1)/52)*2))) + (pi+params[2]))
  if(is.matrix(Mx)) {
    -sum(apply(Mx, 2, function(x) dnorm(x = x, mean = fit, sd = sd, log = T)), na.rm = T)
  } else {
    -sum(dnorm(x= Mx[!is.na(Mx)], mean=fit[!is.na(Mx)], sd=sd, log=TRUE))
  }
}

leastS.cos.5 <- function(params, sd) {
  fit  <- params[1]*cos(pi*((1: length(x1))/(length(x1)/((length(x1)/(52/2))*2))) + (pi+params[2]))
  if(is.matrix(Mx)) {
    -sum(apply(Mx, 2, function(x) dnorm(x = x, mean = fit, sd = sd, log = T)), na.rm = T)
  } else {
    -sum(dnorm(x= Mx[!is.na(Mx)], mean=fit[!is.na(Mx)], sd=sd, log=TRUE))
  }
} 

# outVH <- array(dim = c(dim(r1)[1], dim(r1)[2], 4))
load("outVH.RData")    

for(r in 1:dim(r1)[1]) {
  for(c in 1:dim(r1)[2]) {
    
    if(is.na(outVH[r,c,3])) {
    
    cat("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b")
    cat(r, " - ", c, "\r")

    
    if(!is.na(as.matrix(r1)[r,c]) & as.matrix(r1)[r,c]==1) {
      
      y <- vhA[r, c,]
      if((sum(!is.na(y) & y<0.075)/length(y))>0.1 & as.matrix(snow.r)[r,c]>8)  y[y<0.1] <- min(y, na.rm = T)
      
      outVH[r,c,4] <- diff(range(quantile(y, probs = c(0.975, 0.025), na.rm = T)))
      
      # plot(y, pch = 16, cex = 0.5, type = 'o')
      
      dat <- data.frame(t = 1:dim(vhA)[3], y = y)
      if(quantile(dat$y, prob = 0.95, na.rm = T)>0.12) {
        
        x0 <- loess(y ~ t, data = dat, span = 0.02)
        x1 <- predict(x0, newdata = data.frame(t = 1:nrow(dat)))

        ## normalize x
        Mx <- x1 - mean(x1, na.rm = T)
        
        fit0 <- optim(fn = leastS.cos, par = c(a = 50, b = 0), sd = 0.001)
        if(fit0$par[1]==0){
          Mx <- dat$y
          fit0 <- optim(fn = leastS.cos, par = c(a = 50, b = 0), sd = 0.001)
          Mx <- x1 - mean(x1, na.rm = T)
        }	
        curve.1 <- fit0$par[1]*cos(pi*((1:length(x1))/(length(x1)/((length(x1)/52)*2))) + 
                                     (pi+fit0$par[2])) +  mean(x1, na.rm=T)
        
        spl <- which(diff(curve.1[-length(curve.1)])<0 & diff(curve.1[-1])>0) +1
        
        tmp1 <- split(dat, f = cut(1:nrow(dat), breaks = c(0, spl, nrow(dat))))
        
        k <- sapply(tmp1, function(x) sum(is.na(x)))
        tmp2 <- do.call("rbind", tmp1[-which(k>15)])
          rownames(tmp2) <- 1:nrow(tmp2)
        tmp2$y <- na.approx(tmp2$y, rule = 2)
        # plot(tmp2[,2], type = "o", pch  =16)
        
        
        ### Wavelet analysis
        wt  <- wt(cbind(1:nrow(tmp2),tmp2$y))
        power  <- log2(wt$power.corr)
        
        time   <- wt$t 
        period <- wt$period/52
        
        tmp01.pow <- apply(power, 1, median, na.rm = T) 
        tmp01.sig <- apply(wt$signif, 1, median, na.rm = T)
        tmp01.pow[tmp01.sig<1] <- NA
        
        #plot(period, tmp01.pow, type = "o")
        
        ind <- c(1, which((!is.na(tmp01.pow[-length(tmp01.pow)]) &  is.na(tmp01.pow[-1])) |
                          (is.na(tmp01.pow[-length(tmp01.pow)]) & !is.na(tmp01.pow[-1]))), length(tmp01.pow))
        
        tmp02 <- split(data.frame(period, tmp01.pow, tmp01.sig), f = cut(1:length(tmp01.pow), breaks = unique(ind)))
        
        pck.periods <- function(p) {
          if(sum(!is.na(p[,2]))>1) {
            ind01 <- c(ifelse((p[1, 2]>p[2, 2]), 1, NA), which((p[,2]>c(NA, p[-nrow(p),2]) & c(p[-1,2], NA)<p[,2])),
                       ifelse(p[nrow(p),2]>p[(nrow(p)-1),2], nrow(p), NA))
            p[ind01[!is.na(ind01)],]
          } else p[which.max(p[,2]),]
        }
        
        tmp04 <- do.call(rbind, lapply(tmp02, pck.periods))[,1:3]
        tmp04 <- subset(tmp04, tmp01.pow>=0)
        if(nrow(tmp04)>0) {
        tmp04 <- tmp04[order(tmp04$tmp01.pow, decreasing = T),]
        
        outVH[r,c,1] <- round(tmp04[1,1], 1)
        outVH[r,c,2] <- tmp04[1,2]
        }
    
        if(!is.na(outVH[r,c,1]) & outVH[r,c,1]<5){
        
        ### (Overall) Predictability
        x0 <- loess(y ~ t, data = tmp2, span = 0.02)
        x1 <- predict(x0, newdata = data.frame(t = 1:nrow(tmp2)))

        ## normalize x
        Mx <- x1 - mean(x1, na.rm = T)
        
        if(outVH[r,c,1]>0.25 & outVH[r,c,1]<0.75) {
          fit0 <- optim(fn = leastS.cos.5, par = c(a = 50, b = 0), sd = 0.001)
          if(fit0$par[1]==0){
            Mx <- dat$y
            fit0 <- optim(fn = leastS.cos, par = c(a = 50, b = 0), sd = 0.001)
            Mx <- x1 - mean(x1, na.rm = T)
          }	
          curve.1 <- fit0$par[1]*cos(pi*((1:length(x1))/(length(x1)/((length(x1)/(52/2))*2))) + 
                                       (pi+fit0$par[2])) +  mean(x1, na.rm=T) 
          spl <- which(diff(curve.1[-length(curve.1)])<0 & diff(curve.1[-1])>0) +1
          tmp3 <- split(tmp2, f = cut(1:nrow(tmp2), breaks = c(0, spl[c(TRUE, FALSE)], nrow(tmp2))))
          
        } else {
        fit0 <- optim(fn = leastS.cos, par = c(a = 50, b = 0), sd = 0.001)
        if(fit0$par[1]==0){
          Mx <- dat$y
          fit0 <- optim(fn = leastS.cos, par = c(a = 50, b = 0), sd = 0.001)
          Mx <- x1 - mean(x1, na.rm = T)
        }	
        curve.1 <- fit0$par[1]*cos(pi*((1:length(x1))/(length(x1)/((length(x1)/52)*2))) + 
                                     (pi+fit0$par[2])) +  mean(x1, na.rm=T)
        
        spl <- which(diff(curve.1[-length(curve.1)])<0 & diff(curve.1[-1])>0) +1
        tmp3 <- split(tmp2, f = cut(1:nrow(tmp2), breaks = c(0, spl, nrow(tmp2))))
        }
        
        
        ind2 <- sapply(tmp3, function(x) nrow(x))
        if(ind2[1]<45) tmp3 <- tmp3[-1]
        ind3 <- sapply(tmp3, function(x) nrow(x))
        if(ind3[length(ind3)]<45) tmp3 <- tmp3[-length(ind3)]; ind3 <- ind3[-length(ind3)]
        
        tmp4 <- do.call('rbind', tmp3)
          rownames(tmp4) <- 1:nrow(tmp4)
        
        ind.years <- suppressWarnings(cbind(1:(length(tmp3)-4), 4:(length(tmp3)-1), as.numeric(ind3)[-c(1:4)]))
        
        for(i in 1:nrow(ind.years)) {
          ts <- ts(tmp4$y[((ind.years[i,1]-1)*52):(((ind.years[i,2])*52))], frequency = 52)
          # plot(ts)
          ets <- stlf(ts, method = "arima", h = ind.years[i,3])
          
          if(i==1) {
            outVHR <- cbind(ets$mean, ets$lower, ets$upper)
          } else {
            outVHR <- rbind(outVHR, cbind(ets$mean, ets$lower, ets$upper))
          }
        }
        
        outVHR <- data.frame(y = tmp4$y[((4*52)+1):nrow(tmp4)][1:nrow(outVHR)], as.matrix(outVHR))
        
        
        #plot(NA, xlim = c(1, nrow(outVHR)), ylim = c(range(outVHR, na.rm = T)))
        #polygon(c(1:nrow(outVHR),rev(1:nrow(outVHR))), 
        #        c(outVHR$ets.lower.95., rev(outVHR$ets.upper.95.)),
        #        col = "skyblue", border = NA)
        #polygon(c(1:nrow(outVHR),rev(1:nrow(outVHR))), 
        #        c(outVHR$ets.lower.80., rev(outVHR$ets.upper.80.)),
        #        col = "skyblue3", border = NA)  
        #lines(1:nrow(outVHR), outVHR$ets.mean, col = "darkblue", lwd = 3)
        #lines(1:nrow(outVHR), outVHR$y, col = "violetred1")
        
        RSS <- outVHR$y - mean(outVHR$y, na.rm = T)
        SSE <- sum((outVHR$y-outVHR$ets.mean)^2, na.rm = T)
        outVH[r,c,3] <- (sum(RSS^2, na.rm = T)-SSE)/sum(RSS^2, na.rm = T) 
        }
      }
    }
    
    
    }
  }
if(floor(r/10)==r/10) {
dev.off()
opar <- par(mfrow = c(3,1), mar = c(2,2,1,1))
plot(raster(outVH[,,2], xmn = -180, xmx = 180, ymn = -90, ymx = 90))
  plot(wrld_simpl, add = T)
plot(raster(outVH[,,3], xmn = -180, xmx = 180, ymn = -90, ymx = 90))
  plot(wrld_simpl, add = T)
plot(raster(outVH[,,4], xmn = -180, xmx = 180, ymn = -90, ymx = 90))
  plot(wrld_simpl, add = T)  
par(opar)

# save(outVH, file = "outVH.RData")
}
}
