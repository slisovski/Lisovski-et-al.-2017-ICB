## VH: Vegetation health (seasonality)

library(raster)
library(zoo)
library(biwavelet)
library(forecast)
library(RNetCDF)
library(doParallel)
library(maptools)
  data(wrld_simpl)


load("Results/MODIS/vhWeek_array.RData")
# plot(raster(vhWeek[,,3]))

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
load("temp_data/snowR.RData")


#######################  
######## ocean mask ###
#######################

r.out <- raster(extent(c(-180, 180, -90, 90)), res = c(1,1))  
  proj4string(r.out) <- proj4string(wrld_simpl)
r0 <- rasterize(wrld_simpl, r.out)
t1 <- r.out
t1[] <- vhWeek[,,2]  
t1[] <- ifelse(r0[]>0, t1[], NA)
plot(t1)

ind1 <- apply(vhWeek, 1:2, function(x) !all(is.na(x)))
# ind1[is.na(as.matrix(r0))] <- NA
wrld <- raster(ind1)
  extent(wrld) <- extent(r.out)
plot(wrld)


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

## parallel init ########
cl0 <- detectCores()    #
cl <- makeCluster(cl0-1)#
registerDoParallel(cl)  #
# stopImplicitCluster() #
#########################

vh.s  <- list(Amp = array(dim = c(dim(vhWeek)[1], dim(vhWeek)[2], 35)),
              Pre = array(dim = c(dim(vhWeek)[1], dim(vhWeek)[2], 35)))

# load("Results/MODIS/vh_s.RData")


for(r in 1:dim(vh.s$Amp)[1]) {
  for(c in 1:dim(vh.s$Amp)[2]) {
    
    cat(paste(rep("\b\b\b\b\b\b", 6), collapse = ""))
    cat(r, " - ", c, "\r")
    
    
    if(!is.na(as.matrix(wrld)[r,c]) &  all(is.na(vh.s$Amp[r,c,]))) { 
      
      y <- vhWeek[r, c,]
      # plot(y, type = "o", cex = 0.5, pch =16)
      if(!all(is.na(y))) {
        
      if((sum(!is.na(y) & y<0.075)/length(y))>0.1 & as.matrix(snow.r)[r,c]>8)  y[y<0.1] <- min(y, na.rm = T)
        
      dat <- data.frame(t = 1:dim(vhWeek)[3], 
                        years = rep(1982:2015, each = 52), 
                        y = y)
      if(quantile(dat$y, prob = 0.95, na.rm = T)>0.15) {
        
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
        
        if(length(spl)>0) {
        tmp1 <- split(dat, f = cut(1:nrow(dat), breaks = c(0, spl, nrow(dat))))
          i1 <- lapply(tmp1, function(x) nrow(x))
          tmp1 <- lapply(1:length(i1), function(x) cbind(rep(x, i1[[x]]), tmp1[[x]]))
        
        k <- sapply(tmp1, function(x) sum(is.na(x[,3])))
          i1 <- which(k>=17)
        if(length(i1)>=1) tmp2 <- do.call("rbind", tmp1[-i1]) else tmp2 <- do.call("rbind", tmp1)
        rownames(tmp2) <- 1:nrow(tmp2)
        tmp2$y <- na.approx(tmp2$y, rule = 2)
        # plot(tmp2[,3], type = "o", pch  =16)  
      
        ### Wavelet analysis
        wt  <- wt(cbind(1:nrow(tmp2),tmp2$y))
        power  <- log2(wt$power.corr)
        
        time   <- wt$t 
        period <- wt$period/52
        
        tmp01.pow <- apply(power, 1, function(x) quantile(x, probs = 0.6, na.rm = T)) 
        tmp01.sig <- apply(wt$signif, 1, function(x) quantile(x, probs = 0.6, na.rm = T))
        tmp01.pow[tmp01.sig<1] <- NA
        
        # plot(period, tmp01.pow, type = "o")
        
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
          
          # vh_sum[r,c,2] <- round(tmp04[1,1], 1)
          # vh_sum[r,c,4] <- tmp04[1,2]
        }
        
        ### Power ~ Years
        yrs <- tmp2[,1]
        
        py <- foreach(loop=2:(length(unique(yrs))-1), .combine = rbind) %dopar% {
          
          tmp01.pow <- apply(power[,which(yrs==unique(yrs)[loop])], 1, function(x) quantile(x, probs = 0.6, na.rm = T)) 
          tmp01.sig <- apply(wt$signif[,which(yrs==unique(yrs)[loop])], 1, function(x) quantile(x, probs = 0.6, na.rm = T))
          tmp01.pow[tmp01.sig<1] <- NA
          
          # plot(period, tmp01.pow, type = "o")
          
          ind <- c(1, which((!is.na(tmp01.pow[-length(tmp01.pow)]) &  is.na(tmp01.pow[-1])) |
                              (is.na(tmp01.pow[-length(tmp01.pow)]) & !is.na(tmp01.pow[-1]))), length(tmp01.pow))
          
          tmp02 <- split(data.frame(period, tmp01.pow, tmp01.sig), f = cut(1:length(tmp01.pow), breaks = unique(ind)))
          
          tmp04 <- do.call(rbind, lapply(tmp02, pck.periods))[,1:3]
          tmp04 <- subset(tmp04, tmp01.pow>=0)
          if(nrow(tmp04)>0) {
            tmp04 <- tmp04[order(tmp04$tmp01.pow, decreasing = T),]
          } else tmp04 <- data.frame(period = NA, tmp01.pow = NA, tmp01.sig = NA)
          
          cbind(loop = unique(yrs)[loop], tmp04[1,])
        }
        
        out <- merge(data.frame(loop = 1:max(py[,1])), py, all.x = T)
          out$years = c(1982:2015)[1:nrow(out)]
        # vh_ts$Pow[r,c,1:nrow(out)] <- out[,3]
        
        
        if(sum(!is.na(out$period))>10 & mean(out$period, na.rm = T)<5) {
          
          ### (Overall) Predictability
          x0 <- loess(y ~ t, data = tmp2, span = 0.02)
          x1 <- predict(x0, newdata = data.frame(t = 1:nrow(tmp2)))
          
          ## normalize x
          Mx <- x1 - mean(x1, na.rm = T)
          
          
          if(mean(out$period, na.rm = T)>0.25 & mean(out$period, na.rm = T)<0.75) {
            fit0 <- optim(fn = leastS.cos.5, par = c(a = fit0$par[1], b = fit0$par[2]), sd = 0.001)
            if(fit0$par[1]==0){
              Mx <- dat$y
              fit0 <- optim(fn = leastS.cos, par = c(a = fit0$par[1], b = fit0$par[2]), sd = 0.001)
              Mx <- x1 - mean(x1, na.rm = T)
            }	
            curve.1 <- fit0$par[1]*cos(pi*((1:length(x1))/(length(x1)/((length(x1)/(52/2))*2))) + 
                                         (pi+fit0$par[2])) +  mean(x1, na.rm=T) 
            spl <- which(diff(curve.1[-length(curve.1)])<0 & diff(curve.1[-1])>0) +1
            tmp3 <- split(tmp2, f = cut(1:nrow(tmp2), breaks = c(0, spl[c(TRUE, FALSE)], nrow(tmp2))))
            
          } else {
            fit0 <- optim(fn = leastS.cos, par = c(a = fit0$par[1], b = fit0$par[2]), sd = 0.001)
            if(fit0$par[1]==0){
              Mx <- dat$y
              fit0 <- optim(fn = leastS.cos, par = c(a = fit0$par[1], b = fit0$par[2]), sd = 0.001)
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
          
          fc <- foreach(loop=1:nrow(ind.years), .combine = rbind, .packages = "forecast") %dopar% {
            ts <- ts(tmp4$y[((ind.years[loop,1]-1)*52):(((ind.years[loop,2])*52))], frequency = 52)
            ets <- stlf(ts, method = "arima", h = ind.years[loop,3])
            # plot(ets)
            cbind(c(ets$mean), as.matrix(ets$lower), as.matrix(ets$upper))
          }
          
          outTempR <- data.frame(years = tmp4$years[((4*52)+1):nrow(tmp4)][1:nrow(fc)],
                                 t =     tmp4[((4*52)+1):nrow(tmp4),][1:nrow(fc),1], 
                                 y =     tmp4$y[((4*52)+1):nrow(tmp4)][1:nrow(fc)], as.matrix(fc))
          
          outTempR <- outTempR[outTempR[,1]>4,]
          outTempR <- subset(outTempR, !is.na(y))
          ### Pr ~ years
          fcy <- foreach(loop=unique(outTempR$t), .combine = rbind) %dopar% {
            tmp <- outTempR[outTempR$t==loop,]
            RSS <- tmp$y - mean(tmp$y, na.rm = T)
            SSE <- sum((tmp$y-tmp[,4])^2, na.rm = T)
            
            cbind(yrs  = min(tmp$years),
                  amp  = as.numeric(abs(diff(quantile(tmp$y, probs = c(0.975, 0.025), na.rm = T)))),
                  pred = (sum(RSS^2, na.rm = T)-SSE)/sum(RSS^2, na.rm = T))
          }
          
          fcy <- fcy[!duplicated(fcy[,1], fromLast = T),]
          
          vh.s$Amp[r,c,match(fcy[,1], c(1982:2016))] <- fcy[,2]
          vh.s$Pre[r,c,match(fcy[,1], c(1982:2016))] <- fcy[,3]
          
        } # end predictability
      } # end !all(is.na(x))
    }
  } # is.na(cell)
} # end seasonality (cos != 0)

} # end column
  if(floor(r/5)==r/5) {
    
    png("VH.png", width = 400, height = 700)
    opar <- par(mfrow = c(2,1), mar = c(2,2,1,1))
    plot(raster(apply(vh.s$Amp, 1:2, mean, na.rm = T), xmn = -180, xmx = 180, ymn = -90, ymx = 90))
    plot(wrld_simpl, add = T)
    plot(raster(apply(vh.s$Pre, 1:2, mean, na.rm = T), xmn = -180, xmx = 180, ymn = -90, ymx = 90))
    plot(wrld_simpl, add = T)
    par(opar)
    dev.off()
    
    # save(vh.s, file = "Results/MODIS/vh_s.RData")
    
  }
} # end row
    
    
    