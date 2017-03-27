############################################################################
### Project: Defining the Degree of Seasonality (Lisovski et al 2017 ICB) ##
### Author: Simeon Lisovski                                               ##
### Date: March 2017                                                      ##
############################################################################
### Data manipulation: Net Primary Productivity                           ##
### NPP (original: weelky observations on a 0.1x0.1 degree grid)          ##
### Output: Monthly means on a 1x1 degree grid                            ##
############################################################################

library(raster)
library(rgdal)
library(zoo)
library(forecast)
library(RNetCDF)
library(maptools)
  data(wrld_simpl)

flsN <- list.files("~/Dropbox/Data_Archive/PSN")
   date   <- as.POSIXct(substring(flsN, 15, 24), tz = "GMT")
   
   index <- data.frame(Date = date, Year = as.numeric(format(date, "%Y")), 
                       Yday = as.numeric(format(date, "%j")),
                       NPP  = 1:length(flsN))

  index <- subset(index, Year%in%c(2006:2015))
 
nppA <- array(dim = c(180, 360, nrow(index)))

  r.out <- raster(extent(c(-180, 180, -90, 90)), res = c(1,1))  
           proj4string(r.out) <- proj4string(wrld_simpl)
 
for(i in 1:nrow(index)) {
   
 cat("\b\b\b\b\b\b\b\b\b\b\b\b")
 cat(i, " of ", nrow(index), "\r")
   
 if(!is.na(index$NPP[i])){
   tmp <- raster(paste0("~/Dropbox/Data_Archive/PSN/", flsN[index$NPP[i]]))
     tmp[] <- ifelse(tmp[]>100, NA, tmp[])
     
     tmp1 <- raster:::resample(tmp, r.out, method = "ngb")
 
     nppA[,,i] <-  as.matrix(tmp1)
   }
}

save(nppA, file = "Revision_Analysis/Output/NPP/npp_weekly_array.RData")
load("Revision_Analysis/Output/NPP/npp_weekly_array.RData")


### Monthly means
month.ind <- seq(as.POSIXct("2006-01-15"), as.POSIXct("2015-12-15"), by = "month")

x <- nppA[25,50,]
  plot(index$Date, x, type = "o", pch = 16, cex = 0.5)

nppA_MNT <- array(dim = c(180, 360, length(month.ind)))
  
for(r in 1:dim(nppA_MNT)[1]) {
  for(c in 1:dim(nppA_MNT)[2]) {
    
    cat("\b\b\b\b\b\b\b\b\b\b\b\b")
    cat(r, " - ", c, "\r")
    
    x <- nppA[r,c,]
    if(sum(is.na(x))>(length(x)/2)) {
      nppA_MNT[r,c,] <- rep(NA, length(month.ind))
    } else {
      x <- na.approx(x, rule = 2)
      tmp <- suppressWarnings(cbind(index$Date, predict(loess(x~as.numeric(index$Date), span = 0.025))))
      fnc <- approxfun(tmp)
      nppA_MNT[r,c,] <- fnc(as.numeric(month.ind))
    }
  }
}

save(nppA_MNT, file = "Revision_Analysis/Output/NPP/npp_monthly_array.RData")
load("Revision_Analysis/Output/NPP/npp_monthly_array.RData")


### Monthly means
month.ind <- seq(as.POSIXct("2006-01-15"), as.POSIXct("2015-12-15"), by = "month")

x <- nppA[25,50,]
plot(index$Date, x, type = "o", pch = 16, cex = 0.5)

nppA_MNT <- array(dim = c(180, 360, length(month.ind)))

for(r in 1:dim(nppA_MNT)[1]) {
  for(c in 1:dim(nppA_MNT)[2]) {
    
    cat("\b\b\b\b\b\b\b\b\b\b\b\b")
    cat(r, " - ", c, "\r")
    
    x <- nppA[r,c,]
    if(sum(is.na(x))>(length(x)/2)) {
      nppA_MNT[r,c,] <- rep(NA, length(month.ind))
    } else {
      x <- na.approx(x, rule = 2)
      tmp <- suppressWarnings(cbind(index$Date, predict(loess(x~as.numeric(index$Date), span = 0.025))))
      fnc <- approxfun(tmp)
      nppA_MNT[r,c,] <- fnc(as.numeric(month.ind))
    }
  }
}

save(nppA_MNT, file = "Revision_Analysis/Output/NPP/npp_monthly_array.RData")
load("Revision_Analysis/Output/NPP/npp_monthly_array.RData")










