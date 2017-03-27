############################################################################
### Project: Defining the Degree of Seasonality (Lisovski et al 2017 ICB) ##
### Author: Simeon Lisovski                                               ##
### Date: March 2017                                                      ##
############################################################################
### Data manipulation: Precipitation                                      ##
### GPCP (original: daily observations on a 0.1x0.1 degree grid)          ##
### Output: Monthly means on a 1x1 degree grid                            ##
############################################################################

library(raster)
library(rgdal)
library(zoo)
library(forecast)
library(RNetCDF)
library(maptools)
data(wrld_simpl)


fls <- list.files("~/Dropbox/Data_Archive/GPCP/gpcp_formatted")

tm <- as.POSIXct(seq(as.POSIXct("1997-01-01"), as.POSIXct("2015-10-31"), by = "day"))

precA <- array(dim = c(180, 360, length(tm)))  

for(i in 197:length(fls)) {
  
  cat(paste(rep("\b\b\b\b\b\b", 8), collapse = ""))
  cat(i, " - ", length(fls), "\r")
  
  to.read <- file(paste0("~/Dropbox/Data_Archive/GPCP/gpcp_formatted/", fls[i]), "rb") 
  data <- readBin(to.read, numeric(), size=4, n=31*180*360)
  close(to.read)
  
  l.data <- length(data)
  
  month <- matrix(data, ncol = 180*360, nrow = l.data/360/180, byrow = T)
  
  for(m in 1:(l.data/360/180)) {
    m0  <- t(matrix(month[m,], nrow = 360, ncol = 180))[,c(180:360, 1:179)]
    tmM <- as.Date(as.POSIXct(paste(substring(fls[i],1,4), substring(fls[i], 5,6), m, sep = "-")))
    precA[,,which(tm==tmM)] <- m0
  }
}

save(precA, file = "Revision_Analysis/Output/PRecipitation/prec_day_array.RData")
load("Revision_Analysis/Output/PRecipitation/prec_day_array.RData")



### Monthly means
month.ind <- seq(as.POSIXct("1997-01-15"), as.POSIXct("2015-10-15"), by = "month")

x <- precA[25,50,]
  x <- ifelse(x<0, NA, x)
plot(tm, x, type = "o", pch = 16, cex = 0.5)

precA_MNT <- array(dim = c(180, 360, length(month.ind)))

for(r in 1:dim(precA_MNT)[1]) {
  for(c in 1:dim(precA_MNT)[2]) {
    
    cat("\b\b\b\b\b\b\b\b\b\b\b\b")
    cat(r, " - ", c, "\r")
    
    x <- precA[r,c,]
      x <- ifelse(x<0, NA, x)
    if(sum(is.na(x))>(length(x)/2)) {
      precA_MNT[r,c,] <- rep(NA, length(month.ind))
    } else {
      x <- na.approx(x, rule = 2)
      tmp <- suppressWarnings(cbind(tm, predict(loess(x~as.numeric(tm), span = 0.025))))
      fnc <- approxfun(tmp)
      precA_MNT[r,c,] <- fnc(as.numeric(month.ind))
    }
  }
}

save(precA_MNT, file = "Revision_Analysis/Output/Precipitation/prec_monthly_array.RData")
load("Revision_Analysis/Output/Precipitation/prec_monthly_array.RData")
