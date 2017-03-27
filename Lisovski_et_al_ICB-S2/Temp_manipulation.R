############################################################################
### Project: Defining the Degree of Seasonality (Lisovski et al 2017 ICB) ##
### Author: Simeon Lisovski                                               ##
### Date: March 2017                                                      ##
############################################################################
### Data manipulation: Temperature                                        ##
### NPP (original: monthly observations on a 1x1 degree grid)             ##
### Output: Monthly array                                                 ##
############################################################################

library(raster)
library(rgdal)
library(zoo)
library(forecast)
library(RNetCDF)
library(maptools)
data(wrld_simpl)


### files
fls <- list.files("Data/", pattern = "air.mon.mean.nc")


## Dates
tmp <- open.nc(paste0("Data/", fls))
dat  <- read.nc(tmp)
tm   <- seq(as.POSIXct("1948-01-01"), as.POSIXct("2017-02-01"), by = "month")
lon  <- dat$lon
lat  <- dat$lat
temp <- dat$air
close.nc(tmp)

r.out <- raster(xmn = -180, xmx = 180, ymn = -90, ymx = 90, res = c(1,1))

tempA_MNT <- array(dim = c(180, 360, dim(dat$air)[3]))

for(i in 1:dim(dat$air)[3]) {
  
  cat(paste(rep("\b\b\b\b\b\b\b\b\b", 4), collapse = ""))
  cat(which((1:dim(dat$air)[3])==i), " - ", length(1:dim(dat$air)[3]), "\r")
  
  z   <- temp[,,i]
  tmp <- raster(cbind(t(z)[,((720/2)+1):720], t(z)[,1:(720/2)]))
  extent(tmp) <- c(-179.75, 179.75, -89.75, 89.75)  
  proj4string(tmp) <- proj4string(wrld_simpl)
  out <- resample(tmp, r.out)
  
  tempA_MNT[,,i] <- as.matrix(out)
  
}

save(tempA_MNT, file = "Revision_Analysis/Output/Temperature/temp_monthly_array.RData")
load("Revision_Analysis/Output/Temperature/temp_monthly_array.RData")