## RNCEP II: Temperature (2x2) daily to weekly data manipulation

library(RNetCDF)
library(raster)
library(maptools)
  data("wrld_simpl")

  
### files
fls <- list.files("Data/NCEP.air.2m.gauss", pattern = ".nc")


  ## Dates
  tm <- c()
  for(i in fls) {
    tmp <- open.nc(paste0("Data/NCEP.air.2m.gauss/", i))
    # print.nc(tmp)
    dat<-read.nc(tmp)
    tm <- c(tm, as.POSIXct(dat$time*60*60, origin = "1800-1-1 00:00:0.0"))
    close.nc(tmp)
  }
  tm <- as.POSIXct(tm, origin = "1970-01-01")
  
  date.ind <- data.frame(Date = tm, 
                         Year = as.numeric(format(tm, "%Y")))
  date.ind$Week <- unlist(lapply(split(date.ind, f = date.ind$Year), function(x) rep(1:53, each = 7)[1:nrow(x)]))
  date.ind$Week[date.ind$Week==53] <- 52
  date.ind$Year[1] <- 1982
  date.ind$Week.ind <- apply(date.ind, 1, function(x) paste(x[2], as.numeric(x[3]), sep = "_"))
  
### Weekly array
  tmp <- open.nc(paste0("Data/NCEP.air.2m.gauss/", fls[1]))
  dat <- read.nc(tmp)
  lon <- dat$lon
  lat <- dat$lat
  close.nc(tmp)  
  
tempD <- array(dim = c(90, 180, length(unique(date.ind$Week.ind))))
  dimnames(tempD)[3][[1]] <- unique(date.ind$Week.ind)
  
  r.out <- raster(xmn = -180, xmx = 180, ymn = -90, ymx = 90, res = c(2,2))
  
for(i in fls) {
  
  cat(paste(rep("\b\b\b\b\b\b\b\b\b", 4), collapse = ""))
  cat(which(fls==i), " - ", length(fls), "\r")
  
  
  tmp <- open.nc(paste0("Data/NCEP.air.2m.gauss/", i))
  dat <-read.nc(tmp)
  close.nc(tmp) 
  
  z   <- dat$air
  yr  <- as.numeric(format(as.POSIXct(dat$time*60*60, origin = "1800-1-1 00:00:0.0"), "%Y"))[5]
  
  out.tmp <- array(dim = c(dim(z)[1], dim(z)[2], 52))
  
  for(r in 1:dim(z)[1]) {
    for(c in 1:dim(z)[2]) {
      ind <- as.numeric(substring(date.ind[date.ind$Year==yr, "Week.ind"], 6, 7))
      tmp <- aggregate(z[r,c,], by = list(ind[1:dim(z)[3]]), mean)
      out.tmp[r,c,] <- tmp$x[1:52]
    }
  }

  for(j in 1:dim(out.tmp)[3]) {
    tmp  <- raster(t(out.tmp[,,j])[,c(97:192, 1:96)], xmn = -180, xmx = 180, ymn = -90, ymx = 90)
    tmp2 <- raster:::resample(tmp, r.out)
    tempD[,,which(as.numeric(substring(dimnames(tempD)[3][[1]], 1, 4))==yr)[j]] <- as.matrix(tmp2)
  }
}


tempWeek <- tempD[,,-which(substring(dimnames(tempD)[3][[1]], 1, 4)=="2016")]
save(tempWeek, file = "Results/Temperature/tempWeek_array.RData")  
  