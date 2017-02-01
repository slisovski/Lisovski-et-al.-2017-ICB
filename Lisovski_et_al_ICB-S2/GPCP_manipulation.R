library(raster)
library(maptools)
  data("wrld_simpl")

fls <- list.files("~/Dropbox/Data_Archive/GPCP/gpcp_formatted")

tm <- as.Date(seq(as.POSIXct("1997-01-01"), as.POSIXct("2015-10-31"), by = "day"))

precDay <- array(dim = c(180, 360, length(tm)))  

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
    precDay[,,which(tm==tmM)] <- m0
  }
}
# save(precDay, file = "temp_data/precDay_array.RData")  
load("temp_data/precDay_array.RData")


#### 
## Into weekly data

date.ind <- data.frame(tm = as.Date(seq(as.POSIXct("1997-01-01"), as.POSIXct("2015-12-31"), by = "day")))
  date.ind$Year = as.numeric(format(date.ind$tm, "%Y"))
  date.ind$Week <- unlist(lapply(split(date.ind, f = date.ind$Year), function(x) rep(1:53, each = 7)[1:nrow(x)]))
  date.ind$Week[date.ind$Week==53] <- 52
  date.ind$Week.ind <- apply(date.ind, 1, function(x) paste(x[2], as.numeric(x[3]), sep = "_"))
  

precW <- array(dim = c(180, 360, length(unique(date.ind$Week.ind))))
  dimnames(precW)[3][[1]] <- unique(date.ind$Week.ind)
  
for(r in 1:180) {
  for(c in 1:360) {
    
    cat(paste(rep("\b\b\b\b\b\b", 8), collapse = ""))
    cat(r, " - ", c, "\r")
    
    # opar <- par(mfrow = c(2,1), mar = c(2,2,1,1))
    tmp <- precDay[r,c,]
      tmp[tmp<0] <- NA
      tmp <- c(tmp, rep(NA, 150))[1:nrow(date.ind)]
    # plot(tmp, type = "h")
    out <- aggregate(tmp, by = list(date.ind$Week.ind), function(x) quantile(x, probs = 0.95, na.rm = T))
    out <- merge(data.frame(Group.1 = unique(date.ind$Week.ind)), out, all.x = T, sort = F)
    # plot(out$x, type = "h")
    # par(opar)
    precW[r,c,] <- out$x
  }
}

save(precW, file = "temp_data/precWeek_array.RData")  
