library(raster)
library(maptools); data(wrld_simpl)
library(rgdal)

load("Revision_Analysis/Output/Temperature/temp_monthly_seasonality.RData")

### Map data
wd <- "~/Dropbox"
load(paste(wd, "Science/Projects/Seasonality_Thesis/Analysis/Map/land_l.RData", sep = "/"))
plot(land)
load(paste(wd, "Science/Projects/Seasonality_Thesis/Analysis/Map/water_l.RData", sep = "/"))
plot(water, add = T)

proj = "+proj=loxim +lon_0=11.6"
land.p  <- spTransform(land, proj)
water.p <- spTransform(water, proj)

plot(land.p)
plot(water.p, add  = T)

box <- function(xlim, ylim, raster, proj) {
  xm <- c(seq(xlim[1], xlim[2], length = 100), rep(xlim[2], length = 100), 
          seq(xlim[2], xlim[1], length = 100), rep(xlim[1], length = 100), xlim[1])
  ym<- c(rep(ylim[1], length = length(seq(xlim[1], xlim[2], length = 100))), seq(ylim[1]+0.01, ylim[2], length = 100), 
         rep(ylim[2], length = length(seq(xlim[2], xlim[1], length = 100))),
         seq(ylim[2]-0.01, ylim[1], length = 100), ylim[1])
  
  sp1 <- SpatialPolygons(list(Polygons(list(Polygon(cbind(xm, ym))), 1)), proj4string=CRS("+proj=longlat"))
  box <- spTransform(sp1, CRS(proj))
  rasterize(box, raster)
}

### Out raster
r.out  <- raster(xmn = -180, xmx = 180, ymn = -90, ymx = 90, res = c(1,1))
proj4string(r.out) <- proj4string(wrld_simpl)
###

### Ocean mask
land.r <- rasterize(wrld_simpl, r.out)
land.r[] <- ifelse(land.r[]>0, land.r[], NA)
plot(land.r)

### Color Palette
seas.color <- colorRampPalette(c("aliceblue", "purple4"))
temp.color  <- colorRampPalette(c("white", "darkred"))
    temp.pred.col <- colorRampPalette(c("white", "darkred"))
prec.color  <- colorRampPalette(c("white", "darkblue"))
vh.color  <- colorRampPalette(c("white", "green4"))
change     <- colorRampPalette(c("darkblue", "white", "darkred"))
  
  ### Earth Observatory
  library(adobecolor)
  temp.rgb <- read_aco("Revision_Analysis/ColorScales/modis_lst.aco") 
  temp     <-  colorRampPalette(temp.rgb)
  
  pred.rgb <- read_aco("Revision_Analysis/ColorScales/predictability.aco") 
  pred     <-  colorRampPalette(rev(pred.rgb))

  DoS.rgb  <- read_aco("Revision_Analysis/ColorScales/DoS.aco") 
  DoS     <- colorRampPalette(DoS.rgb)

########################
#### Amplitude #########
########################
temp.y <- temp.s$Amp
temp.a <- apply(temp.y[,,6:10], 1:2, function(x) median(x, na.rm = T))

temp.amp1 <- raster(temp.a, xmn = -180, xmx = 180, ymn = -90, ymx = 90)
  proj4string(temp.amp1) <- proj4string(wrld_simpl)
temp.amp1[] <- ifelse(!is.na(land.r[]), temp.amp1[], NA)

temp.amp <- projectRaster(temp.amp1, crs = CRS(proj))
  brd      <- box(xlim = c(-168.4, 191.6), ylim = c(-60, 89), raster = temp.amp, proj = proj)  
temp.amp[] <- ifelse(!is.na(brd[]), temp.amp[], NA) 

plot(NA, xlim = extent(land.p)[1:2], ylim = extent(land.p)[3:4], bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
plot(land.p, col = "grey50", border = NA, add = T)
  plot(temp.amp, breaks = seq(0, 63, length = 10), col = temp(9), add = T, legend = F)
par(new = T)  
plot(NA, xlim = extent(land.p)[1:2], ylim = extent(land.p)[3:4], bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
plot(land.p, add = T, lwd = 0.5, border = "grey10")
plot(water.p, add = T, col = "white", border = "grey50")


# plot(1:9, rep(0, 9), pch = 15, col = temp.color(9), cex = 10)


########################
#### Predictability ####
########################
temp.y <- temp.s$Pre
temp.p <- apply(temp.y[,,6:10], 1:2, function(x) median(x, na.rm = T))

temp.pred1 <- raster(temp.p, xmn = -180, xmx = 180, ymn = -90, ymx = 90)
  proj4string(temp.pred1) <- proj4string(wrld_simpl)
temp.pred1[] <- ifelse(!is.na(land.r[]), temp.pred1[], NA)

temp.pred <- projectRaster(temp.pred1, crs = CRS(proj))
  brd      <- box(xlim = c(-168.4, 191.6), ylim = c(-60, 89), raster = temp.pred, proj = proj)  
temp.pred[] <- ifelse(!is.na(brd[]), temp.pred[], NA) 

plot(NA, xlim = extent(land.p)[1:2], ylim = extent(land.p)[3:4], bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
plot(land.p, col = "grey50", border = NA, add = T)
plot(temp.pred, breaks = seq(0, 1, length = 10), 
     col = apply(cbind(temp.color(9), seq(0, 1, length = 9)), 1, function(x) adjustcolor(x[1], alpha.f = x[2])), 
     add = T, legend = F)
par(new = T)  
plot(NA, xlim = extent(land.p)[1:2], ylim = extent(land.p)[3:4], bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
plot(land.p, add = T, lwd = 0.5, border = "grey10")
plot(water.p, add = T, col = "white", border = "grey50")




#################################
#### Degree of Seasonality ######
#################################
scale.amp <- function(x) {
  fn <- approxfun(x = seq(0, max(temp.amp[], na.rm = T), length = 10), y = seq(0, 1, length = 10))
  fn(x)
}

temp.seas  <- temp.amp
temp.seas[] <- scale.amp(temp.seas[])
temp.seas[] <- apply(cbind(temp.seas[], temp.pred[]), 1, function(x) ifelse(!any(is.na(x)), median(x), NA))


plot(NA, xlim = extent(land.p)[1:2], ylim = extent(land.p)[3:4], bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
plot(land.p, col = "grey50", border = NA, add = T)
plot(temp.seas, breaks = seq(0, 1, length = 10), col = seas.color(9), add = T, legend = F)
par(new = T)  
plot(NA, xlim = extent(land.p)[1:2], ylim = extent(land.p)[3:4], bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
plot(land.p, add = T, lwd = 0.5, border = "grey10")
plot(water.p, add = T, col = "white", border = "grey50")


# plot(1:9, rep(0, 9), pch = 15, col = temp.color(9), cex = 12, xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
# plot(1:9, rep(0, 9), pch = 15, col = "grey50", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", cex = 12)
# points(1:9, rep(0, 9), pch = 15, col = apply(cbind(temp.color(9), seq(0, 1, length = 9)), 1, function(x) adjustcolor(x[1], alpha.f = x[2])),
#      cex = 12)


save(temp.seas, file = "Revision_Analysis/Output/Temperature/TempSeas.RData")


### Save as Shapefiles
DoSr <- temp.amp1
DoSr[] <- scale.amp(temp.amp1[])

DoSr[] <- round(apply(cbind(DoSr[], temp.pred1[]), 1, mean, na.rm = T),3)
temp.pred1[] <- round(temp.pred1[], 3)
temp.amp1[]  <- round(temp.amp1[], 3)

out_temp <- rasterToPolygons(brick(DoSr, temp.amp1, temp.pred1), dissolve = T)
names(out_temp) <- c("DegreeOfSeasonality", "Amplitude", "Predictability")

#plot(out_npp)
writeOGR(out_temp, ".", "seasonality-Temperature", driver="ESRI Shapefile")









pdf("temp.pdf", height = 10, width = 22)
opar <- par(mfrow = c(2,2), mar = c(0,0,0,0), bty = "n")


plot(NA, xlim = extent(land.p)[1:2], ylim = extent(land.p)[3:4], bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
plot(land.p, col = "grey50", border = NA, add = T)
plot(temp.seas, breaks = seq(0, 1, length = 10), col = seas.color(9), add = T, legend = F)
par(new = T)  
plot(NA, xlim = extent(land.p)[1:2], ylim = extent(land.p)[3:4], bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
plot(land.p, add = T, lwd = 0.5, border = "grey10")
plot(water.p, add = T, col = "white", border = "grey50")

plot(NA, xlim = extent(land.p)[1:2], ylim = extent(land.p)[3:4], bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
plot(land.p, col = "grey50", border = NA, add = T)
plot(temp.amp, breaks = seq(0, 63, length = 10), col = temp.color(9), add = T, legend = F)
par(new = T)  
plot(NA, xlim = extent(land.p)[1:2], ylim = extent(land.p)[3:4], bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
plot(land.p, add = T, lwd = 0.5, border = "grey10")
plot(water.p, add = T, col = "white", border = "grey50")

plot(NA, xlim = extent(land.p)[1:2], ylim = extent(land.p)[3:4], bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
plot(land.p, col = "grey50", border = NA, add = T)
plot(temp.pred, breaks = seq(0, 1, length = 10), 
     col = apply(cbind(temp.color(9), seq(0, 1, length = 9)), 1, function(x) adjustcolor(x[1], alpha.f = x[2])), 
     add = T, legend = F)
par(new = T)  
plot(NA, xlim = extent(land.p)[1:2], ylim = extent(land.p)[3:4], bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
plot(land.p, add = T, lwd = 0.5, border = "grey10")
plot(water.p, add = T, col = "white", border = "grey50")


par(opar)
dev.off()

## Latitude
library(zoo)


matr.temp <- as.matrix(temp.seas)
crds <- coordinates(temp.seas)
lat <- seq(max(crds[,2]), min(crds[,2]), length = nrow(temp.seas))

ax <- data.frame(longlat = c(-75, -50, -25, 0, 25, 50, 75),
                 moll = project(cbind(rep(0, 7), c(-75, -50, -25, 0, 25, 50, 75)), proj = proj4string(temp.seas)))


rmean_50  <- rollapply(matr.temp, 30, median, by.column = F, na.rm = T, fill  = NA)
rmean_sd  <- rollapply(matr.temp, 30, sd, by.column = F, na.rm = T, fill  = NA)

png("temp.lat.png", width = 500, height = 800, res = 100)
plot(rmean_50, lat, type = "l", bty = "n", xaxt = "n", yaxt = "n", lwd = 3, 
     ylab = "", xlab = "", xlim = c(0, 1))
lines(rmean_50-rmean_sd, lat, type = "l", bty = "n", lwd = 2, lty = 3)
lines(rmean_50+rmean_sd, lat, type = "l", bty = "n", lwd = 2, lty = 3)

axis(1, at = seq(0, 1, length = 5), lwd = 1.4, cex.axis = 1.4)
axis(2, at = ax[,3], labels = ax[,1], las = 1, lwd = 1.4, cex.axis = 1.4, line = 1)
dev.off()

#### Ecoregions
load("Metadata/tcn_terr/biome.RData")

col.l <- c("darkblue", "cadetblue4", "darkolivegreen2", "darkmagenta", "burlywood1",
           "darkgreen", "darkorange")


bm <- projectRaster(biome.r, crs = CRS(proj), method = "ngb")
brd      <- box(xlim = c(-168.4, 191.6), ylim = c(-60, 89), raster = bm, proj = proj)  
bm[] <- ifelse(!is.na(brd[]), bm[], NA) 


# pdf("EcoR.pdf", height = 7.5, width = 7.5)
# plot(bm, breaks = 0:7, col = col.l)
# plot(land.p, add = TRUE, 
#      lwd = 0.4, border = "grey30")
# plot(water.p, col = "white", lwd = 0.3, add = T, border = "grey20")
# dev.off()


bm.p <- resample(bm, temp.pred, method = "ngb")

ecor <- cbind(temp.seas[], bm.p[])
ecor[,1] <- ifelse(ecor[,1]>0.05, ecor[,1], NA)
eco.sum <- aggregate(ecor[,1], by = list(ecor[,2]), function(x) quantile(x, probs = c(0.975, 0.6, 0.5, 0.4, 0.025), na.rm = T))
eco.sum$perc <- aggregate(ecor[,1], by = list(ecor[,2]), function(x) sum(!is.na(x) & x>0.05)/length(x))$x

pdf("EcoR_sum_temp.pdf", height = 6, width = 5) 
boxplot(ecor[,1]~ecor[,2], outline = F, names = paste(round(eco.sum$perc,2), "%"), border = col.l, lty = 1,
        col = col.l, lwd = 2, boxwex = rep(0.35, 7),
        ylim = c(0, 1), bty = "n", las = 2, ylab = "", xlab = "")
points(eco.sum$Group.1, eco.sum$x[,3], col = "white", pch = "-", cex = 6)
dev.off()
