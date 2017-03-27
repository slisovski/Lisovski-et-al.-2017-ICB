library(raster)
library(maptools); data(wrld_simpl)
library(rgdal)

load("Revision_Analysis/Output/NPP/npp_monthly_seasonality.RData")

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
prec.color  <- colorRampPalette(c("white", "darkblue"))
vh.color  <- colorRampPalette(c("white", "green4"))
change     <- colorRampPalette(c("darkblue", "white", "darkred"))

### Earth Observatory
library(adobecolor)
  npp.rgb <- read_aco("Revision_Analysis/ColorScales/modis_npp.aco") 
  npp     <-  colorRampPalette(npp.rgb)

  pred.rgb <- read_aco("Revision_Analysis/ColorScales/predictability.aco") 
  pred     <-  colorRampPalette(rev(pred.rgb))

  DoS.rgb  <- read_aco("Revision_Analysis/ColorScales/DoS.aco") 
  DoS     <- colorRampPalette(DoS.rgb)

########################
#### Amplitude #########
########################
npp.y <- npp.s$Amp
npp.a <- apply(npp.y[,,6:10], 1:2, function(x) median(x, na.rm = T))

npp.amp1 <- raster(npp.a, xmn = -180, xmx = 180, ymn = -90, ymx = 90)
proj4string(npp.amp1) <- proj4string(wrld_simpl)
npp.amp1[] <- ifelse(!is.na(land.r[]), npp.amp1[], NA)

npp.amp <- projectRaster(npp.amp1, crs = CRS(proj))
brd      <- box(xlim = c(-168.4, 191.6), ylim = c(-60, 89), raster = npp.amp, proj = proj)  
npp.amp[] <- ifelse(!is.na(brd[]), npp.amp[], NA) 

plot(NA, xlim = extent(land.p)[1:2], ylim = extent(land.p)[3:4], bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
plot(land.p, col = "grey50", border = NA, add = T)
plot(npp.amp, breaks = seq(1, 7, length = 10), col = npp(9), add = T, legend = F)
par(new = T)  
plot(NA, xlim = extent(land.p)[1:2], ylim = extent(land.p)[3:4], bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
plot(land.p, add = T, lwd = 0.5, border = "grey10")
plot(water.p, add = T, col = "white", border = "grey50")


########################
#### Predictability ####
########################
npp.y <- npp.s$Pre
npp.p <- apply(npp.y[,,6:10], 1:2, function(x) median(x, na.rm = T))

npp.pred1 <- raster(npp.p, xmn = -180, xmx = 180, ymn = -90, ymx = 90)
proj4string(npp.pred1) <- proj4string(wrld_simpl)
npp.pred1[] <- ifelse(!is.na(land.r[]) & npp.amp1[]>1, npp.pred1[], NA)

npp.pred <- projectRaster(npp.pred1, crs = CRS(proj))
brd      <- box(xlim = c(-168.4, 191.6), ylim = c(-60, 89), raster = npp.pred, proj = proj)  
npp.pred[] <- ifelse(!is.na(brd[]), npp.pred[], NA) 

plot(NA, xlim = extent(land.p)[1:2], ylim = extent(land.p)[3:4], bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
plot(land.p, col = "grey50", border = NA, add = T)
plot(npp.pred, breaks = seq(0, 1, length = 10), 
     col =  apply(cbind(npp(9), seq(0, 1, length = 9)), 1, function(x) adjustcolor(x[1], alpha.f = x[2])), 
     add = T, legend = F)
par(new = T)  
plot(NA, xlim = extent(land.p)[1:2], ylim = extent(land.p)[3:4], bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
plot(land.p, add = T, lwd = 0.5, border = "grey10")
plot(water.p, add = T, col = "white", border = "grey50")


#################################
#### Degree of Seasonality ######
#################################
scale.amp <- function(x) {
  fn <- approxfun(x = seq(0, max(npp.amp[], na.rm = T), length = 10), y = seq(0, 1, length = 10))
  fn(x)
}

npp.seas  <- npp.amp
npp.seas[] <- scale.amp(npp.seas[])
npp.seas[] <- apply(cbind(npp.seas[], npp.pred[]), 1, function(x) ifelse(!any(is.na(x)), median(x), NA))


plot(NA, xlim = extent(land.p)[1:2], ylim = extent(land.p)[3:4], bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
plot(land.p, col = "grey50", border = NA, add = T)
plot(npp.seas, breaks = seq(0, 1, length = 10), col = seas.color(9), add = T, legend = F)
par(new = T)  
plot(NA, xlim = extent(land.p)[1:2], ylim = extent(land.p)[3:4], bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
plot(land.p, add = T, lwd = 0.5, border = "grey10")
plot(water.p, add = T, col = "white", border = "grey50")


save(npp.seas, file = "Revision_Analysis/Output/NPP/nppSeas.RData")


### Save as Shapefiles
DoSr <- npp.amp1
DoSr[] <- scale.amp(npp.amp1[])

DoSr[] <- round(apply(cbind(DoSr[], npp.pred1[]), 1, mean, na.rm = T),3)
npp.pred1[] <- round(npp.pred1[], 3)
npp.amp1[]  <- round(npp.amp1[], 3)

out_npp <- rasterToPolygons(brick(DoSr, npp.amp1, npp.pred1), dissolve = T)
names(out_npp) <- c("DegreeOfSeasonality", "Amplitude", "Predictability")

#plot(out_npp)
writeOGR(out_npp, ".", "seasonality-NPP", driver="ESRI Shapefile")









pdf("npp.pdf", height = 10, width = 22)
opar <- par(mfrow = c(2,2), mar = c(0,0,0,0), bty = "n")


plot(NA, xlim = extent(land.p)[1:2], ylim = extent(land.p)[3:4], bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
plot(land.p, col = "grey50", border = NA, add = T)
plot(npp.seas, breaks = seq(0, 1, length = 10), col = seas.color(9), add = T, legend = F)
par(new = T)  
plot(NA, xlim = extent(land.p)[1:2], ylim = extent(land.p)[3:4], bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
plot(land.p, add = T, lwd = 0.5, border = "grey10")
plot(water.p, add = T, col = "white", border = "grey50")


plot(NA, xlim = extent(land.p)[1:2], ylim = extent(land.p)[3:4], bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
plot(land.p, col = "grey50", border = NA, add = T)
plot(npp.amp, breaks = seq(1.25, 7, length = 10), col = npp(9), add = T, legend = F)
par(new = T)  
plot(NA, xlim = extent(land.p)[1:2], ylim = extent(land.p)[3:4], bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
plot(land.p, add = T, lwd = 0.5, border = "grey10")
plot(water.p, add = T, col = "white", border = "grey50")

plot(NA, xlim = extent(land.p)[1:2], ylim = extent(land.p)[3:4], bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
plot(land.p, col = "grey50", border = NA, add = T)
plot(npp.pred, breaks = seq(0, 1, length = 10), 
     col =  apply(cbind(npp(9), seq(0, 1, length = 9)), 1, function(x) adjustcolor(x[1], alpha.f = x[2])), 
     add = T, legend = F)
par(new = T)  
plot(NA, xlim = extent(land.p)[1:2], ylim = extent(land.p)[3:4], bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
plot(land.p, add = T, lwd = 0.5, border = "grey10")
plot(water.p, add = T, col = "white", border = "grey50")


par(opar)
dev.off()

plot(1:9, rep(0, 9), pch = 15, col = npp(9), cex = 12, xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
plot(1:9, rep(0, 9), pch = 15, col = "grey50", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", cex = 12)
points(1:9, rep(0, 9), pch = 15, col = apply(cbind( npp(9), seq(0, 1, length = 9)), 1, function(x) adjustcolor(x[1], alpha.f = x[2])),
       cex = 12)




## Latitude
library(zoo)


matr.npp <- as.matrix(npp.seas)
crds <- coordinates(npp.seas)
lat <- seq(max(crds[,2]), min(crds[,2]), length = nrow(npp.seas))

ax <- data.frame(longlat = c(-75, -50, -25, 0, 25, 50, 75),
                 moll = project(cbind(rep(0, 7), c(-75, -50, -25, 0, 25, 50, 75)), proj = proj4string(npp.seas)))


rmean_50  <- rollapply(matr.npp, 30, median, by.column = F, na.rm = T, fill  = NA)
rmean_sd  <- rollapply(matr.npp, 30, sd, by.column = F, na.rm = T, fill  = NA)

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


bm.p <- resample(bm, npp.pred, method = "ngb")

ecor <- cbind(npp.seas[], bm.p[])
ecor[,1] <- ifelse(ecor[,1]>0.05, ecor[,1], NA)
eco.sum <- aggregate(ecor[,1], by = list(ecor[,2]), function(x) quantile(x, probs = c(0.975, 0.6, 0.5, 0.4, 0.025), na.rm = T))
eco.sum$perc <- aggregate(ecor[,1], by = list(ecor[,2]), function(x) sum(!is.na(x) & x>0.05)/length(x))$x

pdf("EcoR_sum_npp.pdf", height = 6, width = 5) 
boxplot(ecor[,1]~ecor[,2], outline = F, names = paste(round(eco.sum$perc,2), "%"), border = col.l, lty = 1,
        col = col.l, lwd = 2, boxwex = rep(0.35, 7),
        ylim = c(0, 1), bty = "n", las = 2, ylab = "", xlab = "")
points(eco.sum$Group.1, eco.sum$x[,3], col = "white", pch = "-", cex = 6)
dev.off()
