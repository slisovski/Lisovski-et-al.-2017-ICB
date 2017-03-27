library(raster)
library(maptools); data(wrld_simpl)
library(rgdal)

load("Revision_Analysis/Output/Temperature/tempSeas.RData")
load("Revision_Analysis/Output/Precipitation/precSeas.RData")
load("Revision_Analysis/Output/NPP/nppSeas.RData")

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


### Landmass
proj = "+proj=loxim +lon_0=11.6"

wrld0  <- rasterize(wrld_simpl, raster(xmn = -180, xmx = 180, ymn = -90, ymx = 90, res = c(0.5, 0.5)))
wrld <- projectRaster(wrld0, crs = CRS(proj))
  brd <- box(xlim = c(-168.4, 191.6), ylim = c(-89.5, 89.5), raster = wrld, proj = proj)  
wrld[] <- ifelse(!is.na(brd[]), wrld[], NA) 
  
plot(wrld)

###

lmass <- data.frame(lat   = seq(85, -85, length = 18),
                    terr  = aggregate(temp.seas[], by = list(cut(project(coordinates(temp.seas), proj, inv = T)[,2], breaks = seq(-90, 90, by = 10), lables = T)), 
                                      function(x) sum(!is.na(x)))$x)  

lmass$temp <- aggregate(temp.seas[], by = list(cut(project(coordinates(temp.seas), proj, inv = T)[,2], breaks = seq(-90, 90, by = 10), lables = T)), 
                        function(x) sum(!is.na(x) & x>0.05))$x


lmass$prec <- aggregate(prec.seas[], by = list(cut(project(coordinates(prec.seas), proj, inv = T)[,2], breaks = seq(-90, 90, by = 10), lables = T)), 
                        function(x) sum(!is.na(x) & x>0.05))$x


lmass$vh <- aggregate(npp.seas[], by = list(cut(project(coordinates(npp.seas), proj, inv = T)[,2], breaks = seq(-90, 90, by = 10), lables = T)), 
                        function(x) sum(!is.na(x) & x>0.05))$x


total <- apply(lmass[,-1], 2, sum)
  perc <- (c(total)/c(total[1]))*100

### Temperature
out  <- data.frame(lat  = seq(-85, 85, length = 18),
                   t.med = aggregate(temp.seas[], by = list(cut(project(coordinates(temp.seas), proj, inv = T)[,2], 
                                                            breaks = seq(-90, 90, by = 10), lables = T)), function(x) mean(x[x>0], na.rm = T))$x)

out$t.sd    <- aggregate(temp.seas[], by = list(cut(project(coordinates(temp.seas), proj, inv = T)[,2], 
                                                breaks = seq(-90, 90, by = 10), lables = T)), function(x) sd(x[x>0], na.rm = T))$x

out$tamp.med <- aggregate(temp.amp[], by = list(cut(project(coordinates(temp.seas), proj, inv = T)[,2], 
                                                     breaks = seq(-90, 90, by = 10), lables = T)), function(x) mean(x[x>0], na.rm = T))$x

out$tpred.med <- aggregate(temp.pred[], by = list(cut(project(coordinates(temp.seas), proj, inv = T)[,2], 
                                                    breaks = seq(-90, 90, by = 10), lables = T)), function(x) mean(x[x>0], na.rm = T))$x


### Precipitation
out$pr.med <- aggregate(prec.seas[], by = list(cut(project(coordinates(temp.seas), proj, inv = T)[,2], 
                                                           breaks = seq(-90, 90, by = 10), lables = T)), function(x) mean(x[x>0], na.rm = T))$x

out$pr.sd  <- aggregate(prec.seas[], by = list(cut(project(coordinates(temp.seas), proj, inv = T)[,2], 
                                                breaks = seq(-90, 90, by = 10), lables = T)), function(x) sd(x[x>0], na.rm = T))$x



out$Pramp.med <- aggregate(prec.amp[], by = list(cut(project(coordinates(temp.seas), proj, inv = T)[,2], 
                                                    breaks = seq(-90, 90, by = 10), lables = T)), function(x) mean(x[x>0], na.rm = T))$x

out$Prpred.med <- aggregate(prec.pred[], by = list(cut(project(coordinates(temp.seas), proj, inv = T)[,2], 
                                                      breaks = seq(-90, 90, by = 10), lables = T)), function(x) mean(x[x>0], na.rm = T))$x


### VH
npp.amp[] <- ifelse(npp.amp[]<1.25, NA, npp.amp[])
npp.pred[] <- ifelse(npp.amp[]<1.25, NA, npp.pred[])

out$vh.med <- aggregate(npp.seas[], by = list(cut(project(coordinates(temp.seas), proj, inv = T)[,2], 
                                                   breaks = seq(-90, 90, by = 10), lables = T)), function(x) mean(x[x>0], na.rm = T))$x

out$vh.sd  <- aggregate(npp.seas[], by = list(cut(project(coordinates(temp.seas), proj, inv = T)[,2], 
                                                   breaks = seq(-90, 90, by = 10), lables = T)), function(x) sd(x[x>0], na.rm = T))$x

out$VHamp.med <- aggregate(npp.amp[], by = list(cut(project(coordinates(temp.seas), proj, inv = T)[,2], 
                                                     breaks = seq(-90, 90, by = 10), lables = T)), function(x) mean(x[x>0], na.rm = T))$x

out$VHpred.med <- aggregate(npp.pred[], by = list(cut(project(coordinates(temp.seas), proj, inv = T)[,2], 
                                                       breaks = seq(-90, 90, by = 10), lables = T)), function(x) mean(x[x>0], na.rm = T))$x




pdf("Fig3.pdf", width = 5, height = 18, useDingbats=FALSE)
opar <- par(mfrow = c(3,1), mar = c(4,4,4,1))
bp <- barplot(height = rbind(lmass$temp*(res(temp.seas)[1]^2)/1000000, 
                             (lmass$terr-lmass$temp)*(res(temp.seas)[1]^2)/1000000), horiz = T, border = NA, col = c("grey50", "grey80"))  
lons <- approxfun(seq(-90, 90, length = length(bp)), bp)
axis(2, at = lons(c(-90, -60, -30, 0, 30, 60, 90)), labels = c(-90, -60, -30, 0, 30, 60, 90), las = 1)

par(new = T)
  plot(NA, xlim = c(0,1), bty = "n", ylim = range(bp), yaxt = "n", xaxt = "n", bty = "n")
  arrows(out$t.med-out$t.sd, lons(out$lat), out$t.med+out$t.sd, lons(out$lat), length = 0, code = 3, col = "darkred", lwd = 1)
  points(out$t.med, lons(out$lat), pch = 21, col = "darkred", bg = "white", type = "b", lwd = 2, cex = 2)
  points(out$tpred.med, lons(out$lat), pch = 22, col = "darkred", bg = "darkred", type = "o", lwd = 2, cex = 1, lty = 2)
axis(3, line = 1)   
par(new = T)
plot(NA, xlim = c(0,60), bty = "n", ylim = range(bp), yaxt = "n", xaxt = "n", bty = "n")
points(out$tamp.med, lons(out$lat), pch = 23, col = "darkred", bg = "darkred", type = "o", lwd = 2, cex = 1, lty = 3)
axis(3, line = 2)   

bp <- barplot(height = rbind(lmass$prec*(res(temp.seas)[1]^2)/1000000, 
                             (lmass$terr-lmass$prec)*(res(temp.seas)[1]^2)/1000000), horiz = T, border = NA, col = c("grey50", "grey80"))  
axis(2, at = lons(c(-90, -60, -30, 0, 30, 60, 90)), labels = c(-90, -60, -30, 0, 30, 60, 90), las = 1)
  
par(new = T)
plot(NA, xlim = c(0,1), bty = "n", ylim = range(bp), yaxt = "n", xaxt = "n", bty = "n")
  arrows(ifelse(out$pr.med-out$pr.sd<0,0,out$pr.med-out$pr.sd), lons(out$lat), out$pr.med+out$pr.sd, lons(out$lat), length = 0, code = 3, col = "darkblue", lwd = 1)
  points(out$pr.med, lons(out$lat), pch = 22, col = "darkblue", bg = "white", type = "b", lwd = 2, cex = 2)
  points(out$Prpred.med, lons(out$lat), pch = 22, col = "darkblue", bg = "darkblue", type = "o", lwd = 2, cex = 1, lty = 2)
axis(3) 
par(new = T)
plot(NA, xlim = c(0,15), bty = "n", ylim = range(bp), yaxt = "n", xaxt = "n", bty = "n")
points(out$Pramp.med, lons(out$lat), pch = 23, col = "darkblue", bg = "darkblue", type = "o", lwd = 2, cex = 1, lty = 3)
axis(3, line = 2)  

bp <- barplot(height = rbind(lmass$vh*(res(temp.seas)[1]^2)/1000000, 
                             (lmass$terr-lmass$vh)*(res(temp.seas)[1]^2)/1000000), horiz = T, border = NA, col = c("grey50", "grey80"))  
axis(2, at = lons(c(-90, -60, -30, 0, 30, 60, 90)), labels = c(-90, -60, -30, 0, 30, 60, 90), las = 1)

par(new = T)
plot(NA, xlim = c(0,1), bty = "n", ylim = range(bp), yaxt = "n", xaxt = "n", bty = "n")
arrows(ifelse(out$vh.med-out$vh.sd<0,0,out$vh.med-out$vh.sd), lons(out$lat), out$vh.med+out$vh.sd, lons(out$lat), length = 0, code = 3, col = "darkgreen", lwd = 1)
points(out$vh.med, lons(out$lat), pch = 24, col = "darkgreen", bg = "white", type = "b", lwd = 2, cex = 2)
points(out$VHpred.med, lons(out$lat), pch = 22, col = "darkgreen", bg = "darkgreen", type = "o", lwd = 2, cex = 1, lty = 2)
axis(3) 
par(new = T)
plot(NA, xlim = c(0, 7), bty = "n", ylim = range(bp), yaxt = "n", xaxt = "n", bty = "n")
points(out$VHamp.med, lons(out$lat), pch = 23, col = "darkgreen", bg = "darkgreen", type = "o", lwd = 2, cex = 1, lty = 3)
axis(3, line = 2)  
#
  
par(opar)
dev.off()
