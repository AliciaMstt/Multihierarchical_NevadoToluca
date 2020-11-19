

library(raster)

################ plots Resistance rasters #####################

png(filename="../figures/PlotRasters.png", width=920 , height=1536, units="px") # set size of the file to plot 
par(mfrow=c(8,4)) #the number of rows and columns the figure would have


# use different grey colors depending on values for conduductaces for altitude A and B
# 0 = white
# .1 ="grey90"
# .2 = "grey80"
# .3 = "grey70"
# .5 = "grey50"
# .7 = "grey30"
# .9 = "grey10"
# 1 ="black·

######Elevation
# crude altitdes conductances .2, 1
for (i in c(2000, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700)) {
  x<-raster(paste0("../spatial/Elevation/NevTol_Alt_", i, "_reclass.asc"))
  plot(x, col=c("grey80","black"), legend=FALSE, xaxt='n', yaxt='n')
  text(label=as.character(i), x=-100.09, y=19.03, adj=0, cex=1.2)
  }

# altitude A conductaces 0, .1, .2, .3, 1 
x<-raster("../spatial/Elevation/rcl_Alt_A.asc")
plot(x, col=c("white", "grey90", "grey80", "grey70", "black"), 
     breaks=c(0, .09, .19, .29, .99, 1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')
text(label="Alt. A", x=-100.09, y=19.03, adj=0, cex=1.2)

# altitude B conductaces 0, .2, .3, 1
x<-raster("../spatial/Elevation/rcl_Alt_B.asc")
plot(x, col=c("white", "grey80", "grey70" , "black"), 
     breaks=c(0, .19, .29, .99, 1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')
text(label="Alt. B", x=-100.09, y=19.03, adj=0, cex=1.2)

# altitude C conductaces 0, .2, 1
x<-raster("../spatial/Elevation/rcl_Alt_C.asc")
plot(x, col=c("white", "grey80", "black"), 
     breaks=c(0, .19, .99, 1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')
text(label="Alt. C", x=-100.09, y=19.03, adj=0, cex=1.2)


## Slope
# Slope A .2, .5, 1
x<-raster("../spatial/Slope/rcl_S_A.asc")
plot(x, col=c("grey80", "grey50", "black"), 
     breaks=c(.19, .49, .99, 1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')
text(label="Slo. A", x=-100.09, y=19.03, adj=0, cex=1.2)

# Slope B .2, .7, 1
x<-raster("../spatial/Slope/rcl_S_B.asc")
plot(x, col=c("grey80", "grey30", "black"), 
     breaks=c(.19, .69, .99, 1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')
text(label="Slo. B", x=-100.09, y=19.03, adj=0, cex=1.2)

# Slope C and D .2, 1
for (i in c("C", "D")){
  x<-raster(paste0("../spatial/Slope/rcl_S_", i, ".asc"))
  plot(x, col=c("gray80", "black"), legend=FALSE, xaxt='n', yaxt='n')
  text(label=paste("Slo.", i), x=-100.09, y=19.03, adj=0, cex=1.2)
}

 
# Vegetation type
# A  0, .1, .7, 1
x<-raster("../spatial/VegetationType/rcl_A_G_3.asc")
plot(x,  col=c("white", "grey90", "grey30", "black"), 
     breaks=c(0, .09, .69, .99, 1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')
text(label="Veg. A", x=-100.09, y=19.03, adj=0, cex=1.2)

# B 0, .1, .5, 1
x<-raster("../spatial/VegetationType/rcl_B_D_3.asc")
plot(x,  col=c("white", "grey90", "grey50", "black"), 
     breaks=c(0, .09, .49, .99, 1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')
text(label="Veg. B", x=-100.09, y=19.03, adj=0, cex=1.2)

# C 0, .1, .2, 1
x<-raster("../spatial/VegetationType/rcl_C_C_3.asc")
plot(x,  col=c("white", "grey90", "grey80", "black"), 
     breaks=c(0, .09, .19, .99, 1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')
text(label="Veg. C", x=-100.09, y=19.03, adj=0, cex=1.2)

# D 0, .1, .2, .9, 1 
x<-raster("../spatial/VegetationType/rcl_D_A_3.asc")
plot(x,  col=c("white", "grey90", "grey80", "grey10", "black"), 
     breaks=c(0, .09, .19, .89, .99, 1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')
text(label="Veg. D", x=-100.09, y=19.03, adj=0, cex=1.2)

# E 0, .1, .2, .7, 1
x<-raster("../spatial/VegetationType/rcl_E_B_3.asc")
plot(x,  col=c("white", "grey90", "grey80", "grey30", "black"), 
     breaks=c(0, .09, .19, .69, .99, 1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')
text(label="Veg. E", x=-100.09, y=19.03, adj=0, cex=1.2)

# F 0, .1, .2, .5, 1
x<-raster("../spatial/VegetationType/rcl_F_E_3.asc")
plot(x,  col=c("white", "grey90", "grey80", "grey50", "black"), 
     breaks=c(0, .09, .19, .49, .99, 1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')
text(label="Veg. F", x=-100.09, y=19.03, adj=0, cex=1.2)

# G 0, .1, .2, 1
x<-raster("../spatial/VegetationType/rcl_G_F_3.asc")
plot(x,  col=c("white", "grey90", "grey80", "black"), 
     breaks=c(0, .09, .19, .99, 1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')
text(label="Veg. G", x=-100.09, y=19.03, adj=0, cex=1.2)

# H 0, .1, .2, .7, 1
x<-raster("../spatial/VegetationType/rcl_H_H_3.asc")
plot(x,  col=c("white", "grey90", "grey80", "grey30", "black"), 
     breaks=c(0, .09, .19, .69, .99, 1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')
text(label="Veg. H", x=-100.09, y=19.03, adj=0, cex=1.2)

# I 0, .1, .2, 1
x<-raster("../spatial/VegetationType/rcl_I_I_3.asc")
plot(x, col=c("white", "grey90", "grey80", "black"), 
     breaks=c(0, .09, .19, .99, 1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')
text(label="Veg. I", x=-100.09, y=19.03, adj=0, cex=1.2)

# J 0, .1, .5, .7, 1 
x<-raster("../spatial/VegetationType/rcl_I_I_3.asc")
plot(x,  col=c("white", "grey90", "grey50", "grey30", "black"), 
     breaks=c(0, .09, .49, .69, .99, 1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')
 text(label="Veg. J", x=-100.09, y=19.03, adj=0, cex=1.2)

## Flat
x<-raster("../spatial/Elevation/NevTol_Alt_flat.asc")
plot(x, col=c("black"),  legend=FALSE, xaxt='n', yaxt='n')
text(label="Flat", x=-100.09, y=19.03, adj=0, cex=1.2)

##  Sampling points
# get data points
Ppoints<-read.delim("../spatial/surveyed_mountainNevadoToluca.csv", header = TRUE, sep = ",")

#plot
plot(x, col=c("grey"),  legend=FALSE, xaxt='n', yaxt='n')
points(Ppoints[,5], Ppoints[,6], col="darkgreen", cex=0.3)
text(label="Samp.", x=-100.09, y=19.06, adj=0, cex=1.2)
text(label="points", x=-100.09, y=19.03, adj=0, cex=1.2)

## pseudo legend
x<- c(rep(1,4), rep(1.7,4))
y<- c(seq(.1,4), seq(.1,4))
plot(x,y, pch=22, cex=3,
     bg= c("grey70","grey80", "grey90", "white", "black", "grey10", "grey30", "grey50"),
     xaxt='n', yaxt='n', xlab='', ylab='', bty='n', # delete axis and box stuff
     xlim=c(-.1,3), ylim=c(-2,6)) # this makes points to be closer to each other
text(x=x+.3, y=y, cex=1.2,
     labels=c("0.3", "0.2", "0.1", "0", "1", "0.9", "0.7", "0.5"))
text(x=.8, y=4.5, labels="Conductances", adj=0, cex=1.5)

## close plotting device
dev.off() #needed because we used png in line 7

