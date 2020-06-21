

library(raster)

########### plots Resistance rasters #########

par(mfrow=c(10,3)) #the number of rows and columns the figure would have

# use different grey colors depending on values for conduductaces for altitude A and B
# 0 = white
# .1 ="grey90"
# .2 = "grey80"
# .3 = "grey70"
# .5 = "grey50"
# .7 = "grey30"
# .9 = "grey10"

## Elevation
# crude altitdes conductances .2, 1
for (i in c(2000, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700)) {
  x<-raster(paste0("../spatial/Elevation/NevTol_Alt_", i, "_reclass.asc"))
  plot(x, col=c("grey80","black"), legend=FALSE, xaxt='n', yaxt='n')
  text(label=as.character(i), x=-100, y=19.01)
  }

# altitude A conductaces 0, .1, .2, .3, 1 ### ALICIA PONER BREAKS COMO AQUÍ REVISAR LOS DE ABAJO
x<-raster("../spatial/Elevation/rcl_Alt_A.asc")
plot(x, col=c("white", "grey90", "grey80", "grey70" , "black"), 
     breaks=c(0, .1, .2, .3, .99, 1), # add breakpoints so colors correspond to conductances
     legend=TRUE, xaxt='n', yaxt='n')

# altitude B conductaces 0, .2, .3, 1
x<-raster("../spatial/Elevation/rcl_Alt_B.asc")
plot(x, col=c("white", "grey80", "grey70" , "black"), 
     breaks=c(0, .2, .3, .31, 1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')


## Slope
# Slope A .2, .5. 1
x<-raster("../spatial/Slope/rcl_S_A.asc")
plot(x, col=c("grey80", "grey50", "black"), 
     breaks=c(.2, .5, .51,  1), # add breakpoints so colors correspond to conductances
     legend=TRUE, xaxt='n', yaxt='n')

# Slope B .2, .7, 1
x<-raster("../spatial/Slope/rcl_S_B.asc")
plot(x, col=c("grey80", "grey30", "black"), 
     breaks=c(.2, .7, .71,  1), # add breakpoints so colors correspond to conductances
     legend=TRUE, xaxt='n', yaxt='n')

# Slope C and D
for (i in c("C", "D")){
  x<-raster(paste0("../spatial/Slope/rcl_S_", i, ".asc"))
  plot(x, col=c("gray80", "black"), legend=FALSE, xaxt='n', yaxt='n')
}

 
# Vegetation type
# A  0 .1 .7 1
x<-raster("../spatial/VegetationType/rcl_A_G_3.asc")
plot(x,  col=c("white", "grey90", "grey30", "black"), 
     breaks=c(0, .09, .1, .7, 1), # add breakpoints so colors correspond to conductances
     legend=TRUE, xaxt='n', yaxt='n')

# B 0 .1 .7 1
x<-raster("../spatial/VegetationType/rcl_B_D_3.asc")
plot(x,  col=c(), 
     breaks=c(0), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')

#
x<-raster("../spatial/VegetationType/rcl_A_3.asc")
plot(x,  col=c(), 
     breaks=c(), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')

#
x<-raster("../spatial/VegetationType/rcl_A_3.asc")
plot(x,  col=c(), 
     breaks=c(), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')

#
x<-raster("../spatial/VegetationType/rcl_A_3.asc")
plot(x,  col=c(), 
     breaks=c(), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')

#
x<-raster("../spatial/VegetationType/rcl_A_3.asc")
plot(x,  col=c(), 
     breaks=c(), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')

#
x<-raster("../spatial/VegetationType/rcl_A_3.asc")
plot(x,  col=c(), 
     breaks=c(), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')

#
x<-raster("../spatial/VegetationType/rcl_A_3.asc")
plot(x,  col=c(), 
     breaks=c(), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')

#
x<-raster("../spatial/VegetationType/rcl_A_3.asc")
plot(x,  col=c(), 
     breaks=c(), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')

for (i in c("A_G", "B_D", "C_C", "D_A", "E_B", "F_E", "G_F", "H_H", "I_I")){
  x<-raster(paste0("../spatial/VegetationType/rcl_", i, "_3.asc"))
  plot(x, col=gray.colors(2, start=1, end=.2), legend=FALSE, xaxt='n', yaxt='n')
}


#  Sampling points and Background
i="flat"
x<-raster(paste0("../spatial/Elevation/NevTol_Alt_", i, ".asc"))
plot(x, col=c("grey"),  legend=FALSE, xaxt='n', yaxt='n')
text(label=as.character("S. points"), x=-100, y=20)
points(Ppoints[,5], Ppoints[,6], col="black", cex=0.2)
