

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

# altitude A conductaces 0, .1, .2, .3, 1 
x<-raster("../spatial/Elevation/rcl_Alt_A.asc")
plot(x, col=c("white", "grey90", "grey80", "grey70" , "black"), 
     breaks=c(0, .09, .19, .29, .99, 1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')

# altitude B conductaces 0, .2, .3, 1
x<-raster("../spatial/Elevation/rcl_Alt_B.asc")
plot(x, col=c("white", "grey80", "grey70" , "black"), 
     breaks=c(0, .19, .29, .99, 1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')


## Slope
# Slope A .2, .5. 1
x<-raster("../spatial/Slope/rcl_S_A.asc")
plot(x, col=c("grey80", "grey50", "black"), 
     breaks=c(.19, .49, .99,  1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')

# Slope B .2, .7, 1
x<-raster("../spatial/Slope/rcl_S_B.asc")
plot(x, col=c("grey80", "grey30", "black"), 
     breaks=c(.19, .69, .99,  1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')

# Slope C and D .2, 1
for (i in c("C", "D")){
  x<-raster(paste0("../spatial/Slope/rcl_S_", i, ".asc"))
  plot(x, col=c("gray80", "black"), legend=FALSE, xaxt='n', yaxt='n')
}

 
# Vegetation type
# A  0 .1 .7 1
x<-raster("../spatial/VegetationType/rcl_A_G_3.asc")
plot(x,  col=c("white", "grey90", "grey30", "black"), 
     breaks=c(0, .09, .69, .99, 1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')

# B 0 .1 .5 1
x<-raster("../spatial/VegetationType/rcl_B_D_3.asc")
plot(x,  col=c("white", "grey90", "grey50", "black"), 
     breaks=c(0, .09, .49, .99, 1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')

# C 0 .1 .2 1
x<-raster("../spatial/VegetationType/rcl_C_C_3.asc")
plot(x,  col=c("white", "grey90", "grey80", "black"), 
     breaks=c(0, .09, .19, .99, 1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')

# D 0 .1 .2 .9 1 
x<-raster("../spatial/VegetationType/rcl_D_A_3.asc")
plot(x,  col=c("white", "grey90", "grey80", "grey10", "black"), 
     breaks=c(0, .09, .19, .89, .99, 1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')

# E 0 .1 .2 .7 1
x<-raster("../spatial/VegetationType/rcl_E_B_3.asc")
plot(x,  col=c("white", "grey90", "grey80", "grey30", "black"), 
     breaks=c(0, .09, .19, .69, .99, 1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')

# F 0 .1 .2 .5 1
x<-raster("../spatial/VegetationType/rcl_F_E_3.asc")
plot(x,  col=c("white", "grey90", "grey80", "grey50", "black"), 
     breaks=c(0, .09, .19, .49, .99, 1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')

# G 0 .1 .2 1
x<-raster("../spatial/VegetationType/rcl_G_F_3.asc")
plot(x,  col=c("white", "grey90", "grey80", "black"), 
     breaks=c(0, .09, .19, .99, 1), # add breakpoints so colors correspond to conductances
     legend=TRUE, xaxt='n', yaxt='n')

# H 0 .1 .2 .7 1
x<-raster("../spatial/VegetationType/rcl_H_H_3.asc")
plot(x,  col=c("white", "grey90", "grey80", "grey30", "black"), 
     breaks=c(0, .09, .19, .69, .99, 1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')

# I 0 .1 .2 1
x<-raster("../spatial/VegetationType/rcl_I_I_3.asc")
plot(x,  col=c("white", "grey90", "grey80", "black"), 
     breaks=c(0, .09, .19, .99, 1), # add breakpoints so colors correspond to conductances
     legend=FALSE, xaxt='n', yaxt='n')

# J 0 .1 .5 .7 1 ### Raster is missing, complete file name and comment when raster is ready
#x<-raster("../spatial/VegetationType/rcl_XXXXXX_3.asc")
#plot(x,  col=c("white", "grey90", "grey50", "grey30", "black"), 
#     breaks=c(0, .09, .49, .69, .99, 1), # add breakpoints so colors correspond to conductances
#     legend=FALSE, xaxt='n', yaxt='n')


## Flat
x<-raster("../spatial/Elevation/NevTol_Alt_flat.asc")
plot(x, col=c("black"),  legend=FALSE, xaxt='n', yaxt='n')


##  Sampling points
# get data points
Ppoints<-read.delim("../spatial/surveyed_mountainNevadoToluca.csv", header = TRUE, sep = ",")

#plot
plot(x, col=c("grey"),  legend=FALSE, xaxt='n', yaxt='n')
points(Ppoints[,5], Ppoints[,6], col="black", cex=0.3)



