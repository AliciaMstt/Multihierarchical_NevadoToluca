## script to reclassify altitude output rasters to desired values and to create a flat lanscape
# plot rasters afterwards

rm(list = ls())
library(sp)
library(raster)
library(rgdal)


##### Create "flat" landscape #####
source("../bin/create_Flat.R")

create_Flat(inraster="../spatial/Elevation/NevTol_Alt.tif",
            outname="../spatial/Elevation/NevTol_Alt_flat.asc")

##### Reclassify models  #####
# so that non suitable habitat is set to NotS and suitable habitat to 1
source("../bin/reclass4circuit.R")
datafolder="../spatial/SDM/out"


##### Reclassify elevation models  #####
# so that not suitable elevation is set to NotS and suitable to 1

# path to elevation raw data
datafolder="../spatial/Elevation"

# define in data and function
inraster=paste0(datafolder, "/NevTol_Alt.tif")
source("../bin/reclassElev4circuit.R")

## Do reclassification for several elevations

#for (i in c(2000, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700)) {
#    print(paste("for min elevation:", i))
#    outraster=paste0(datafolder, "/NevTol_Alt_", i, "_reclass.asc")
#    # reclassify
#    reclassElev4circuit(inraster, minElev=i, outraster, NotS=0.2)
#    }


##### Raster Elevation Hipotesis "A"######### 
myraster<-raster("../spatial/Elevation/NevTol_Alt.tif")
myraster
plot(myraster)

# read reclasification matrix (make sure no header)
# reclassify the values into three groups 
rcl_Alt_A_reclass <- read.table("../spatial/Elevation/RasterAltitude_A.txt", sep = ",", header=F)
dim(rcl_Alt_A_reclass)
class(rcl_Alt_A_reclass)
rcl_Alt_A_reclass<- as.matrix(rcl_Alt_A_reclass) 
class(rcl_Alt_A_reclass)
rcl_Alt_A_reclass

# reclasify input raster 
xa<-reclassify(myraster,  rcl=rcl_Alt_A_reclass, include.lowest=FALSE, right=FALSE)
xa
plot(xa)
plot(xa, col=c("grey", "black"), legend=FALSE, xaxt='n', yaxt='n')

writeRaster(xa, filename="../spatial/Elevation/rcl_Alt_A.tif", format="GTiff", overwrite=TRUE)
writeRaster(xa, filename="../spatial/Elevation/rcl_Alt_A.asc", format="ascii", overwrite=TRUE)


###### Raster Altitude Hipotesis "B"############ 
# reclassify the values into three groups 

rcl_Alt_B_reclass <- read.table("../spatial/Elevation/RasterAltitude_B.txt",sep = ",", header=F)
dim(rcl_Alt_B_reclass)
class(rcl_Alt_B_reclass)
rcl_Alt_B_reclass<- as.matrix(rcl_Alt_B_reclass) 
class(rcl_Alt_B_reclass)
rcl_Alt_B_reclass

# reclasify input raster 
xb<-reclassify(myraster,  rcl=rcl_Alt_B_reclass, include.lowest=FALSE, right=FALSE)
xb
plot(xb)
plot(xb, col=c("grey", "black"), legend=FALSE, xaxt='n', yaxt='n')


writeRaster(xb, filename="../spatial/Elevation/rcl_Alt_B.tif", format="GTiff", overwrite=TRUE)
writeRaster(xb, filename="../spatial/Elevation/rcl_Alt_B.asc", format="ascii", overwrite=TRUE)


##################SLOPE#############################


#Raster Slope Nevado de Toluca
myrasterS<-raster("../spatial/Slope/NevTol_Pen.tif")
myrasterS
plot(myrasterS)

#######Raster Slope hipotesis "A"#####################
rcl_S_A <- read.table("../spatial/Slope/RasterValue_SlopeA.txt",sep = ",", header=F)
dim(rcl_S_A)
class(rcl_S_A)
rcl_S_A<- as.matrix(rcl_S_A) 
class(rcl_S_A)
rcl_S_A

# reclasify input raster Slope NT
x_SA<-reclassify(myrasterS,  rcl=rcl_S_A, include.lowest=FALSE, right=FALSE)
x_SA
plot(x_SA)
writeRaster(x_SA, filename="../spatial/Slope/rcl_S_A.tif", format="GTiff", overwrite=TRUE)
writeRaster(x_SA, filename="../spatial/Slope/rcl_S_A.asc", format="ascii", overwrite=TRUE)


#############Raster Slope Hipotesis "B"###################
rcl_S_B <- read.table("../spatial/Slope/RasterValue_SlopeB.txt",sep = ",", header=F)
dim(rcl_S_B)
class(rcl_S_B)
rcl_S_A<- as.matrix(rcl_S_B) 
class(rcl_S_B)
rcl_S_B

# reclasify input raster Slope NT
x_SB<-reclassify(myrasterS,  rcl=rcl_S_B, include.lowest=FALSE, right=FALSE)
x_SB
plot(x_SB)
writeRaster(x_SB, filename="../spatial/Slope/rcl_S_B.tif", format="GTiff", overwrite=TRUE)
writeRaster(x_SB, filename="../spatial/Slope/rcl_S_B.asc", format="ascii", overwrite=TRUE)

#############Raster Slope Hipotesis "C"######################
rcl_S_C <- read.table("../spatial/Slope/RasterValue_SlopeC.txt",sep = ",", header=F)
dim(rcl_S_C)
class(rcl_S_C)
rcl_S_C<- as.matrix(rcl_S_C) 
class(rcl_S_C)
rcl_S_C

# reclasify input raster Slope NT
x_SC<-reclassify(myrasterS,  rcl=rcl_S_C, include.lowest=FALSE, right=FALSE)
x_SC
plot(x_SC)
writeRaster(x_SC, filename="../spatial/Slope/rcl_S_C.tif", format="GTiff", overwrite=TRUE)
writeRaster(x_SC, filename="../spatial/Slope/rcl_S_C.asc", format="ascii", overwrite=TRUE)

##############3Raster Slope Hipotesis "D"#######################
rcl_S_D <- read.table("../spatial/Slope/RasterValue_SlopeD.txt",sep = ",", header=F)
dim(rcl_S_D)
class(rcl_S_D)
rcl_S_C<- as.matrix(rcl_S_D) 
class(rcl_S_D)
rcl_S_D

# reclasify input raster Slope NT
x_SD<-reclassify(myrasterS,  rcl=rcl_S_D, include.lowest=FALSE, right=FALSE)
x_SD
plot(x_SD)
writeRaster(x_SD, filename="../spatial/Slope/rcl_S_D.tif", format="GTiff", overwrite=TRUE)
writeRaster(x_SD, filename="../spatial/Slope/rcl_S_D.asc", format="ascii", overwrite=TRUE)

######################################################################


###################Vegetation Type###############################


#Raster Vegetation Type in Nevado de Toluca
myrasterVT<-raster("../spatial/VegetationType/nevado_f.tif")
myrasterVT
plot(myrasterVT)

###########Vegetation Hipotesis "A"##########
#read reclasification matrix
rclA <- read.table("../spatial/VegetationType/RasterValue_A_G.txt",sep = ",", header=F)
dim(rclA)
class(rclA)
rclA<- as.matrix(rclA) 
class(rclA)
rclA

# reclasify input raster 
xVA<-reclassify(myrasterVT,  rcl=rclA, include.lowest=FALSE, right=FALSE)
xVA
plot(xVA)
writeRaster(xVA, filename="../spatial/VegetationType/rcl_A_G_3.tif", format="GTiff", overwrite=TRUE)
writeRaster(xVA, filename="../spatial/VegetationType/rcl_A_G_3.asc", format="ascii", overwrite=TRUE)


###########Vegetation Hipotesis "B"##########
#read reclasification matrix
rclB <- read.table("../spatial/VegetationType/RasterValue_B_D.txt",sep = ",", header=F)
dim(rclB)
class(rclB)
rclB<- as.matrix(rclB) 
class(rclB)
rclB

# reclasify input raster 
xVB<-reclassify(myrasterVT,  rcl=rclB, include.lowest=FALSE, right=FALSE)
xVB
plot(xVB)
writeRaster(xVB, filename="../spatial/VegetationType/rcl_B_D_3.tif", format="GTiff", overwrite=TRUE)
writeRaster(xVB, filename="../spatial/VegetationType/rcl_B_D_3.asc", format="ascii", overwrite=TRUE)

###########Vegetation Hipotesis "C"##########
#read reclasification matrix

rclC <- read.table("../spatial/VegetationType/RasterValue_C_C.txt",sep = ",", header=F)
dim(rclC)
class(rclC)
rclC<- as.matrix(rclC) 
class(rclC)
rclC

# reclasify input raster 
xVC<-reclassify(myrasterVT,  rcl=rclC, include.lowest=FALSE, right=FALSE)
xVC
plot(xVC)
writeRaster(xVC, filename="../spatial/VegetationType/rcl_C_C_3.tif", format="GTiff", overwrite=TRUE)
writeRaster(xVC, filename="../spatial/VegetationType/rcl_C_C_3.asc", format="ascii", overwrite=TRUE)


###########Vegetation Hipotesis "D"##########
#read reclasification matrix

rclD <- read.table("../spatial/VegetationType/RasterValue_D_A.txt",sep = ",", header=F) 
dim(rclD)
class(rclD)
rclD<- as.matrix(rclD) 
class(rclD)
rclD

# reclasify input raster 
xVD<-reclassify(myrasterVT,  rcl=rclD, include.lowest=FALSE, right=FALSE)
xVD
plot(xVD)
writeRaster(xVD, filename="../spatial/VegetationType/rcl_D_A_3.tif", format="GTiff", overwrite=TRUE)
writeRaster(xVD, filename="../spatial/VegetationType/rcl_D_A_3.asc", format="ascii", overwrite=TRUE)


###########Vegetation Hipotesis "E"##########
#read reclasification matrix

rclE <- read.table("../spatial/VegetationType/RasterValue_E_B.txt",sep = ",", header=F)
dim(rclE)
class(rclE)
rclE<- as.matrix(rclE) 
class(rclE)
rclE

# reclasify input raster 
xVE<-reclassify(myrasterVT,  rcl=rclE, include.lowest=FALSE, right=FALSE)
xVE
plot(xVE)
writeRaster(xVE, filename="../spatial/VegetationType/rcl_E_B_3.tif", format="GTiff", overwrite=TRUE)
writeRaster(xVE, filename="../spatial/VegetationType/rcl_E_B_3.asc", format="ascii", overwrite=TRUE)


###########Vegetation Hipotesis "F"##########
#read reclasification matrix
rclF <- read.table("../spatial/VegetationType/RasterValue_F_E.txt",sep = ",", header=F)
dim(rclF)
class(rclF)
rclF<- as.matrix(rclF) 
class(rclF)
rclF

# reclasify input raster 
xVF<-reclassify(myrasterVT,  rcl=rclF, include.lowest=FALSE, right=FALSE)
xVF
plot(xVF)
writeRaster(xVF, filename="../spatial/VegetationType/rcl_F_E_3.tif", format="GTiff", overwrite=TRUE)
writeRaster(xVF, filename="../spatial/VegetationType/rcl_F_E_3.asc", format="ascii", overwrite=TRUE)

###########Vegetation Hipotesis "G"##########
#read reclasification matrix

rclG <- read.table("../spatial/VegetationType/RasterValue_G_F.txt",sep = ",", header=F)
dim(rclG)
class(rclG)
rclG<- as.matrix(rclG) 
class(rclG)
rclG

# reclasify input raster 
xVG<-reclassify(myrasterVT,  rcl=rclG, include.lowest=FALSE, right=FALSE)
xVG
plot(xVG)
writeRaster(xVG, filename="../spatial/VegetationType/rcl_G_F_3.tif", format="GTiff", overwrite=TRUE)
writeRaster(xVG, filename="../spatial/VegetationType/rcl_G_F_3.asc", format="ascii", overwrite=TRUE)


###########Vegetation Hipotesis "H"##########
#read reclasification matrix

rclH <- read.table("../spatial/VegetationType/RasterValue_H_H.txt",sep = ",", header=F)
dim(rclH)
class(rclH)
rclH<- as.matrix(rclH) 
class(rclH)
rclH

# reclasify input raster 
xVH<-reclassify(myrasterVT,  rcl=rclH, include.lowest=FALSE, right=FALSE)
xVH
plot(xVH)
writeRaster(xVH, filename="../spatial/VegetationType/rcl_H_H_3.tif", format="GTiff", overwrite=TRUE)
writeRaster(xVH, filename="../spatial/VegetationType/rcl_H_H_3.asc", format="ascii", overwrite=TRUE)


###########Vegetation Hipotesis "I"##########
#read reclasification matrix

rclI <- read.table("../spatial/VegetationType/RasterValue_I_I.txt",sep = ",", header=F)
dim(rclI)
class(rclI)
rclI<- as.matrix(rclI) 
class(rclI)
rclI

# reclasify input raster 
xVI<-reclassify(myrasterVT,  rcl=rclI, include.lowest=FALSE, right=FALSE)
xVI
plot(xVI)
writeRaster(xVI, filename="../spatial/VegetationType/rcl_I_I_3.tif", format="GTiff", overwrite=TRUE)
writeRaster(xVI, filename="../spatial/VegetationType/rcl_I_I_3.asc", format="ascii", overwrite=TRUE)




########### plots Resistance rasters #########

par(mfrow=c(5,3))

##Elevation
# path to elevation raw data
datafolder="../spatial/Elevation"
for (i in c(2000, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700)) {
  x<-raster(paste0(datafolder, "/NevTol_Alt_", i, "_reclass.asc"))
  plot(x, col=c("grey", "black"), legend=FALSE, xaxt='n', yaxt='n')
  text(label=as.character(i), x=-100, y=20)}


datafolder="../spatial/Elevation"

plot(xa, col=c("grey", "black"), legend=FALSE, xaxt='n', yaxt='n')
plot(xb, col=c("grey", "black"), legend=FALSE, xaxt='n', yaxt='n')


datafolder="../spatial/Slope"

plot(x_SA, col=c("grey", "black"), legend=FALSE, xaxt='n', yaxt='n')
plot(x_SB, col=c("grey", "black"), legend=FALSE, xaxt='n', yaxt='n')
plot(x_SC, col=c("grey", "black"), legend=FALSE, xaxt='n', yaxt='n')
plot(x_SD, col=c("grey", "black"), legend=FALSE, xaxt='n', yaxt='n')


datafolder="../spatial/VegetationType"

plot(xVA, col=c("grey", "black"), legend=FALSE, xaxt='n', yaxt='n')
plot(xVB, col=c("grey", "black"), legend=FALSE, xaxt='n', yaxt='n')
plot(xVC, col=c("grey", "black"), legend=FALSE, xaxt='n', yaxt='n')
plot(xVD, col=c("grey", "black"), legend=FALSE, xaxt='n', yaxt='n')
plot(xVE, col=c("grey", "black"), legend=FALSE, xaxt='n', yaxt='n')
plot(xVF, col=c("grey", "black"), legend=FALSE, xaxt='n', yaxt='n')
plot(xVG, col=c("grey", "black"), legend=FALSE, xaxt='n', yaxt='n')
plot(xVH, col=c("grey", "black"), legend=FALSE, xaxt='n', yaxt='n')
plot(xVI, col=c("grey", "black"), legend=FALSE, xaxt='n', yaxt='n')

datafolder="../spatial/Elevation"
i="flat"
x<-raster(paste0(datafolder, "/NevTol_Alt_", i, ".asc"))
plot(x, col=c("grey", "black"),  legend=FALSE, xaxt='n', yaxt='n')
text(label=as.character(i), x=-100, y=20)


### data points
Ppoints<-read.delim("../spatial/surveyed_mountainNevadoToluca.csv", header = TRUE, sep = ",")
FP <- c("TLC", "SJH", "ASB", "AAB")
FP <- Ppoints$Key %in% FP 
Ppoint<-Ppoints[FP,]


# Background
i="flat"
x<-raster(paste0(datafolder, "/NevTol_Alt_", i, ".asc"))
plot(x, col=c("grey"),  legend=FALSE, xaxt='n', yaxt='n')
text(label=as.character("S. points"), x=-100, y=20)
points(Ppoints[,5], Ppoints[,6], col="black", cex=0.2)


