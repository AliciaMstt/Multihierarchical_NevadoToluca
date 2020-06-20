reclass4circuit<-function(inraster, Ppoints, outraster, NotS){
  # Function to  Reclassify SDM 
  # so that non suitable habitat is set to NotS and suitable habitat to 1
   library(raster)
  inraster=inraster #input raster
  Ppoints=Ppoints # 2 cols dataframe with longitude and latitude data, in that order
  NotS=NotS # value to give to cells below the min value found for the presence points
  outraster=outraster #Path to the desired output file
  
  ### do the stuff
  x<-raster(inraster)
  print("original raster")
  print(x)
  
  # extract raster values at sampling points
  extract(x, Ppoints)
  minv<-min(extract(x, Ppoints), na.rm=TRUE)
  print("min Ppoint value")
  print(minv)
  
  # reclasify so that values below min found tat Ppoints is 0.1 and above 1.
  m<-rbind(c(0, minv , NotS), c(minv, 1, 1))
  y<-reclassify(x, rcl=m,
    filename=outraster,
    include.lowest=TRUE, overwrite=TRUE, right=FALSE) 
  print("output raster")
  
  par(mfrow=c(1,2))
  plot(x)
  points(Ppoints[,1], Ppoints[,2], col="black", cex=0.2)
  plot(y)
  points(Ppoints[,1], Ppoints[,2], col="black", cex=0.2)
  par(mfrow=c(1,1))
}