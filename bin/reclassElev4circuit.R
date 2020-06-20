reclassElev4circuit<-function(inraster, minElev, outraster, NotS){
  # Function to  Reclassify an elevation  raster
  # so that non suitable elevation is set to NotS and suitable to 1
  library(raster)
  inraster=inraster #input raster
  minElev=minElev # numeric value with desired minElevation to be reclasified to 1
  NotS=NotS # value to give to cells below the min value found for the presence points
  outraster=outraster #Path to the desired output file
  
  ### do the stuff
  x<-raster(inraster)
  print("original raster")
  print(x)
  
  
  # reclasify so that values below min found tat Ppoints is 0.2 and above 1.
  m<-rbind(c(cellStats(x, stat="min"), minElev, NotS), c(minElev, cellStats(x, stat="max"), 1))
  y<-reclassify(x, rcl=m,
    filename=outraster,
    include.lowest=TRUE, overwrite=TRUE, right=FALSE) 
  print("output raster")
  print(y)
  
  par(mfrow=c(1,2))
  plot(x)
  plot(y)
  par(mfrow=c(1,1))
}