create_Flat<-function(inraster, outname){
  # this function creates a ‘flat’ landscape surface (e.g. Lee-Yaw et al. 2009), 
  # in which all grid cells had the same value
  
  # inraster = path to the raster file to convert
  # outname = filename for output
  
  inraster=inraster
  outname=outname
  
  # load raster
  library(raster)
  
  # get data
  x<-raster(inraster)
  print("original raster")
  print(x)
  
  # reclasify so that: Old values: 0-1 →  New values: 1.
  m<-c(cellStats(x, stat="min"), cellStats(x, stat="max"), 1)
  y<-reclassify(x, rcl=m,
    filename=outname,
    include.lowest=TRUE, overwrite=TRUE) 
  print("output raster")
  return(y)  
  
}

