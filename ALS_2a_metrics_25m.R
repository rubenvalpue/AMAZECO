

## Amazeco: 
#Covering the Amazon with an Ecosystem Structure EBV product combining 
#satellite and airborne lidar

#========================================================================#
#------------------------------------------------------------------------#
# Title: ALS-derived metrics with resolution of GEDI (25m)               #    
# Author: Danilo R. A. de Almeida                                        #
# contact: daniloflorestas@gmail.com                                     #
#                                                                        #
# Products:                                                              #
#           - DTM - mean and roughness                                   #
#           - CHM - mean, max, roughness, and cover                      #
#           - LAI, LAI.under, LAHV, and LAWV
#                                                                        #
#------------------------------------------------------------------------#

#========================================================================#



# remove all objects
rm(list = ls())

# load packages
pacman::p_load(raster, sf, stars)

# paths inputs (dtm, chm and lad)
path.dtm = "1_output/1_DTM"
path.chm = "1_output/2_CHM"
path.lad= "1_output/4_LAD"

##files
files.dtm = list.files(path.dtm, pattern = ".tif")
if(length(-c(grep("aux.", files.dtm)))>1){
  files.dtm = files.dtm[-c(grep("aux.", files.dtm))]
}

files.chm = list.files(path.chm, pattern = ".tif")
if(length(-c(grep("aux.", files.chm)))>1){
  files.chm = files.chm[-c(grep("aux.", files.chm))]
}

files.lad = list.files(path.lad, pattern = ".tif")
files.lad = files.lad[-c(grep("aux.|LadStack", files.lad))]

# ouputs path
path.out = "1_output/7_Layers_25m" 
dir.create(path.out, showWarnings = FALSE)

## get sites to run in the looping
sites = unique(unlist(lapply(strsplit(files.dtm, split = "_"), function(x) (x[2]))))

# function cover (5 and 10 m threshold)
cover5 = function(x, na.rm = T){
  return(sum(as.numeric(na.omit(x))>=5)/length(as.numeric(na.omit(x))))
}

cover10 = function(x, na.rm = T){
  return(sum(as.numeric(na.omit(x))>=10)/length(as.numeric(na.omit(x))))
}

# LOOPING
for(i in 1:length(sites)){
  
  
  #--------------------------------------------------------------------------#
  ## DTM
  
  #get site i files
  t.files = grep(sites[i], files.dtm, value = T)
  
  # load files
  t.layer = lapply(as.list(t.files), function(x) (raster(file.path(path.dtm, x))))
  t.layer = do.call(merge, t.layer) #merge all
  
  # Visual Check
  plot(t.layer)
  
  ## DTM
  # mean
  t.mean = aggregate(t.layer, fact=25, fun=mean, expand=F, na.rm = T)
  #plot(t.mean)
  
  # sd
  t.sd = aggregate(t.layer, fact=25, fun=sd, expand=F, na.rm = T)
  #plot(t.sd)
  
  #roughness
  t.roug = t.sd/t.mean
  # plot(t.roug)
  
  ## write
  raster::writeRaster(t.mean, paste0(path.out, "/", sites[i], "_DTM_MEAN25m.tif"))
  raster::writeRaster(t.roug, paste0(path.out, "/", sites[i], "_DTM_ROUG25m.tif"))
  
  rm(t.files, t.layer, t.mean, t.sd, t.roug)
  
  #print progress
  print(paste(i, "/", length(sites), "- DTM done"))
  
  #END DTM
  #--------------------------------------------------------------------------#
  
  
  
  #--------------------------------------------------------------------------#
  ## CHM
  
  #get site i files
  t.files = grep(sites[i], files.chm, value = T)
  
  # load files
  t.layer = lapply(as.list(t.files), function(x) (raster(file.path(path.chm, x))))
  t.layer = do.call(merge, t.layer) #merge all
  
  #Visual check
  # plot(t.layer)
  
  # mean
  t.mean = aggregate(t.layer, fact=25, fun=mean, expand=F, na.rm = T)
  #plot(t.mean)
  
  # max
  t.max = aggregate(t.layer, fact=25, fun=max, expand=F, na.rm = T)
  #plot(t.max)
  
  # sd
  t.sd = aggregate(t.layer, fact=25, fun=sd, expand=F, na.rm = T)
  #plot(t.sd)
  
  #roughness
  t.roug = t.sd/t.mean
  # plot(t.roug)
  
  #cover
  t.cover5 = aggregate(t.layer, fact=25, fun=cover5, expand=F, na.rm = T)
  t.cover10 = aggregate(t.layer, fact=25, fun=cover10, expand=F, na.rm = T)
  
  # plot(t.cover10)
  
  ## write
  raster::writeRaster(t.mean, paste0(path.out, "/", sites[i], "_CHM_MEAN_25m.tif"))
  raster::writeRaster(t.max, paste0(path.out, "/", sites[i], "_CHM_MAX_25m.tif"))
  raster::writeRaster(t.roug, paste0(path.out, "/", sites[i], "_CHM_ROUG_25m.tif"))
  raster::writeRaster(t.cover5, paste0(path.out, "/", sites[i], "_CHM_COVER5_25m.tif"))
  raster::writeRaster(t.cover10, paste0(path.out, "/", sites[i], "_CHM_COVER10_25m.tif"))
  
  #grid ref for lad layers (using resample)
  grid.ref = t.mean
  rm(t.files, t.layer, t.mean, t.max ,t.sd, t.roug, t.cover5, t.cover10)
  
  #print progress
  print(paste(i, "/", length(sites), "- CHM done"))
  
  #END CHM
  #--------------------------------------------------------------------------#
  
  
  #--------------------------------------------------------------------------#
  ## LAI
  
  #get site i files
  t.files = grep(sites[i], files.lad, value = T)
  t.files = grep("lai.tif", t.files, value = T)
  
  # load files
  t.layer = lapply(as.list(t.files), function(x) (raster(file.path(path.lad, x))))
  
  #resample 25m
  t.layer.25m = list()
  for(j in 1:length(t.layer)){
    t.layer.25m[[j]] = raster::resample(t.layer[[j]], grid.ref, method="bilinear")
  }
  
  if(length(t.layer.25m)>1){
    t.mean = do.call(merge, t.layer.25m)
  }else{
    t.mean = t.layer.25m[[1]]
  } #merge all
  
  # plot(t.mean)

  #object to calculate lwhv
  t.lai = t.mean
  
  ## write
  raster::writeRaster(t.mean, paste0(path.out, "/", sites[i], "_LAI_MEAN_25m.tif"))
  rm(t.files, t.layer, t.layer.25m, t.mean)
  
  #print progress
  print(paste(i, "/", length(sites), "- LAI done"))
  
  #END LAI
  #--------------------------------------------------------------------------#
  
  
  #--------------------------------------------------------------------------#
  ## LAI.UNDER
  
  #get site i files
  t.files = grep(sites[i], files.lad, value = T)
  t.files = grep("understory", t.files, value = T)
  
  # load files
  t.layer = lapply(as.list(t.files), function(x) (raster(file.path(path.lad, x))))
  
  #resample 25m
  t.layer.25m = list()
  
  for(j in 1:length(t.layer)){
    t.layer.25m[[j]] = raster::resample(t.layer[[j]], grid.ref, method="bilinear")
  }
  
  # t.mean
  if(length(t.layer.25m)>1){
    t.mean = do.call(merge, t.layer.25m)
  }else{
    t.mean = t.layer.25m[[1]]
  } #merge all
  
  # plot(t.mean)
  
  
  ## write
  raster::writeRaster(t.mean, paste0(path.out, "/", sites[i], "_LAIunder_MEAN_25m.tif"))
  
  rm(t.files, t.layer, t.layer.25m, t.mean)
  
  #print progress
  print(paste(i, "/", length(sites), "- LAIunder done"))
  
  #END LAI.under
  #--------------------------------------------------------------------------#
  
  
  
  #--------------------------------------------------------------------------#
  ## LAHV
  
  #get site i files
  t.files = grep(sites[i], files.lad, value = T)
  t.files = grep("lahv", t.files, value = T)
  
  # load files
  t.layer = lapply(as.list(t.files), function(x) (raster(file.path(path.lad, x))))
  
  #resample 25m
  t.layer.25m = list()
  for(j in 1:length(t.layer)){
    t.layer.25m[[j]] = raster::resample(t.layer[[j]], grid.ref, method="bilinear")
  }
  
  # t.mean
  if(length(t.layer.25m)>1){
    t.mean = do.call(merge, t.layer.25m)
  }else{
    t.mean = t.layer.25m[[1]]
  } #merge all
  
  # plot(t.mean)
  
  # object to calculate lawv
  t.lahv = t.mean
  
  ## write
  raster::writeRaster(t.mean, paste0(path.out, "/", sites[i], "_LAHV_MEAN_25m.tif"))
  
  rm(t.files, t.layer, t.layer.25m, t.mean)
  
  #print progress
  print(paste(i, "/", length(sites), "- LAHV done"))
  
  #END LAHV
  #--------------------------------------------------------------------------#
  
  
  
  #--------------------------------------------------------------------------#
  ## LAWV
  t.lawv = t.lahv/t.lai
  
  ## write
  raster::writeRaster(t.lawv, paste0(path.out, "/", sites[i], "_LAWV_MEAN_25m.tif"))
  
  rm(t.lai, t.lahv, lawv)
  
  #print progress
  print(paste(i, "/", length(sites), "- LAWV done"))
  
  #END LAWV
  #--------------------------------------------------------------------------#
  
  print(paste(i, "/", length(sites), "- ALL DONE!"))
  
  rm(grid.ref)
  
};rm(i)


#### END


