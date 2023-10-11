
## Amazeco: 
#Covering the Amazon with an Ecosystem Structure EBV product combining 
#satellite and airborne lidar

#========================================================================#
#------------------------------------------------------------------------#
# Title: Rescale ALS 25m maps to 1km                                     #    
# Author: Danilo R. A. de Almeida                                        #
# contact: daniloflorestas@gmail.com                                     #
#------------------------------------------------------------------------#

# remove all objects
rm(list = ls())

# load packages
pacman::p_load(raster, sf, stars, sp, fasterize, e1071)

# read reference 1km grid (from GEDI maps)
grid.ref = stars::read_stars("1_output/GEDI_maps_1km/_rh80_M1.tif")
# plot(grid.ref)

# paths inputs 25m layers (dtm, chm and lad)
path.in = "1_output/7_Layers_25m"

# input 25m files
files.in = list.files(path.in, pattern = ".tif")
if(length(grep("aux.|LadStack", files.in))>0){
  files.in= files.in[-c(grep("aux.|LadStack", files.in))]
}

# outputs path
path.out = "1_output/8_Layers_1km" 
dir.create(path.out, showWarnings = FALSE)

#moments functions
M2 = function(x){e1071::moment(x, order = 2, center = TRUE, na.rm = T) * length(x)}
M3 = function(x){e1071::moment(x, order = 3, center = TRUE, na.rm = T) * length(x)}
M4 = function(x){e1071::moment(x, order = 4, center = TRUE, na.rm = T) * length(x)}


####
# LOOPING 
for(i in 1:length(files.in)){
  
  #read 25m t.layer and set crs
  t.layer = raster(file.path(path.in, files.in[i]))
  
  # temp reproject 25m-layer to latlong wgs84
  t.25m = st_transform(stars::st_as_stars(t.layer), sf::st_crs(grid.ref))
  # plot(t, axes = T)
  
  
  ### resample ~1km
  
  #n
  t = aggregate(t.25m, t.ref, function(x) (sum(!is.na(x))))
  df.t = as.data.frame(grid.ref)
  df.t[,3] = as.data.frame(t)[,2]
  t = stars::st_as_stars(df.t)
  st_crs(t) = st_crs(grid.ref)
  stars::write_stars(t, 
              paste0(path.out, "/", gsub("25m.tif", "1km", files.in[i]), "_N.tif"))
  rm(t)
  
  
  #mean
  t = aggregate(t.25m, grid.ref, mean, na.rm = T)
  df.t = as.data.frame(grid.ref)
  df.t[,3] = as.data.frame(t)[,2]
  t = stars::st_as_stars(df.t)
  st_crs(t) = st_crs(grid.ref)
  stars::write_stars(t, 
                     paste0(path.out, "/", gsub("25m.tif", "1km", files.in[i]), "_MEAN.tif"))
  rm(t)
  
  #max
  t = aggregate(t.25m, grid.ref, max, na.rm = T)
  df.t = as.data.frame(grid.ref)
  df.t[,3] = as.data.frame(t)[,2]
  t = stars::st_as_stars(df.t)
  st_crs(t) = st_crs(grid.ref)
  stars::write_stars(t, 
                     paste0(path.out, "/", gsub("25m.tif", "1km", files.in[i]), "_MAX.tif"))
  rm(t)
  
  #min
  t = aggregate(t.25m, grid.ref, min, na.rm = T)
  df.t = as.data.frame(grid.ref)
  df.t[,3] = as.data.frame(t)[,2]
  t = stars::st_as_stars(df.t)
  st_crs(t) = st_crs(grid.ref)
  stars::write_stars(t, 
                     paste0(path.out, "/", gsub("25m.tif", "1km", files.in[i]), "_MIN.tif"))
  rm(t)
  
  # M2
  t = aggregate(t.25m, grid.ref, M2)
  df.t = as.data.frame(grid.ref)
  df.t[,3] = as.data.frame(t)[,2]
  t = stars::st_as_stars(df.t)
  st_crs(t) = st_crs(grid.ref)
  stars::write_stars(t, 
                     paste0(path.out, "/", gsub("25m.tif", "1km", files.in[i]), "_M2.tif"))
  rm(t)
  
  # M3
  t = aggregate(t.25m, grid.ref, M3)
  df.t = as.data.frame(grid.ref)
  df.t[,3] = as.data.frame(t)[,2]
  t = stars::st_as_stars(df.t)
  st_crs(t) = st_crs(grid.ref)
  stars::write_stars(t, 
                     paste0(path.out, "/", gsub("25m.tif", "1km", files.in[i]), "_M3.tif"))
  rm(t)
  
  # M4
  t = aggregate(t.25m, grid.ref, M4)
  df.t = as.data.frame(grid.ref)
  df.t[,3] = as.data.frame(t)[,2]
  t = stars::st_as_stars(df.t)
  st_crs(t) = st_crs(grid.ref)
  stars::write_stars(t, 
                     paste0(path.out, "/", gsub("25m.tif", "1km", files.in[i]), "_M4.tif"))
  rm(t)
  
  
  #rm temp objects
  rm(t.layer,  t.25m)
  
  #print process
  print(paste(i, "/", length(files.in)))
  
}

