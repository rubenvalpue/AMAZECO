

## Amazeco: 
#Covering the Amazon with an Ecosystem Structure EBV product combining 
#satellite and airborne lidar

#========================================================================#
#------------------------------------------------------------------------#
# Title: Rescale GEDI.like metrics to 1km resolution                     #    
# Author: Danilo R. A. de Almeida                                        #
# contact: daniloflorestas@gmail.com                                     #
#------------------------------------------------------------------------#

# remove all objects
rm(list = ls())

# load packages
pacman::p_load(raster, sf, stars, sp, fasterize, e1071)

# read reference 1km grid
grid.ref = stars::read_stars("1_output/GEDI_maps_1km/_rh80_M1.tif")

# path and files of gedilike data (25m res)
path.in = "1_output/9_GediLike/"
gedilike.files = list.files(path.in, pattern = ".csv", full.names = T)

# get coord references and crs from 25m layers files
path.25m = "1_output/7_Layers_25m"
files.25m = list.files(path.25m, pattern = ".tif")
if(length(grep("aux.", files.25m))>0){
  files.25m = files.25m[-c(grep("aux.", files.25m))]
}
files.25m = grep("CHM_MEAN_25m", files.25m, value = T)

# outputs path
path.out = "1_output/10_Layers_1km_Gedilike" 
dir.create(path.out, showWarnings = FALSE)

#set sites (to get crs from reference crs datum)
sites = substr(basename(gedilike.files), 1, 3)

# vector of variables
v.var = c("cover", 
          "rhGauss.50", "rhGauss.70", "rhGauss.75", "rhGauss.80", "rhGauss.85",
          "rhGauss.90", "rhGauss.95", "rhGauss.100", 
          "gLAI0t10", "gLAI10t20", "gLAI20t30", "gLAI30t40", 
          "FHD", "FHDcanGauss")

#moments functions
M2 = function(x){e1071::moment(x, order = 2, center = TRUE, na.rm = T) * length(x)}
M3 = function(x){e1071::moment(x, order = 3, center = TRUE, na.rm = T) * length(x)}
M4 = function(x){e1071::moment(x, order = 4, center = TRUE, na.rm = T) * length(x)}

####
# LOOPING will start here

for(i in 1:length(sites)){
  
  #read 25m gedi.like data from site i 
  gedilike = read.csv(gedilike.files[grep(sites[i], basename(gedilike.files))])
  t.25m = grep(sites[i], files.25m, value = T)          
  t.raster = raster::raster(file.path(path.25m, t.25m[j]))
  
  # looping to get each variable
  for(j in 1:length(v.var)){
    
    #read 25m t.layer
    t.layer = t.raster
    values(t.layer)[!is.na(values(t.raster))] = gedilike[, v.var[j]]
    t.25m = st_transform(stars::st_as_stars(t.layer), sf::st_crs(grid.ref))
    
    # Visual check
    # par(mfrow = c(1, 2))
    # plot(t.raster)
    # plot(t.layer)
     
    
    ### resample ~1km
    
    #n
    t = aggregate(t.25m, t.ref, function(x) (sum(!is.na(x))))
    df.t = as.data.frame(grid.ref)
    df.t[,3] = as.data.frame(t)[,2]
    t = stars::st_as_stars(df.t)
    st_crs(t) = st_crs(grid.ref)
    # plot(t, axes = T)
    stars::write_stars(t, 
                       paste0(path.out, "/", sites[i], "_GediLike_", v.var[j], "_N.tif"))
    rm(t)
    
    
    #mean
    t = aggregate(t.25m, grid.ref, mean, na.rm = T)
    df.t = as.data.frame(grid.ref)
    df.t[,3] = as.data.frame(t)[,2]
    t = stars::st_as_stars(df.t)
    st_crs(t) = st_crs(grid.ref)
    stars::write_stars(t, 
                       paste0(path.out, "/", sites[i], "_GediLike_", v.var[j], "_MEAN.tif"))
    rm(t)
    
    #max
    t = aggregate(t.25m, grid.ref, max, na.rm = T)
    df.t = as.data.frame(grid.ref)
    df.t[,3] = as.data.frame(t)[,2]
    t = stars::st_as_stars(df.t)
    st_crs(t) = st_crs(grid.ref)
    stars::write_stars(t, 
                       paste0(path.out, "/", sites[i], "_GediLike_", v.var[j], "_MAX.tif"))
    rm(t)
    
    #min
    t = aggregate(t.25m, grid.ref, min, na.rm = T)
    df.t = as.data.frame(grid.ref)
    df.t[,3] = as.data.frame(t)[,2]
    t = stars::st_as_stars(df.t)
    st_crs(t) = st_crs(grid.ref)
    stars::write_stars(t, 
                       paste0(path.out, "/", sites[i], "_GediLike_", v.var[j], "_MIN.tif"))
    rm(t)
    
    # M2
    t = aggregate(t.25m, grid.ref, M2)
    df.t = as.data.frame(grid.ref)
    df.t[,3] = as.data.frame(t)[,2]
    t = stars::st_as_stars(df.t)
    st_crs(t) = st_crs(grid.ref)
    stars::write_stars(t, 
                       paste0(path.out, "/", sites[i], "_GediLike_", v.var[j], "_M2.tif"))
    rm(t)
    
    # M3
    t = aggregate(t.25m, grid.ref, M3)
    df.t = as.data.frame(grid.ref)
    df.t[,3] = as.data.frame(t)[,2]
    t = stars::st_as_stars(df.t)
    st_crs(t) = st_crs(grid.ref)
    stars::write_stars(t, 
                       paste0(path.out, "/", sites[i], "_GediLike_", v.var[j], "_M3.tif"))
    rm(t)
    
    # M4
    t = aggregate(t.25m, grid.ref, M4)
    df.t = as.data.frame(grid.ref)
    df.t[,3] = as.data.frame(t)[,2]
    t = stars::st_as_stars(df.t)
    st_crs(t) = st_crs(grid.ref)
    stars::write_stars(t, 
                       paste0(path.out, "/", sites[i], "_GediLike_", v.var[j], "_M4.tif"))
    rm(t)
    
    
    #rm temp objects
    rm(t.layer,  t.25m)
    
    #print process
    print(paste0("i: ", i, "/", length(sites), " - j: ", j, "/", length(v.var)))
    
  }# end j
  
  rm(gedilike, t.raster)
  
}# end i

