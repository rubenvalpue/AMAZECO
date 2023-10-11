


## Amazeco: 
#Covering the Amazon with an Ecosystem Structure EBV product combining 
#satellite and airborne lidar

#========================================================================#
#------------------------------------------------------------------------#
# Title: Creating shape points from 1km ALS and GEDI variables           #    
# Author: Danilo R. A. de Almeida                                        #
# contact: daniloflorestas@gmail.com                                     #
#------------------------------------------------------------------------#

################################################################################
################################################################################
################################################################################

#####
# 1 - create points from all 1km ALS rasters 

#remove all objects
rm(list = ls())

pacman::p_load(raster, sf, stars)

# path of 1km maps
path.1km = "1_output/8_Layers_1km"
files.1km = list.files(path.1km, pattern = "tif")

#get sites
site = unique(substr(files.1km, 1, 3))
if(length(grep("aux.", files.1km))>0){
  files.1km = files.1km[-c(grep("aux.", files.1km))]
}

#output path
path.out = "1_output/8_Layers_1km/ALL_Points/"
dir.create(path.out, showWarnings = FALSE)

## Looping i do shape points for each site
# points use coordinates of GEDI 1km grid reference

T.LIST.ALL = c() # List to fill with data.frames per site

for(i in 1:length(site)){

  # set site o
  t.files = grep(site[i], files.1km, value = T)

  
  #looping to get all raster values and coordinates
  T.LIST = list()
  for(j in 1:length(t.files)){

    # get 1km ref crop and coordinates

    #load raster j
    t.raster = stars::read_stars(file.path(path.1km, t.files[j]))
    
    
    # get points from coordinates
    if(j == 1){

      t.ref = st_as_sf(x = st_coordinates(t.raster),
                       coords = c("x", "y"),
                       crs = sf::st_crs(t.raster))

      #plot(t.ref)
      
      t.coords = st_coordinates(t.ref)

    }# end if

    # get var name, extract values and fill final list
    t.var.name = gsub(".tif", "", paste(strsplit(t.files[j], split = "_")[[1]][2:5], collapse = "_"))
    t.values = aggregate(t.raster, t.ref, mean)
    T.LIST[[t.var.name]] =  as.data.frame(t.values)[,2]

    #remove temp j files
    rm(t.raster, t.var.name)

    #print process
    print(paste(i, "/", length(site), "-", j, "/", length(t.files)))

  }; rm(j)# end j

  # df with all variables
  t.df = do.call(cbind, T.LIST)

  if(i == 1){write.csv(names(data.frame(t.df)), "1km_ALS_varnames.csv", row.names = F)}

  T.LIST.ALL[[i]] = data.frame(t.coords, site = site[i], t.df)
  
  # create ref.crs to create shape file after looping
  if(i == length(site)){
    ref.crs = sf::st_crs(t.ref)
  }
  
  #remove temp objects
  #rm(t.files, T.LIST, t.ref, t.coords, t.df, t.sf)
  rm(t.files, T.LIST, t.ref, t.coords, t.df)

  #print
  print(paste(i, "/", length(site), "- DONE"))

}; rm(i)# end looping i

# join all sites
all = do.call(rbind, T.LIST.ALL)
names(all)

#fix some names
names(all) = gsub("_NA", "", names(all))

# The option above is to save shape file (I disabled it and chose to save just the data.frame)
# convert data frame to sf object
t.sf = st_as_sf(x = all,
                coords = c("X", "Y"),
                crs = ref.crs)
t.sf = as(t.sf, "Spatial")
# plot(t.sf, axes = T)
#write final shape points
shapefile(t.sf, paste0(path.out, "ALL_ALS_1km.shp"), overwrite = TRUE)

# write final data.frame as csv
write.csv(all, paste0(path.out, "ALL_ALS_1km.csv"), row.names = F)


################################################################################
################################################################################
################################################################################



################################################################################
################################################################################
################################################################################
# 3 - Join GEDI variables

# rm all objects
rm(list = ls())

#load packages
pacman::p_load(sf, raster, stars)

# read shp all ALS points
all = read.csv("1_output/8_Layers_1km/ALL_Points/ALL_ALS_1km.csv")

# ref coordinates
ref.crs = st_crs(stars::read_stars("1_output/GEDI_maps_1km/_rh80_M1.tif"))
coords = st_as_sf(x = all, 
                 coords = c("X", "Y"),
                 crs = ref.crs)
rm(ref.crs)

# path of 1km GEDI maps
path.1km = "1_output/GEDI_maps_1km"
files.1km = list.files(path.1km, pattern = "tif")
if(length(grep("aux.", files.1km))>0){
  files.1km = files.1km[-c(grep("aux.", files.1km))]
}

#looping t extract GEDI variables from ALS_1km points
T.LIST = list()
for(i in 1:length(files.1km)){
  
  # get 1km ref cropand coordinates 
  
  #load raster j
  t.raster = stars::read_stars(file.path(path.1km, files.1km[i]))
  
  # get var name, extract values and fill final list
  t.var.name = gsub(".tif", "", paste(strsplit(files.1km[i], split = "_")[[1]][2:3], collapse = "_"))
  t.values = aggregate(t.raster, coords, mean)
  T.LIST[[t.var.name]] =  as.data.frame(t.values)[,2]
  
  #remove temp j files
  rm(t.raster, t.var.name)
  
  #print process
  print(paste(i, "/", length(files.1km)))
  
}# end i

# df with all variables
t.df = do.call(cbind, T.LIST)
t.df = data.frame(t.df)

# rename variables and join ALS and GEDI table
names(all)[-c(1:3)] = paste0("ALS_", names(all)[-c(1:3)])
names(t.df) = paste0("GED_", names(t.df))
t.df = cbind(all, t.df)

names(t.df)

#write table in csv
write.csv(t.df, "1_output/8_Layers_1km/ALL_Points/ALL_ALSandGEDI_1km.csv",row.names = F)

## write spatial file (SHP)
# convert data frame to sf object
t.sf = st_as_sf(x = t.df,
                coords = c("X", "Y"),
                crs = st_crs(coords))

t.sf = as(t.sf, "Spatial")
plot(t.sf, axes = T)

#write
shapefile(t.sf, "1_output/8_Layers_1km/ALL_Points/ALL_ALSandGEDI_1km.shp", overwrite = T)


## DONE!
#------------------------------------------------------------------------------#

