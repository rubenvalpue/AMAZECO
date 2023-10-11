

## Amazeco: 
#Covering the Amazon with an Ecosystem Structure EBV product combining 
#satellite and airborne lidar

#========================================================================#
#------------------------------------------------------------------------#
# Title: Complete dataframe of 1km GEDI.like metrics                     #    
# Author: Danilo R. A. de Almeida                                        #
# contact: daniloflorestas@gmail.com                                     #
#------------------------------------------------------------------------#


################################################################################
################################################################################
################################################################################

#####
# 1 - create a one df by points (pixel centroids) for all 1km gedi.like rasters 

#remove all objects
rm(list = ls())

pacman::p_load(raster, sf, stars)

# path of 1km gedilike maps
path.1km = "1_output/10_Layers_1km_Gedilike"
files.1km = list.files(path.1km, pattern = "tif")
if(length(grep("aux.", files.1km))>0){
  files.1km = files.1km[-c(grep("aux.", files.1km))]
}

site = unique(substr(files.1km, 1, 3))

#output path
path.out = "1_output/10_Layers_1km_Gedilike/ALL_Points/"
dir.create(path.out, showWarnings = FALSE)

#ref coordinates
ref.1km = read.csv("1_output/8_Layers_1km/ALL_Points/ALL_ALSandGEDI_1km.csv")

# get wgs64 crs
ref.crs = stars::read_stars("1_output/GEDI_maps_1km/_rh80_M1.tif")
ref.crs = st_crs(ref.crs)

## Looping i do shape points for each site
# points use coordinates of GEDI 1km grid reference

T.LIST.ALL = list()
for(i in 1:length(site)){

  # set site o
  t.files = grep(site[i], files.1km, value = T)

  T.LIST = list()
  for(j in 1:length(t.files)){

    #load and crs transform raster j
    t.raster = stars::read_stars(file.path(path.1km, t.files[j]))
    
    # get coordinates
    if(j == 1){t.coords = st_coordinates(t.raster)}
    
    # get var name, and values and fill final list
    t.var.name = gsub(".tif", "", paste(strsplit(t.files[j], split = "_")[[1]][2:4], collapse = "_"))
    T.LIST[[t.var.name]] =  as.data.frame(t.raster)[,3]

    #remove temp j files
    rm(t.raster, t.var.name)

    #print process
    print(paste(i, "/", length(site), "-", j, "/", length(t.files)))

  }# end j


  # df with all variables
  t.df = do.call(cbind, T.LIST)

  if(i == 1){write.csv(names(data.frame(t.df)), "1km_GediLike_varnames.csv", row.names = F)}
  
  # write final csv
  T.LIST.ALL[[site[i]]] = data.frame(site = site[i], t.coords, t.df)
  
  #remove temp objects
  rm(t.files, T.LIST, t.coords, t.df)
  
  #print
  print(paste(i, "/", length(site), "- DONE"))

}# end looping i


all = do.call(rbind, T.LIST.ALL)

write.csv(all, file.path(path.out, "ALL_GediLike_1km.csv"), row.names = F)

# convert data frame to sf object
t.sf = st_as_sf(x = all,
                 coords = c("x", "y"),
                 crs = ref.crs)
t.sf = as(t.sf, "Spatial")
# plot(t.sf, axes = T)

#write
shapefile(t.sf, file.path(path.out, "ALL_GediLike_1km.shp"))


# END
################################################################################
################################################################################
################################################################################




