
## Amazeco: 
#Covering the Amazon with an Ecosystem Structure EBV product combining 
#satellite and airborne lidar

#========================================================================#
#------------------------------------------------------------------------#
# Title: Simulating 25m resolution GEDI metrics from ALS data (GEDI.like)#    
# Author: Danilo R. A. de Almeida                                        #
# contact: daniloflorestas@gmail.com                                     #
#------------------------------------------------------------------------#

# remove all objects
rm(list = ls())

# load packages
pacman::p_load(lidR, rGEDI, raster, stars, sf, stars)

# las files
cloud.path = "C:/Lidar_EBA/RProcessing_Tutorial/ALS_Data"
ctg = lidR::catalog(cloud.path, recursive = T)
site = unique(substr(basename(ctg@data$filename), 1, 3))

# get coord references from 25m layers files
path.25m = "1_output/7_Layers_25m"
files.25m = list.files(path.25m, pattern = ".tif")
if(length(grep("aux.", files.25m))>0){
  files.25m = files.25m[-c(grep("aux.", files.25m))]
}
files.25m = grep("CHM_MEAN_25m", files.25m, value = T)

# outputs path
path.out = "1_output/9_GediLike"
dir.create(path.out, showWarnings = FALSE)


##Looping
# this looping get each site (i) sample "n" coordinates, get them gedilike
#metrics and save a final data frame with "n" lines and 148 columns (metrics)

for(i in 1:length(site)){

  # get site i files (las and 25m)
  t.las = ctg@data$filename[grep(site[i], basename(ctg@data$filename), value = F)]
  t.25m = grep(site[i], files.25m, value = T)

  #get all coords of site i
  T.25M = list()
  for(j in 1:length(t.25m)){

    t.raster = raster::raster(file.path(path.25m, t.25m[j]))
    T.25M[[j]] = data.frame(coordinates(t.raster)[!is.na(values(t.raster)),])
    rm(t.raster)
    #T.25M[[j]] = st_coordinates(read_stars(file.path(path.25m, t.25m[j])))

    #points = sf::st_as_sf(t.layer, as_points = T)
  }; rm(j)
  coords = do.call(rbind, T.25M)

  # write temp txt coords
  write.table(coords, "2_TempProcessingTiles/listCoord.txt", row.names = F, col.names = F)

  #simulate gedi data from ALS (GEDI.like)
  gedilike = rGEDI::gediWFSimulator(input = t.las,
                                    output = paste0(path.out, "/", site[i], "_gedilike.h5"),
                                    listCoord = "2_TempProcessingTiles/listCoord.txt")
  # remove temp txt files
  file.remove("2_TempProcessingTiles/listCoord.txt")


  # #get values
  gedilike.metrics = gediWFMetrics(input = gedilike,
                                    outRoot = paste0(path.out, "/", site[i]),
                                    linkNoise= c(3.0103,0.95),
                                    maxDN= 4096,
                                    sWidth= 0.5,
                                    varScale= 3)


  # add coordinates
  gedilike.metrics$GeoX = coords[,1]
  gedilike.metrics$GeoY = coords[,2]

  write.csv(gedilike.metrics, paste0(path.out, "/", site[i], "_GediLikeMetrics.csv"))

  # close
  close(gedilike)

  #remove temp objects
  rm(t.las, t.25m, T.25M, coords, gedilike)

  #print process
  print(paste(i, "/", length(site)))

}# end looping

rm(i)




