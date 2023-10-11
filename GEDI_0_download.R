
## Amazeco: 
#Covering the Amazon with an Ecosystem Structure EBV product combining 
#satellite and airborne lidar

#========================================================================#
#------------------------------------------------------------------------#
# Title: Download GEDI data using rGEDI package                          #    
# Author: Danilo R. A. de Almeida, Carlos A. Silva and Caio Hamamura     #
# contact: daniloflorestas@gmail.com                                     #
#                                                                        #
#------------------------------------------------------------------------#


# remove all objects
rm(list = ls())

# loading packages
if(!require("pacman")) install.packages("pacman") 
pacman::p_load("rnaturalearth", "rnaturalearthdata", "rGEDI", "raster")


## boundary box (bbox) for the Amazon biome
aoi = raster::shapefile("Shapes/DUC_A01_2020_LiDAR_lasboundary.shp")
ul_lat = aoi@bbox[2, 2] # Upper left latitude for the bounding box
ul_lon = aoi@bbox[1, 1] # Upper left longitude for the bounding box
lr_lat = aoi@bbox[2, 1] # Lower right latitude for the bounding box
lr_lon = aoi@bbox[1, 2] # Lower right longitude for the bounding box

# directory for saving GEDI data
outdir_2a = "C:/Lidar_EBA/RProcessing_Tutorial/GEDI_Data/L2A"
outdir_2b = "C:/Lidar_EBA/RProcessing_Tutorial/GEDI_Data/L2B"

# create world map to plot boxes projections
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# Check the bbox in the world map
windows()
plot(world[1]$geometry, axes=T, xlim=c(ul_lon, lr_lon), ylim=c(lr_lat,ul_lat))
points(c(ul_lon,lr_lon),c(ul_lat,lr_lat), col="red", cex=2, pch=16)
polygon(c(ul_lon,lr_lon,lr_lon,ul_lon),c(ul_lat,ul_lat,lr_lat,lr_lat))
grid()


# Gedifinder function (finds the GEDI data for a given region of interest and date range)
#Level1B - Geolocation Waveforms
#Level2A - Ground elevation, canopy top height, relative height (RH) metrics
#Level2B - Canopy cover fraction (CCF), CCF profile, LAI and LAI profile

g2a = gedifinder(product="GEDI02_A",ul_lat, ul_lon, lr_lat, lr_lon,  version = "002")
g2b = gedifinder(product="GEDI02_B",ul_lat, ul_lon, lr_lat, lr_lon,  version = "002")

# Downloading GEDI data
# add your NASA's Earth Data login and password info
gediDownload(g2a, outdir_2a)
gediDownload(g2b, outdir_2b)



### END
#------------------------------------------------------------------------------#
