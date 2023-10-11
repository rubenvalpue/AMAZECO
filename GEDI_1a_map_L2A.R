

## Amazeco: 
#Covering the Amazon with an Ecosystem Structure EBV product combining 
#satellite and airborne lidar

#========================================================================#
#------------------------------------------------------------------------#
# Title: Mapping GEDI L2A variables                                      #    
# Author: Danilo R. A. de Almeida, Carlos A. Silva and Caio Hamamura     #
# contact: daniloflorestas@gmail.com                                     #
#------------------------------------------------------------------------#

# Install package
# require(devtools)
# devtools::install_github("carlos-alberto-silva/rGEDI", "rasterizeHDF")

## requires to install GIT:
# https://gitforwindows.org/

# remove all objects
rm(list = ls())

# load packages
pacman::p_load(
  rGEDI,
  rgdal
)

source("GEDI_function_level2aRasterizeStats2.R")

# The directory paths where the H5 GEDI files are stored
gedi_dir = "C:/Lidar_EBA/RProcessing_Tutorial/GEDI_Data/L2A"

# A vector of metrics available from Level2B, as in the getLevel2BVPM documentation
metrics = c("rh80", "rh98")

# The root name for the raster output files, the pattern is out_rootmetriccount/m1/m2/m3/m4.tif. This should include the full path for the file.
out_root = "1_output/GEDI_maps_1km/temp_johan"
dir.create(out_root, showWarnings = FALSE)

#Area of interest
aoi = raster::shapefile("Shapes/DUC_A01_2020_LiDAR_lasboundary.shp")
ul_lat = aoi@bbox[2, 2] # Upper left latitude for the bounding box
ul_lon = aoi@bbox[1, 1] # Upper left longitude for the bounding box
lr_lat = aoi@bbox[2, 1] # Lower right latitude for the bounding box
lr_lon = aoi@bbox[1, 2] # Lower right longitude for the bounding box

# Resolution lon lat for the output raster in coordinates decimal degrees
degree.to.km.factor <- 1 / 111.32

default_finalizer = list()


###
agg.function <- ~data.table::data.table(
  n = length(x),
  M1 = mean(x,na.rm = T),
  M2 = e1071::moment(x, order = 2, center = TRUE, na.rm = T) * length(x)
)

default_agg_join <- function(x1, x2) {
  combined = data.table::data.table()
  x1$n[is.na(x1$n)] = 0
  x1$M1[is.na(x1$M1)] = 0
  x1$M2[is.na(x1$M2)] = 0
  
  combined$n = x1$n + x2$n
  
  delta = x2$M1 - x1$M1
  delta2 = delta * delta
  
  combined$M1 = (x1$n * x1$M1 + x2$n * x2$M1) / combined$n
  
  combined$M2 = x1$M2 + x2$M2 +
    delta2 * x1$n * x2$n / combined$n
  
  return(combined)
}


projstring = 'GEOGCS["WGS 84",
    DATUM["WGS_1984",
        SPHEROID["WGS 84",6378137,298.257223563,
            AUTHORITY["EPSG","7030"]],
        AUTHORITY["EPSG","6326"]],
    PRIMEM["Greenwich",0,
        AUTHORITY["EPSG","8901"]],
    UNIT["degree",0.01745329251994328,
        AUTHORITY["EPSG","9122"]],
    AUTHORITY["EPSG","4326"]]'

def_co = c("COMPRESS=DEFLATE",
           "BIGTIFF=YES",
           "TILED=YES",
           "BLOCKXSIZE=512",
           "BLOCKYSIZE=512"
)

level2aRasterizeStats2(
  l2aDir = gedi_dir,
  metrics = metrics,
  out_root = out_root,
  ul_lat = ul_lat,
  ul_lon = ul_lon,
  lr_lat = lr_lat,
  lr_lon = lr_lon,
  res = c(degree.to.km.factor, -degree.to.km.factor),
  creation_options = def_co,
  agg_function = agg.function,
  agg_join = default_agg_join,
  finalizer = default_finalizer,
  polygon_spdf=aoi
)

## END
#------------------------------------------------------------------------------#

