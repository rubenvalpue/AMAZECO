
## Amazeco: 
#Covering the Amazon with an Ecosystem Structure EBV product combining 
#satellite and airborne lidar

#========================================================================#
#------------------------------------------------------------------------#
# Title: Leaf area density modeling (25m)                                #    
# Author: Danilo R. A. de Almeida                                        #
# contact: daniloflorestas@gmail.com                                     #
#                                                                        #
# Products:                                                              #
#           - Leaf Area Density (LAD) and Leaf Area Index (LAI) modeling #
#           - LAI, LAI.under, LAHV and LAWV                              #
#                                                                        #
#------------------------------------------------------------------------#

#========================================================================#


#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
### ENVIRONMENT SETTINGS ----

# Clean environment
rm(list = ls(all = T))

# Load packages
pacman::p_load(units,
               rgdal, 
               rgeos,
               sf,
               lidR,
               raster,
               ForestGapR,
               leafR,
               tictoc,
               future, 
               sf, 
               stars, 
               devtools,
               mapview)

#------------------------------------------------------------------------------#
## Set paths
# IMPORTANT: To use this script it is necessary just to set:
# - a folder.path: where everything will be saved
# - a cloud.path with lidar files (las or laz).

# IMPORTANT: All clouds (las or laz) files need to be named with the site name 
#at the beginning using "_" as separator. 
#We recommend using site name as a code with tree letters.
# Example of file names for "Ducke" site: 
#"DUC_cloud1.laz", "DUC_cloud2.laz", "DUC_cloud3.laz", ...
# This tagged name is important for the next Scripts.
#If you are not going to use the next Scripts, the file names do not matter.

# Set main folder
folder.path = "C:/Users/danil/OneDrive - usp.br/POSDOC/DOCUMENTOS/Bangor_POSDOC/RProcessing_Tutorial"

# set normalize cloud path
cloud.path = "1_output/5_NormClouds"

## End of Set paths 
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
## Set parameters

# lad modeling parameters
lad.res = 25
angle_filter = 10 # filter pulses from max angle to lad modeling
homogenize.density = 4

#end of Set parameters
#------------------------------------------------------------------------------#


### NO NEED TO CHANGE ANYTHING ELSE FROM HERE.
### NO NEED TO CHANGE ANYTHING ELSE FROM HERE.
# create output paths
output.path = paste0(folder.path,"/1_output/")
rm(folder.path)

# Cloud files and sites 
ctg = lidR::catalog(cloud.path, recursive = T)
sites = unique(unlist(lapply(strsplit(basename(ctg@data$filename), split = "_"), function(x) (x[1]))))

# Visual plot using "mapview" package
#plot(ctg, mapview = TRUE, map.type = "Esri.WorldImagery")
#plot(ctg)

# load function LADVOX2 function
source("ALS_function_LADVOX2.R") # read function "lad.voxels2()"
#This function replaces "lad.voxels()" function from leafR package. 
#With this new function there is no need to write a temp las file. 

# END setting environment
#----------------------------------------------------------#





#----------------------------------------------------------#

### LOOPING AND PROCESSING ----
for(i in 1:length(sites)){ 
  
  #tic() and toc() functions are used to print the time-lapse after the looping
  tic()
  
  t.files = grep(sites[i], basename(ctg@data$filename))
  t.las = readMSLAS(catalog(ctg@data$filename[t.files]), 
                    select="xyz", 
                    filter = "-drop_abs_scan_angle_above 10 -keep_number_of_returns 1")
  rm(t.files)
  
  #----------------------------------------------------------------------------#
  ### LAD and LAI ----
  print(paste0(i,"/", length(sites) ," ", sites[i], "_LAD processing"))
  
  ####LAD MODELING
  
  ### LAD and LAI ----
  
  # if there is returns above 1 m 
  if(length(t.las@data$Z) > 1 & max(t.las@data$Z) > 1){
    
    #Visual check of the previous density
    #density = lidR::grid_density(t.las, lad.res)
    #plot(density)
    
    # homogeneize point cloud
    t.las = decimate_points(t.las, homogenize(density = homogenize.density, res = lad.res))
    
    # Visual check
    # density = lidR::grid_density(t.las, lad.res)
    # plot(density)
    # rm(density)
    
    print(paste0(i,"/",length(sites), " ", sites[i],"_VOXELSLAD processing"))
    
    # Voxels lad
    VOXELS_LAD = lad.voxels2(t.las, grain.size = lad.res, k=1)
    
    
    print(paste0(i,"/",length(ctg),"_LAD_VOXELSLAD_ok"))
    
    # LAD profile
    lad_profile = lad.profile(VOXELS_LAD)
    # plot(lad_profile$height ~ lad_profile$lad, type = 'l')
    
    ## Create stack LADs and LAI maps
    df = data.frame(VOXELS_LAD$LAD)
    df = df[,ncol(df):1]
    
    
    ## If: eliminating LAD raster processing for very small clouds (small areas got ERROR)
    if(nrow(df) > 10){# if for too small areas
      
      LAD.RASTERS = list()
      for(j in 1:ncol(df)){
        LAD.RASTERS[[j]] = rasterFromXYZ(data.frame(x = VOXELS_LAD$coordenates[1], 
                                                    y = VOXELS_LAD$coordenates[2], 
                                                    z = df[,j]))
      }; rm(j)# end i
      
      lad.raster = do.call(stack, LAD.RASTERS) # stack created
      rm(LAD.RASTERS)
      
      # create NA maks 
      NA.mask = apply(df, 1, function(x) (sum(is.na(x)) == ncol(df)))
      NA.mask = c(NA.mask == FALSE) #invert true and false
      NA.mask = rasterFromXYZ(data.frame(x = VOXELS_LAD$coordenates[1], 
                                         y = VOXELS_LAD$coordenates[2], 
                                         z = NA.mask))
      NA.mask[NA.mask == 0] = NA
      # plot(NA.mask)
      
      # LAI raster
      lai_raster = sum(lad.raster, na.rm = T)
      lai_raster = lai_raster*NA.mask
      
      # LAI_under raster
      understory.limit.height = 5
      lai_under = sum(lad.raster[[which(lad_profile$height < understory.limit.height)]], na.rm = T)
      lai_under = lai_under*NA.mask
      rm(understory.limit.height)
      
      # LAHV
      df.h = matrix(rep(1:ncol(df), each = nrow(df)), 
                    ncol = ncol(df), nrow = nrow(df))
      
      lahv = df*df.h
      rm(df, df.h)
      lahv = apply(lahv, 1, sum, na.rm = T)
      
      lahv = rasterFromXYZ(data.frame(x = VOXELS_LAD$coordenates[1], 
                                      y = VOXELS_LAD$coordenates[2], 
                                      z = lahv))
      lahv = lahv*NA.mask
      # plot(lahv)
      
      lawv = lahv/lai_raster
      # plot(lawv)
      
      # Visual check
      # t = stack(lai_raster, lai_under)
      # names(t) = c("LAI", "LAI_under")
      # plot(t, zlim = c(0, 5.5))
      # rm(t)
  
      crs(lad.raster) = crs(lai_raster) = crs(lai_under) = crs(lahv) = crs(lawv) = crs(t.las)
      rm(t.las)
      
      # write raster LADs-stack, lai, and lai_under
      writeRaster(lad.raster, paste0(output.path,'/7_Layers_25m/',sites[i],'_LadStack_25mv2.tif'), overwrite=TRUE)
      writeRaster(lai_raster, paste0(output.path,'/7_Layers_25m/',sites[i],'_lai_25mv2.tif'), overwrite=TRUE)
      writeRaster(lai_under, paste0(output.path,'/7_Layers_25m/',sites[i],'_lai_understory_25mv2.tif'), overwrite=TRUE)
      writeRaster(lahv, paste0(output.path,'/7_Layers_25m/',sites[i],'_lahv_25mv2.tif'), overwrite=TRUE)
      writeRaster(lawv, paste0(output.path,'/7_Layers_25m/',sites[i],'_lawv_25mv2.tif'), overwrite=TRUE)
    
      rm(lad.raster, lai_raster, lai_under, lahv, lawv, NA.mask)
      }# end if for to small areas
    ## END IF SMALL AREAS
    
    print(paste0(i,"/",length(ctg),"_LAD_Layers_ok"))
    
  }# end if there is returns
  
  # END LAD and LAI
  #----------------------------------------------------------------------------#
  
  toc()
  
  
}; rm(i)### End loop  ----

## LOOPING END
#----------------------------------------------------------#



