
## Amazeco: 
#Covering the Amazon with an Ecosystem Structure EBV product combining 
#satellite and airborne lidar

#========================================================================#
#------------------------------------------------------------------------#
# Title: Airborne lidar data processing                                  #    
# Author: Danilo R. A. de Almeida                                        #
# contact: daniloflorestas@gmail.com                                     #
#                                                                        #
# Products:                                                              #
#           - Normalized denoised clouds                                 #
#           - Digital Terrain Model (DTM)                                #
#           - Canopy Height Model (CHM)                                  #
#           - Polygon layers of the area covered by the lidar data       #
#------------------------------------------------------------------------#

#========================================================================#


#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
### ENVIRONMENT SETTINGS ----

# Clean environment
rm(list = ls(all = T))

# Load packages
if(!require("pacman")) install.packages("pacman") 
pacman::p_load(units,
               rgdal, 
               rgeos,
               sf,
               lidR,
               raster,
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

# set cloud paths
cloud.path = "C:/Lidar_EBA/RProcessing_Tutorial/ALS_Data/DUC_A01_2020"

## End of Set paths 
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
## Set parameters

# set lad resolution (grain.size)
cloud.classified = FALSE # If the cloud points is already classified set it as "TRUE". If you don't know, leave it as "FALSE" 
break.tiles.size = NULL # set the break tiles size if cloud files is too big.

# lad modeling parameters
angle_filter = 10 # filter pulses from max angle to lad modeling
homogenize.density = 4

#end of Set parameters
#------------------------------------------------------------------------------#


### NO NEED TO CHANGE ANYTHING ELSE FROM HERE.
### NO NEED TO CHANGE ANYTHING ELSE FROM HERE.
# create output paths
output.path = paste0(folder.path,"/1_output/")
dir.create(output.path, showWarnings = F)
dir.create(paste0(output.path, "/1_DTM/"), showWarnings = F) # Digita terrain models
dir.create(paste0(output.path, "/2_CHM/"), showWarnings = F) # Canopy hight model
dir.create(paste0(output.path, "/3_GAP/"), showWarnings = F) # Gap layers
dir.create(paste0(output.path, "/4_LAD/"), showWarnings = F) # LAD layers 
dir.create(paste0(output.path, "/5_NormClouds/"), showWarnings = F) # Normalized and denoised cloud
dir.create(paste0(output.path, "/6_Shapes/"), showWarnings = F) # Polygon layers of the cloud

# create "TempProcessingTiles" path (temporary files using during processing)
temp.path = paste0(folder.path, "/2_TempProcessingTiles/")
dir.create(temp.path, showWarnings = F)
rm(folder.path)

# Cloud files' list
ctg = lidR::catalog(cloud.path, recursive = T)

# Visual plot using "mapview" package
#plot(ctg, mapview = TRUE, map.type = "Esri.WorldImagery")
#plot(ctg)


# Lasnoise function used into looping
lasnoise = function(cluster, res = 5, n = 6){
  las = readLAS(cluster)
  if (is.empty(las)) return(NULL)
  las = classify_noise(las, ivf(res = res, n = n))
  las = filter_poi(las, Classification != 18 & Z < 100 & Z > -3)
  
  return(las)
}

# lad grain resolution
lad.res = 10


# load function LADVOX2 function
source("ALS_function_LADVOX2.R") # read function "lad.voxels2()"
#This function replaces "lad.voxels()" function from leafR package. 
#With this new function there is no need to write a temp las file. 

# END setting environment
#----------------------------------------------------------#





#----------------------------------------------------------#

### LOOPING AND PROCESSING ----
for(i in 1:length(ctg)){ 
  
  #tic() and toc() functions are used to print the time-lapse after the looping
  tic()
  
  # Clouds pre-processing and Digital Models ----
  
  # selecting cloud "i"
  t.ctg = ctg[i, ]
  t.cloudname = basename(t.ctg@data$filename)
  t.cloudname = gsub(".las|.laz","", t.cloudname)
  
  # temporary i cloud path
  t.cloudinprocess = paste0(temp.path, t.cloudname)
  dir.create(t.cloudinprocess, showWarnings = F) # Polygon layers of the cloud
  
  #----------------------------------------------------------------------------#
  # Begining las processing
  
  if(cloud.classified == FALSE){
    
    ### CLASSIFY ----
    print(paste0("Classifying cloud: ", i,"_", t.cloudname))
    
    # opt_chunk_buffer(t.ctg) = 10 # set tile buffer. Default is 30m
    if(!is.null(break.tiles.size)){
      opt_chunk_size(t.ctg) = break.tiles.size # set tile size
    }
    
    opt_output_files(t.ctg) = paste0(t.cloudinprocess, "/1_classifed/{ID}_classified") # set ouput path
    # plot(t.ctg, chunk = T)
    
    classify_ground(t.ctg, csf(), last_returns = F)
    
    t.ctg <- readLAScatalog(paste0(t.cloudinprocess, "/1_classifed/"))
    
  }# end if cloud.classified
  
  #----------------------------------------------------------------------------#
  
  
  #----------------------------------------------------------------------------#
  ### DTM ----
  print(paste0(i,"_",t.cloudname,"_DTM"))
  
  # opt_chunk_buffer(t.ctg) = 10 # set tile buffer. Default is 30m
  if(!is.null(break.tiles.size)){
    opt_chunk_size(t.ctg) = break.tiles.size # set tile size
  }
  opt_output_files(t.ctg) = paste0(t.cloudinprocess, "/2_DTM/{ID}_dtm")
  
  # plot(t.ctg, chunk = T)
  
  # run DTM
  dtm = grid_terrain(t.ctg, tin(), res = 1)
  
  # crs(dtm) # checking crs
  dtm[dtm < 0] = NA # filtering negative values
  # plot(dtm) # visual check
  # Writing raster
  raster::writeRaster(dtm, paste0(output.path, "/1_DTM/DTM_", t.cloudname, ".tif"), overwrite=TRUE)
  
  # Create shape file polygon
  t.shape = dtm
  raster::values(t.shape) = raster::values(t.shape)/raster::values(t.shape)
  t.shape = sf::st_as_sf(stars::st_as_stars(t.shape),as_points = FALSE, merge = TRUE)
  t.shape$Cloud = t.cloudname
  t.shape = t.shape[,c(2,3)]
  
  # Visual Check
  #plot(dtm); plot(as_Spatial(t.shape), add = T) #visual check
  
  #write shapefile 
  st_write(t.shape,
           paste0(output.path, "/6_Shapes/SHP_", t.cloudname, ".shp"), 
           driver = "ESRI Shapefile", append=TRUE)
  rm(t.shape)
  # END DTM
  #----------------------------------------------------------------------------#
  
  
  #----------------------------------------------------------------------------#
  ### NORMALIZE ----
  print(paste0(i,"_",t.cloudname,"_normalize"))
  
  opt_output_files(t.ctg) = paste0(t.cloudinprocess, "/3_Norm/{ID}_norm")
  normalize_height(t.ctg, tin(), res = 1, na.rm = T)
  
  #END Normalize
  #----------------------------------------------------------------------------#
  
  Sys.sleep(2) # 2s wait period to avoid ERROR
  
  #----------------------------------------------------------------------------#
  ## DENOISE ----
  print(paste0(i,"_",t.cloudname,"_denoise"))
  
  t.ctg = readLAScatalog(paste0(t.cloudinprocess, "/3_Norm/"))
  opt_output_files(t.ctg) = paste0(t.cloudinprocess, "/4_NormDenoised/{ID}_normdenoise")
  
  opt    <- list(need_buffer = F,   # catalog_apply will throw an error if buffer = 0
                 automerge   = T)   # catalog_apply will merge the outputs into a single object
  
  output = catalog_sapply(t.ctg, lasnoise, .options = opt)
  
  Sys.sleep(2) # 1s wait period to avoid ERROR
  
  # read all denoised norm tiles and save it as one laz file
  print(paste0(i,"_",t.cloudname,"_readMSLAS"))
  
  t.las = readMSLAS(output, select="xyzarc")
  lidR::writeLAS(t.las, file = paste0(output.path, "/5_NormClouds/", t.cloudname, "_norm.laz"))
  
  
  # t.las = lidR::readLAS(paste0(output.path,'/5_NormClouds/',t.cloudname,'_norm.laz'), select="xyzar")
  # END Denoise
  #----------------------------------------------------------------------------#
  
  
  #----------------------------------------------------------------------------#
  ### CHM ----
  t.ctg = readLAScatalog(paste0(t.cloudinprocess, "/4_NormDenoised/"))
  opt_output_files(t.ctg) = paste0(t.cloudinprocess, "/5_CHM/{ID}_chm")
  # plot(t.ctg, chunk = T)
  
  chm = grid_canopy(t.ctg, res = 1, p2r())
  
  #crs(dtm) # checking crs
  chm[chm < 0] = 0 # filtering negative values
  
  #Visual check
  # x11(); plot(chm, col=height.colors(40))
  
  # Writing raster
  raster::writeRaster(chm, paste0(output.path, "/2_CHM/CHM_", t.cloudname, ".tif"), overwrite=TRUE)
  # chm =  raster(paste0(output.path, "/2_CHM/CHM_", t.cloudname, ".tif"))
  
  #----------------------------------------------------------------------------#
  
  
  #----------------------------------------------------------------------------#
  ### LAD and LAI ----
  print(paste0(i,"/", length(ctg) ," ", t.cloudname, "_LAD processing"))
  
  ####LAD MODELING
 
  ### LAD and LAI ----
  
  # filtering just 10 degrees from nadir
  scan_angle = grep("Scan",names(t.las@data),value=T)
  if (scan_angle == "ScanAngle"){
    # filtering just 10 degrees from nadir
    las_filtered = filter_poi(t.las, ScanAngle < angle_filter+1 & ScanAngle > (angle_filter+1)*-1)
  } else {
    # filtering just 10 degrees from nadir
    las_filtered = filter_poi(t.las, ScanAngleRank < angle_filter+1 & ScanAngleRank > (angle_filter+1)*-1)
  }# end ifelse
  rm(scan_angle, t.las)
  
  # filtering just first returns
  las_filtered = filter_poi(las_filtered, ReturnNumber == 1)
  
  # Visual check
  # plot(las_filtered)
  
  # if there is returns
  if(length(las_filtered@data$Z) > 1 & max(las_filtered@data$Z) > 1){
    
    #Visual check of the previous density
    # density = lidR::grid_density(las_filtered, lad.res)
    # plot(density)
    
    # homogeneize point cloud
    las_filtered = decimate_points(las_filtered, homogenize(density = homogenize.density, res = lad.res))
    
    # Visual check
    # density = lidR::grid_density(las_filtered, lad.res)
    # plot(density)
    # rm(density)
    
    print(paste0(i,"/",length(ctg), " ", t.cloudname,"_VOXELSLAD processing"))
    
    # Voxels lad
    VOXELS_LAD = lad.voxels2(las_filtered, grain.size = lad.res, k=1)
    
    # Save VOXELS_LAD as RData
    # save(VOXELS_LAD, file = paste0(output.path,'/4_LAD/', t.cloudname,'_LAD.Rdata'))
    
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
      
      crs(lad.raster) = crs(lai_raster) = crs(lai_under) = crs(lahv) = crs(lawv) = crs(las_filtered)
      
      # write raster LADs-stack, lai, and lai_under
      writeRaster(lad.raster, paste0(output.path,'/4_LAD/',t.cloudname,'_LadStack.tif'), overwrite=TRUE)
      writeRaster(lai_raster, paste0(output.path,'/4_LAD/',t.cloudname,'_lai.tif'), overwrite=TRUE)
      writeRaster(lai_under, paste0(output.path,'/4_LAD/',t.cloudname,'_lai_understory.tif'), overwrite=TRUE)
      writeRaster(lahv, paste0(output.path,'/4_LAD/',t.cloudname,'_lahv.tif'), overwrite=TRUE)
      writeRaster(lawv, paste0(output.path,'/4_LAD/',t.cloudname,'_lawv.tif'), overwrite=TRUE)
    }# end if for to small areas
    ## END IF SMALL AREAS
    
    # save LAD profile jpeg  
    lad.mean = lad_profile[,2]
    #lad.mean[lad.mean < 0.003] = NA
    
    jpeg(paste0(output.path,'/4_LAD/',t.cloudname,'_LADProfile.jpg'),unit = "in",h = 9, w = 9,res = 300)
    plot(NA, xlim=c(0,0.35),  
         ylim=c(1, 80),
         xlab= expression(paste("LAD (", m^2*m^-3, ")", sep = "")),
         ylab="Canopy height (m)", cex.axis=1, tck = 0, frame=F,
         yaxt = "n", cex.lab = 1.2, cex.axis=1.2)
    axis(2, at = seq(0,80, 5), las = 2, cex.axis = 1.3)
    lines(lad.mean[which(!is.na(lad.mean))], 1:length(which(!is.na(lad.mean))), lty = 1, lwd = 3, col="green3")
    dev.off()
    rm(lad.mean, lad_profile)
    
  }# end if there is returns
  
  # END LAD and LAI
  #----------------------------------------------------------------------------#
  
  toc()
  
  # remove temp folder
  unlink(t.cloudinprocess, recursive = TRUE)
  
  # remove temp objects
  rm(t.cloudname, t.cloudinprocess, output, opt)
  
}; rm(i)### End loop  ----

## LOOPING END
#----------------------------------------------------------#



