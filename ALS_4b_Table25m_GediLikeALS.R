
## Amazeco: 
#Covering the Amazon with an Ecosystem Structure EBV product combining 
#satellite and airborne lidar

#========================================================================#
#------------------------------------------------------------------------#
# Title: Create data frame of 25m with Gedilike and ALS variables        #    
# Author: Danilo R. A. de Almeida                                        #
# contact: daniloflorestas@gmail.com                                     #
#------------------------------------------------------------------------#

# remove all objects
rm(list = ls())

#load packages
pacman::p_load(sf)

# read csv all, with gedilike variables and coordinates
files.gedilike = list.files("1_output/9_GediLike/", pattern = ".csv", full.names = T)

# path 25m ALS
path.als = "1_output/7_Layers_25m"
files.als = list.files(path.als, pattern = ".tif", full.names = T)
if(length(grep("aux.|LadStack", files.als ))>0){
  files.als = files.als [-c(grep("aux.|LadStack", files.als))]
}

#sites vector to run into looping
site = substr(basename(files.gedilike), 1, 3)

# output path
path.out = "1_output/11_Table_Gedilike_ALS_25m/"
dir.create(path.out, showWarnings = FALSE)

### FINAL object that will be filled in the looping
FINAL.LIST = list()
final.coord = data.frame(matrix(NA, nrow = 0, ncol = 2))

#looping
for(i in 1:length(site)){
  
  # tif files of the site i
  t.files = grep(site[i], files.als, value = T)
  
  # get site coordinates
  t.gedilike = read.csv(files.gedilike[grep(site[i], basename(files.gedilike))])
  xy = t.gedilike[, c("GeoX", "GeoY")]; names(xy) = c("X", "Y")
  
  #create sf objects
  t.coords = st_as_sf(x = xy, 
                  coords = c("X", "Y"),
                  crs = st_crs(stars::read_stars(t.files[1])))
  
  final.coord = rbind(final.coord, st_coordinates(t.coords))
  
  # looping getting var values of each tif file
  T.LIST = list()
  
  for(j in 1:length(t.files)){
    
    #load raster j
    t.raster = stars::read_stars(t.files[j])
    
    # get var name, extract values and fill final list
    t.var.name = gsub(".tif", "", paste(strsplit(basename(t.files[j]), split = "_")[[1]][2:3], collapse = "_"))
    t.values = aggregate(t.raster, t.coords, mean)
    T.LIST[[t.var.name]] =  as.data.frame(t.values)[,2]
    
    #remove temp j files
    rm(t.raster, t.var.name)
    
    #print process
    print(paste("j:",j, "/", length(t.files), "- i:", i, "/", length(site)))
    
  };rm(j)# end i
  
  lapply(T.LIST, function(x) (length(x)))
  
  # fill final list
  t.df = data.frame(do.call(cbind, T.LIST))
  
  # fill final LIST
  names(t.df) = paste0("ALS", "_", names(t.df))
  t.gedilike = t.gedilike[-c(which(names(t.gedilike) == "X"))]
  names(t.gedilike) = paste0("GL", "_", names(t.gedilike))
  FINAL.LIST[[site[i]]] = data.frame(site = site[i], t.df, t.gedilike)
  
  # remove temp objects
  rm(t.files, xy, t.coords, T.LIST)
  
}; rm(i)

# rbind all
final = do.call(rbind, FINAL.LIST)
names(final)

#write
write.csv(final, file.path(path.out, "Table25m_Gedilike_ALS.csv"), row.names = F)



### DONE
################################################################################
#### curious analysis

