

## Amazeco: 
#Covering the Amazon with an Ecosystem Structure EBV product combining 
#satellite and airborne lidar

#========================================================================#
#------------------------------------------------------------------------#
# Title: Create data frame of 1km with Gedi, Gedilike and ALS variables  #    
# Author: Danilo R. A. de Almeida                                        #
# contact: daniloflorestas@gmail.com                                     #
#------------------------------------------------------------------------#


# rm all objects
rm(list = ls())

#load packages
pacman::p_load(raster)

# load ALS and GEDI 1 km data
gedi.als = read.csv("1_output/8_Layers_1km/ALL_Points/ALL_ALSandGEDI_1km.csv")
sites = unique(gedi.als$site)

# load gedilike 1km
gedilike = read.csv("1_output/10_Layers_1km_Gedilike/ALL_Points/ALL_GediLike_1km.csv")

# check columns
names(gedilike)
names(gedi.als)

#reorder dfs by sites
gedi.als = gedi.als[order(gedi.als$site),]
gedilike = gedilike[order(gedilike$site),]

#visual check of match coordiaates
t = data.frame(gedi.als[,c("X", "Y")], gedilike[,c("x", "y")])
nrow(t)
sum(t[,1] == t[,3])
sum(t[,2] == t[,4])
# All above sould have the same value
names(gedilike)
final = cbind(gedi.als, gedilike[,-c(1:3)])

# write final table
path.output = "1_output/12_Table_Gedilike_Gedi_ALS_1km"
dir.create(path.output, showWarnings = F)
write.csv(final, file.path(path.output, "ALL_AlsGediGedilike_1km.csv"), row.names = F)



## END
#------------------------------------------------------------------------------#