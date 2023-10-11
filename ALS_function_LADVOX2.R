# Functions
myMetrics2 = function(Z, maxZ){
  heightSlices = floor(Z) # Round down
  zSlice = data.table::data.table(Z=Z, heightSlices=heightSlices) # Create a data.table (Z, slices)
  sliceCount = zSlice[, as.double(length(Z)), by = .(heightSlices)] # Count number of returns by slice
  
  ##############################################
  # Add columns to equalize number of columns
  ##############################################
  colRange = 0:maxZ
  addToList = colRange[!(colRange %in% sliceCount$heightSlices)]
  if(length(addToList) != 0){ 
    bindDt = data.frame(heightSlices = addToList, V1=0)
  }else{
    bindDt = data.frame()}
  sliceCount = rbind(sliceCount, bindDt)
  
  sliceCount = sliceCount[order(sliceCount$heightSlices)] # Order by height
  
  colNames = as.character(sliceCount$heightSlices)
  colNames[1] = "ground_0_1m"
  colNames[-1] = paste0("pulses_", colNames[-1], "_", sliceCount$heightSlices[-1]+1, "m")
  metrics = list()
  metrics[colNames] = sliceCount$V1
  
  return(metrics)
  
} #end function myMetrics

## LAD.VOXEL2S functions
lad.voxels2 = function(.las, grain.size = 1, k=1){
  
  #empty list object that will be fueling with binneds data.frames
  LAD_VOXELS = list()
  Z = NA
  
  #load normalized las cloud
  #.las = lidR::readLAS(normlas.file)
  
  .las@data$Z[.las@data$Z < 0] = 0
  
  maxZ = floor(max(.las@data$Z))
  lazyFunc = formula(paste0("~myMetrics2(Z, ", maxZ, ")"))
  t.binneds = lidR::grid_metrics(.las, lazyFunc, res = grain.size,
                                 start = c(min(.las@data$X), max(.las@data$Y)))
  t.binneds = data.frame(sp::coordinates(t.binneds), raster::values(t.binneds))
  names(t.binneds)[1:2] = c("X", "Y")
  
  
  #getting the coordinates X and Y
  #t.binneds$X = coordinates(t.binneds)[,1]
  #t.binneds$Y = coordinates(t.binneds)[,2]
  #t.binneds = as.data.frame(t.binneds) #transforming in a data.frame
  
  #clip product by las files limits
  #t.binneds = t.binneds[t.binneds$X < xmax(.las) &
  #                        t.binneds$X > xmin(.las) &
  #                        t.binneds$Y > ymin(.las) &
  #                        t.binneds$Y < ymax(.las),]
  
  
  #select ground returns
  ground.returns = t.binneds[, grep("ground", names(t.binneds))]
  
  #select columns vegetation above 1m:
  if(nrow(t.binneds) != 1){ #this if is necessary when grain size is the whole plot
    pulses.profile.dz1 = t.binneds[, c(grep("pulses", names(t.binneds)))]
  }else{
    pulses.profile.dz1 = data.frame(matrix(as.numeric(as.character(t.binneds[, c(grep("pulses", names(t.binneds)))])), ncol = length(grep("pulses", names(t.binneds)))))
    names(pulses.profile.dz1) = names(t.binneds)[c(grep("pulses", names(t.binneds)))]
  }
  
  #invert data.frames for the sky be first
  pulses.profile.dz1 = pulses.profile.dz1[,length(pulses.profile.dz1):1] #invert columns
  
  #add grounds returns (0-1m)
  pulses.profile.dz1 = cbind(pulses.profile.dz1, ground.returns)
  rm(ground.returns)
  
  ### total matriz and cumsum.matrix:
  total.pulses.matrix.dz1 = matrix(apply(pulses.profile.dz1, 1, sum), ncol = length(pulses.profile.dz1), nrow = nrow(pulses.profile.dz1))
  cumsum.matrix.dz1 = matrix(apply(pulses.profile.dz1, 1, cumsum), ncol = length(pulses.profile.dz1), nrow = nrow(pulses.profile.dz1), byrow = T)
  
  rm(pulses.profile.dz1)
  
  #Pulses out for each voxel
  pulse.out.dz1 = total.pulses.matrix.dz1 - cumsum.matrix.dz1
  
  #The pulses.out of voxel 1 is the pulses.in of voxel 2 and so on...
  #Therefore, pulse.in is pulse.out without the last line and adding in the
  #first line the total pulses:
  if(nrow(t.binneds) != 1){ #if used when grain size of the whole plot
    pulse.in.dz1 <- cbind(total.pulses.matrix.dz1[,1], pulse.out.dz1[,-c(ncol(pulse.out.dz1))])
  }else{
    pulse.in.dz1 <- c(total.pulses.matrix.dz1[,1], pulse.out.dz1[,-c(ncol(pulse.out.dz1))])
  } #enf if
  
  rm(total.pulses.matrix.dz1, cumsum.matrix.dz1)
  
  # MacArthur-Horn eqquation
  # LAD = ln(S_bottom/S_top)*(1/(dz*K))
  #k value for LAD equation
  dz = 1
  
  LAD.dz1 = log(pulse.in.dz1/pulse.out.dz1) * 1/k * 1/dz
  
  rm(pulse.in.dz1, pulse.out.dz1)
  
  # Remove infinite and NaN values
  #Inf ocorre qndo pulses.out eh zero
  #NaN ocorre qndo pulses.in eh zero
  LAD.dz1[is.infinite(LAD.dz1)] <- NA; LAD.dz1[is.nan(LAD.dz1)] <- NA;
  
  #remove the first 1 meter close to the ground (and the ground too)
  LAD.dz1 = LAD.dz1[, -c(ncol(LAD.dz1))]
  
  #fuel list object
  LAD_VOXELS[["LAD"]] = LAD.dz1
  LAD_VOXELS[["coordenates"]] = t.binneds[,c("X", "Y")]
  
  rm(LAD.dz1, t.binneds)
  
  return(LAD_VOXELS)
}#End function

