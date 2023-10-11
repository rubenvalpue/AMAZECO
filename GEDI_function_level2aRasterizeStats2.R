# pacman::p_load(
#   rGEDI,
#   rgdal
# )


# l2aDir = gedi_dir;
# metrics = metrics;
# out_root = out_root;
# ul_lat = ul_lat;
# ul_lon = ul_lon;
# lr_lat = lr_lat;
# lr_lon = lr_lon;
# res = c(degree.to.km.factor, -degree.to.km.factor);
# creation_options = def_co;
# agg_function = agg.function;
# agg_join = default_agg_join;
# finalizer = default_finalizer;
# polygon_spdf=aoi

level2aRasterizeStats2<-function (l2aDir, metrics, out_root, ul_lat, ul_lon, lr_lat,
                                  lr_lon, res, creation_options = def_co, agg_function = default_agg_function,
                                  agg_join = default_agg_join, finalizer = default_finalizer, polygon_spdf)
{
  block_inds = block_xind = block_yind = inds = quality_flag = lat_lowestmode = lon_lowestmode = x_block = y_block = x_ind = y_ind = NULL
  projstring = "GEOGCS[\"WGS 84\",\n    DATUM[\"WGS_1984\",\n        SPHEROID[\"WGS 84\",6378137,298.257223563,\n            AUTHORITY[\"EPSG\",\"7030\"]],\n        AUTHORITY[\"EPSG\",\"6326\"]],\n    PRIMEM[\"Greenwich\",0,\n        AUTHORITY[\"EPSG\",\"8901\"]],\n    UNIT[\"degree\",0.01745329251994328,\n        AUTHORITY[\"EPSG\",\"9122\"]],\n    AUTHORITY[\"EPSG\",\"4326\"]]"
  l2a_list = sapply(l2aDir, function(search_path) list.files(search_path,
                                                             pattern = "GEDI02_A.*h5", recursive = TRUE, full.name = TRUE))
  total_files = length(l2a_list)
  xres = res[1]
  yres = res[2]
  cols.coord = c("lat_lowestmode", "lon_lowestmode",
                 "quality_flag")
  metricCounter = 0
  nMetrics = length(metrics)
  func = lazyeval::f_interp(agg_function)
  call = lazyeval::as_call(func)
  x = 1
  stats = eval(call)
  classes = lapply(stats, class)
  stats = names(stats)
  for (metric in metrics) {
    metricCounter = metricCounter + 1
    message(sprintf("Metric %s (%d/%d)", metric, metricCounter,
                    nMetrics), appendLF = T)
    cols = c(cols.coord, metric)
    rast_paths = sprintf("%s_%s_%s.tif", out_root,
                         metric, stats)
    rasts = list()
    for (stat_ind in 1:length(stats)) {
      datatype = GDALDataType$GDT_Float64
      nodata = -9999
      if (classes[[stat_ind]] == "integer") {
        datatype = GDALDataType$GDT_Int32
        nodata = 0
      }
      rasts[[stats[[stat_ind]]]] = createDataset(raster_path = rast_paths[[stat_ind]],
                                                 nbands = 1, datatype = datatype, projstring = projstring,
                                                 lr_lat = lr_lat, ul_lat = ul_lat, ul_lon = ul_lon,
                                                 lr_lon = lr_lon, res = c(xres, yres), nodata = nodata,
                                                 co = creation_options)
    }
    xsize = rasts[[1]]$GetRasterXSize()
    ysize = rasts[[1]]$GetRasterYSize()
    bands = lapply(rasts, function(x) x[[1]])
    block_x_size = bands[[1]]$GetBlockXSize()
    block_y_size = bands[[1]]$GetBlockYSize()
    file_index = 0
    for (l2a_path in l2a_list) {
      file_index = file_index + 1
      message(sprintf("Reading file %s (%d/%d)",
                      basename(l2a_path), file_index, total_files),
              appendLF = T)
      l2a = readLevel2A(l2a_path)
      vals = getLevel2AM(l2a)
      vals= vals[,  .SD, .SDcols=cols]
      
      vals = clipLevel2AM(vals, ul_lon, lr_lon, lr_lat,
                          ul_lat)
      
      #vals = clipLevel2AMGeometry(vals,polygon_spdf, split_by= NULL)
      
      if(nrow(vals) > 0){
        cols_without_quality = c(setdiff(cols, "quality_flag"))
        vals = vals[quality_flag == 1, cols_without_quality,
                    with = FALSE]
        if (nrow(vals) == 0)
          next
        vals[, `:=`(x_ind, as.integer(vals[, floor((lon_lowestmode -
                                                      ul_lon)/xres)]))]
        vals[, `:=`(y_ind, as.integer(vals[, floor((lat_lowestmode -
                                                      ul_lat)/yres)]))]
        vals[, `:=`(inds, 1 + x_ind + y_ind * block_x_size)]
        names(vals) = gsub(metric, "x", names(vals))
        aggs = vals[, eval(call), by = list(inds, x_ind,
                                            y_ind)]
        aggs[, `:=`(c("x_block", "y_block"),
                    lapply(.SD, function(x) as.integer(floor(x/block_x_size)))),
             .SDcols = c("x_ind", "y_ind")]
        aggs[, `:=`(block_xind = x_ind - x_block *
                      block_x_size, block_yind = y_ind - y_block *
                      block_y_size)]
        aggs[, `:=`(block_inds, 1 + block_xind + block_yind *
                      block_x_size)]
        blocks = aggs[, list(vals = list(.SD)), by = list(x_block,
                                                          y_block), .SDcols = c(stats, "block_inds")]
        thisEnv = new.env()
        assign("ii", 0, thisEnv)
        total_rows = nrow(blocks)
        invisible(apply(blocks, 1, function(row) {
          ii = get("ii", thisEnv) + 1
          assign("ii", ii, thisEnv)
          message(sprintf("\rProcessing blocks...%.2f%%",
                          (100 * ii)/total_rows), appendLF = F)
          agg1 = data.table::as.data.table(lapply(bands,
                                                  function(x) x[[row$x_block, row$y_block]]))
          agg1[row$vals$block_inds] = agg_join(agg1[row$vals$block_inds],
                                               row$vals[, 1:(ncol(row$vals) - 1)])
          lapply(stats, function(x) bands[[x]][[row$x_block,
                                                row$y_block]] = agg1[[x]])
        }))
        message()
        rm(list = ls(envir = thisEnv), envir = thisEnv)
        rm(thisEnv)
        close(l2a)
      }# enf if of does not overlap
      
    }
    finalize_rasts = lapply(names(finalizer), function(x) {
      rast_name = sprintf("%s_%s_%s.tif", out_root,
                          metric, x)
      message(sprintf("Writing raster: %s", rast_name))
      rast = createDataset(raster_path = rast_name, nbands = 1,
                           datatype = GDALDataType$GDT_Float64, projstring = projstring,
                           lr_lat = lr_lat, ul_lat = ul_lat, ul_lon = ul_lon,
                           lr_lon = lr_lon, res = c(xres, yres), nodata = -9999,
                           co = creation_options)
      band = rast[[1]]
      formula = finalizer[[x]]
      formulaCalculate(formula, bands, band)
      rast$Close()
    })
    lapply(rasts, function(x) x$Close())
  }
}

#
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


agg.sum <- ~data.table::data.table(
  sum = sum(x, na.rm = T),
  n = length(x[!is.na(x)])
)


agg.join.sum <- function(x1, x2) {
  combined = data.table::data.table()
  x1$n[is.na(x1$n)] = 0
  x1$sum[is.na(x1$sum)] = 0
  
  combined$n = x1$n + x2$n
  combined$sum = x1$sum + x2$sum
  
  return (combined)
}

# default_finalizer = list()
# gedi_dir <- "Z:\\GEDI_2A_v2\\02_maps"
# outdir <- 'K:\\01_Manuscripts\\12_SIF_GEDI\\02_data\\02_GEDI'
# ## GET BBOX
# aoi <- rgdal::readOGR("K:\\01_Manuscripts\\12_SIF_GEDI\\02_data\\01_SHP\\continents.shp")
# 
# degree.to.km.factor <- 1 / 111.32
# ul_lat = aoi@bbox[2, 2]
# ul_lon = aoi@bbox[1, 1]
# lr_lat = aoi@bbox[2, 1]
# lr_lon = aoi@bbox[1, 2]
# 
# plot(aoi, axes=T)
# 
# #world <- ne_countries(scale = "medium", returnclass = "sf")
# #class(world)
# #plot(world[1]$geometry, axes=T)#, xlim=c(ul_lon,lr_lon), ylim=c(lr_lat,ul_lat))
# points(c(ul_lon,lr_lon),c(ul_lat,lr_lat), col="red", cex=2, pch=16)
# polygon(c(ul_lon,lr_lon,lr_lon,ul_lon),c(ul_lat,ul_lat,lr_lat,lr_lat))
# grid()
# 
# projstring = 'GEOGCS["WGS 84",
#     DATUM["WGS_1984",
#         SPHEROID["WGS 84",6378137,298.257223563,
#             AUTHORITY["EPSG","7030"]],
#         AUTHORITY["EPSG","6326"]],
#     PRIMEM["Greenwich",0,
#         AUTHORITY["EPSG","8901"]],
#     UNIT["degree",0.01745329251994328,
#         AUTHORITY["EPSG","9122"]],
#     AUTHORITY["EPSG","4326"]]'
# 
# def_co = c("COMPRESS=DEFLATE",
#            "BIGTIFF=YES",
#            "TILED=YES",
#            "BLOCKXSIZE=512",
#            "BLOCKYSIZE=512"
# )
# level2aRasterizeStats2(
#   l2aDir = gedi_dir,
#   metrics = c("rh75","rh80","rh85","rh90","rh95","rh96","rh97","rh98","rh99"),
#   out_root = file.path(outdir,'map_continents'),
#   ul_lat = ul_lat,
#   ul_lon = ul_lon,
#   lr_lat = lr_lat,
#   lr_lon = lr_lon,
#   res = c(0.05, -0.05),
#   creation_options = def_co,
#   #res = c(degree.to.km.factor, -degree.to.km.factor),
#   agg_function = agg.function,
#   agg_join = default_agg_join,
#   finalizer = default_finalizer,
#   polygon_spdf=aoi
# )
# 
# 
# 
# 
# 
