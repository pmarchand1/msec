# Prepare mask layers for filtering and interpolation 
#  of Coast Watch Net Primary Productivity data

library(raster)
library(rgdal)
library(rgeos)
source("utils.R")

# Data directories
prod_dir <- "{{Insert path to Coast Watch NPP data folder}}"
mars_dir <- "{{Insert path to MARSPEC data folder}}"
reef_dir <- "{{Insert path to reef location data folder}}"


#### Pre-process original datasets ####

# Load a NPP layer
prod_lyr <- raster(file.path(prod_dir, "prod_20101207.nc"))

# Remove rightmost column of prod_lyr (all NA)
prod_ext <- extent(prod_lyr)
xmax(prod_ext) <- xmax(prod_ext) - res(prod_lyr)[1]
prod_lyr <- crop(prod_lyr, prod_ext)

# Load 30-arcsecond bathymetry and distance to shore layers from MARSPEC
bathymetry <- raster(file.path(mars_dir, "bathymetry_30s/bathy_30s.tif"))
shoredist <- raster(file.path(mars_dir, "biogeo01_07_30s/shoredist_30s.tif"))

# Load the coral reef buffer shapefile
# NOTE: This file was generated in ArcGIS by creating a 9.3km buffer around the
# the Reefs at Risk maps of reef locations, then reprojecteing to long/lat coords. 
# centered on the 180th meridian
reef_buffer_shp <- readOGR(file.path(reef_dir, "reef_500_buffer_longlat180.shp"),
                           layer = "reef_500_buffer_longlat180")

# Convert reef buffer to raster at same resolution as shoredist
reef_buffer <- raster(nrows = nrow(shoredist), ncols = ncol(shoredist), 
                      ext = extent(shoredist))
projection(reef_buffer) <- projection(shoredist)
reef_buffer <- rasterize(reef_buffer_shp, reef_buffer, field = 1, background = 0, 
                         filename = "npp/intermediate/reef_buffer.tif")
extent(reef_buffer) <- extent(0, 360, -90, 90)

# Rotate MARSPEC layers from (-180, 180) to (0, 360) longitude range
bathymetry <- inv_rotate(bathymetry)
shoredist <- inv_rotate(shoredist)


#### Derive mask layers ####

# depth30: is depth under 30m? (1=YES, 0=NO, NA=land)
# interpolation layer: depth > 30m and distance to shore or reef <= 9.3km
depth30 <- bathymetry > -30
writeRaster(depth30, "npp/intermediate/depth30.grd")
interp_lyr <- depth30 == 0 & ((shoredist <= 9.3) | reef_buffer)
writeRaster(interp_lyr, "npp/intermediate/interp_lyr.grd")


#### Convert depth mask to prod_lyr extent and resolution ####

#   1) Double cell resolution (from about 5x to 10x that of prod_lyr)
depth30_disag <- disaggregate(depth30, fact = 2, 
                              filename = "npp/intermediate/depth30_disag.grd")

#   2) 'Wrap' 5 cells from right edge to left edge to match prod_lyr extent
extent_left <- extent(360 - 5 * res(depth30_disag)[1], 360, -90, 90)
extent_right <- extent(0, 360 - 5 * res(depth30_disag)[1], -90, 90)
depth30_left <- crop(depth30_disag, extent_left)
xmin(depth30_left) <- xmin(depth30_left) - 360
xmax(depth30_left) <- 0
depth30_right <- crop(depth30_disag, extent_right)
depth30_rot <- merge(depth30_left, depth30_right,
                     filename = "npp/intermediate/depth30_rot.grd")

#   3) resample to exactly 10x prod. layer resolution
depth30_resamp <- raster(nrows = nrow(prod_lyr) * 10, 
                         ncol = ncol(prod_lyr) * 10, ext = extent(prod_lyr))
projection(depth30_resamp) <- projection(prod_lyr)
depth30_resamp <- resample(depth30_rot, depth30_resamp, method = "ngb", 
                           filename = "npp/intermediate/depth30_resamp.grd")

#   4) aggregate with max (large cell is '1' if there's any '1' in small cells)
#       and na.rm=TRUE (large cell is only NA if all small cells are NA)
depth30_agg <- aggregate(depth30_resamp, fact = 10, fun = max, na.rm = TRUE,
                         filename = "npp/intermediate/depth30_agg.grd")

#   5) Apply more precise land mask

# Load land_final, recode to NA over land, 1 over water
land_final <- raster("reeflandarea/land_final.grd")
land_final[land_final == 1] <- NA
land_final[land_final == 0] <- 1 

# Merge with depth30_agg to replace extra land cells (NA in depth30_agg,
#  but not in land_final) with shallow ocean (depth30 = 1)
depth30_final <- merge(depth30_agg, land_final,
                       filename = "npp/masks/depth30_final.grd")

#### Repeat steps for interpolation mask ####

interp_disag <- disaggregate(interp_lyr, fact = 2,
                             filename = "npp/intermediate/interp_disag.grd")
extent_left <- extent(360 - 5 * res(interp_disag)[1], 360, -90, 90)
extent_right <- extent(0, 360 - 5 * res(interp_disag)[1], -90, 90)
interp_left <- crop(interp_disag, extent_left)
xmin(interp_left) <- xmin(interp_left) - 360
xmax(interp_left) <- 0
interp_right <- crop(interp_disag, extent_right)
interp_rot <- merge(interp_left, interp_right, 
                    filename =  "npp/intermediate/interp_rot.grd")
interp_resamp <- raster(nrows = nrow(prod_lyr) * 10, 
                        ncol = ncol(prod_lyr) * 10, ext = extent(prod_lyr))
projection(interp_resamp) <- projection(prod_lyr)
interp_resamp <- resample(interp_rot, interp_resamp, method = "ngb", 
                          filename = "npp/intermediate/interp_resamp.grd")
interp_final <- aggregate(interp_resamp, fact=10, fun=max, na.rm=TRUE,
                          filename = "npp/intermediate/interp_final.grd")

