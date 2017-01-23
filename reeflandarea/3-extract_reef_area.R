# Produce rasters of reef area within 15km and 200km of marine grid cells

library(raster)
library(rgeos)
source("utils.R")

reef_dir <- "{{Insert path to Reefs at Risk raster}}"

# Load coral reefs layer 
reefs <- raster(file.path(reef_dir, "reef_500"))
projection(reefs) <- "+proj=cea +lon_0=-160 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

# Load land mask
land <- raster("land_final.grd")

# Function to calculate reef area within dist of point (long, lat)
reef_area <- function(long, lat, dist) {
    tryCatch({
        # Create a circular buffer in a equidistant projection centered on point
        pt <- SpatialPoints(cbind(0, 0), 
                            proj = CRS(paste0("+proj=aeqd +lon_0=", long, 
                                              " +lat_0=", lat, " +unit=m")))
        buf <- gBuffer(pt, width = dist, quadsegs = 20)
        buf <- spTransform(buf, projection(reefs))
        # Compute reef area in buffer (# cells x 0.25km^2 per cell)
        reef_crop <- crop(reefs, buf, snap = "out")
        cell_area <- 0.25
        rarea <- suppressWarnings(
            extract(reef_crop, buf, fun = sum, na.rm = TRUE) * cell_area)
        c(long = long, lat = lat, reef_area = rarea)
    }, error = function(e) {
        print(c(long, lat, e))
        c(long = long, lat = lat, reef_area = NA)
    })
}


#### Compute reef area within 15km radius ####

# 20km buffer around reef areas (pre-computed in ArcGIS)
reefs_buf20 <- raster("reeflandarea/buffers/reef_20km_rast.tif")
reefs_buf20 <- inv_rotate(reefs_buf20)

# Resample reef buffer to final grid
land_crop <- crop(land, reefs_buf20, snap = "out")
reefs_buf20 <- resample(reefs_buf20, land_crop, method = "ngb")

# Remove points over land
reefs_buf20 <- mask(reefs_buf20, land_crop, maskvalue = 1)

# Get grid points for reef area calculation
grid_pts <- rasterToPoints(reefs_buf20)
grid_pts <- as.data.frame(grid_pts)
colnames(grid_pts) <- c("long", "lat", "dist")
grid_pts$dist <- 15000

# Calculate reef area within 15km of each point
# NOTE: This calculation (and the one for 200km below) was processed in parallel on a HPC cluster.
res <- Map(reef_area, grid_pts$long, grid_pts$lat, grid_pts$dist)

# Convert result to SpatialPointsDataFrame and rasterize
res <- as.data.frame(do.call(rbind, res))
coordinates(res) <- ~long + lat
res_rast <- rasterize(res, land, field = "reef_area", background = 0,
                      filename = "reef_area_15km.grd")
res_mask <- mask(res_rast, land, maskvalue = 1, 
                 filename = "reef_area_15km_masked.grd")

#### Repeat for 200km radius ####

# Load a 205km buffer around reef areas
reefs_buf200 <- raster("reeflandarea/buffers/reef_205kmbuff_rast.tif")
reefs_buf200 <- inv_rotate(reefs_buf200)

land_crop <- crop(land, reefs_buf200, snap = "out")
reefs_buf200 <- resample(reefs_buf200, land_crop, method = "ngb")
reefs_buf200 <- mask(reefs_buf200, land_crop, maskvalue = 1)

grid_pts <- rasterToPoints(reefs_buf200)
grid_pts <- as.data.frame(grid_pts)
colnames(grid_pts) <- c("long", "lat", "dist")
grid_pts$dist <- 200000

res <- Map(reef_area, grid_pts$long, grid_pts$lat, grid_pts$dist)

# Convert result to SpatialPointsDataFrame and rasterize
res <- as.data.frame(do.call(rbind, res))
coordinates(res) <- ~long + lat
res_rast <- rasterize(res, land, field = "reef_area", background = 0,
                      filename = "reef_area_200km.grd")
res_mask <- mask(res_rast, land, maskvalue = 1, 
                 filename = "reef_area_200km_masked.grd")
