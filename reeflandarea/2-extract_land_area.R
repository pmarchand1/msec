# Produce rasters of land area within 15km and 50km of marine grid cells

library(raster)
library(rgeos)
source("utils.R")

gshhg_dir = "{{Insert path to GSHHG data}}"

#### Load input data ####

# High-resolution (15 arc-sec.) land rasters
#  (longitude in (-180, 180) for gshhs and (0, 360) for gshhs_rotate)
gshhs <- raster("reeflandarea/gshhs_rast.grd")
gshhs_rotate <- raster("reeflandarea/gshhs_rotated.grd")

# land_mask is at final grid resolution (2.5 arc-min.)
land_mask <- raster("reeflandarea/land_final.grd")

# Distance to shore raster from GSHHG
dist_land <- raster(file.path(gshhg_dir, "dist_to_GSHHG_v2.3.4_1m.nc"))


# Function to calculate land area within dist of point (long, lat)
land_area <- function(long, lat, dist) {
    tryCatch({
        # Create a circular buffer in a equidistant projection centered on point
        pt <- SpatialPoints(cbind(0, 0), 
                            proj = CRS(paste0("+proj=aeqd +lon_0=", long, 
                                              " +lat_0=", lat, " +unit=m")))
        buf <- gBuffer(pt, width = dist, quadsegs = 20)
        buf <- spTransform(buf, projection(gshhs))
        
        # For points close to 180 meridian, used rotated (0, 360) longitude
        if (abs(long - 180) < 3) {
            buf <- rotate_poly(buf)
            land_crop <- crop(gshhs_rotate, buf, snap = "out")
        } else {
            land_crop <- crop(gshhs, buf, snap = "out")
        }
        # Calculate cell areas within cropped section, mask those not on land,
        #  and get total area within buffer
        cell_areas <- area(land_crop)
        area_mask <- mask(cell_areas, land_crop, maskvalue = 0)
        larea <- suppressWarnings(
            extract(area_mask, buf, fun = sum, na.rm = TRUE)
        )
        c(long = long, lat = lat, land_area = larea)
    }, error = function(e) {
        print(c(long, lat, e))
        c(long = long, lat = lat, land_area = NA)
    })
}


#### Compute land area within 15km radius ####

# Only keep points within a 20km buffer around land
land_buf20 <- calc(dist_land, function(x) x > -20, datatype = "INT1U",
                   filename = "reeflandarea/buffers/land_buf20.grd")
land_buf20 <- resample(land_buf20, land_mask, method = "ngb", datatype = "INT1U",
                       filename = "reeflandarea/buffers/land_buf20_resamp.grd")

# Remove points over land
land_buf20 <- mask(land_buf20, land_mask, maskvalue = 1)

# Get grid points for land area calculation
grid_pts <- rasterToPoints(land_buf20, fun = function(x) {x == 1})
grid_pts <- as.data.frame(grid_pts)
colnames(grid_pts) <- c("long", "lat", "dist")
grid_pts$dist <- 15000

# Calculate land area within 15km of each point
# NOTE: This calculation (and the one for 50km below) was processed in parallel on a HPC cluster.
res <- Map(land_area, grid_pts$long, grid_pts$lat, grid_pts$dist)

# Convert result to SpatialPointsDataFrame and rasterize
res <- as.data.frame(do.call(rbind, res))
coordinates(res) <- ~long + lat
res_rast <- rasterize(res, land_mask, field = "land_area", background = 0,
                      filename = "reeflandarea/land_area_15km.grd")
res_mask <- mask(res_rast, land_mask, maskvalue = 1, 
                 filename = "reeflandarea/land_area_15km_masked.grd")


#### Repeat for 50km radius ####

# Keep points within a 55km buffer
land_buf55 <- calc(dist_land, function(x) x > -55, datatype = "INT1U",
                   filename = "reeflandarea/buffers/land_buf55.grd")
land_buf55 <- resample(land_buf55, land_mask, method = "ngb", datatype = "INT1U",
                       filename = "reeflandarea/buffers/land_buf55_resamp.grd")

# Remove points over land
land_buf55 <- mask(land_buf55, land_mask, maskvalue = 1)

# Get grid points for land area calculation
grid_pts <- rasterToPoints(land_buf55, fun = function(x) {x == 1})
grid_pts <- as.data.frame(grid_pts)
colnames(grid_pts) <- c("long", "lat", "dist")
grid_pts$dist <- 50000

# Compute land area for all points, then rasterize
res <- Map(land_area, grid_pts$long, grid_pts$lat, grid_pts$dist)
res <- as.data.frame(do.call(rbind, res))
coordinates(res) <- ~long + lat

res_rast <- rasterize(res, land_mask, field = "land_area", background = 0,
                      filename = "reeflandarea/land_area_50km.grd")
res_mask <- mask(res_rast, land_mask, maskvalue = 1, 
                 filename = "reeflandarea/land_area_50km_masked.grd")

