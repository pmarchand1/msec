# Compute MSEC human population rasters

library(raster)
library(rgeos)
source("utils.R")

#### Load input data ####

gshhg_dir = "{{Insert path to GSHHG data}}"

# SEDAC raster bricks (1990-1995 and 2000-2020)
#  'rotated' versions have (0,360) longitude range instead of (-180,180)
SEDAC_90_brick <- brick("humanpop/SEDAC_90_brick.grd")
SEDAC_90_rot_brick <- brick("humanpop/SEDAC_90_rot_brick.grd")
SEDAC_00_brick <- brick("humanpop/SEDAC_00_brick.grd")
SEDAC_00_rot_brick <- brick("humanpop/SEDAC_00_rot_brick.grd")

# Distance to land raster from GSHHG
dist_land <- raster(file.path(gshhg_dir, "dist_to_GSHHG_v2.3.4_1m.nc"))

# land_mask is at final grid resolution (2.5 arc-min.)
land_mask <- raster("reeflandarea/land_final.grd")


# Function to calculate human population count within dist of point (long, lat)
#  version is GPW version (3 or 4)
extract_pop <- function(long, lat, dist, version) {
    tryCatch({
        if (version == 3) {
            SEDAC <- SEDAC_90_brick
            SEDAC_rot <- SEDAC_90_rot_brick
            var_names = c("human_pop_90", "human_pop_95")
        } else {
            SEDAC <- SEDAC_00_brick
            SEDAC_rot <- SEDAC_00_rot_brick
            var_names = c("human_pop_00", "human_pop_05", "human_pop_10", 
                          "human_pop_15", "human_pop_20")
        }
        pt <- SpatialPoints(cbind(0, 0), 
                            proj = CRS(paste0("+proj=aeqd +lon_0=", long, 
                                              " +lat_0=", lat, " +unit=m")))
        buf <- gBuffer(pt, width = dist, quadsegs = 20)
        buf <- spTransform(buf, projection(SEDAC))
        
        # For points close to 180 meridian, used rotated (0-360) longitude
        if (abs(long - 180) < 3) {
            buf <- rotate_poly(buf)
            SEDAC_crop <- crop(SEDAC_rot, buf, snap = "out")
        } else {
            SEDAC_crop <- crop(SEDAC, buf, snap = "out")
        }
        humanpop <- suppressWarnings(
            extract(SEDAC_crop, buf, fun = sum, na.rm = TRUE)
        )
        setNames(c(long, lat, humanpop), c("long", "lat", var_names))
    }, error = function(e) {
        print(c(long, lat, e))
        setNames(c(long, lat, rep(NA, length(var_names))),
                 c("long", "lat", var_names))
    })
}


#### Compute human population within 20km radius ####

# Only keep points that could be within 20km of land (25km buffer)
land_buf25 <- calc(dist_land, function(x) x > -25,
                   filename = "humanpop/land_buf25.grd", datatype = "INT1U")

land_buf25 <- resample(land_buf25, land_mask, method = "ngb",
                       filename = "humanpop/land_buf25_resamp.grd", datatype = "INT1U")

# Remove points over land
land_buf25 <- mask(x = land_buf25, mask = land_mask, maskvalue = 1)

# Get grid cell midpoints within 25km buffer
grid_pts <- rasterToPoints(land_buf25, fun = function(x) {x == 1})
grid_pts <- as.data.frame(grid_pts)
colnames(grid_pts) <- c("long", "lat", "dist")
grid_pts$dist <- 20000

# Calculate population within 20km of each point for 1990-1995 (v3)
# NOTE: This calculation (as similar ones below) was processed in parallel on a HPC cluster.
res <- lapply(1:nrow(grid_pts), function(i) {
    extract_pop(grid_pts$long[i], grid_pts$lat[i], grid_pts$dist[i], version = 3)
})

# Convert result to SpatialPointsDataFrame and rasterize
res <- as.data.frame(do.call(rbind, res))
coordinates(res) <- ~long + lat
res <- res[!is.na(res$human_pop_90), ] # remove NA
res@data <- round(res@data, digits = 0) # round to nearest person
res_rast_90 <- rasterize(res, land_mask, field = colnames(res@data), 
                         background = 0, filename = "humanpop/humanpop_1990s_20km.grd")
# Apply land mask
res_mask_90 <- mask(res_rast_90, land_mask, maskvalue = 1, 
                    filename = "humanpop/humanpop_1990s_20km_masked.grd")

# Repeat for 2000-2020 (v4)
res <- lapply(1:nrow(grid_pts), function(i) {
    extract_pop(grid_pts$long[i], grid_pts$lat[i], grid_pts$dist[i], version = 4)
})
res <- as.data.frame(do.call(rbind, res))
coordinates(res) <- ~long + lat
res <- res[!is.na(res$human_pop_00), ]
res@data <- round(res@data, digits = 0)
res_rast_00 <- rasterize(res, land_mask, field = colnames(res@data), 
                         background = 0, filename = "humanpop/humanpop_2000s_20km.grd")
res_mask_00 <- mask(res_rast_00, land_mask, maskvalue = 1, 
                    filename = "humanpop/humanpop_2000s_20km_masked.grd")


#### Compute human population within 50km radius ####

# Create a 55km distance to land buffer
land_buf55 <- calc(dist_land, function(x) x > -55,
                   filename = "humanpop/land_buf55.grd", datatype = "INT1U")
land_buf55 <- resample(land_buf55, land_mask, method = "ngb",
                       filename = "humanpop/land_buf55_resamp.grd", datatype = "INT1U")

# Remove points over land
land_buf55 <- mask(land_buf55, land_mask, maskvalue = 1)

# Get grid points for land area calculation
grid_pts <- rasterToPoints(land_buf55, fun = function(x) {x == 1})
grid_pts <- as.data.frame(grid_pts)
colnames(grid_pts) <- c("long", "lat", "dist")
grid_pts$dist <- 50000

# Extract population for GPWv3 and v4 as above

res <- lapply(1:nrow(grid_pts), function(i) {
    extract_pop(grid_pts$long[i], grid_pts$lat[i], grid_pts$dist[i], version = 3)
})
res <- as.data.frame(do.call(rbind, res))
coordinates(res) <- ~long + lat
res <- res[!is.na(res$human_pop_90), ] # remove NA
res@data <- round(res@data, digits = 0) # round to nearest person
res_rast_90 <- rasterize(res, land_mask, field = colnames(res@data), 
                         background = 0, filename = "humanpop/humanpop_1990s_50km.grd")
res_mask_90 <- mask(res_rast_90, land_mask, maskvalue = 1, 
                    filename = "humanpop/humanpop_1990s_50km_masked.grd")

res <- lapply(1:nrow(grid_pts), function(i) {
    extract_pop(grid_pts$long[i], grid_pts$lat[i], grid_pts$dist[i], version = 4)
})
res <- as.data.frame(do.call(rbind, res))
coordinates(res) <- ~long + lat
res <- res[!is.na(res$human_pop_00), ]
res@data <- round(res@data, digits = 0)
res_rast_00 <- rasterize(res, land_mask, field = colnames(res@data), 
                         background = 0, filename = "humanpop/humanpop_2000s_50km.grd")
res_mask_00 <- mask(res_rast_00, land_mask, maskvalue = 1, 
                    filename = "humanpop/humanpop_2000s_50km_masked.grd")

