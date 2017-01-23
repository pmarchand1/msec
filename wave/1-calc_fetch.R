# Calculate fetch for all grid cells points within 50km from a coastline

library(rgdal)
library(rgeos)
library(raster)
library(waver)
source("utils.R")

gshhs_dir <- "{{Insert path to GSHHS shapefiles}}"

# Combine shorelines from GSHHS L1 (World except Antartica) and L5 (Antarctica sea-ice border)
gshhs1 <- readOGR(file.path(gshhs_dir, "f"), "GSHHS_f_L1")
gshhs5 <- readOGR(file.path(gshhs_dir, "h"), "GSHHS_h_L5")
gshhs_full <- rbind(gshhs1, gshhs5, makeUniqueIDs = TRUE)
gshhs_full <- as(gshhs_full, "SpatialPolygons")
rm(gshhs1, gshhs5)

# Get grid points to calculate fetch at

# 55km land buffer from land area calculation
land_buf55 <- raster("reeflandarea/buffers/land_buf55_resamp.grd")
# Remove land cells
land_mask <- raster("reeflandarea/land_final.grd")
land_buf55 <- mask(land_buf55, land_mask, maskvalue = 1)
grid_pts <- rasterToPoints(land_buf55, fun = function(x) {x == 1})
grid_pts <- as.data.frame(grid_pts)
grid_pts$layer <- NULL
# Need to rotate coordinates to (-180, 180) longitude range
grid_pts$x[grid_pts$x > 180] <- grid_pts$x[grid_pts$x > 180] - 360
# Convert to SpatialPointsDataFrame
coordinates(grid_pts) <- ~x + y
proj4string(grid_pts) <- CRS(proj4string(gshhs_full))


# Parameters for fetch calculation
bearings <- seq(0, 337.5, 22.5)
spread <- seq(-10, 10, 2.5)
dmax <- 50000 # Only calculate fetch up to 50km

# Find bounding box intersections between 50km rectangle around points and
#  shoreline polygons. (To speed up later calculation)
rects <- do.call(rbind, c(lapply(1:length(grid_pts),
                   function(i) get_clip_rect(grid_pts[i], dmax, FALSE)
                 ), makeUniqueIDs = TRUE))
btree <- gBinarySTRtreeQuery(gshhs_full, rects)
rm(rects)

# Function to calculate fetch for point at index "i"
#  first subsetting the shoreline layer based on btree, to save processing time
# Returns a vector with names corresponding to bearings
fetch_i <- function(i) {
    if (is.null(btree[[i]])) {
        # If no shoreline polygons around, put max fetch
        setNames(rep(dmax, length(bearings)), bearings)
    } else {
        tryCatch(
            fetch_len(grid_pts[i], bearings, gshhs_full[btree[[i]]], 
                      dmax, spread), 
            error = function(e) {
                    print(paste("Error at", i, ":", e))
                    setNames(rep(NA, length(bearings)), bearings)
            }
        )
    }
}

# NOTE: This calculation was parallelized on a HPC cluster
fetch_res <- lapply(1:length(grid_pts), fetch_i)
fetch_res <- do.call(rbind, fetch_res) # Forms n_points x n_bearings matrix

# Merge coordinates and fetch data, only keep points with at least one fetch < 15km
fetch_res <- SpatialPointsDataFrame(grid_pts, fetch_res)
fetch_res <- fetch_res[which(apply(fetch_res@data, 1, function(x) any(x < 50000))), ]
# Rotate coordinates back to (0, 360) latitude range, and save
fetch_res <- rotate_pts(fetch_res)
saveRDS(fetch_res, "wave/fetch_res.RData")
