# This script corrects an aggregate NPP layer
#  by interpolating shallow cells (depth < 30m) based on values
#  of nearest cells within a given distance, except that cells within 9.3km
#  of coast or reef are preferably interpolated using other coastal/reef cells.

# Output includes a final data raster and a 'flags' raster 
#   (0 = non-interpolated, 1 = interpolated with coastal/reef cells,
#    2 = interpolated with other cells)

library(raster)
source("utils.R")

stat <- "mean" # which aggregate layer to process 
               # options: mean, min, max, sd, interann_sd
interp_d <- 125000 # max. distance for interpolation (in meters)
ninterp_max <- 3 # max. number of cells to use for interpolation

# Read in the depth and interpolation (coast/reef) masks
depth30 <- raster("npp/masks/depth30_final.grd")
interp <- raster("npp/masks/interp_final.grd")

# Read productivity layer and crop out last (duplicate) column
prod_lyr <- raster(paste0("npp/aggregates/prod_", stat, ".grd"))
prod_lyr <- crop(prod_lyr, depth30)

# Filter out shallow cells and land cells
prod_lyr[depth30 == 1 | is.na(depth30)] <- NA
# We want to interpolate missing values that aren't land
to_interp <- is.na(prod_lyr) & !is.na(depth30)

# Crop out extreme north (>84.5, since all cells with lat >83.22 are NA)
to_interp <- crop(to_interp, extent(c(xmin(to_interp), xmax(to_interp), 
                                      ymin(to_interp), 84.5)))

# Get points to interpolate in 2-column matrix format
pts_interp <- rasterToPoints(to_interp, function(x) x == 1)
pts_interp <- pts_interp[,1:2]
coast_reef <- as.logical(raster::extract(interp, pts_interp))
pts_interp <- cbind(as.data.frame(pts_interp), coast_reef)

# Now remove shallow cells from interp
interp <- interp & (depth30 == 0)

# For a given point (x, y), this function returns
#  a vector (interpolated value, flag)
#  coast_reef is TRUE if point is near coast/reef
get_interp <- function(x, y, coast_reef) {
    # Clip prod_lyr to rectangular buffer of at least interp_d around p
    p <- c(x, y)
    prod_crop <- clip_raster_pt(prod_lyr, p, interp_d)
    # Distance of each raster cell to p
    dists <- distanceFromPoints(prod_crop, p)
    # Extract values as vectors
    prod_val <- getValues(prod_crop)
    dist_val <- getValues(dists)
    # Cells further than interp_d are not useable
    prod_val[dist_val > interp_d] <- NA
    
    # No data available within interp_d
    if (all(is.na(prod_val))) return(c(NA, NA))
    
    # If point near coast/reef, prioritize using cells in interpolation mask, 
    # but if none, use any cell within interp_d 
    if (coast_reef) {
        interp_crop <- clip_raster_pt(interp, p, interp_d)
        interp_val <- getValues(interp_crop)
        interp_idx <- which(interp_val == 1 & !is.na(prod_val))
        flag <- 1
    }
    if (!coast_reef || length(interp_idx) == 0) {
        interp_idx <- which(!is.na(prod_val))
        flag <- 2 # interpolated with other cells
    } 
    # If more than ninterp_max cells available, take closest
    if (length(interp_idx) > ninterp_max) {
        interp_idx <- interp_idx[order(dist_val[interp_idx])][1:ninterp_max]
    } 
    
    return(c(mean(prod_val[interp_idx]), flag))
}

# Get interpolated values (column 1) and flags (column 2) for all points
# NOTE: This was computed in parallel on a HPC cluster
interp_values <- Map(get_interp, pts_interp$x, pts_interp$y,
                     pts_interp$coast_reef)
interp_values <- do.call(rbind, interp_values)

# Filter out NA values
has_val <- which(!is.na(interp_values[, 1]))
pts_interp <- pts_interp[has_val, ]
interp_values <- interp_values[has_val, ]

# Rasterize interpolated values and merge back into prod_lyr
prod_interp <- rasterize(pts_interp[, 1:2], prod_lyr, field = interp_values[, 1], 
    update = TRUE, filename = paste0("npp/aggregates/prod_", stat, "_final.grd"))

# Rasterize flags, set 0 for non-interpolated
prod_flags <- rasterize(pts_interp[, 1:2], prod_lyr, 
                        field = interp_values[, 2], background = 0)
# Flag is NA for all cells where data is NA
prod_flags[is.na(prod_interp)] <- NA

writeRaster(prod_flags, paste0("npp/aggregates/prod_", stat, "_flags.grd"))
