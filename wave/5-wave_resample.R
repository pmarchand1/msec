# Resample wave energy statistics from each WAVEWATCH III (WW3) grid to
# our final grid, and pick the highest-resolution value for each point

library(raster)
source("utils.R")

# WW3 grids from high to low resolution
grid_names <- c("ak_4m", "ecg_4m", "nsb_4m", "oz_4m", "wc_4m", 
                "ak_10m", "ecg_10m", "med_10m", "nsb_10m", "nwio_10m", 
                "oz_10m", "pi_10m", "wc_10m", "glo_30m")

# Resolution codes (1 = 30m, 2 = 10m, 3 = 4m)
grid_res <- c(rep(3, 5), rep(2, 8), 1)

# Land mask (and template for resampling)
land_mask <- raster("reeflandarea/land_final.grd")

# Combine grids for a certain statistic (mean, sd, interann_sd)
#  First resample all to final grid, then merge to keep first non-NA for each cell
combine_grids <- function(stat) {
    rasts <- lapply(paste0("wave/we_agg/", grid_names, "_", stat, ".grd"), raster)
    # Some grids have min. longitude < 0, convert to (0, 360) longitude range
    to_rotate <- which(vapply(rasts, xmin, 0) < -0.25)
    rasts[to_rotate] <- lapply(rasts[to_rotate], extend_inv_rotate)
    rasts_resamp <- lapply(rasts, function(r) resample(r, land_mask))
    merged_rast <- do.call(merge, rasts_resamp)
    mask(merged_rast, land_mask, maskvalue = 1,
         filename = paste0("wave/we_resamp/", stat, "_resamp.grd"))
    # For mean, also save original WW3 resolution for each cell
    if (stat == "mean") {
        # has_vals will be equal to resolution code (1,2,3) for cells with values, 0 otherwise
        has_vals <- lapply(rasts_resamp, function(r) !is.na(r))
        for (i in seq_along(grid_res)) {
            has_vals[[i]] <- has_vals[[i]] * grid_res[[i]]
        }
        # Aggregate with maximum to get highest resolution at grid cell
        ww3_res <- do.call(mosaic, c(has_vals, fun = max))
        ww3_res <- mask(ww3_res, land_mask, maskvalue = 1, datatype = "INT1U", 
                        filename = "ww3_resolution_final.grd")
    }
}

# Apply function above to each summary statistic
mean_resamp <- combine_grids("mean")
sd_resamp <- combine_grids("sd")
interann_sd_resamp <- combine_grids("interann_sd")


#### Incorporate sheltered points ####

# Add statistics to SpatialPoints layer, remove points with no data
fetch_pts <- readRDS("wave/fetch_sheltered.RData")
fetch_pts$mean <- readRDS("wave/we_agg/fetch_mean.RData")
fetch_pts$sd <- readRDS("wave/we_agg/fetch_sd.RData")
fetch_pts$interann_sd <- readRDS("wave/we_agg/fetch_interann_sd.RData")
fetch_pts <- fetch_pts[!is.na(fetch_pts$mean), ]

# Create mask indicating at which points wave energy was calculated from wind and fetch
fetch_mask <- rasterize(fetch_pts, mean_resamp, field = 1, background = 0)

# Input calculated values into global rasters
mean_final <- rasterize(fetch_pts, mean_resamp, field = "mean", update = TRUE)
sd_final <- rasterize(fetch_pts, sd_resamp, field = "sd", update = TRUE)
interann_sd_final <- rasterize(fetch_pts, interann_sd_resamp, 
                               field = "interann_sd", update = TRUE)

# Apply land mask to all final layers and save
fetch_mask <- mask(fetch_mask, land_mask, maskvalue = 1,
                   filename = "wave/wind_fetch_mask.grd")
mean_final <- mask(mean_final, land_mask, maskvalue = 1, 
                   filename = "wave/wave_mean_final.grd")
sd_final <- mask(sd_final, land_mask, maskvalue = 1, 
                 filename = "wave/wave_sd_final.grd")
interann_sd_final <- mask(interann_sd_final, land_mask, maskvalue = 1, 
                          filename = "wave/wave_interann_sd_final.grd")
