# Process WAVEWATCH III (WW3) reanalysis data to calculate daily values of wave energy

library(raster)
library(dplyr)
library(tidyr)
library(lubridate)
library(waver)

# Days per month (leaving out Feb.29 for consistency across years)
month_days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# Time resolution of original data
steps_by_day <- 8 

# Get data for sheltered points
fetch_pts <- readRDS("wave/fetch_sheltered.RData")
bearings <- as.numeric(setdiff(colnames(fetch_pts@data), c("grid", "depth")))
# Duplicate minimum bearing (usually 0) to make sure values near 360 match to it
bearings <- c(bearings, 360 + min(bearings))
fetch_mat <- as.matrix(fetch_pts@data[, colnames(fetch_pts@data) %in% bearings])
fetch_grid <- fetch_pts$grid
fetch_depth <- fetch_pts$depth

# This function takes a year-month string of the format "MMMMYY", downloads the
#  corresponding WW3 files and computes the daily wave energy for each grid
#  as well as the sheltered points, with output rasters saved to "we_rast" folder
get_we_month <- function(month) {
    # Download a month of data from NOAA
    tmp_dir <- file.path("wave/tmp", month)
    ndays <- month_days[as.numeric(substr(month, 5, 6))]
    src_dir <- paste0("http://polar.ncep.noaa.gov/waves/nopp-phase1/", month, "/grib")
    system(paste0("wget -q -r -nd -np -P '", tmp_dir,"' -A '*grb2.gz' ", src_dir))
    system(paste0("gunzip ", tmp_dir, "/*.gz"))
    
    # Parse source file names in a table, keep grid code, variable code and month
    src_files <- data.frame(name = dir(tmp_dir, pattern = "grb2", full.names = TRUE),
                            stringsAsFactors = FALSE)
    src_files <- separate(src_files, name, c("type", "grid", "var", "month", "ext"), 
                          sep = "\\.", remove = FALSE) %>%
        dplyr::select(-type, -month, -ext)
    grid_codes <- unique(src_files$grid)
    
    # Initialize wave energy matrix for sheltered points
    we_fetch <- matrix(0, nrow = nrow(fetch_pts), ncol = ndays)

    for (grid_cd in grid_codes) {
        # Get subset of files matching grid
        grid_files <- filter(src_files, grid == grid_cd)

        # Part 1: Oceanic wave energy

        # Load height and period raster bricks
        height <- brick(grid_files$name[grid_files$var == "hs"])
        period <- brick(grid_files$name[grid_files$var == "tp"])
        # Check for errors (values < 0 or > 50) and if any, remove whole time step
        layer_errors <- which(cellStats(height, min) < 0 | cellStats(height, max) > 50 |
                              cellStats(period, min) < 0 | cellStats(period, max) > 50)
        if (length(layer_errors) > 0) 
            print(paste0(month, "_", grid_cd, ": error at ", 
                         paste(layer_errors, collapse = ",")))
        for (i in seq_along(layer_errors)) {
            height[[i]] <- NA
            period[[i]] <- NA
        }
        # Calculate wave energy raster brick and average by day
        #  (subset to discard layers above ndays * 8)
        we_tmp <- overlay(height, period, fun = wave_energy, 
                    filename = paste0(tmp_dir, "/", grid_cd, "_", month, "_tmp.grd"))
        we_rast <- stackApply(subset(we_tmp, 1:(ndays * steps_by_day)), 
                    rep(1:ndays, each = steps_by_day), fun = mean, na.rm = TRUE,
                    filename = paste0("wave/we_rast/", grid_cd, "_", month, ".grd"))
        
        # Part 2: Local wave energy for sheltered points (based on wind, fetch)

        wind <- brick(grid_files$name[grid_files$var == "wind"])
        # Extract u and v components for points of interest
        fetch_idx <- which(fetch_grid == grid_cd)
        u_idx <- 1:(ndays * steps_by_day)
        v_idx <- nlayers(wind)/2 + 1:(ndays * steps_by_day) 
        wu <- raster::extract(subset(wind, u_idx), fetch_pts[fetch_idx, ], 
                              method = "bilinear")
        wv <- raster::extract(subset(wind, v_idx), fetch_pts[fetch_idx, ], 
                              method = "bilinear")    
        # Compute wind speed, and bearing where wind is *from* (in degrees, 0 is N, 90 is E)
        wspd <- sqrt(wu^2 + wv^2)
        wdir <- atan2(-wu, -wv) * 180/pi
        wdir[which(wdir < 0)] <- wdir[which(wdir < 0)] + 360
        rm(wu, wv)
        
        # Calculate wind energy for all points and time steps, then average by day
        bearing_idx <- matrix(NA_integer_, nrow = nrow(wdir), ncol = ncol(wdir))
        # Find fetch bearing nearest to each wind direction and get corresponding fetch
        for (i in 1:nrow(wdir)) {
            for (j in 1:ncol(wdir)) {
                if (!is.na(wdir[i, j])) {
                    bearing_idx[i, j] <- which.min(abs(wdir[i, j] - bearings))
                }
            }
        }
        bearing_idx[bearing_idx == length(bearings)] <- which.min(bearings)
        fetch <- t(vapply(seq_along(fetch_idx), function(i) {
            fetch_mat[fetch_idx[i], bearing_idx[i, ]]
        }, rep(0, ncol(bearing_idx))))
        rm(bearing_idx)
        # Calculate wave energy and average per day
        we <- wave_energy(wind = wspd, fetch = fetch, depth = fetch_depth[fetch_idx])
        we_fetch[fetch_idx, ] <- t(apply(we, 1, function(v) {
            tapply(v, rep(1:ndays, each = steps_by_day), mean, na.rm = TRUE)    
        }))
    }
    # Save output of all sheltered points for month
    saveRDS(we_fetch, paste0("wave/we_rast/fetch_", month, ".RData"))
    unlink(paste0(tmp_dir, "/*"))
}


# Apply function above to all 372 months in dataset
#  Note: This computation was parallelized on a HPC cluster
year_list <- 1979:2009
month_list <- c("01", "02", "03", "04", "05", "06",
                "07", "08", "09", "10", "11", "12")
month_list <- as.vector(outer(year_list, month_list, paste0))

lapply(month_list, get_we_month)
