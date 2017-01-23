# Compute summary statistics from daily wave energy values

library(raster)
library(dplyr)
library(tidyr)

# Get filelist and parse into (grid, year, month) table
we_files <- data.frame(name = dir("wave/we_rast", pattern = "grd"), 
                       stringsAsFactors = FALSE) %>%
    tidyr::extract(name, c("grid", "year", "month"), 
                   "(.*)_([[:digit:]]{4})([[:digit:]]{2})", remove = FALSE)

# Function to get daily means (across years) and annual means (full years only)
#  Save output rasters in we_agg folder
calc_day_means <- function(grid, month) {
    filelist <- dir("wave/we_rast", full.names = TRUE, 
                    pattern = paste0(grid, ".*", month, "\\.grd"))
    month_stack <- stack(filelist)
    nyear <- length(filelist)
    day_means <- stackApply(month_stack, 
                    indices = rep(1:(nlayers(month_stack)/nyear), nyear),
                    fun = mean, na.rm = TRUE,
                    filename = paste0("wave/we_agg/", grid, "_", month, ".grd"))
    return(0)
}

calc_year_means <- function(grid, year) {
    filelist = dir("wave/we_rast", full.names = TRUE,
                   pattern = paste0(grid, "_", year, ".*grd"))
    if (length(filelist) < 12) return(NA) # incomplete year
    year_stack <- stack(filelist)
    year_mean <- calc(year_stack, fun = mean, na.rm = TRUE,
                      filename = paste0("wave/we_agg/", grid, "_", year, ".grd"))
    return(0)
}

# Apply functions above across grid-month and grid-year combinations
#  Note: This was performed in parallel on a HPC cluster
grid_months <- unique(we_files[, c("grid", "month")])
Map(calc_day_means, we_files$grid, we_files$month)

grid_years <- unique(we_files[, c("grid", "year")])
Map(calc_year_means, we_files$grid, we_files$year)


# Get daily and yearly means for fetch-limited (sheltered) points

for (month in unique(we_files$month)) {
    filelist <- dir("wave/we_rast", full.names = TRUE, 
                    pattern = paste0("fetch.*", month, "\\.RData"))
    fetch <- do.call(cbind, lapply(filelist, readRDS))
    # Reshape into 3D (point, day, year) array and take average by point and day
    dim(fetch) <- c(nrow(fetch), ncol(fetch)/length(filelist), length(filelist))
    day_means <- apply(fetch, 1:2, mean, na.rm = TRUE)
    saveRDS(day_means, paste0("wave/we_agg/fetch_", month, ".RData"))
}

for (year in unique(we_files$year)) {
    filelist <- dir("wave/we_rast", full.names = TRUE, 
                    pattern = paste0("fetch_", year))
    fetch <- do.call(cbind, lapply(filelist, readRDS))
    year_mean <- rowMeans(fetch, na.rm = FALSE)
    saveRDS(year_mean, paste0("wave/we_agg/fetch_", year, ".RData"))
}


# Calculate global mean, intra-annual and inter-annual std. dev.

calc_global_stats <- function(grid) {
    months <- unique(we_files$month[we_files$grid == grid])
    month_stack <- stack(paste0("wave/we_agg/", grid, "_", months, ".grd"))
    mean_rast <- calc(month_stack, fun = mean, na.rm = TRUE,
                      filename = paste0("wave/we_agg/", grid, "_mean.grd"))
    sd_rast <- calc(month_stack, fun = sd, na.rm = TRUE,
                    filename = paste0("wave/we_agg/", grid, "_sd.grd"))
    year_files <- dir("wave/we_agg", full.names = TRUE,
                      pattern = paste0(grid, "_[[:digit:]]{4}\\.grd"))
    year_stack <- stack(year_files)
    interann_sd_rast <- calc(year_stack, fun = sd, na.rm = TRUE,
                filename = paste0("wave/we_agg/", grid, "_interann_sd.grd"))
    return(0)
}

lapply(unique(we_files$grid), calc_global_stats)

# Calculate same statistics for fetch-limited points

fetch_month_files <- dir("wave/we_agg", full.names = TRUE,
                         pattern = "fetch_[[:digit:]]{2}.RData")
fetch_months <- do.call(cbind, lapply(fetch_month_files, readRDS))
saveRDS(rowMeans(fetch_months, na.rm = TRUE), "wave/we_agg/fetch_mean.RData")
saveRDS(apply(fetch_months, 1, sd, na.rm = TRUE), "wave/we_agg/fetch_sd.RData")
fetch_year_files <- dir("wave/we_agg", full.names = TRUE,
                        pattern = "fetch_[[:digit:]]{4}.RData")
fetch_years <- do.call(cbind, lapply(fetch_year_files, readRDS))
saveRDS(apply(fetch_years, 1, sd, na.rm = TRUE), "wave/we_agg/fetch_interann_sd.RData")


