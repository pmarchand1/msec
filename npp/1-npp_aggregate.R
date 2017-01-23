# Calculate temporal aggregates of Net Primary Productivity
#  based on the 8-day layers downloaded from Coast Watch (July 2002 to Oct. 2013)

library(raster)
library(lubridate)
library(stringr)

prod_dir <- "{{Insert path to Coast Watch Net Primary Producitivity layers}}"

# Set rasterOptions to limit amount of cells loaded into memory at once
rasterOptions(chunksize = 1E6, maxmemory = 1E6)

# Combine all individual rasters into one RasterBrick (506 layers)
prod_layers <- as.list(dir(prod_dir, pattern = "prod_20*", full.names = TRUE))
prod_stack <- stack(prod_layers)
prod_brick <- brick(prod_stack, bandorder = "BIP",
                    filename = file.path(prod_dir, "prod_brick.grd"))

# Get vector of dates corresponding to each layer
prod_times <- sapply(prod_layers, function(x) getZ(raster(x)))
prod_dates <- as.Date(as.POSIXct(prod_times, origin = "1970-01-01"))


# Calculate annual means and inter-annual standard deviation for each cell
#  Note: limit to layers 1:495 to have exactly 11 years (July to June) of data 
prod_brick_sub <- subset(prod_brick, 1:495)
prod_year_means <- stackApply(prod_brick_sub, rep(1:11, each = 45),
                              fun = mean, na.rm = TRUE,
                              filename = "npp/aggregates/prod_year_means.grd")
prod_inter_sd <- calc(prod_year_means, function(x) sd(x, na.rm = TRUE),
                      filename = "npp/aggregates/prod_interann_sd.grd")


# Turns dates into intra-annual sampling point from 1 (Jan.5) to 45 (Dec.22-23)
day_idx <- as.integer(as.factor(yday(prod_dates)))

# Get mean NPP (over years) and count (non-NA) data values for each cell
#  and sampling day (i.e. 45-layer bricks)
prod_day_means <- stackApply(prod_brick, day_idx, fun = mean, na.rm = TRUE,
                             filename = "npp/aggregates/prod_day_means.grd")

prod_day_counts <- stackApply(prod_brick, day_idx,
                              fun = function(x, ...) sum(!is.na(x)),
                              filename = "npp/aggregates/prod_day_counts.grd")

rasterOptions(maxmemory = 1E7)

# Compute intra-annual mean, min, max and standard deviation per cell
prod_mean <- calc(prod_day_means, function(x) mean(x, na.rm = TRUE),
                  filename = "npp/aggregates/prod_mean.grd")

prod_min <- calc(prod_day_means, function(x) min(x, na.rm = TRUE),
                 filename = "npp/aggregates/prod_min.grd")

prod_max <- calc(prod_day_means, function(x) max(x, na.rm = TRUE),
                 filename = "npp/aggregates/prod_max.grd")

prod_sd <- calc(prod_day_means, function(x) sd(x, na.rm = TRUE),
                filename = "npp/aggregates/prod_sd.grd")
