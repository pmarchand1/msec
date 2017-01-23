# Determine which points are sheltered based on fetch values 
#  and WAVEWATCH III (WW3) grid resolution

library(raster)
source("utils.R")

mars_dir <- "{{Insert path to MARSPEC data}}"

# Download first month of WW3 reanalysis data
month <- "197901"
src_dir <- paste0("http://polar.ncep.noaa.gov/waves/nopp-phase1/", month, "/grib")
system(paste("wget -q -r -nd -np -P 'wave/tmp' -A '*grb2.gz'", src_dir))
system("gunzip wave/tmp/*.gz")

# List of WW3 grid codes from high to low resolution
grid_codes <- c("ak_4m", "ecg_4m", "nsb_4m", "oz_4m", "wc_4m",
                "ak_10m", "ecg_10m", "med_10m", "nsb_10m", "nwio_10m", "oz_10m",
                "pi_10m", "wc_10m", "glo_30m")

# Resolution codes matching each grid (1 = 30 arcmin, 2 = 10 arcmin, 3 = 4 arcmin)
grid_res <- c(rep(3, 5), rep(2, 8), 1)


# Load land mask (also serves as template for final grid)
land_mask <- raster("reeflandarea/land_final.grd")

# Create masks for area covered by each grid
has_vals <- lapply(grid_codes, function(g) {
    calc(brick(paste0("wave/tmp/multi_reanal.", g, ".wind.", month, ".grb2")), 
         function(x) any(!is.na(x)))
})

# Convert grids with negative longitude values to (0, 360) longitude range
to_rotate <- which(vapply(has_vals, xmin, 0) < -0.25)
has_vals[to_rotate] <- lapply(has_vals[to_rotate], extend_inv_rotate)
# Resample each WW3 grid mask to final grid
has_vals <- lapply(has_vals, function(r) resample(r, land_mask) > 0)

# Convert "1" values in mask to resolution codes (1, 2 or 3)
for (i in seq_along(grid_res)) {
    has_vals[[i]] <- has_vals[[i]] * grid_res[[i]]
}

# Determine highest resolution code for each cell
ww3_res <- do.call(mosaic, c(has_vals, fun = max))
ww3_res <- mask(ww3_res, land_mask, maskvalue = 1, datatype = "INT1U")

# Using fetch data, find 'sheltered' points where we need to
#  calculate local wave energy
#  (i.e. where more than half of fetch values are smaller than size of WW3 cell)
fetch <- readRDS("wave/fetch_res.RData")
fetch <- extract(ww3_res, fetch, sp = TRUE)
colnames(fetch@data)[ncol(fetch@data)] <- "res_code"
# Cells with fetch code == 0 are not covered by any grid
fetch$res_code[fetch$res_code == 0] <- NA

thresh <- c(50000, 18400, 7400) # threshold for each resolution code in meters
sheltered <- rowSums(sweep(fetch@data[, 1:16], 1, thresh[fetch$res_code], "<")) > 8
fetch_sheltered <- fetch[which(sheltered), 1:16]

# Find highest resolution grid with values at each sheltered point 
has_vals <- stack(has_vals)
shelt_vals <- extract(has_vals, fetch_sheltered)
# Since grids were ordered from high to low resolution, 
# we keep the first grid with values at point
grid_num <- apply(shelt_vals, 1, function(x) which(x > 0)[1])
fetch_sheltered$grid <- grid_codes[grid_num]

# Get depth for sheltered points from MARSPEC bathymetry layer 
#  (note: MARSPEC uses negative value for depth)
bathy <- raster(file.path(mars_dir, "bathymetry_30s/bathy_30s.tif"))
bathy <- inv_rotate(bathy)
# Use bilinear interpolation to reduce NAs
fetch_sheltered$depth <- -extract(bathy, fetch_sheltered, method = "bilinear")
fetch_sheltered <- fetch_sheltered[!is.na(fetch_sheltered$depth), ]

saveRDS(fetch_sheltered, "wave/fetch_sheltered.RData")

