# Prepare raster bricks from SEDAC data

library(raster)
source("utils.R")

sedac_dir = "{{Insert path to SEDAC data directory}}"

# Load SEDAC human population count layers for years 1990, 1995, 2000, 2005, 2010, 2015, 2020
SEDAC_90 <- raster(file.path(sedac_dir, "GPWv3/PopGrid/glcount90/glp90ag/"))
SEDAC_95 <- raster(file.path(sedac_dir, "GPWv3/PopGrid/glcount95/glp95ag/"))
SEDAC_00 <- raster(file.path(sedac_dir, "GPWv4",
    "gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-2000", 
    "gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals_2000.tif"))
SEDAC_05 <- raster(file.path(sedac_dir, "GPWv4",
    "gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-2005", 
    "gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals_2005.tif"))
SEDAC_10 <- raster(file.path(sedac_dir, "GPWv4",
    "gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-2010", 
    "gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals_2010.tif"))
SEDAC_15 <- raster(file.path(sedac_dir, "GPWv4",
    "gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-2015", 
    "gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals_2015.tif"))
SEDAC_20 <- raster(file.path(sedac_dir, "GPWv4",
    "gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-2020", 
    "gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals_2020.tif"))

# Rotate the SEDAC layers to a (0, 360) longitude range
SEDAC_90_rot <- inv_rotate(SEDAC_90)
SEDAC_95_rot <- inv_rotate(SEDAC_95)
SEDAC_00_rot <- inv_rotate(SEDAC_00)
SEDAC_05_rot <- inv_rotate(SEDAC_05)
SEDAC_10_rot <- inv_rotate(SEDAC_10)
SEDAC_15_rot <- inv_rotate(SEDAC_15)
SEDAC_20_rot <- inv_rotate(SEDAC_20)

# Create raster brick of GPWv4 layers
#  Use bandorder = "BIP" (storage by pixel) to facilitate spatial slicing
SEDAC_00_stack <- stack(SEDAC_00, SEDAC_05, SEDAC_10, SEDAC_15, SEDAC_20)
SEDAC_00_brick <- brick(SEDAC_00_stack, bandorder = "BIP",
                        filename = "humanpop/SEDAC_00_brick.grd")
SEDAC_00_rot_stack <- stack(SEDAC_00_rot, SEDAC_05_rot, SEDAC_10_rot, SEDAC_15_rot, SEDAC_20_rot)
SEDAC_00_rot_brick <- brick(SEDAC_00_rot_stack, bandorder = "BIP",
                            filename = "humanpop/SEDAC_00_rot_brick.grd")

# Create raster brick of GPWv3 layers
SEDAC_90_stack <- stack(SEDAC_90, SEDAC_95)
SEDAC_90_brick <- brick(SEDAC_90_stack, bandorder = "BIP",
                        filename = "humanpop/SEDAC_90_brick.grd")
SEDAC_90_rot_stack <- stack(SEDAC_90_rot, SEDAC_95_rot)
SEDAC_90_rot_brick <- brick(SEDAC_90_rot_stack, bandorder = "BIP",
                            filename = "humanpop/SEDAC_90_rot_brick.grd")
