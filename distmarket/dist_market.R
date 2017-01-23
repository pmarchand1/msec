# Calculate distance from each grid cell's midpoint to the nearest provincial captial

library(rgdal)
library(geosphere)
library(raster)

data_dir = "{{insert path to Provincial_capitals shapefile}}"

# Read in capitals shapefile
caps <- readOGR(dsn = data_dir, layer = "Provincial_capitals")

# Get land mask (also serves as a grid template)
land <- raster("reeflandarea/land_final.grd")

# Get midpoint coordinates of all points not on land
ocean_cells <- rasterToPoints(land, fun = function(x) {x == 0})
ocean_cells <- ocean_cells[, -3] # no need for data column

# Rotate longitude form (0, 360) to (-180, 180) range
ocean_cells[, 1] <- ocean_cells[, 1] - (ocean_cells[, 1] > 180) * 360
# Ensure maximum latitude is 90
ocean_cells[ocean_cells[, 2] > 90, 2] <- 90

# Calculate distance in km from row "i" of ocean cells to nearest capital
dist_capital <- function(i) {
    min(distGeo(ocean_cells[i, ], caps)) / 1000
}

# This calculation was processed in parallel on a HPC cluster.
res <- lapply(1:nrow(ocean_cells), dist_capital)


# Add distances to ocean_cells and transform longtudes back to (0, 360) degrees
ocean_cells <- data.frame(cbind(ocean_cells, do.call(rbind, res)))
colnames(ocean_cells) <- c("long", "lat", "dist")
ocean_cells$long <- ocean_cells$long + (ocean_cells$long < 0) * 360

# Convert to SpatialPointsDataFrame and rasterize to grid given by 'land' layer
coordinates(ocean_cells) <- ~long + lat 

dist_rast <- rasterize(ocean_cells, land, field = "dist",
                       background = NA, filename = "dist_market.grd")




