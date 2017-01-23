# Create high-resolution land mask from GSHHS shoreline shapefile

library(rgdal)
library(raster)
source("utils.R")

gshhs_dir = "{{Insert path to GSHHS shapefiles}}"

# Combine full resolution ('f') shoreline for world except Antarctica (L1)
#  and high resolution ('h') sea-ice line for Antarctica (L5)
#  (not using 'f' for L5 due to a missing vertex at (180, -90))
gshhs1 <- readOGR(file.path(gshhs_dir, "f"), "GSHHS_f_L1")
gshhs5 <- readOGR(file.path(gshhs_dir, "h"), "GSHHS_h_L5")
gshhs_full <- rbind(gshhs1, gshhs5, makeUniqueIDs = TRUE)

# Rasterize shoreline polygons to 10x our final grid resolution
#  1 = land, 0 = sea
new_rast <- raster(nrows = 43190, ncols = 86390, 
                   crs = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
gshhs_rast <- rasterize(gshhs_full, new_rast, field = 1, background = 0,
                        filename = "reeflandarea/gshhs_rast.grd", datatype = "INT1U")

# Rotate to (0, 360) longitude range
gshhs_rotate <- inv_rotate(gshhs_rast)
writeRaster(gshhs_rotate, filename = "reeflandarea/gshhs_rotated.grd", datatype = "INT1U")

# 'Wrap' 5 cells from right edge to left edge, and add 5 cells to top and bottom
#  to match extent of final layer
extent_left <- extent(360 - 5*res(gshhs_rotate)[1], 360, -90, 90)
extent_right <- extent(0, 360 - 5*res(gshhs_rotate)[1], -90, 90)
gshhs_left <- crop(gshhs_rotate, extent_left)
xmin(gshhs_left) <- xmin(gshhs_left) - 360
xmax(gshhs_left) <- 0
gshhs_right <- crop(gshhs_rotate, extent_right)
gshhs_top <- raster(nrows = 5, ncols = ncol(gshhs_rotate), vals = 0,
                    ext = extent(xmin(gshhs_left), xmax(gshhs_right), 
                                 90, 90 + 5*res(gshhs_rotate)[2]))
gshhs_bottom <- raster(nrows = 5, ncols = ncol(gshhs_rotate), vals = 1,
                       ext = extent(xmin(gshhs_left), xmax(gshhs_right), 
                                    -90 - 5*res(gshhs_rotate)[2], -90))

gshhs_wrap <- merge(gshhs_left, gshhs_right, gshhs_top, gshhs_bottom,
                     filename = "reeflandarea/gshhs_wrap.grd", datatype = "INT1U")

# Aggregate by factor of 10 to get resolution of final products
#  fun = min so that cell is only land if all underlying values are land
land_final <- aggregate(gshhs_wrap, fact=10, fun = min, datatype = "INT1U",
                        filename = "reeflandarea/land_final.grd")


