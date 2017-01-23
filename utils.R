# Utility functions for the computation of the
#  Marine Socio-Environmental Covariates data layers

# Clip a raster r to a distance of at least dmax (in meters) 
#  around a point vector p (long, lat)
clip_raster_pt <- function(r, p, dmax) {
    long <- p[1]
    lat <- p[2]
    lat_dist <- 111600 # approx. distance (in m) between degrees of latitude
    ybuf = dmax / lat_dist
    xbuf = ybuf / cospi(abs(lat) / 180)
    if (long - xbuf < xmin(r)) {
        west_r <- crop(r, extent(xmin(r), long + xbuf, 
                                 lat - ybuf, lat + ybuf), snap = "out")
        east_r <- crop(r, extent(long - xbuf + 360, xmax(r), 
                                 lat - ybuf, lat + ybuf), snap = "out")
        extent(east_r) <- extent(east_r) - c(360, 360, 0, 0)
        crop_r <- merge(west_r, east_r)
    } else if(long + xbuf > xmax(r)) {
        west_r <- crop(r, extent(xmin(r), long + xbuf - 360, 
                                 lat - ybuf, lat + ybuf), snap = "out")
        east_r <- crop(r, extent(long - xbuf,  xmax(r), 
                                 lat - ybuf,lat + ybuf), snap = "out")
        extent(west_r) <- extent(west_r) + c(360, 360, 0, 0)
        crop_r <- merge(west_r, east_r)
    } else {
        ext <- extent(long - xbuf, long + xbuf, lat - ybuf, lat + ybuf)
        crop_r <- crop(r, ext, snap = "out")
    }
    crop_r
}

# Create clipping rectangle around point p (longlat coordinates)
#  to guarantee at least dmax (in meters) on each side
get_clip_rect <- function(p, dmax, projected) {
    lat_dist <- 111600 # approx. distance (in m) between degrees of latitude
    long <- coordinates(p)[1]
    lat <- coordinates(p)[2]
    ybuf = dmax / lat_dist
    xbuf = ybuf / cospi(abs(lat) / 180)
    # Split clip_rect in two if it would overlap 180th meridian
    if (long - xbuf < -180) {
        westr <- poly_rect(-180, lat - ybuf, long + xbuf, lat + ybuf)
        eastr <- poly_rect(long - xbuf + 360, lat - ybuf, 180, lat + ybuf)
        clip_rect <- SpatialPolygons(list(Polygons(list(westr, eastr), ID = 1)),
                                     proj4string = CRS(proj4string(p)))
    } else if(long + xbuf > 180) {
        westr <- poly_rect(-180, lat - ybuf, long + xbuf - 360, lat + ybuf)
        eastr <- poly_rect(long - xbuf, lat - ybuf, 180, lat + ybuf)
        clip_rect <- SpatialPolygons(list(Polygons(list(westr, eastr), ID = 1)),
                                     proj4string = CRS(proj4string(p)))
    } else {
        r1 <- poly_rect(long - xbuf, lat - ybuf, long + xbuf, lat + ybuf)
        clip_rect <- SpatialPolygons(list(Polygons(list(r1), ID = 1)),
                                     proj4string = CRS(proj4string(p)))
    }
    clip_rect
}

# Transforms an input raster with x-extent from (-180, 180) to (0, 360)
# i.e. inverse of raster::rotate
inv_rotate <- function(r) {
    s <- r
    xmin(s) <- 0
    xmax(s) <- 360
    s <- rotate(s)
    xmin(s) <- 0
    xmax(s) <- 360
    s
}

# Same as previous, but first increase raster extent to full world
extend_inv_rotate <- function(r) {
    r <- extend(r, extent(-180, 180, -90, 90))
    xmin(r) <- 0
    xmax(r) <- 360
    r <- rotate(r)
    xmin(r) <- 0
    xmax(r) <- 360
    r
}

# Create a Polygon object corresponding to rectangle with given coords
poly_rect <- function(xmin, ymin, xmax, ymax) {
    Polygon(cbind(c(rep(xmin, 2), rep(xmax, 2), xmin),
                  c(ymin, rep(ymax, 2), rep(ymin, 2))))
}

# Function to change polygon from (-180, 180) to (0, 360) longitude range
#  Assumes input is a SpatialPolygons with a single Polygon
rotate_poly <- function(poly) {
    poly
    coords <- poly@polygons[[1]]@Polygons[[1]]@coords
    coords[coords[, 1] < 0, 1] <- coords[coords[, 1] < 0, 1] + 360
    SpatialPolygons(
        list(Polygons(list(Polygon(coords, hole = FALSE)), 1)),
        proj4string = CRS(proj4string(poly))
    )
}

# Rotate SpatialPointsDataFrame from (-180, 180) to (0, 360) longitude
rotate_pts <- function(sp_df, tol = 1E-15) {
    west_pts <- sp_df@coords[, 1] < 0 - tol
    sp_df@coords[west_pts, 1] <- sp_df@coords[west_pts, 1] + 360
    sp_df@bbox["x", "min"] <- min(sp_df@coords[, 1])
    sp_df@bbox["x", "max"] <- max(sp_df@coords[, 1])
    sp_df
}