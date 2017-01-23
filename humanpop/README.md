# Marine Socio-Environmental Covariates - Human Population

The R scripts in this directory compute the human population for two radii (20km and 50km)
around each marine grid cell, based on input datasets of gridded population data.

**1-prep_SEDAC_bricks.R**

Combines the individual year SEDAC rasters into two RasterBricks, for the GPWv3 and GPWv4 grids.
Also creates a "rotated" (0 to 360 degrees longitude range) version of the RasterBricks to 
extract values for points close to the 180th meridian.

**2-human_pop_count.R**

Using the distance to shoreline raster, extracts the MSEC grid cell midpoints within 25km
(resp., 55km) of land for calculation of the population within a 20km (resp., 50km) radius.
Calculates the human population around each point by creating a circular buffer of the
appropriate radius, projecting it to geographical coordinates and aggregating the population
of cells covered by that buffer via the `raster::extract` function. Rasterizes the values back
to the MSEC grid, assigning a values of 0 to marine cells beyond the specified distance from land.


## Data Input Requirements

* The Socioeconomic Data and Applications Center (SEDAC) Gridded Population of the World rasters:
    - [GPWv3](http://sedac.ciesin.columbia.edu/data/collection/gpw-v3) for 1990 and 1995;
    - [GPWv4](http://sedac.ciesin.columbia.edu/data/collection/gpw-v4) in 5-year intervals
    from 2000 to 2020.

* The [distance to shoreline raster](ftp://ftp.soest.hawaii.edu/gshhg/dist_to_GSHHG_v2.3.4_1m.nc) from the Global Self-consistent, Hierarchical, High-resolution Geography (GSHHG) Database.





