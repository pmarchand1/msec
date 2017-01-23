# Marine Socio-Environmental Covariates - Reef and Land Area

The R scripts in this directory compute, for each marine grid cell, the reef area
within 15km and 200km, and the land area within 15km and 50km.

**1-create_land_mask.R**

Creates a high-resolution raster version of the GSHHS shoreline and produces a 
land mask that will be applied to all MSEC layers.

**2-land_area_global_layers.R**

Using the distance to shoreline raster, extracts the MSEC grid cell midpoints within
20km (resp., 55km) of land for calculation of the land area within a 15km (resp., 50km)
radius. For each point, creates a circular buffer of the appropriate radius, projects
it to geographical coordinates and aggregates the areas of land cells covered by the buffer,
using to a high-resolution rasterized version of the GSHHS shoreline. Rasterizes the values
back to the MSEC grid, assigning a value of 0 to marine cells beyond the specified
distance from land.

**3-reef_area_global_layers.R**

Using pre-computed buffers, extracts the MSEC grid cell midpoints within 20km
(resp., 205km) of coral reef locations (Reefs at Risk Revisited dataset),
For each point, creates a 15km (resp., 200km) radius buffer, projects it to the
Reefs at Risk 500-meter equal-area grid, and computes the reef area covered by
the buffer. Rasterizes the values backto the MSEC grid, as above.


## Data Input Requirements

* Coral reef locations from the [Reefs at Risk Revisited](http://www.wri.org/publication/reefs-risk-revisited)
dataset (equal-area 500-meter global raster). We also created 20km and 205km buffers
around those reef areas in ArcGIS, rather than in the included R scripts. 

* The [coastline shapefiles](ftp://ftp.soest.hawaii.edu/gshhg/gshhg-shp-2.3.6.zip)
and [distance to shoreline raster](ftp://ftp.soest.hawaii.edu/gshhg/dist_to_GSHHG_v2.3.4_1m.nc)
from the Global Self-consistent, Hierarchical, High-resolution Geography (GSHHG) Database.





