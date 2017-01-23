# Marine Socio-Environmental Covariates - Net Primary Productivity

The R scripts in this directory produce aggregate statistics of marine net primary productivity
from NOAA Coast Watch 8-day layers (2.5 arc-minute grid), and correct the values
of shallow areas based on value at neighboring cells.

**1-npp_aggregate.R**

From the 8-day layers, computes annual averages and averages for the same sampling
day across years. The min, mean, max and standard deviation (sd) of the day-of-year
means are saved as separate raster layers, as is the standard deviation of
annual means (interann_sd).

**2-prep_npp_masks.R**

Produces two binary spatial layers (masks) that are used in the correction step below.
The depth mask identifies all grid cells with depth < 30m, whereas the interpolation
mask identifies cells that are within 9.3km of a coast or reef.

**3-npp_correct.R**

For one of the 5 aggregate NPP layers (mean, min, max, sd and interann_sd),
this script replaces the NPP values at shallow cells (< 30m depth), and any
other missing values, with an average of nearby cells, preferably reef/coastal 
cells (interpolation mask). It outputs the corrected raster as well as a "flag"
raster that indicates which cells where interpolated. The script must be run 
for all 5 aggregates.

*Note on NPP flags*

The "flag" layer produced at the last step can take the following values:

* 0 - not interpolated
* 1 - interpolated from nearby reef/coastal cells
* 2 - interpolated from other nearby cells

The *npp_correct.R* script produces a separate flag layer for each summary statistic.
The flags for mean, min and max are the same, but those for the standard 
deviations may differ for <1% of the cells (almost all near the poles, where data
is missing on most days of the year). To produce the single flag raster distributed
with MSEC, we took the maximum of the flag values over the 5 statistics at each point. 


## Data Input Requirements

* Primary Productivity rasters from [NOAA Coast Watch](http://coastwatch.pfeg.noaa.gov/erddap/griddap/erdPPbfp28day.html) 
(506 files, 8-day intervals from July 2002 to October 2013).

* Bathymetry and distance to shore rasters (30 arc-second resolution) from 
[MARSPEC](http://www.marspec.org).

* Shapefile representing a 9.3km buffer around reef locations from the 
[Reefs at Risk Revisited](http://www.wri.org/publication/reefs-risk-revisited) dataset. 
We created this buffer in ArcGIS rather than in the included R scripts. 





