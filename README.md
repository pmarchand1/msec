# Marine Socio-Environmental Covariates

This repository contains the R source code to create the following summary 
statistics of environmental and anthropogenic variables on a common
2.5 arc-minute grid, for use in marine science research:

- Net Primary Productivity (with correction applied to shallow locations)
- Reef Area (within 15km and 200km) and Land Area (within 15km and 50km)
- Wave Energy
- Human Population within 20km and 55km for 7 years between 1990 and 2020
- Distance to Market

The details for each variable's calculation can be found in the `README.md` file
of the corresponding subfolder.

## Citation

These data products were published in the following study:

* Yeager, L.A., Marchand, P., Gill, D.A., Baum, J.K., and McPherson, J.M. In review. 
Queryable global layers of environmental and anthropogenic variables for marine ecosystem studies. 

Please cite this article for any work that re-uses this code.

## Web Application

A web application to extract these variables for specific points of interest can
be found at: [http://shiny.sesync.org/apps/msec](http://shiny.sesync.org/apps/msec).

## Input Data Sources

The data products were created from the following publicly available data sources:

- primary productivity rasters from [NOAA Coast Watch](http://coastwatch.pfeg.noaa.gov/erddap/griddap/erdPPbfp28day.html);

- bathymetry and distance to shore rasters (30 arc-second resolution) from 
[MARSPEC](http://www.marspec.org);

- coral reef locations from the 
[Reefs at Risk Revisited](http://www.wri.org/publication/reefs-risk-revisited) dataset;

- high-resolution coastlines from the Global Self-consistent, Hierarchical, 
High-resolution Geography ([GSHHG](https://www.soest.hawaii.edu/pwessel/gshhg/)) database;

- wind and wave data from the NOAA [WAVEWATCH III hindcast reanalysis](http://polar.ncep.noaa.gov/waves/nopp-phase1/);

- the Socioeconomic Data and Applications Center (SEDAC) 
[Gridded Population of the World](http://sedac.ciesin.columbia.edu/data/collection/gpw-v4)
rasters;

- provincial capital locations from the ESRI
[World Cities map layer](https://www.arcgis.com/home/item.html?id=dfab3b294ab24961899b2a98e9e8cd3d).


## R Packages Required

- The following [tidyverse](https://blog.rstudio.org/2016/09/15/tidyverse-1-0-0/) packages: 
*dplyr*, *tidyr*, *stringr* and *lubridate*;
- the *sp*, *rgdal*, *rgeos* and *raster* packages to manipulate spatial objects;
- the *geosphere* package to calculate geodesic distances; and
- the *waver* package to calculate fetch and wave energy.

## Acknowledgement

This work was supported by the National Socio-Environmental Synthesis Center (SESYNC)
under funding received from the National Science Foundation DBI-1052875.

## Miscellaneous notes

- The `utils.R` script in the main folder contains a few functions to manipulate
spatial objects, e.g. creating a rectangular buffer of minimum distance, 
converting vector and raster layers from a (-180, 180) to a (0, 360) longitude range.
These functions are re-used in multiple parts of the project.

- The land mask raster (`land_final.grd`), which serves as a grid template for each 
final layer, should be generated first via the script found in the `reeflandarea` subfolder.

- As indicated in the various scripts, the most computationally-intensive steps were
performed in parallel on a HPC cluster. More specifically, these calculations were
run on a 20-node SLURM cluster at SESYNC, each node having 8 cores and 60 Gb of RAM.

- All processing steps from the original data sources to the final products were
performed in R, with the exception of the distance buffers around the Reefs at Risk
coral reef map, which were computed in ArcGIS. This is due to the absence of 
efficient R functions to compute distance buffers in geographic (unprojected)
coordinates.

- While all raster calculations use the native `.grd` format from the R raster package,
the final products were saved as NetCDF files with a compression level of 5, 
using the ncdf4 package.


