# Marine Socio-Environmental Covariates - Distance to Market

The R script in this directory produces a 2.5 arc-minute raster layer with values
corresponding to the distance (in km) from each raster cell's midpoint to the nearest
provincial capital. Geodesic distances (WGS84 ellipsoid) are calculated with the
`geosphere::distGeo` function.

## Data Input Requirements

A shapefile of point locations for provincial capitals, such as the ESRI
[World Cities map layer](https://www.arcgis.com/home/item.html?id=dfab3b294ab24961899b2a98e9e8cd3d).





