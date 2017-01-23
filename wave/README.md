# Marine Socio-Environmental Covariates - Wave Energy

The R scripts in this directory produce aggregate statistics of wave energy
based on hindcast reanalysis data from the NOAA WAVEWATCH III (WW3) wave model. 
For points that are sheltered from oceanic waves at a scale smaller than the
WW3 grid, we calculate the energy of locally-generated waves from
wind velocity, fetch length and depth.

**1-calc_fetch.R**

Using the `fetch_len` function from the waver R package, calculates the
fetch (up to 50km) for all grid cells with midpoints less than 55km from a coast,
for 16 bearings (every 22.5 degree).

**2-get_sheltered_pts.R**

Determines which points in our final grid are "sheltered", i.e. for more than 50% of
wind bearings, the fetch (distance to shore) is less than the cell size for 
the highest resolution WW3 grid with values at that point.

**3-calc_wave_energy.R**

For each grid and each of the 372 months of WW3 data, calculates the wave
energy (based on significant height and peak period rasters) at each 3-hour time
step, then computes daily averages. For sheltered points, calculates the 
wave energy at each time point based on wind velocity, fetch (for bearing
closest to wind direction) and depth.

**4-wave_aggregate.R**

For each grid, calculates wave energy means by year and by day of the year,
then computes the global mean, intra-annual and inter-annual standard deviations.
Performs the same aggregation for the values at sheltered points.

**5-wave_resample.R**

For each of our three summary statistics, resamples from each WW3 grid to our
final (2.5 arc-minute) grid, keeping the value from the highest-resolution
grid available at each point. Overrides the values for sheltered points based
on those calculated from wind and fetch. In addition to the three summary
statistics, also produces flag rasters indicating the resolution at each point,
as well as whether or not the point was fetch-limited (sheltered).


## Data Input Requirements

* NOAA [WAVEWATCH III hindcast reanalysis](http://polar.ncep.noaa.gov/waves/nopp-phase1/) 
data (31 years from 1979-2009). This shouldn't be downloaded in advance, as the
`3-calc_wave_energy.R` script automatically downloads the data for a month into a
temporary folder, computes daily wave energy rasters, then deletes the temporary
files before repeating for the next month.

* [Coastline shapefiles](ftp://ftp.soest.hawaii.edu/gshhg/gshhg-shp-2.3.6.zip)
from the Global Self-consistent, Hierarchical, High-resolution Geography (GSHHG) Database.

* Bathymetry raster (30 arc-second resolution) from [MARSPEC](http://www.marspec.org).






