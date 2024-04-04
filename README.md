
<!-- README.md is generated from README.Rmd. Please edit that file -->

# metnodi

<!-- badges: start -->
<!-- badges: end -->

The goal of metnodi is to provide an interface to the met Nordic data
(metnodi stands for met Nordic data interface). For a detailed
description of this data and links to the original source please see the
[met Nordic github page](https://github.com/metno/NWPdocs/wiki).
Currently, the package merely simplifies the download of targeted
weather data from thredds.met.no, and saves this data in convenient
netcdf format. It supports downloading for target locations (finding the
nearest neighbor gridcells of the MET Nordic data), or a rectangular
area (in lon/lat). Moreover, downloads are pooled in time into either
daily-, monthly-, yearly- or a single file, depending on the size of the
downloaded area. This format is more convenient to work with than the
original format, which is a single netcdf per hour. The package also has
a function to plot the downloaded data.

## Installation

You can install the development version of metnodi from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("NorskRegnesentral/metnodi")
```

## Example

This is a basic example which shows you how to download data:

``` r
library(metnodi)

### downloading a rectangular area ###

# area around Oslo
lon_min = 8
lon_max = 12
lat_min = 58
lat_max = 62

time_start = "2024-01-01"
time_end = "2024-01-03"

# download all available weather variables in this area and store it in a temporary file:

download_metno_analysis_rect(file_out = 'oslo_weather_202401.nc', 
                             out_dir = tempdir(),
                             bbox = c(lon_min,lat_min,lon_max,lat_max),
                             time_start = time_start,
                             time_end = time_end)

### downloading data for specific locations ###

# Here, we only download precipitation. Use weather_variables() to see the names of all available weather variables 

# some locations in the area around Oslo:
coords = data.table(lon = runif(50,min = lon_min, max = lon_max),
                    lat = runif(50,min = lat_min, max = lat_max))

download_metno_analysis_nn(file_out = 'oslo_weather_202401.nc', 
                           out_dir = tempdir(),
                           coords = coords,
                           time_start = time_start,
                           time_end = time_end,
                           vars = 'precipitation_amount')


### plotting ###

st = read_stars_metno(paste0(tempdir(),'/oslo_weather_202401.nc'),var = 'precipitation_amount')
pp = gmaps_plot(st) + ggplot2::scale_fill_gradient(name = 'mm',low = 'white',high = 'blue')
plot(pp)
```

![](https://github.com/NorskRegnesentral/metnodi/blob/master/example_plot.png)
