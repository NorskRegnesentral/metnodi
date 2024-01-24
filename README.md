
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
netcdf format. It supports downloading and combining data for a time
range (original format is one separate file per hour), and for target
locations or a rectangular area.

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

# area around Oslo
lon_min = 8
lon_max = 12
lat_min = 58
lat_max = 62

time_start = "2024-01-01"
time_end = "2024-01-05"

# download all available data in this area and store it in a temporary file:

download_metno_analysis_rect(file_out = 'oslo_weather_202401.nc', 
                             out_dir = tempdir(),
                             bbox = c(lon_min,lat_min,lon_max,lat_max),
                             time_start = time_start,
                             time_end = time_end)

# download precipitation for specific locations:

coords = data.table(lon = runif(50,min = lon_min, max = lon_max),
                    lat = runif(50,min = lat_min, max = lat_max))

download_metno_analysis_nn(file_out = 'oslo_weather_202401.nc', 
                           out_dir = tempdir(),
                           coords = coords,
                           time_start = time_start,
                           time_end = time_end,
                           vars = 'precipitation_amount')

# use weather_variables() to see the names of all weather variables provided by met Nordic
```
