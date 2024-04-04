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



nc = nc_open(paste0(tempdir(),'/oslo_weather_202401.nc'))

### plotting ###

#' Function for reading downloaded netcdfs as stars
#'
#' Corrects the projection, which is not correctly identified by [stars::read_stars()] and [stars::read_ncdf()]
#' @param ... passed on to [stars::read_ncdf()].
#' @export

read_stars_metno = function(...)
{
  st = stars::read_ncdf(...)
  p4s = '+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +pm=0 +no_defs +units=m'
  sf::st_crs(st) = sf::st_crs(p4s)
  return(st)
}


