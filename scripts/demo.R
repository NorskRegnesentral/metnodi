# install.packages("devtools")
# devtools::install_github("NorskRegnesentral/metnodi")
# rm(list = ls())

devtools::load_all()

#library(metnodi)

### downloading a rectangular area ###

# area around Oslo
lon_min = 8
lon_max = 12
lat_min = 58
lat_max = 62

time_start = "2024-01-01"
time_end = "2024-01-01"

# download all available weather variables in this area and store it in a temporary file:

download_metno_analysis_rect(file_out = 'oslo_weather_202401.nc',
                             out_dir = tempdir(),
                             bbox = c(lon_min,lat_min,lon_max,lat_max),
                             time_start = time_start,
                             time_end = time_end)

nc = ncdf4::nc_open(paste0(tempdir(),'/oslo_weather_202401.nc'))
nc

### downloading data for specific locations ###

# Here, we only download precipitation. Use weather_variables() to see the names of all available weather variables

# some locations in the area around Oslo:
coords = data.table(lon = runif(50,min = lon_min, max = lon_max),
                    lat = runif(50,min = lat_min, max = lat_max))

download_metno_analysis_nn(file_out = 'oslo_weather2_202401.nc',
                           out_dir = tempdir(),
                           coords = coords,
                           time_start = time_start,
                           time_end = time_end,
                           vars = 'precipitation_amount')


nc2 = ncdf4::nc_open(paste0(tempdir(),'/oslo_weather2_202401.nc'))
nc2


### plotting ###

# here we plot the first time-slice of precipitation downloaded above:

st = read_stars_metno(paste0(tempdir(),'/oslo_weather_202401.nc'),var = 'precipitation_amount')
pp = gmaps_plot(st) + ggplot2::scale_fill_gradient(name = 'mm',low = 'white',high = 'blue')
plot(pp)
