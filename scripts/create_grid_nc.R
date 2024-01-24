rm(list = ls())

devtools::load_all()

# create the spatial_grid.nc file in inst/extdata, which contains the spatial grid information.

# one random netcdf file. They should all have the same grid
nc = nc_open('https://thredds.met.no/thredds/dodsC/metpparchivev3/2022/12/31/met_analysis_1_0km_nordic_20221231T23Z.nc')


# variables:
nc_lat = nc$var$latitude
nc_lat$name = 'lat'
lats = ncvar_get(nc,'latitude')

nc_lon = nc$var$longitude
nc_lon$name = 'lon'
lons = ncvar_get(nc,'longitude')

nc_proj = nc$var$projection_lcc
proj = ncvar_get(nc,'projection_lcc')


nc_alt = nc$var$altitude
alts = ncvar_get(nc,'altitude')

nc_land_area_fraction = nc$var$land_area_fraction
lafs = ncvar_get(nc,'land_area_fraction')

global_atts = ncatt_get(nc,0)
# get rid of time
global_atts = global_atts[-c(15,10,9)]


# create netcdf:
grid_nc = nc_create(nc_out,vars = list(nc_lat,nc_lon,nc_alt,nc_land_area_fraction,nc_proj))

# put values:
ncvar_put(nc = grid_nc,varid = 'lat',vals = lats)
ncvar_put(nc = grid_nc,varid = 'lon',vals = lons)
ncvar_put(nc = grid_nc,varid = 'altitude',vals = alts)
ncvar_put(nc = grid_nc,varid = 'land_area_fraction',vals = lafs)
ncvar_put(nc = grid_nc,varid = 'projection_lcc',vals = proj)

for(i in seq_along(global_atts))
{
  ncatt_put(grid_nc,0,attname = names(global_atts)[i],attval = global_atts[[i]])
}
nc_close(grid_nc)

test = SeaVal::netcdf_to_dt(nc_out)
