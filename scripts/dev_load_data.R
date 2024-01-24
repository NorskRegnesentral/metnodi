rm(list = ls())

devtools::load_all()


vars = weather_variables()

# if(any(c('precipitation_amount_gt','precipitation_quantiles') %in% vars)){
#   stop('precipitation_amount_gt and precipitation_quantiles not yet implemented (they have an additional dimension)')
# }

# get data table of locations:

N = 1000

nc = nc_open('https://thredds.met.no/thredds/dodsC/metpparchivev3/2022/12/31/met_analysis_1_0km_nordic_20221231T23Z.nc')

# variables:
nc_lat = nc$var$latitude
nc_lat$name = 'lat'
lats = ncvar_get(nc,'latitude')

nc_lon = nc$var$longitude
nc_lon$name = 'lon'
lons = ncvar_get(nc,'longitude')


bbox = c(10,50,15,60)

time_start = '2022-01-01'
time_end = '2022-01-03'
use_rerun = TRUE

download_metno_analysis_rect(file_out = 'test.nc',out_dir = '/nr/samba/PostClimDataNoBackup/KMdata/Met_Nordic/',bbox = bbox,time_start = time_start,time_end = time_end)

test = SeaVal::netcdf_to_dt( paste0(out_dir,file_out))


non_spatial_length = prod(dim(var_vals)[which(!dim_names %in% c('x','y'))])

# subset










