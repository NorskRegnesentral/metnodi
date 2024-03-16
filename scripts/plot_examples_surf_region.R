devtools::load_all()
library(stars)


data_dir = '/nr/samba/PostClimDataNoBackup/sensKlimGJ/MetNo/MetNo_for_SURF_region/'

yy = 2016
mm = 8


fn = paste0(yy,'-',ifelse(mm<10,0,''),mm,'.nc')

st = read_ncdf(paste0(data_dir,fn))

times_new = time(st) + 2*3600# tried force_tz and with_tz, but they didn't work as expected.
times_new = lubridate::with_tz(times_new,tz = 'CEST')

st_new = aggregate(st,by = 'days', FUN = sum)

st_new = aperm(st_new,c('x','y','time'))

#get map:

map = gmaps_plot(st_new[,,,1],fix_crs = 'metno',only_map = TRUE)

pp = gmaps_plot(st_new[,,,1],fix_crs = 'metno',map = map)
pp


# get surf data for plotting:


dir0 = '/nr/samba/PostClimDataNoBackup/sensKlimGJ/SURF/Oslo_study_new/'

d1 = paste0(yy,'-',ifelse(mm<10,0,''),mm,'-01') |> as.Date()
d2 = paste0(yy,'-',ifelse(mm+1<10,0,''),mm+1,'-01') |> as.Date() - 1

dates = surf_dates(dates = seq(d1,d2,by = 1))

# get SURF spatial indices:

nc = nc_open(paste0(dir0, '2016_8.nc'))

lons = ncvar_get(nc,'lon')
lats = ncvar_get(nc,'lat')

bbox = c(range(lons),range(lats))[c(1,3,2,4)]


times = ncvar_get(nc,'time')
times_new = as.POSIXct(60*times-7200,origin = '1900-01-01')
times_new[781] # selected time
dt1 = SeaVal::netcdf_to_dt(paste0(dir0, '2016_8.nc'),vars = c('xa'),subset_list = list(time = times[781]))
dt2 = SeaVal::netcdf_to_dt(paste0(dir0, '2016_8.nc'),vars = c('lon','lat'))
dt = merge(dt1,dt2,by = 'location_index')
