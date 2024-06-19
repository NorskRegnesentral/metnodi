
# For R CMD check:
.datatable.aware = TRUE

lat = lat_metno = lon = lon_metno = vector_ind = x_ind = y_ind = . = NULL


#'
#' Convert to POSIXct with the right timezone
#'
#' @param x object that can be converted to POSIXct.
#' @examples
#' as.GMT('2024-01-01')
#' as.GMT('2024-01-01 04:00:30')
#' as.GMT(Sys.time())
#'
#'
#' @import data.table
#' @import ncdf4
#' @export

as.GMT = function(x) {return(as.POSIXct(x, tz = "GMT", origin = "1970-01-01") |> lubridate::with_tz(tzone = "GMT"))}


#' Function for getting all hourly timestamps contained in a vector of date
#'
#' Midnight is assigned to the date that is just beginning. Returned times are GMT (as is met Nordic).
#'
#' @param dates vector of dates
#' @export
#'
#' @examples
#' date_hours('2024-01-01')

date_hours = function(dates)
{
  ret_val = c()
  for(dd in dates){
    t0 = as.GMT(as.Date(dd,origin = '1970-01-01'))
    ret_val = c(ret_val, seq(t0,t0 + 23.5*60*60,by = '1 hour'))
  }

  #loops always screw up POSIXct:
  ret_val = as.GMT(ret_val)
  return(ret_val)
}


#' Auxiliary for estimating the size of the written netcdf file
#'
#' Returns an estimate in Megabytes.
#' @param nloc number of loaded gridcells
#' @param ntime number of loaded time-slizes
#' @param nvars number of loaded variables


est_fileout_size = function(nloc,ntime,nvars){
  fs_in_byte = 15967.3 + 8 * nloc + 8 * ntime + 4*nloc*ntime*nvars # This was derived on rectangle-data, but it seems to work good enough for nn as well
  return(fs_in_byte/10^6)
}




#' In which data streams does a time stamp exist
#'
#' @description The analyses datasets are organized in three data streams (latest, operational archive, rerun archive), see [met Nordic documentation](https://github.com/metno/NWPdocs/wiki/MET-Nordic-dataset) for details.
#' These functions check which of the provided timestamps exist in one of the three data streams.
#'
#' @param time A vector of times
#' @return logical of the same length as `time`
#'
#' @examples
#' exists_ts_latest(Sys.Date())
#' exists_ts_operational_archive(Sys.Date())
#' exists_ts_rerun_archive(Sys.Date())
#'
#' exists_ts_latest('2023-05-01')
#' exists_ts_operational_archive('2023-05-01')
#' exists_ts_rerun_archive('2023-05-01')
#'
#' @export

exists_ts_latest = function(time, forecast = FALSE){
  version_char = ifelse(forecast, yes = 'forecast', no = 'analysis')
  fns = which_files_exist_latest()
  analysis_time_stamps = fns[grepl(version_char,fns)] |> tail(-1) # first element is always "latest", which we don't deal with right now
  analysis_time_stamps = gsub(paste0("metpplatest/met_",version_char,"_1_0km_nordic_"),"",analysis_time_stamps) |> gsub(pattern = '.nc',replacement = '')
  analysis_time_stamps = as.POSIXct(analysis_time_stamps, tz = 'GMT',format = "%Y%m%dT%HZ")
  return(time %in% analysis_time_stamps)
}

#' @rdname exists_ts_latest
#' @export

exists_ts_operational_archive = function(time){
  time = as.GMT(time)
  min_time = as.GMT("2018-03-01")
  max_time = as.GMT(Sys.Date()-2)
  return(time %between% c(min_time,max_time) | exists_ts_latest(time)) # We assume that, for time stamps 2 days ago or less from today, the same timestamps are available in latest and operational archive
                                                                       # (seems to be the case).
}

#' @rdname exists_ts_latest
#' @export

exists_ts_rerun_archive = function(time){
  time = as.GMT(time)
  min_time = as.GMT("2012-09-01")
  max_time = as.GMT("2023-01-31 23:00:00")
  return(time %between% c(min_time,max_time))
}


#' Get filenames of Met Nordic netcdfs
#'
#' Given a time-range or a vector of times, this function returns a vector of all met Nordic filenames/paths corresponding to those times.
#' Times are converted to POSIXct using \code{\link{as.GMT}}. In particular, character strings such as `"2024-01-01 00:00:00"` are interpreted to be in GMT(=UTC), which is the timezone used by met Nordic.
#' This differs from `as.POSIXct` which by default uses the timezone of the machine it is running on.
#'
#'@param time_start,time_end the time stamps you want to extract data for. If `time_end` is provided, an hourly sequence from (the lowest entry of) `time_start` to `time_end` is considered. Alternatively,
#' a vector of hourly time stamps can be provided as `time_start` (if `time_end == NULL`). If a vector of dates is provided, all hourly time stamps contained in these dates are used (using UTC timezone, see \code{\link{date_hours}})
#' @param use_rerun logical. Decides whether the operational archive or the rerun archive should be preferred (for timestamps for which both is available),
#' see the [met Nordic data documentation](https://github.com/metno/NWPdocs/wiki/MET-Nordic-dataset) for details.
#' If `TRUE` (the default), the rerun archive is used (which uses one consistent analysis), else the operational archive is used.
#'
#'@examples
#'# time range:
#' fns1 = get_analysis_fns(time_start = "2023-09-01 01:23:45", time_end =  "2023-09-30")
#'
#'# selected times:
#' times = c(as.GMT('2024-01-03 01:23:45'), as.GMT('2024-02-03 00:12:34'))
#' fns2 = get_analysis_fns(time_start = times)
#'
#'# selected dates:
#'
#' dates = seq(as.Date('2024-01-01'), Sys.Date(), by = 1)
#'fns = get_analysis_fns(time_start = dates) # ALL time stamps for these dates are loaded
#'                                           # (not only the midnight time stamps,
#'                                           # which would be the case if the dates
#'                                           # are simply converted to times).
#'
#'@export


get_analysis_fns = function(time_start, time_end = NULL, use_rerun = TRUE)
{
  # we have to deal with dates differently than with times, but both can be provided as character.
  # That means we have to recognize whether a character is date-character or time-character:
  is_date = function(d){
    return((lubridate::is.Date(d[1]) | (is.character(d[1]) & nchar(d[1]) <= 10)))
  }

  # dates are sometimes given as characters
  # get an hourly time-vector:

  # for dates get all time-stamp during those dates:
  if(is.null(time_end) & is_date(time_start)){
    times = date_hours(time_start)
  } else if(is.null(time_end)){
    times = as.GMT(time_start)
    if(! all(c(second(times),minute(times)) == 0)){
      warning("The times contain non-zero minutes and/or seconds, but the data is hourly. Times are rounded to the full hour.")
      times = round(times,"hour")
    }
  } else {
    time_start = as.GMT(time_start)[1]
    if(is_date(time_end)){
      time_end = tail(date_hours(time_end),1)
    } else{
      time_end = as.GMT(time_end)
      }
    times = seq.POSIXt(time_start, time_end, by = "1 hour")
    if(! all(c(second(times),minute(times),second(time_end),minute(time_end)) == 0)){
      warning("The times contain non-zero minutes and/or seconds, but the data is hourly. Times are rounded down to the next-lower full hour.")
    }
  }

  times = sort(unique(times))

  latest = exists_ts_latest(times)
  operational_archive = exists_ts_operational_archive(times)
  rerun_archive = exists_ts_rerun_archive(times)

  operational_archive[which(latest)] = FALSE # prefer 'latest' over 'operational archive'

  # should we use rerun archive or operational archive if both are available:
  if(use_rerun) operational_archive[which(rerun_archive)] = FALSE else rerun_archive[which(operational_archive)] = FALSE

  # generate loading message:
  mstr = "The following datastreams are used:\n"
  if(!all(latest | operational_archive | rerun_archive)){
    mstr = paste0(mstr,"< 2012-09-01: No data available\n")
  }
  if(any(rerun_archive)){
    if(use_rerun) mstr = paste0(mstr,"2012-09-01 to 2023-01-31: rerun archive version 3\n")
    if(!use_rerun) mstr = paste0(mstr,"2012-09-01 to 2018-02-28: rerun archive version 3\n")
  }
  if(any(operational_archive)){
    if(use_rerun) mstr = paste0(mstr,"2023-02-01 to ",Sys.Date()-3,": operational archive\n")
    if(!use_rerun) mstr = paste0(mstr,"2018-03-01 to ",Sys.Date()-3,": operational archive\n")
  }
  if(any(latest)) mstr = paste0(mstr,Sys.Date()-2," to ",Sys.Date(),": operational real-time\n")

  message(mstr)

  if(any(!latest & !operational_archive & !rerun_archive)) {
    warning("For some of the provided timestamps I could not find Data. I am downloading what's there.")
  }

  # for correct string-format of months,days,hours:
  number_str = function(numbers)
  {
    ret_val = as.character(numbers)
    ret_val[numbers < 10] = paste0(0,ret_val[numbers < 10])
    return(ret_val)
  }

  fns = c()
  fn0 =  "https://thredds.met.no/thredds/dodsC/"
  for(i in seq_along(times))
  {
    if(!latest[i] & !operational_archive[i] & !rerun_archive[i]) next
    timestampstr = paste0(year(times[i]),number_str(month(times[i])),number_str(mday(times[i])),'T',number_str(hour(times[i])),'Z')
    if(latest[i]) fn = paste0(fn0,"metpplatest/met_analysis_1_0km_nordic_",timestampstr,".nc")
    if(operational_archive[i]) fn = paste0(fn0,"metpparchive/",year(times[i]),"/",number_str(month(times[i])),"/",number_str(mday(times[i])),"/met_analysis_1_0km_nordic_",timestampstr,".nc")
    if(rerun_archive[i]) fn = paste0(fn0,"metpparchivev3/",year(times[i]),"/",number_str(month(times[i])),"/",number_str(mday(times[i])),"/met_analysis_1_0km_nordic_",timestampstr,".nc")

    fns = c(fns,fn)
  }

  return(fns)
}


#' Function for matching coordinates to the grid used in the met Nordic data
#'
#' This function takes a set of coordinates and finds the corresponding met-Nordic gridcells.
#' For this we use approximate nearest neighbor using fast kd-tree search. This search only works with Euclidean distance,
#' which introduces error because the x-y coordinates in metno are only given in lon/lat (and d lon is not the same as d lat, especially off the equator).
#' This is counteracted by transforming longitudes to cos(lat) * lon (both for the given coords and the netcdf coordinates), which gives a first order correction
#' and should be pretty close to exact. Everything is back-transformed at the end.
#' The function returns a data table providing the original coordinates and the coordinates of the met Nordic grid cell, as well
#' as the indices (both vector indices and x- and y-indices) of the met Nordic grid cell in the met Nordic netcdfs
#'
#' @param coords A data table with columns `'lon'`, `'lat'`
#' @param max_dist A maximum distance (in km) to be away from the center of the met-no gridpoint. In theory, the grid is 1km x 1km, so no location within the covered area
#' should be more than 1/sqrt(2) = 0.7 km from its grid center. But since distance calculation is approximate, we put the default to 1.5 for some tolerance. This has been loosely verified
#' to be a good distance by plotting, so only change this if you know what you're doing.
#'
#'
#' @return A data table with columns `'lon'`, `'lat'` (same as in `coords`), `'lon_metno'`, `'lat_metno'` (center coordinates of corresponding met Nordic gridcell),
#' `'x_ind'`, `'y_ind'`, `'vector_ind'` (indices of the gridcell within the metno netcdfs).
#'
#' @examples
#' # met gridcell for Oslo city hall:
#' get_metno_gridcells_nn(data.table(lat = 59.912,lon = 10.733))
#'
#' @export

get_metno_gridcells_nn = function(coords,max_dist = 1.5)
{
  coords = unique(coords[,.(lon,lat)]) |> setkey(lon,lat)
  sg_nc = ncdf4::nc_open(system.file("extdata", "spatial_grid.nc", package="metnodi"))

  coords_transformed = copy(coords)
  coords_transformed[,lon := sin(2*pi*lat/90)*lon]

  lons_nc = ncdf4::ncvar_get(sg_nc,varid = 'lon')
  lats_nc = ncdf4::ncvar_get(sg_nc,varid = 'lat')

  lons_nc_transformed = sin(2*pi*lats_nc/90)*lons_nc

  ll_nc = cbind(as.vector(lons_nc_transformed), as.vector(lats_nc))

  nns =RANN::nn2(data = ll_nc, query = coords_transformed[, .(lon,lat)], k = 1)
  vector_indices = nns$nn.idx
  dists = nns$nn.dists * 111.1 # approximate distance in km

  matrix_indices = arrayInd(ind = vector_indices, .dim = c(sg_nc$dim$x$len, sg_nc$dim$y$len)) |> data.table() |> setnames(old = c('x_ind','y_ind'))

  ret_dt = coords[,.(lon,lat)] |> cbind(matrix_indices)
  ret_dt[,vector_ind := vector_indices]
  ret_dt[,lon_metno := lons_nc[vector_ind]]
  ret_dt[,lat_metno := lats_nc[vector_ind]]

  if(any(dists > max_dist)){
    warning(paste0(sum(dists > max_dist)," of the provided coordinates are farther than ",max_dist," km from the nearest gridcell-center. They are not matched to a coordinate."))
    ret_dt[(as.numeric(dists) > max_dist),c('x_ind','y_ind','vector_ind','lon_metno','lat_metno') := NA]
  }

  return(ret_dt[])
}

#' Function for getting met Nordic gridcells in a rectangle area
#'
#' Given a bounding box of the format `lon_min`, `lat_min`, `lon_max`, `lat_max`, find all met Nordic gridcells between those boundaries.
#'
#' @param bbox numeric vector of length four giving the boundaries in the order `lon_min`, `lat_min`, `lon_max`, `lat_max`.
#'
#'
#' @return A data table with columns `'lon_metno'`, `'lat_metno'` (center coordinates of corresponding met Nordic gridcell),
#' `'x_ind'`, `'y_ind'`, `'vector_ind'` (indices of the gridcell within the metno netcdfs).
#'
#' @seealso get_metno_gridcells_nn
#' @examples
#' get_metno_gridcells_rect(c(8,58,12,62))
#'
#' @export
get_metno_gridcells_rect = function(bbox)
{
  sg_nc = ncdf4::nc_open(system.file("extdata", "spatial_grid.nc", package="metnodi"))

  lons_nc = ncvar_get(sg_nc,'lon')
  lats_nc = ncvar_get(sg_nc,'lat')

  vector_indices = which(lons_nc %between% bbox[c(1,3)] & lats_nc %between% bbox[c(2,4)])
  if(length(vector_indices) == 0){
    warning("No MET Nordic gridcells contained in area.")
    return(data.table())
  }

  matrix_indices = arrayInd(ind = vector_indices, .dim = c(sg_nc$dim$x$len, sg_nc$dim$y$len)) |> data.table() |> setnames(old = c('x_ind','y_ind'))

  ret_dt = matrix_indices
  ret_dt[,vector_ind := vector_indices]
  ret_dt[,lon_metno := lons_nc[vector_ind]]
  ret_dt[,lat_metno := lats_nc[vector_ind]]

  return(ret_dt[])
}


#' Function to attempt opening a NetCDF file with exponential backoff
#' @param url URL of netcdf file
#' @param max_retries Number of retries. We wait 2^(n-1) secs after the nth try.
#' @export

open_netcdf_from_url_with_backoff = function(url, max_retries = 10) {
  attempt = 1
  nc = NULL
  while (is.null(nc) & attempt <= max_retries) {
    nc = tryCatch({
      # Attempt to open the NetCDF file
      nc = nc_open(url)
      # If successful, return the file handle
      return(nc)
    }, error = function(e) {
      if (attempt == max_retries) {
        # If maximum retries reached, stop with an error
        stop("Failed to open NetCDF file after ", max_retries, " attempts.")
      } else {
        # On failure, print the error and wait for the next attempt
        message(paste("Attempt", attempt, "failed. Error:", e$message))
        delay = 2^(attempt - 1) # Calculate delay with exponential backoff
        message(paste("Waiting", delay, "seconds before retrying..."))
        Sys.sleep(delay) # Wait before retrying
        return(NULL)
      }
    })
    if(is.null(nc)) attempt = attempt + 1
  }
}



#' Separates filenames into named lists
#'
#' Auxiliary function used in the download in order to group downloads by year, month, day, or hour.
#' Takes download filenames and a level to group by, returns a named list.
#'
#' @param fns filenames of download files
#' @param level Either `'year'`, `'month'`, `'day'` or `'hour'`.

separate_fns = function(fns,level){

  ret_list = list()
  if(level == 'year') {
    present_ys = unique(year(time_from_fn(fns)))
    for(i in seq_along(present_ys)){
      yy = present_ys[i]
      fns_temp = fns[year(time_from_fn(fns)) == yy]
      append_list = list(fns_temp)
      names(append_list) = yy
      ret_list = c(ret_list,append_list)
    }
  }
  if(level == 'month') {
    ym_dt = unique(data.table(year = year(time_from_fn(fns)),month = month(time_from_fn(fns))))
    for(i in 1:ym_dt[,.N]){
      yy = ym_dt[i,year]
      mm = ym_dt[i,month]
      fns_temp = fns[year(time_from_fn(fns)) == yy &
                       month(time_from_fn(fns)) == mm]
      append_list = list(fns_temp)
      names(append_list) = paste0(yy,'-',ifelse(mm<10,yes = paste0('0',mm),no = mm))
      ret_list = c(ret_list,append_list)
    }
  }

  if(level == 'day') {
    ymd_dt = unique(data.table(year = year(time_from_fn(fns)),
                               month = month(time_from_fn(fns)),
                               day = mday(time_from_fn(fns))))
    for(i in 1:ymd_dt[,.N]){
      yy = ymd_dt[i,year]
      mm = ymd_dt[i,month]
      dd = ymd_dt[i,day]
      fns_temp = fns[year(time_from_fn(fns)) == yy &
                       month(time_from_fn(fns)) == mm &
                       mday(time_from_fn(fns)) == dd]
      append_list = list(fns_temp)
      names(append_list) = paste0(yy,'-',
                                  ifelse(mm<10,yes = paste0('0',mm),no = mm),'-',
                                  ifelse(dd<10,yes = paste0('0',dd),no = dd))
      ret_list = c(ret_list,append_list)
    }
  }

  if(level == 'hour') {
    ymdh_dt = unique(data.table(year = year(time_from_fn(fns)),
                                month = month(time_from_fn(fns)),
                                day = mday(time_from_fn(fns)),
                                hour = hour(time_from_fn(fns))))
    for(i in 1:ymdh_dt[,.N]){
      yy = ymdh_dt[i,year]
      mm = ymdh_dt[i,month]
      dd = ymdh_dt[i,day]
      hh = ymdh_dt[i,hour]
      fns_temp = fns[year(time_from_fn(fns)) == yy &
                       month(time_from_fn(fns)) == mm &
                       mday(time_from_fn(fns)) == dd &
                       hour(time_from_fn(fns)) == hh ]
      append_list = list(fns_temp)
      names(append_list) = paste0(yy,'-',
                                  ifelse(mm<10,yes = paste0('0',mm),no = mm),'-',
                                  ifelse(dd<10,yes = paste0('0',dd),no = dd),'-',
                                  ifelse(hh<10,yes = paste0('0',hh),no = hh))
      ret_list = c(ret_list,append_list)
    }
  }

  return(ret_list)
}



#' Get a POSIXct from a download-filename
#'
#' @param fns Vector of file names of the online netcdfs
#' @export

time_from_fn = function(fns){
  pos_start = nchar(fns) - 14
  pos_stop = nchar(fns) - 3
  time_str = substr(fns,start = pos_start,stop = pos_stop)
  times = as.POSIXct(time_str,format = "%Y%m%dT%HZ",tz = "GMT")
  return(times)
}


#' Names of available weather variables
#'
#' "precipitation_amount_quantiles" and "precipitation_amount_gt" are excluded - they seem to only be available for some of the analysis products
#' @examples
#' weather_variables()
#' @export

weather_variables = function() return(c('air_pressure_at_sea_level',
                                        'air_temperature_2m',
                                        'cloud_area_fraction',
                                        'integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time',
                                        'precipitation_amount',
                                        'relative_humidity_2m',
                                        'wind_speed_10m',
                                        'wind_direction_10m'))

#' List files in the METPP-latest catalogue
#'
#' used for checking which files can already be downloaded.
#'
#' @import XML
#' @export

which_files_exist_latest = function()
{
  catalogue = xmlParse(httr::GET("https://thredds.met.no/thredds/catalog/metpplatest/catalog.xml"))
  xmllength = xmlSize(xmlRoot(catalogue)[[2]])

  xmlnames = c()
  for(i in 1:xmllength) {
    xmlnames = c(xmlnames, XML::xmlName(XML::xmlRoot(catalogue)[[2]][[i]]))
  }

  fns = c()
  for(i in 1:xmllength){
    if(xmlnames[i] == 'dataset'){
      fns = c(fns,xmlAttrs(xmlRoot(catalogue)[[2]][[i]])['ID'])
    }
  }

  return(unname(fns))
}





#' Auxiliary function
#'
#' Fetches the forecast file name and the time-indices for netcdf subsetting.
#' @param time_start As in [download_metno_analysis_nn()], but you can also put `'all'` for downloading all available future timestamps.
#' @param time_end As in [download_metno_analysis_nn()].
#' @param initialization_time Default is `'latest'`, which will download the most recent forecast. You can specify downloading older forecasts,
#' by giving their initialization time (which is rounded to the full hour and passed to [as.GMT()]).

get_fc_fn_and_times = function(time_start,
                               time_end,
                               initialization_time){

  if(identical(initialization_time, 'latest')) {
    fn = 'https://thredds.met.no/thredds/dodsC/metpplatest/met_forecast_1_0km_nordic_latest.nc'
  } else {
    number_str = function(numbers)
    {
      ret_val = as.character(numbers)
      ret_val[numbers < 10] = paste0(0,ret_val[numbers < 10])
      return(ret_val)
    }

    time = as.GMT(initialization_time)
    time = round(time,"hours")-2
    timestampstr = paste0(year(time),number_str(month(time)),number_str(mday(time)),'T',number_str(hour(time)),'Z')
    fn = paste0("https://thredds.met.no/thredds/dodsC/metpplatest/met_forecast_1_0km_nordic_",timestampstr,".nc")
  }

  nc = nc_open(fn)
  on.exit(nc_close(nc))
  message('Loading the forecast with the following specifications:\n',
          'history: ',ncatt_get(nc,0,"history")$value,'\n',
          'meps_forecast_reference_time: ',ncatt_get(nc,0,'meps_forecast_reference_time')$value |> as.numeric() |> as.GMT(), ' UTC' )


  # Get timestamps:
  if(identical(time_start,'all')) {
    times = ncvar_get(nc,'time')
  } else if(is.null(time_end) & is_date(time_start)){  # for dates get all time-stamp during those dates:
    times = date_hours(time_start)
  } else if(is.null(time_end)){
    times = as.GMT(time_start)
    if(! all(c(second(times),minute(times)) == 0)){
      warning("The times contain non-zero minutes and/or seconds, but the data is hourly. Times are rounded to the full hour.")
      times = round(times,"hour")
    }
  } else {
    time_start = as.GMT(time_start)[1]
    if(is_date(time_end)){
      time_end = tail(date_hours(time_end),1)
    } else{
      time_end = as.GMT(time_end)
    }
    times = seq.POSIXt(time_start, time_end, by = "1 hour")
    if(! all(c(second(times),minute(times),second(time_end),minute(time_end)) == 0)){
      warning("The times contain non-zero minutes and/or seconds, but the data is hourly. Times are rounded down to the next-lower full hour.")
    }
  }

  times = sort(unique(times))
  if('POSIXct' %in% is(times)) times = as.numeric(times, origin = '1970-01-01')

  time_inds = match(times,ncvar_get(nc,'time'))
  if(any(is.na(time_inds))) {
    warning('Some requested times are not part of the forecast file')
    time_inds = time_inds[!is.na(time_inds)]
  }
  if(length(time_inds) == 0) stop('Nothing to load')

  time_start = min(time_inds)
  time_length = max(time_inds) - min(time_inds) + 1
  time_which = time_inds - min(time_inds) + 1


  return(list(fn = fn,
              time_start = time_start,
              time_length = time_length,
              time_which = time_which))

}

