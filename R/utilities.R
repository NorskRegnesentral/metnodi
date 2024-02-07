
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

exists_ts_latest = function(time){
  time = as.GMT(time)
  min_time = as.GMT(Sys.Date() - 2)
  return(time >= min_time)
}

#' @rdname exists_ts_latest
#' @export

exists_ts_operational_archive = function(time){
  time = as.GMT(time)
  min_time = as.GMT("2018-03-01")
  return(time >= min_time)
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


get_analysis_fns = function(time_start,time_end = NULL,use_rerun = TRUE)
{
  # get an hourly time-vector:

  # for dates get all time-stamp during those dates:
  if(is.null(time_end) & lubridate::is.Date(time_start)){
    times = date_hours(time_start)
  } else if(is.null(time_end)){
    times = as.GMT(time_start)
    if(! all(c(second(times),minute(times)) == 0)){
      warning("The times contain non-zero minutes and/or seconds, but the data is hourly. Times are rounded to the full hour.")
      times = round(times,"hour")
    }
  } else {
    time_start = as.GMT(time_start)[1]
    time_end = as.GMT(time_end)
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
    warning("No met Nordic gridcells contained in area.")
    return(data.table())
  }

  matrix_indices = arrayInd(ind = vector_indices, .dim = c(sg_nc$dim$x$len, sg_nc$dim$y$len)) |> data.table() |> setnames(old = c('x_ind','y_ind'))

  ret_dt = matrix_indices
  ret_dt[,vector_ind := vector_indices]
  ret_dt[,lon_metno := lons_nc[vector_ind]]
  ret_dt[,lat_metno := lats_nc[vector_ind]]

  return(ret_dt[])
}


#' Names of available weather variables
#'
#' "precipitation_amount_quantiles" and "precipitation_amount_gt" are excluded - they seem to only be available for some of the analysis products
#' @examples
#' weather_variables()
#' @export

weather_variables = function() return(c('air_pressure_at_sea_level','air_temperature_2m','cloud_area_fraction','integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time','integral_of_surface_downwelling_longwave_flux_in_air_wrt_time',
                                        'precipitation_amount','relative_humidity_2m','wind_speed_10m','wind_direction_10m'))

#' Names of all variables we can potentially get from netcdfs
#' @examples
#' all_variables()
#' @export

all_variables = function() return(sort(c(weather_variables(),"projection_lcc","forecast_reference_time","precipitation_amount_gt","precipitation_amount_quantiles","altitude","land_area_fraction")))


