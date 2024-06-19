devtools::load_all()
library(ncdf4)

wfel = which_files_exist_latest()

nc_ana_latest_ts = nc_open("https://thredds.met.no/thredds/dodsC/metpplatest/met_analysis_1_0km_nordic_20240618T10Z.nc")
nc_ana_latest = nc_open("https://thredds.met.no/thredds/dodsC/metpplatest/met_analysis_1_0km_nordic_latest.nc")
nc_fc_latest_ts = nc_open("https://thredds.met.no/thredds/dodsC/metpplatest/met_forecast_1_0km_nordic_20240618T10Z.nc")
nc_fc_latest = nc_open("https://thredds.met.no/thredds/dodsC/metpplatest/met_forecast_1_0km_nordic_latest.nc")
ts_analysis = exists_ts_latest()

fc_times = ncvar_get(nc_fc_latest_ts,'time') |> as.POSIXct(origin = "1970-01-01", tz = "UTC")
ana_time = ncvar_get(nc_ana_latest_ts,'time') |> as.POSIXct(origin = "1970-01-01", tz = "UTC")


### latest available analysis timestamp ###

latest_analysis_ts = function(){
  nc_ana_latest = nc_open("https://thredds.met.no/thredds/dodsC/metpplatest/met_analysis_1_0km_nordic_latest.nc")
  on.exit(nc_close(nc_ana_latest))
  ana_time = ncvar_get(nc_ana_latest,'time') |> as.GMT()
  return(ana_time)
}


latest_forecast_ts = function(){
  nc_fc_latest = nc_open("https://thredds.met.no/thredds/dodsC/metpplatest/met_forecast_1_0km_nordic_latest.nc")
  on.exit(nc_close(nc_fc_latest))
  fc_time = ncvar_get(nc_fc_latest,'time') |> as.GMT()
  return(fc_time)
}



#' Auxiliary function
#'
#' Fetches the forecast file name and the time-indices for netcdf subsetting.
#' @param time_start As in [download_metno_analysis_nn()], but you can also put `'all'` for downloading all possible future timestamps.
#' @param time_end As in [download_metno_analysis_nn()].
#' @param initialization_time Default is `'latest'`, which will download the most recent forecast. You can specify downloading older forecasts.

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


download_metno_analysis_rect('Oslo.nc',out_dir = tempdir(),bbox = c(10,59,11,60),time_start = Sys.Date())

get_analysis_fns(seq(Sys.Date()-3,Sys.Date(),1))

# latest forecast




#'@description
#' These functions download the MET Nordic forecast.
#'
#' The function `download_metno_fc_nn` takes a data table of coordinates and downloads the analysis data for the nearest-neighbor gridcells.
#' The resulting netcdf contains the loaded weather variables and as dimension index a `coord_index` which identifies the coordinate.
#' It moreover has variables `lon`,`lat` (indexed by `coord_index`) for retrieving the actual coordinate, as well as
#' `lon_metno`, `lat_metno`, the center coordinate for the associate met Nordic gridcell.
#' Nearest neighbor matching is performed by [get_metno_gridcells_nn()].
#'
#' `download_metno_fc_rect` downloads the data for a rectangular region specified by `bbox`.
#'
#'
#' @param file_out,out_dir file name and directory for saving the results.
#' @inheritParams get_metno_gridcells_nn
#' @inheritParams get_metno_gridcells_rect
#' @inheritParams get_fc_fn_and_times
#' @param vars character vector of (weather) variables to load, see [weather_variables()] to see all available variables.
#' @param discard_missing_locations logical. If `TRUE` (the default), then coordinates not falling within the met Nordic area are discarded.
#' If you put this to `FALSE`, the function will find a nearest neighbor for all gridcells, even if they are way outside the covered region.
#'
#'
#' @details
#' Unlike for the analysis-download functions, the loaded data is always saved as single file.
#'
#'@examples \dontrun{
#' # Oslo city hall:
#' coords = data.table(lat = 59.912,lon = 10.733)
#'
#' download_metno_fc_nn(file_out = 'oslo_ch_prec_fc.nc',
#'                      out_dir = tempdir(),
#'                      coords = coords,
#'                      vars = 'precipitation_amount')
#'
#'download_metno_fc_rect(file_out = 'oslo_prec_fc.nc',
#'                       out_dir = tempdir(),
#'                       bbox = c(8,58,12,62)
#'                       vars = 'precipitation_amount')
#'
#'}
#'
#' @export

download_metno_fc_nn = function(file_out,
                                out_dir,
                                coords,
                                time_start = 'all',
                                time_end = NULL,
                                initialization_time = 'latest',
                                vars = weather_variables(),
                                discard_missing_locations = TRUE)

{

  file_out = paste0(gsub('.nc','',file_out),'.nc') # for resilience, if filename does not end on .nc
  dir.create(out_dir,showWarnings = FALSE,recursive = TRUE)

  message('getting nearest neighbors...')
  gridcells = get_metno_gridcells_nn(coords = coords)
  if(discard_missing_locations) gridcells = gridcells[!is.na(x_ind)]

  # spatial dimvar and information
  nc_coord_index = ncdim_def(name = 'coord_index',
                             units = 'index of spatial coordinate, see vars lon/lat for actual coordinate.',
                             vals = 1:gridcells[,.N])

  nc_lon = ncvar_def(name = 'lon',units = 'degree longitude',longname = 'longitude of the provided coordinate',dim = nc_coord_index)
  nc_lat = ncvar_def(name = 'lat',units = 'degree latitude',longname = 'latitude of the provided coordinate',dim = nc_coord_index)

  nc_lonmetno = ncvar_def(name = 'lonmetno',units = 'degree longitude',longname = 'longitude of the center of the associated met Nordic gridcell',dim = nc_coord_index)
  nc_latmetno = ncvar_def(name = 'latmetno',units = 'degree latitude',longname = 'latitude of the center of the associated met Nordic gridcell',dim = nc_coord_index)


  # for extracting weather variables at the corresponding coordinates:

  x_start = gridcells[,min(x_ind)]
  x_length = gridcells[,max(x_ind) - min(x_ind) + 1]
  x_which = gridcells[, x_ind - min(x_ind) + 1]

  y_start = gridcells[,min(y_ind)]
  y_length = gridcells[,max(y_ind) - min(y_ind) + 1]
  y_which = gridcells[, y_ind - min(y_ind) + 1]


  fc_fn_and_times = get_fc_fn_and_times(time_start = time_start,time_end = time_end,initialization_time = initialization_time)

  fn = fc_fn_and_times$fn

  #nc = open_netcdf_from_url_with_backoff(fn)
  nc = nc_open(fn)
  on.exit(nc_close(nc))

  ### fix dimvars from first file ###


  # get all required dimension variables
  all_dimvars = names(nc$dim)

  # which ones do we need?
  req_dimvars = c()
  for(vv in vars){
    ldv = length(nc$var[[vv]]$dim)
    for(j in 1:ldv){
      req_dimvars = c(req_dimvars,nc$var[[vv]]$dim[[j]]$name)
    }
  }
  req_dimvars = unique(req_dimvars)


  # get dimvars directly from netcdf:

  for(dv in req_dimvars){
    assign(paste0('nc_',dv),nc$dim[[dv]])
  }

  # truncate space accordingly:
  if("x" %in% req_dimvars){
    nc_x$vals = nc_x$vals[x_start:gridcells[,max(x_ind)]]
    nc_x$len = length(nc_x$vals)

    nc_y$vals = nc_y$vals[y_start:gridcells[,max(y_ind)]]
    nc_y$len = length(nc_y$vals)
  }

  if("time" %in% req_dimvars){
    nc_time$vals = nc_time$vals[fc_fn_and_times$time_start:(fc_fn_and_times$time_start + fc_fn_and_times$time_length) - 1]
    nc_time$len = length(nc_time$vals)

  }

  ### which variables to write? ###

  variable_indices = match(vars,names(nc$var))
  if(any(is.na(variable_indices))){
    na_vars = which(is.na(variable_indices))
    stop(paste(paste0("The variable names '",paste(vars[na_vars],collapse = "', '"),"' do not exist in the forecast file."),
               "Run all_variables() to see all available variable names.",sep = '\n'))
    #vars = vars[-na_vars]
    #variable_indices = variable_indices[-na_vars]
  }

  ### code for getting variable values: ###

  for(var_ind in variable_indices)
  {
    dim_names = c()
    for(ind in seq_along(nc$var[[var_ind]]$dim)){
      dim_names = c(dim_names,nc$var[[var_ind]]$dim[[ind]]$name)
    }

    start = rep(1,length(dim_names))
    start[dim_names == "x"] = x_start
    start[dim_names == "y"] = y_start
    start[dim_names == "time"] = fc_fn_and_times$time_start

    count = rep(-1,length(dim_names))
    count[dim_names == "x"] = x_length
    count[dim_names == "y"] = y_length
    count[dim_names == "time"] = fc_fn_and_times$time_length

    var_vals = ncvar_get(nc, varid = names(nc$var)[var_ind],start = start,count = count,collapse_degen = FALSE)

    time_which = fc_fn_and_times$time_which
    index_mat = matrix(c(rep(x_which,length(time_which)),rep(y_which,length(time_which)),rep(time_which, each = length(x_which))),ncol = 3)

    # read all other dimensions out completely:
    if(length(dim_names) > 3)
    {
      for(ind in 3:length(dim_names)){
        for(j in 1:dim(var_vals)[ind]){
          if(j == 1) index_mat_temp = cbind(index_mat,j)
          if(j>1) index_mat_temp = rbind(index_mat_temp,cbind(index_mat,j))
        }
        index_mat = index_mat_temp
      }
    }


    assign(paste0("var_vals_",vars[match(var_ind,variable_indices)]),value = var_vals[index_mat]) # initialization
    #initialize ncvar:

    dim_names_new = c('coord_index',setdiff(dim_names,c("x","y")))

    dimlist = lapply(X = dim_names_new,FUN = function(x) get(paste0("nc_",x)))

    #define variable and put values:
    assign(paste0('ncvar_',vars[match(var_ind,variable_indices)]),
           value = ncvar_def(name = nc$var[[var_ind]]$name,
                             units = nc$var[[var_ind]]$units,
                             dim = dimlist,
                             missval = nc$var[[var_ind]]$missval))



  }

  # create netcdf:

  message('writing file ',file.path(out_dir,file_out))

  var_list = c(lapply(vars,function(x) get(paste0('ncvar_',x))),list(nc_lon,nc_lat,nc_lonmetno,nc_latmetno))
  nc_out = nc_create(filename = file.path(out_dir,file_out),vars = var_list)

  # put location information:

  ncvar_put(nc = nc_out,
            varid = nc_lon,
            vals = gridcells[,lon])
  ncvar_put(nc = nc_out,
            varid = nc_lat,
            vals = gridcells[,lat])
  ncvar_put(nc = nc_out,
            varid = nc_lonmetno,
            vals = gridcells[,lon_metno])
  ncvar_put(nc = nc_out,
            varid = nc_latmetno,
            vals = gridcells[,lat_metno])

  for(vv in vars){

    ncvar_put(nc = nc_out,
              varid = get(paste0('ncvar_',vv)),
              vals = get(paste0("var_vals_",vv)))
  }

  nc_close(nc_out)
}


#'@rdname download_metno_fc_nn
#'@export

download_metno_fc_rect = function(file_out,
                                  out_dir,
                                  bbox,
                                  time_start = 'all',
                                  time_end = NULL,
                                  initialization_time = 'latest',
                                  vars = weather_variables())

{
  file_out = paste0(gsub('.nc','',file_out),'.nc') # for resilience, if filename does not end on .nc
  dir.create(out_dir,showWarnings = FALSE,recursive = TRUE)

  gridcells = get_metno_gridcells_rect(bbox = bbox)


  # for extracting weather variables at the corresponding coordinates:

  x_start = gridcells[,min(x_ind)]
  x_length = gridcells[,max(x_ind) - min(x_ind) + 1]

  y_start = gridcells[,min(y_ind)]
  y_length = gridcells[,max(y_ind) - min(y_ind) + 1]


  fc_fn_and_times = get_fc_fn_and_times(time_start = time_start,time_end = time_end,initialization_time = initialization_time)

  fn = fc_fn_and_times$fn

  nc = nc_open(fn)
  on.exit(nc_close(nc))
  ### fix dimvars ###

  # get all required dimension variables
  all_dimvars = names(nc$dim)

  # which ones do we need?
  req_dimvars = c()
  for(vv in vars){
    ldv = length(nc$var[[vv]]$dim)
    for(j in 1:ldv){
      req_dimvars = c(req_dimvars,nc$var[[vv]]$dim[[j]]$name)
    }
  }
  req_dimvars = unique(req_dimvars)


  # get dimvars directly from netcdf:

  for(dv in req_dimvars){
    assign(paste0('nc_',dv),nc$dim[[dv]])
  }

  # truncate space accordingly:
  # truncate space accordingly:
  if("x" %in% req_dimvars){
    nc_x$vals = nc_x$vals[x_start:gridcells[,max(x_ind)]]
    nc_x$len = length(nc_x$vals)

    nc_y$vals = nc_y$vals[y_start:gridcells[,max(y_ind)]]
    nc_y$len = length(nc_y$vals)
  }

  if("time" %in% req_dimvars){
    nc_time$vals = nc_time$vals[fc_fn_and_times$time_start:(fc_fn_and_times$time_start + fc_fn_and_times$time_length) - 1]
    nc_time$len = length(nc_time$vals)

  }

  ### which variables to write? ###

  variable_indices = match(vars,names(nc$var))
  if(any(is.na(variable_indices))){
    na_vars = which(is.na(variable_indices))
    stop(paste(paste0("The variable names '",paste(vars[na_vars],collapse = "', '"),"' do not exist in the forecast file."),
               "Run all_variables() to see all available variable names.",sep = '\n'))
  }
  ### code for getting variable values: ###

  # spatial dimvar and information

  nc_lonmetno = ncvar_def(name = 'lonmetno',units = 'degree longitude',
                          longname = 'longitude of the center of the associated met Nordic gridcell',
                          dim = list(nc_x,nc_y))
  nc_latmetno = ncvar_def(name = 'latmetno',units = 'degree latitude',
                          longname = 'latitude of the center of the associated met Nordic gridcell',
                          dim = list(nc_x,nc_y))

  for(var_ind in variable_indices)
  {
    dim_names = c()
    for(ind in seq_along(nc$var[[var_ind]]$dim)){
      dim_names = c(dim_names,nc$var[[var_ind]]$dim[[ind]]$name)
    }

    start = rep(1,length(dim_names))
    start[dim_names == "x"] = x_start
    start[dim_names == "y"] = y_start
    start[dim_names == "time"] = fc_fn_and_times$time_start

    count = rep(-1,length(dim_names))
    count[dim_names == "x"] = x_length
    count[dim_names == "y"] = y_length
    count[dim_names == "time"] = fc_fn_and_times$time_length

    var_vals = ncvar_get(nc, varid = names(nc$var)[var_ind],start = start,count = count,collapse_degen = FALSE)

    index_mat = c(rep(1:x_length, times = y_length*length(time_which)),
                  rep(1:y_length, each = x_length) |> rep(times = length(time_which)),
                  rep(time_which, each = x_length * y_length)) |> matrix(ncol = 3)

    # read all other dimensions out completely:
    if(length(dim_names) > 3)
    {
      for(ind in 3:length(dim_names)){
        for(j in 1:dim(var_vals)[ind]){
          if(j == 1) index_mat_temp = cbind(index_mat,j)
          if(j>1) index_mat_temp = rbind(index_mat_temp,cbind(index_mat,j))
        }
        index_mat = index_mat_temp
      }
    }


    assign(paste0("var_vals_",vars[match(var_ind,variable_indices)]),value = var_vals[index_mat]) # initialization
    #initialize ncvar:

    dimlist = lapply(X = dim_names,FUN = function(x) get(paste0("nc_",x)))

    #define variable and put values:
    assign(paste0('ncvar_',vars[match(var_ind,variable_indices)]),
           value = ncvar_def(name = nc$var[[var_ind]]$name,
                             units = nc$var[[var_ind]]$units,
                             dim = dimlist,
                             missval = nc$var[[var_ind]]$missval))
  }

  # extract coordinates once:
  lonmetnos = ncvar_get(nc = nc,varid = 'longitude',start = c(x_start,y_start),count = c(x_length,y_length))
  latmetnos = ncvar_get(nc = nc,varid = 'latitude',start = c(x_start,y_start),count = c(x_length,y_length))

    # create netcdf:

  message('writing file ',file.path(out_dir,file_out))

  var_list = c(lapply(vars,function(x) get(paste0('ncvar_',x))),list(nc_lonmetno,nc_latmetno))
  nc_out = nc_create(filename = file.path(out_dir,file_out),vars = var_list)

  # put location information:

  ncvar_put(nc = nc_out,
            varid = nc_lonmetno,
            vals = lonmetnos) # use latest open netcdf, should be fine because they should all be the same
  ncvar_put(nc = nc_out,
            varid = nc_latmetno,
            vals = latmetnos)

  for(vv in vars){

    ncvar_put(nc = nc_out,
              varid = get(paste0('ncvar_',vv)),
              vals = get(paste0("var_vals_",vv)))
  }

  ncatt_put(nc = nc_out, 0, attname = 'proj4str',
            attval = '+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +pm=0 +no_defs +units=m')

  nc_close(nc_out)
}

download_metno_fc_rect('temp',tempdir(),c(9,59,10,60),vars = 'precipitation_amount')


nc = nc_open('/tmp/Rtmp8CYubu/temp.nc')

test = ncvar_get(nc, 'precipitation_amount')

