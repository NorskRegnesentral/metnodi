
#' Download Met Nordic data for selected coordinates
#'
#' This function takes a data table of coordinates and downloads the met nordic hourly analysis for these locations and a specified timerange and saves it as netcdf.
#' The resulting netcdf contains the loaded weather variables and as dimension index a `coord_index` which identifies the coordinate.
#' It moreover has variables `lon`,`lat` (indexed by `coord_index`) for retrieving the actual coordinate, as well as
#' `lon_metno`, `lat_metno`, the center coordinate for the associate met Nordic gridcell.
#' Nearest neighbor matching is performed by \code{\link{get_metno_gridcells}}.
#'
#' @param file_out,out_dir file name and directory for saving the results. `file_out` needs to end on `'.nc'`.
#' @inheritParams get_analysis_fns
#' @inheritParams get_metno_gridcells_nn
#' @param vars character vector of (weather) variables to load, see \code{\link{weather_variables}} to see all available variables.
#' @param discard_missing_locations logical. If TRUE (the default), then coordinates not falling within the met Nordic area are discarded.
#'
#'@examples \dontrun{
#' # Oslo city hall:
#' data.table(lat = 59.912,lon = 10.733)
#'
#' download_metno_analysis_nn(file_out = 'oslo_ch_prec_20240101.nc',
#'                            out_dir = tempdir(),
#'                            coords = coords,
#'                            time_start = 2024-01-01,
#'                            vars = 'precipitation_amount')
#'
#'}
#'
#' @export

download_metno_analysis_nn = function(file_out,
                                      out_dir,
                                      coords,
                                      time_start,
                                      time_end = NULL,
                                      use_rerun = TRUE,
                                      vars = weather_variables(),
                                      discard_missing_locations = TRUE)

{

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

  # file names of met Nordic files we have to look at:

  fns = get_analysis_fns(time_start = time_start,
                         time_end = time_end,
                         use_rerun = use_rerun)


  # warn when in interactive mode and download will take long:
  if(interactive() & (length(fns) * length(vars) > 500)){
    choices = menu(choices = c('continue','abort'),title = 'The download will probably take longer than 5 minutes. Do you want to continue?')
    if(choices != 1) stop('download aborted.')
  }


  message('loading files:')
  for(i in seq_along(fns))
  {
    cat(paste0('\r',i,'/',length(fns)))

    fn = fns[i]

    nc = nc_open(fn)

    if(i == 1){
      ### which variables to write? ###

      variable_indices = match(vars,names(nc$var))
      if(any(is.na(variable_indices))){
        na_vars = which(is.na(variable_indices))
        warning(paste0("The variable names '",paste(vars[na_vars],collapse = "', '"),"' do not exist.\n Run all_variables() to see all variable names."))
        vars = vars[-na_vars]
        variable_indices = variable_indices[-na_vars]
      }

      # get all required dimension variables
      all_dimvars = names(nc$dim)

      # which ones do we need?
      req_dimvars = c()
      for(ind in variable_indices){
        ldv = length(nc$var[[ind]]$dim)
        for(j in 1:ldv){
          req_dimvars = c(req_dimvars,nc$var[[ind]]$dim[[j]]$name)
        }
      }
      req_dimvars = unique(req_dimvars)

      new_dimvars = c('coord_index',setdiff(req_dimvars,c('x','y')))


      # get other dimvars directly from netcdf:

      for(dv in setdiff(req_dimvars,c('x','y'))){
        assign(paste0('nc_',dv),nc$dim[[dv]])
      }

      if("time" %in% req_dimvars) {
        nc_time$len = length(fns)
        time_vals = nc_time$vals
      }

    }

    # append times:
    if(i > 1 & "time" %in% req_dimvars)
    {
      time_vals = c(time_vals,nc$dim$time$vals)
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

      count = rep(-1,length(dim_names))
      count[dim_names == "x"] = x_length
      count[dim_names == "y"] = y_length

      var_vals = ncvar_get(nc, varid = names(nc$var)[var_ind],start = start,count = count,collapse_degen = FALSE)

      index_mat = matrix(c(x_which,y_which),ncol = 2)

      # read all other dimensions out completely:
      if(length(dim_names) > 2)
      {
        for(ind in 3:length(dim_names)){
          for(j in 1:dim(var_vals)[ind]){
            if(j == 1) index_mat_temp = cbind(index_mat,j)
            if(j>1) index_mat_temp = rbind(index_mat_temp,cbind(index_mat,j))
          }
          index_mat = index_mat_temp
        }
      }

      if(i == 1) {# initialize:
        assign(paste0("var_vals_",var_ind),value = var_vals[index_mat]) # initialization
        #initialize ncvar:

        dim_names_new = c('coord_index',setdiff(dim_names,c("x","y")))

        dimlist = lapply(X = dim_names_new,FUN = function(x) get(paste0("nc_",x)))

        #define variable and put values:
        assign(paste0('ncvar_',var_ind), value = ncvar_def(name = nc$var[[var_ind]]$name,
                                                           units = nc$var[[var_ind]]$units,
                                                           dim = dimlist,
                                                           missval = nc$var[[var_ind]]$missval))


      }
      if(i > 1) {
        #only append variable values:
        assign(paste0("var_vals_",var_ind),value = c(get(paste0("var_vals_",var_ind)),var_vals[index_mat]))
      }

    }

  }



  # fix time dimension variable:

  if("time" %in% req_dimvars){
    nc_time$vals = time_vals
  }


  for(var_ind in variable_indices){

    #fix time dimension variable
    if("time" %in% req_dimvars){
      temp = get(paste0('ncvar_',var_ind))
      # which dim is time
      time_ind = 0
      for(j in seq_along(temp$dim)){
        if(temp$dim[[j]]$name == 'time') time_ind = j
      }

      # overwrite
      temp$dim[[time_ind]]$vals = time_vals
      assign(paste0('ncvar_',var_ind),temp)
    }
  }

  # create netcdf:

  message('writing file ',file.path(out_dir,file_out))

  var_list = c(lapply(variable_indices,function(x) get(paste0('ncvar_',x))),list(nc_lon,nc_lat,nc_lonmetno,nc_latmetno))
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

  for(var_ind in variable_indices){

    ncvar_put(nc = nc_out,
              varid = get(paste0('ncvar_',var_ind)),
              vals = get(paste0("var_vals_",var_ind)))
  }

  nc_close(nc_out)
}




#' Download Met Nordic data for rectangle region
#'
#' This function takes a data table of coordinates and downloads the met nordic hourly analysis for these locations and a specified timerange and saves it as netcdf.
#' The resulting netcdf contains the loaded weather variables and as dimension index a `coord_index` which identifies the coordinate.
#' It moreover has variables `lon`,`lat` (indexed by `coord_index`) for retrieving the actual coordinate, as well as
#' `lon_metno`, `lat_metno`, the center coordinate for the associate met Nordic gridcell.
#' Nearest neighbor matching is performed by \code{\link{get_metno_gridcells}}.
#'
#' @param file_out,out_dir file name and directory for saving the results. `file_out` needs to end on `'.nc'`.
#' @inheritParams get_analysis_fns
#' @inheritParams get_metno_gridcells_rect
#' @param vars character vector of (weather) variables to load, see \code{\link{weather_variables}} to see all available variables.
#'
#'@examples \dontrun{
#'
#' download_metno_analysis_rect(file_out = 'oslo_prec_20240101.nc',
#'                               out_dir = tempdir(),
#'                               bbox = c(8,58,12,62),
#'                               time_start = 2024-01-01,
#'                               vars = 'precipitation_amount')
#'
#'}
#'
#' @export

download_metno_analysis_rect = function(file_out,
                                        out_dir,
                                        bbox,
                                        time_start,
                                        time_end = NULL,
                                        use_rerun = TRUE,
                                        vars = weather_variables())

{
  gridcells = get_metno_gridcells_rect(bbox = bbox)


  # for extracting weather variables at the corresponding coordinates:

  x_start = gridcells[,min(x_ind)]
  x_length = gridcells[,max(x_ind) - min(x_ind) + 1]

  y_start = gridcells[,min(y_ind)]
  y_length = gridcells[,max(y_ind) - min(y_ind) + 1]


  # file names of met Nordic files we have to look at:

  fns = get_analysis_fns(time_start = time_start,
                         time_end = time_end,
                         use_rerun = use_rerun)


  # warn when in interactive mode and download will take long:
  if(interactive() & (length(fns) * length(vars) > 500)){
    choices = menu(choices = c('continue','abort'),title = 'The download will probably take longer than 5 minutes. Do you want to continue?')
    if(choices != 1) stop('download aborted.')
  }

  message('loading files:')
  for(i in seq_along(fns))
  {
    cat(paste0('\r',i,'/',length(fns)))

    fn = fns[i]

    nc = nc_open(fn)

    if(i == 1){
      ### which variables to write? ###

      variable_indices = match(vars,names(nc$var))
      if(any(is.na(variable_indices))){
        na_vars = which(is.na(variable_indices))
        warning(paste0("The variable names '",paste(vars[na_vars],collapse = "', '"),"' do not exist.\n Run all_variables() to see all variable names."))
        vars = vars[-na_vars]
        variable_indices = variable_indices[-na_vars]
      }

      # get all required dimension variables
      all_dimvars = names(nc$dim)

      # which ones do we need?
      req_dimvars = c()
      for(ind in variable_indices){
        ldv = length(nc$var[[ind]]$dim)
        for(j in 1:ldv){
          req_dimvars = c(req_dimvars,nc$var[[ind]]$dim[[j]]$name)
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

      if("time" %in% req_dimvars) {
        nc_time$len = length(fns)
        time_vals = nc_time$vals
      }

    }

    # append times:
    if(i > 1 & "time" %in% req_dimvars)
    {
      time_vals = c(time_vals,nc$dim$time$vals)
    }


    ### code for getting variable values: ###

    # spatial dimvar and information

    nc_lonmetno = ncvar_def(name = 'lonmetno',units = 'degree longitude',longname = 'longitude of the center of the associated met Nordic gridcell',dim = list(nc_x,nc_y))
    nc_latmetno = ncvar_def(name = 'latmetno',units = 'degree latitude',longname = 'latitude of the center of the associated met Nordic gridcell',dim = list(nc_x,nc_y))

    for(var_ind in variable_indices)
    {
      dim_names = c()
      for(ind in seq_along(nc$var[[var_ind]]$dim)){
        dim_names = c(dim_names,nc$var[[var_ind]]$dim[[ind]]$name)
      }

      start = rep(1,length(dim_names))
      start[dim_names == "x"] = x_start
      start[dim_names == "y"] = y_start

      count = rep(-1,length(dim_names))
      count[dim_names == "x"] = x_length
      count[dim_names == "y"] = y_length

      var_vals = ncvar_get(nc, varid = names(nc$var)[var_ind],start = start,count = count,collapse_degen = FALSE)


      if(i == 1) {# initialize:
        assign(paste0("var_vals_",var_ind),value = var_vals) # initialization
        #initialize ncvar:

        dimlist = lapply(X = dim_names,FUN = function(x) get(paste0("nc_",x)))

        #define variable:
        assign(paste0('ncvar_',var_ind), value = ncvar_def(name = nc$var[[var_ind]]$name,
                                                           units = nc$var[[var_ind]]$units,
                                                           dim = dimlist,
                                                           missval = nc$var[[var_ind]]$missval))


      }
      if(i > 1) {
        #only append variable values:
        assign(paste0("var_vals_",var_ind),value = c(get(paste0("var_vals_",var_ind)),var_vals))
      }

    }

  }

  # fix time dimension variable:

  if("time" %in% req_dimvars){
    nc_time$vals = time_vals
  }


  for(var_ind in variable_indices){

    #fix time dimension variable (stored as part of the var)
    if("time" %in% req_dimvars){
      temp = get(paste0('ncvar_',var_ind))
      # which dim is time
      time_ind = 0
      for(j in seq_along(temp$dim)){
        if(temp$dim[[j]]$name == 'time') time_ind = j
      }

      # overwrite
      temp$dim[[time_ind]]$vals = time_vals
      assign(paste0('ncvar_',var_ind),temp)
    }

  }

  # create netcdf:

  message('writing file ',file.path(out_dir,file_out))

  var_list = c(lapply(variable_indices,function(x) get(paste0('ncvar_',x))),list(nc_lonmetno,nc_latmetno))
  nc_out = nc_create(filename = file.path(out_dir,file_out),vars = var_list)

  # put location information:

  ncvar_put(nc = nc_out,
            varid = nc_lonmetno,
            vals = ncvar_get(nc = nc,varid = 'longitude',start = c(x_start,y_start),count = c(x_length,y_length))) # use latest open netcdf, should be fine because they should all be the same
  ncvar_put(nc = nc_out,
            varid = nc_latmetno,
            vals = ncvar_get(nc = nc,varid = 'latitude',start = c(x_start,y_start),count = c(x_length,y_length)))

  for(var_ind in variable_indices){

    ncvar_put(nc = nc_out,
              varid = get(paste0('ncvar_',var_ind)),
              vals = get(paste0("var_vals_",var_ind)))
  }

  nc_close(nc_out)
}

