
#' Download Met Nordic analysis data
#'
#'@description
#' These functions download the MET Nordic analysis data.
#'
#' The function `download_metno_analysis_nn` takes a data table of coordinates and downloads the analysis data for the nearest-neighbor gridcells.
#' The resulting netcdf contains the loaded weather variables and as dimension index a `coord_index` which identifies the coordinate.
#' It moreover has variables `lon`,`lat` (indexed by `coord_index`) for retrieving the actual coordinate, as well as
#' `lon_metno`, `lat_metno`, the center coordinate for the associate met Nordic gridcell.
#' Nearest neighbor matching is performed by [get_metno_gridcells_nn()].
#'
#' `download_metno_analysis_rect` downloads the data for a rectangular region specified by `bbox`.
#'
#'
#' @param file_out,out_dir file name and directory for saving the results.
#' If everything is saved as single file (see details), a .nc file is created in `out_dir`.
#' Else, if multiple output files need to be created, a new directory `file_out` (without .nc) is created in `out_dir`,
#' and the files are saved inside.
#' @inheritParams get_analysis_fns
#' @inheritParams get_metno_gridcells_nn
#' @inheritParams get_metno_gridcells_rect
#' @param vars character vector of (weather) variables to load, see [weather_variables()] to see all available variables.
#' @param discard_missing_locations logical. If TRUE (the default), then coordinates not falling within the met Nordic area are discarded.
#' @param save_resolution Either `"auto"`, `"year"`, `"month"`, `"day"`, or `"hour"`, see details.
#' @param max_filesize_in_MB File size limit (in MB) for the saved netcdfs. Only used if `save_resolution == 'auto'`.
#' @param fns,nfiles,start_at For internal use only. Keep at default.
#'
#'
#' @details
#' The loaded data is saved either as single file or as yearly-, monthly-, daily-, or hourly files.
#' The save format is determined by `save_resolution`. This can be either `year`, `month`, `day`, or `hour`
#' (forcing the corresponding format), or `auto` (the default). When `auto` is selected, the files are saved in
#' the largest format, such that single files remain (approximately) below the max size specified by `max_filesize_in_MB`.
#' In particular, if you want to force to save everything in a single file, you can set `save_resolution = 'auto'` and
#' `max_filesize_in_MB = Inf` (but of course this may crash).
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
#'  download_metno_analysis_rect(file_out = 'oslo_prec_20240101.nc',
#'                               out_dir = tempdir(),
#'                               bbox = c(8,58,12,62),
#'                               time_start = 2024-01-01,
#'                               vars = 'precipitation_amount')
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
                                      discard_missing_locations = TRUE,
                                      save_resolution = 'auto',
                                      max_filesize_in_MB = 1000,
                                      fns = NULL, nfiles = NULL, start_at = 0)

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

  # file names of met Nordic files we have to look at:
  is_null_fns = is.null(fns) # we need that for later
  if(is_null_fns){
    fns = get_analysis_fns(time_start = time_start,
                           time_end = time_end,
                           use_rerun = use_rerun)
    # for progress indication
    nfiles = length(fns)
    start_at = 0
    # warn when in interactive mode and download will take long:
    if(interactive() & (length(fns) * length(vars) > 500)){
      choices = utils::menu(choices = c('continue','abort'),title = 'The download will probably take longer than 5 minutes. Do you want to continue?')
      if(choices != 1) stop('download aborted.')
    }

  }




  # automatic save-resolution selection:
  if(!(tolower(save_resolution) %in% c('year','month','day','hour'))){
    nlocs = gridcells[,.N]
    ntimes = length(fns)
    efs = est_fileout_size(nlocs,ntimes,nvars = length(vars))
    if(efs <= max_filesize_in_MB) {is_null_fns = FALSE} # we keep save_resolution as is, e.g. 'auto'. Thus we will skip the next if-part.
    # Setting is_null_fns to FALSE ensures that we enter the actual download part of the function.
    if(efs > max_filesize_in_MB){  #else
      fs_ymdh = est_fileout_size(nlocs,ntime = c(366*24,31*24,24,1),nvars = length(vars))
      index = which(fs_ymdh < max_filesize_in_MB)[1]
      save_resolution = c('year','month','day','hour')[index]
      message(paste0('saving ',save_resolution,'ly files, else files become too large.'))
    }
  }

  if(tolower(save_resolution) %in% c('year','month','day','hour')){
    fn_list = separate_fns(fns,level = save_resolution)
    new_out_dir = paste0(out_dir,gsub('.nc','',file_out),'/')
    dir.create(new_out_dir,showWarnings = FALSE,recursive = TRUE)
    for(i in seq_along(fn_list)){
      # Only for correct progress-indication
      start_at = 0
      if(i > 1){
        for(j in 1:(i-1)){
          start_at = start_at + length(fn_list[[j]])
        }
      }
      suppressMessages(
        download_metno_analysis_nn(file_out = names(fn_list)[i],
                                   out_dir = new_out_dir,
                                   coords = coords,
                                   time_start = time_start,# does not matter because we provide fns
                                   time_end = NULL,# does not matter because we provide fns
                                   use_rerun = TRUE,# does not matter because we provide fns
                                   vars = vars,
                                   discard_missing_locations = discard_missing_locations,
                                   save_resolution = 'auto', # The configuration of save_resolution and max_filesize_in_MB ensures that
                                   max_filesize_in_MB = Inf, # in this function call we end up in the same if-statement and start cascading down
                                   fns = fn_list[[i]], nfiles = nfiles, start_at = start_at) #In particular, is_null_fns = FALSE
      )
    }
  }


  # rest of the function is only executed if(!is_null_fns)
  if(!is_null_fns){
    message('loading files:')
    for(i in seq_along(fns)){
      cat(paste0('\r',i + start_at,'/',nfiles))

      fn = fns[i]

      nc = open_netcdf_from_url_with_backoff(fn)
      #nc = nc_open(fn)
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
      nc_close(nc) # important for not having too many open files
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
}





#'@rdname download_metno_analysis_nn
#'@export
download_metno_analysis_rect = function(file_out,
                                        out_dir,
                                        bbox,
                                        time_start,
                                        time_end = NULL,
                                        use_rerun = TRUE,
                                        vars = weather_variables(),
                                        save_resolution = 'auto',
                                        max_filesize_in_MB = 1000,
                                        fns = NULL, nfiles = NULL, start_at = 0)

{
  file_out = paste0(gsub('.nc','',file_out),'.nc') # for resilience, if filename does not end on .nc
  dir.create(out_dir,showWarnings = FALSE,recursive = TRUE)

  gridcells = get_metno_gridcells_rect(bbox = bbox)


  # for extracting weather variables at the corresponding coordinates:

  x_start = gridcells[,min(x_ind)]
  x_length = gridcells[,max(x_ind) - min(x_ind) + 1]

  y_start = gridcells[,min(y_ind)]
  y_length = gridcells[,max(y_ind) - min(y_ind) + 1]


  # file names of met Nordic files we have to look at:
  # file names of met Nordic files we have to look at:
  is_null_fns = is.null(fns) # we need that for later
  if(is_null_fns){
    fns = get_analysis_fns(time_start = time_start,
                           time_end = time_end,
                           use_rerun = use_rerun)
    # for progress indication
    nfiles = length(fns)
    start_at = 0
    # warn when in interactive mode and download will take long:
    if(interactive() & (length(fns) * length(vars) > 500)){
      choices = utils::menu(choices = c('continue','abort'),title = 'The download will probably take longer than 5 minutes. Do you want to continue?')
      if(choices != 1) stop('download aborted.')
    }

  }

  # automatic save-resolution selection:
  if(!(tolower(save_resolution) %in% c('year','month','day','hour'))){
    nlocs = gridcells[,.N]
    ntimes = length(fns)
    efs = est_fileout_size(nlocs,ntimes,nvars = length(vars))
    if(efs <= max_filesize_in_MB) {is_null_fns = FALSE} # we keep save_resolution as is, e.g. 'auto'. Thus we will skip the next if-part.
    # Setting is_null_fns to FALSE ensures that we enter the actual download part of the function.
    if(efs > max_filesize_in_MB){  #else
      fs_ymdh = est_fileout_size(nlocs,ntime = c(366*24,31*24,24,1),nvars = length(vars))
      index = which(fs_ymdh < max_filesize_in_MB)[1]
      save_resolution = c('year','month','day','hour')[index]
      message(paste0('saving ',save_resolution,'ly files, else files become too large.'))
    }
  }

  if(tolower(save_resolution) %in% c('year','month','day','hour')){
    fn_list = separate_fns(fns,level = save_resolution)
    new_out_dir = paste0(out_dir,gsub('.nc','',file_out),'/')
    dir.create(new_out_dir,showWarnings = FALSE,recursive = TRUE)
    for(i in seq_along(fn_list)){
      # Only for correct progress-indication
      start_at = 0
      if(i > 1){
        for(j in 1:(i-1)){
          start_at = start_at + length(fn_list[[j]])
        }
      }
      suppressMessages(
        download_metno_analysis_rect(file_out = names(fn_list)[i],
                                     out_dir = new_out_dir,
                                     bbox = bbox,
                                     time_start = time_start,# does not matter because we provide fns
                                     time_end = NULL,# does not matter because we provide fns
                                     use_rerun = TRUE,# does not matter because we provide fns
                                     vars = vars,
                                     save_resolution = 'auto', # The configuration of save_resolution and max_filesize_in_MB ensures that
                                     max_filesize_in_MB = Inf, # in this function call we end up in the same if-statement and start cascading down
                                     fns = fn_list[[i]], nfiles = nfiles, start_at = start_at) #In particular, is_null_fns = FALSE
      )
    }
  }

  # rest of the function is only executed if(!is_null_fns)
  if(!is_null_fns){
    message('loading files:')
    for(i in seq_along(fns)){
      cat(paste0('\r',i + start_at,'/',nfiles))

      fn = fns[i]

      nc = open_netcdf_from_url_with_backoff(fn)
      #nc = nc_open(fn)

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

      # extract coordinates once:
      if(i == 1)
      {
        lonmetnos = ncvar_get(nc = nc,varid = 'longitude',start = c(x_start,y_start),count = c(x_length,y_length))
        latmetnos = ncvar_get(nc = nc,varid = 'latitude',start = c(x_start,y_start),count = c(x_length,y_length))
      }
      nc_close(nc)
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
              vals = lonmetnos) # use latest open netcdf, should be fine because they should all be the same
    ncvar_put(nc = nc_out,
              varid = nc_latmetno,
              vals = latmetnos)

    for(var_ind in variable_indices){

      ncvar_put(nc = nc_out,
                varid = get(paste0('ncvar_',var_ind)),
                vals = get(paste0("var_vals_",var_ind)))
    }

    ncatt_put(nc = nc_out, 0, attname = 'proj4str',
              attval = '+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +pm=0 +no_defs +units=m')

    nc_close(nc_out)
  }
}

