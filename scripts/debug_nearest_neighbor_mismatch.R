rm(list = ls())
devtools::load_all()

# the nearest neighbor function does not load the right data - why's that?
# we realized that for hans date, so let's check that:



out_dir = '/nr/samba/PostClimDataNoBackup/sensKlimGJ/MetNo/'
out_fns = 'MetNo_Hans_test'

hans_date = '2023-08-08'

N = 1e6
no_lims = c(4,57,32,72)

coords = data.table(lon = runif(N,min = no_lims[1],max = no_lims[3]),
                    lat = runif(N,min = no_lims[2],max = no_lims[4]))


file_out = out_fns
out_dir = out_dir
coords = coords
time_start = hans_date
vars = 'precipitation_amount'



time_end = NULL
use_rerun = TRUE
discard_missing_locations = TRUE
save_resolution = 'auto'
max_filesize_in_MB = 1000
fns = NULL
nfiles = NULL
start_at = 0


file_out = paste0(gsub('.nc','',file_out),'.nc') # for resilience, if filename does not end on .nc
dir.create(out_dir,showWarnings = FALSE,recursive = TRUE)

message('getting nearest neighbors...')
gridcells = get_metno_gridcells_nn(coords = coords)

######## check distances #########
dists = geosphere::distHaversine(p1 = as.matrix(gridcells[,.(lon,lat)]),p2 = as.matrix(gridcells[,.(lon_metno,lat_metno)]))
# distribution of distances makes sense, but we get numbers higher than 1000/sqrt(2)
gridcells[,dists:=dists]
hist(gridcells[,dists])

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















##############################################
dt = SeaVal::netcdf_to_dt('/nr/samba/PostClimDataNoBackup/sensKlimGJ/MetNo/MetNo_Hans_test.nc')
dt[,c('lon','lat'):=NULL]
setnames(dt,c('lonmetno','latmetno'),c('lon','lat'))

plot_dt = unique(dt[,.(time,precipitation_amount,lon,lat)])[ ,.(prec = sum(precipitation_amount)),by = .(lon,lat)]

