devtools::load_all()


coords = data.table(lon = 10,lat = 60)


out_dir = '/nr/samba/PostClimDataNoBackup/sensKlimGJ/MetNo/'
out_fn = 'MetNo_for_SURF_region.nc'

dates = seq(as.Date('2023-01-01'),as.Date('2024-01-01'),by = 1)

fns = get_analysis_fns(time_start = dates)

for(i in 1:2000){
  cat(paste0('\r',i))
  fn = fns[i]
  nc = nc_open(fn,return_on_error = TRUE)
  nc_close(nc)
}

metnodi::download_metno_analysis_nn(file_out = out_fn,
                                    out_dir = out_dir,
                                    coords = coords,
                                    time_start = dates,
                                    vars = 'precipitation_amount')
