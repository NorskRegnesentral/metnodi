rm(list = ls())
library(XML)
devtools::load_all()

todays_fns = get_analysis_fns(Sys.Date())

for(i in seq_along(todays_fns)){
fn = todays_fns[i]

#nc = open_netcdf_from_url_with_backoff(fn)
nc = nc_open(fn)
nc_close(nc)
}

which_files_exist_latest = list.files("https://thredds.met.no/thredds/dodsC/metpplatest/")


library(XML)



df <- httr::GET("https://thredds.met.no/thredds/dodsC/metpplatest/")

which_files_exist_latest = function()
{
  catalogue = XML::xmlParse(httr::GET("https://thredds.met.no/thredds/catalog/metpplatest/catalog.xml"))
  xmllength = XML::xmlSize(XML::xmlRoot(catalogue)[[2]])

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


fn = "https://thredds.met.no/thredds/dodsC/metpplatest/met_forecast_1_0km_nordic_latest.nc"

nc = nc_open(fn)
as.GMT(ncvar_get(nc,'time'))

exists_ts_latest = function(time){
  fns = which_files_exist_latest()
  analysis_time_stamps = fns[grepl('analysis',fns)] |> tail(-1) # first element is always "latest", which we don't deal with right now
  analysis_time_stamps = gsub("metpplatest/met_analysis_1_0km_nordic_","",analysis_time_stamps) |> gsub(pattern = '.nc',replacement = '')
  analysis_time_stamps = as.POSIXct(analysis_time_stamps, tz = 'GMT',format = "%Y%m%dT%HZ")
  return(time %in% analysis_time_stamps)
}

catalogue = XML::xmlParse(httr::GET("https://thredds.met.no/thredds/catalog/metpparchive/catalog.xml"))


which_files_exist_latest()

test = XML::xmlToList(httr::GET("https://thredds.met.no/thredds/catalog/metpplatest/catalog.xml"))
is(xmlRoot(test))
xmlSize(test)
xmlName(xmlRoot(test)[[1]])
xmlSize(xmlRoot(test)[[2]])
xmlRoot(test)[[2]][[3]]
xmlRoot(test)[[2]][["dataset"]]


xmllength = xmlSize(xmlRoot(test)[[2]])
xmlnames = c()

for(i in 1:xmllength) {
  xmlnames = c(xmlnames, xmlName(xmlRoot(test)[[2]][[i]]))
}
fns = c()
for(i in 1:xmllength) {
  if(xmlnames[i] == 'dataset'){
    fns = c(fns,xmlSize(xmlRoot(test)[[2]][[i]])[["ID"]])
  }
}

xmlValue(xmlRoot(test)[[2]][[2]])
test

ids = c()
for(i in 1:xmllength){
  if(xmlnames[i] == 'dataset'){
  ids = c(ids,xmlAttrs(xmlRoot(test)[[2]][[i]])['ID'])
  }
}

ids = unname(ids)


### explore forecasts: ###


fn = "https://thredds.met.no/thredds/dodsC/metpplatest/met_forecast_1_0km_nordic_latest.nc"

nc = nc_open(fn)
times = as.GMT(ncvar_get(nc,'time'))
times = lubridate::with_tz(times,'CET')
(lubridate::with_tz(as.GMT(ncvar_get(nc,'time')),'CET') |> tail(1) - Sys.time())*24

fn2 = "https://thredds.met.no/thredds/dodsC/metpplatest/met_forecast_1_0km_nordic_20240423T04Z.nc"
nc2 = nc_open(fn2)

fn3 = "https://thredds.met.no/thredds/dodsC/metpplatest/met_analysis_1_0km_nordic_latest.nc"
nc3 = nc_open(fn3)
as.GMT(ncvar_get(nc3,'time'))


