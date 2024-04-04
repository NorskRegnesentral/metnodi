
#' Function for reading downloaded netcdfs as stars
#'
#' Corrects the projection, which is not correctly identified by [stars::read_stars()] and [stars::read_ncdf()]
#' @param ... passed on to [stars::read_ncdf()].
#' @export

read_stars_metno = function(...)
{
  st = stars::read_ncdf(...)
  p4s = '+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +pm=0 +no_defs +units=m'
  sf::st_crs(st) = sf::st_crs(p4s)
  return(st)
}



#' Function for plotting a stars object overlaying google maps
#'
#' This function take a `stars`-object and creates a single map-plot of the data, overlaying a googlemap.
#' Requires google maps API access. If the `stars`-object data for more than 1 time, the first time-slice is plotted.
#' Returns a gg-object, adjust the color-scale by adding a `scale_fill` to it. When calling this many times with the same bbox,
#' you should provide the map, see parameters `map` and `only_map`.
#'
#' @param st Data as `star`.
#' @param bbox Optional bounding box in the form `c(lon_min, lat_min, lon_max, lat_max)` (always in lon/lat, no matter the CRS of `st`).
#' If no bounding box is provided, it is derived from `st`. Values are approximate, see `tol`.
#' @param var_name Name of the variable/attribute in `st` that you want to plot. Default is to plot the first one.
#' @param downsample Should we downsample the data for faster plotting? Default is to downsample to get less than 100.000 pixels.
#' You can either provide a number (0 means no downsampling, 1 means every second value is kept etc.) or a string of the form
#' `'auto<some number n>'` which will downsample to plot fewer than `n`*1000 points, so downsampling rate depends on the size of `st`.
#' @param alpha Alpha-value for the data. Higher values make the data-plot less opaque (making it harder to see the googlemaps).
#' @param tol tolerance value (in lon/lat) for reducing the spatial bounding-box. Only used if `bbox = NULL`. Extend of the googlemap is approximate,
#' so plots frequently look awkward if the full bbox of the data is used, because the map includes areas for which no data is present.
#' @param fix_crs This is used to overwrite/set the CRS for netcdfs with wrongly specified CRS: Set this to `'metno'`
#' for plotting MET Nordic data, and to `'surf'` for plotting SURF data. Else, `st` needs to have the correct CRS
#' @param map Optional. You can provide a `ggmap` object, then it does not need to be retrieved from googlemaps, since you pay per API-ping.
#' You can get the map by calling this function with `only_map = TRUE`.
#' @param only_map Logical. If `TRUE`, only the googlemap is returned.
#' @param ... Passed on to `ggmap::get_map()`
#'
#' @export


gmaps_plot = function(st,
                      bbox = NULL,
                      var_name = NULL,
                      downsample = 'auto100',
                      alpha = 0.7,
                      tol = 0,
                      fix_crs = FALSE,
                      map = NULL,
                      only_map = FALSE,...){
  if(tolower(fix_crs) == 'metno'){
    p4s = '+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +pm=0 +no_defs +units=m'
    suppressWarnings({sf::st_crs(st) = sf::st_crs(p4s)})
  }
  st = sf::st_transform(st,crs = 4326)

  if(is.null(bbox)) bbox = sf::st_bbox(st) + tol*c(1,1,-1,-1)
  names(bbox) = c('left','bottom','right','top')

  if(is.null(map)) map  = ggmap::get_map(bbox,color = 'bw',darken = c(0.5, "white"),...)

  if(only_map){
    pp = ggmap::ggmap(map)
    if(interactive()) plot(pp)
    return(map)
  }

  plot_st = st[,,,1] # first time slice
  if(is.null(var_name)) var_name = names(st)[1]

  if(is.character(downsample)){
    max_pix = readr::parse_number(downsample) * 1000
    downsample = ceiling(sqrt(dim(st)[1]*dim(st)[2]/max_pix))-1
    if(downsample > 0) message('Downsampling = ',downsample)
  }

  pp = suppressMessages({ggmap::ggmap(map) +
    ggplot2::coord_sf(crs = 4326) +
    stars::geom_stars(data = plot_st,
               ggplot2::aes(fill = get(var_name)),
               inherit.aes = FALSE,
               alpha = alpha,
               downsample = c(downsample,downsample,0)) +
      ggplot2::scale_fill_continuous(name = var_name) +
      ggplot2::ggtitle(time(plot_st))})
  if(interactive()) plot(pp)
  return(pp)
}



#' Plot of MET Nordic analysis
#'
#' Generate a plot of MET Nordic analysis for a specific timeslice and a rectangular region.
#' The
#' @param time Time of the slice you want to plot. See [as.GMT()] for details.
#' @param bbox Bounding box of the form `c(lon_min,lat_min,lon_max,lat_max)` for plotting.
#' @param var Name of the variable to plot. See [weather_variables()] for what's available.
#' @param ... passed to [gmaps_plot()]
#'
#' @export
plot_analysis_gmap = function(time,
                             bbox = c(-180,-90,180,90),
                             var = 'precipitation_amount',
                             ...){
  time_start = time
  if(length(time_start) > 1){
    warning(paste0('time_start has length > 1. Only values for ',time_start[1],' are plotted.'))
    time_start = time_start[1]
  }

  # avoid downloading a full date:
  time_start = as.GMT(time_start)
  out_dir = tempdir()

  download_metno_analysis_rect(file_out = 'temp.nc',out_dir = out_dir,bbox = bbox,time_start = time_start,vars = var)
  st = suppressWarnings(stars::read_ncdf(paste0(out_dir,'/temp.nc'),var = var))

  pp = gmaps_plot(st,fix_crs = 'metno',...)

  return(pp)
}

