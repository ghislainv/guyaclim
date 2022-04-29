##===
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# license         :GPLv3
##===

# Libraries
library(here)
library(sf)
library(glue)
get_extent <- function(ISO_country_code, EPSG, area_borders=NULL, verbose=TRUE, write=TRUE, latlon=FALSE){
  
  if(is.null(area_borders)){
    URL <- paste0("https://geodata.ucdavis.edu/gadm/gadm3.6/gpkg/gadm36_", ISO_country_code, "_gpkg.zip")
    # Download file
    # The coordinate reference system is longitude/latitude and the WGS84 datum.
    download.file(URL, quiet = !verbose, 
                  here("data_raw", "fao_gaul", paste0("gpkg_gadm36_", ISO_country_code, ".zip")))
    # Unzip 
    unzip(here("data_raw", "fao_gaul", paste0("gpkg_gadm36_", ISO_country_code, ".zip")),
          exdir=here("data_raw", "fao_gaul"), overwrite=TRUE)
    # Read vector (level 0 for country borders)
    borders <- sf::st_read(here("data_raw", "fao_gaul", paste0("gadm36_", ISO_country_code, ".gpkg")),
                           layer=paste0("gadm36_", ISO_country_code, "_0"), quiet=!verbose)
  } else{
    borders <- area_borders  
  }
  crs_borders <- st_crs(borders)
  bb_ll <- st_bbox(borders, crs=crs_borders)
  bb_utm <- bb_ll %>%
    st_as_sfc() %>%
    st_transform(crs=EPSG) %>%
    st_bbox()
  latmin = bb_ll[2]
  latmax = bb_ll[1]
  longmin = bb_ll[3]
  longmax = bb_ll[4]
  # Bounding box for New Caledonia
  xmin <- as.integer((floor(bb_utm$xmin/1000)-5)*1000)
  xmax <- as.integer((ceiling(bb_utm$xmax/1000)+5)*1000)
  ymin <- as.integer((floor(bb_utm$ymin/1000)-5)*1000)
  ymax <- as.integer((ceiling(bb_utm$ymax/1000)+5)*1000)
  
  # Print extent with description (xmin, ymin, xmax, ymax)
  msg <- paste0("Extent of ", ISO_country_code," in EPSG:", EPSG)
  extent <- glue("xmin: {xmin}, ymin: {ymin}, xmax: {xmax}, ymax: {ymax}")
  extent_ll <- glue("latmin: {latmin}, longmin: {longmin}, latmax: {latmax}, longmax: {longmax}")
  if(write){
    writeLines(c(msg, extent), here("output", "extent.txt"))
    writeLines(c(msg,extent_ll), here("output", "extent_ll.txt"))
    }
  
  # Print extent (xmin, ymin, xmax, ymax)
  extent <- glue("{xmin} {ymin} {xmax} {ymax}")
  bb_ll <- glue("{latmin} {longmin} {latmax} {longmax}")
  if(write){
    writeLines(extent, here("output", "extent_short.txt"))
    writeLines(bb_ll, here("output", "extent_short_latlong.txt"))
    }
  if (latlon) { return(bb_ll)}
  else{ return(c(xmin, ymin, xmax, ymax))}
}
# End