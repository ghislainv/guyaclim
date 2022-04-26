#!/usr/bin/python

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

# ==== Using google earth engine ======
# library(reticulate)
# https://earthengine.google.com/
# # install miniconda
# # Path to conda binary
# #conda <- "~/.pyenv/versions/miniconda3-latest/bin/conda"
# conda <- "~/miniconda3/bin/conda"
# # Install the python virtual environment
# if (!("conda-gee" %in% conda_list(conda)$name)) {
#     # Vector of packages
#     conda_pkg <- c("python=3.8", "pip", "earthengine-api")
#     # Create conda virtual environment and install packages
#     conda_create("conda-gee", packages=conda_pkg,
#                  forge=TRUE, conda=conda)
# }
# # python=3.8 because ee.Initialize() issue with python=3.10 
# # Specify the Python environment to use
# use_condaenv("conda-gee", required=TRUE, conda=conda)
# # Check python configuration
# py_config()
# 
# # Run python script
# # terminal : 
# # sudo snap install google-cloud-cli --classic
# # earthengine authenticate
# py_run_file(here("Python", "FAO_GAUL_GUY.py"))
# print(py$URL)
# https://earthengine.googleapis.com/v1alpha/projects/earthengine-legacy/tables/b41df7168964496b4375fb8f4b8d3631-3f6b51d518de231ca4fac96a1eca37e6:getFeatures
# Download file
#download.file(py$URL, here("data_raw", "fao_gaul", "FAO_GAUL_GUY.kml"))

#===== Using GADM version 3.6 ========
# https://gadm.org/download_country36.html
#ISO_country_code="GUF" for French Guyana and EPSG=2972
get_extent <- function(ISO_country_code, EPSG, area_borders=NULL, verbose=TRUE, write=TRUE){
  
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
  
  # Bounding box for French Guiana
  xmin <- as.integer((floor(bb_utm$xmin/1000)-5)*1000)
  xmax <- as.integer((ceiling(bb_utm$xmax/1000)+5)*1000)
  ymin <- as.integer((floor(bb_utm$ymin/1000)-5)*1000)
  ymax <- as.integer((ceiling(bb_utm$ymax/1000)+5)*1000)
  
  # Print extent with description (xmin, ymin, xmax, ymax)
  msg <- paste0("Extent of ", ISO_country_code," in EPSG:", EPSG)
  extent <- glue("xmin: {xmin}, ymin: {ymin}, xmax: {xmax}, ymax: {ymax}")
  if(write){writeLines(c(msg, extent), here("output", "extent.txt"))}
  
  # Print extent (xmin, ymin, xmax, ymax)
  extent <- glue("{xmin} {ymin} {xmax} {ymax}")
  if(write){writeLines(extent, here("output", "extent_short.txt"))}
  
return(c(xmin, ymin, xmax, ymax))
}
# End