#!/usr/bin/python

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# license         :GPLv3
# ==============================================================================

# Libraries
library(reticulate)
library(here)
library(sf)
library(glue)

# Path to conda binary
conda <- "~/.pyenv/versions/miniconda3-latest/bin/conda"

# Install the python virtual environment
if (!("conda-gee" %in% conda_list(conda)$name)) {
    # Vector of packages
    conda_pkg <- c("python", "pip", "earthengine-api")
    # Create conda virtual environment and install packages
    conda_create("conda-gee", packages=conda_pkg, forge=TRUE, conda=conda)
}

# Specify the Python environment to use
use_condaenv("conda-gee", required=TRUE, conda=conda)

# Check python configuration
py_config()

# Run python script
py_run_file(here("Python", "FAO_GAUL_GUY.py"))
print(py$URL)
# https://earthengine.googleapis.com/v1alpha/projects/earthengine-legacy/tables/b41df7168964496b4375fb8f4b8d3631-3f6b51d518de231ca4fac96a1eca37e6:getFeatures

# Download file
download.file(py$URL, here("data_raw", "fao_gaul", "FAO_GAUL_GUY.kml"))

# Read vector
guy <- sf::st_read(here("data_raw", "fao_gaul", "FAO_GAUL_GUY.kml"))
bb_ll <- st_bbox(guy, crs=st_crs(4326))
bb_utm <- bb_ll %>%
  st_as_sfc() %>%
  st_transform(crs=32622) %>%
  st_bbox()

# Bounding box for French Guiana
xmin <- (floor(bb_utm$xmin/1000)-5)*1000
xmax <- (ceiling(bb_utm$xmax/1000)+5)*1000
ymin <- (floor(bb_utm$ymin/1000)-5)*1000
ymax <- (ceiling(bb_utm$ymax/1000)+5)*1000

# Print extent (xmin, xmax, ymin, ymax)
msg <- "Extent of French Guyana in UTM32N (epsg: 32622):"
extent <- glue("xmin: {xmin}, xmax: {xmax}, ymin: {ymin}, ymax: {ymax}")
writeLines(c(msg, extent), here("output", "extent.txt"))

# End