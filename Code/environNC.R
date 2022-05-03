##=====================================================
##
## Environmental data for New Caledonia
##
## Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
## Jeanne Clément <jeanne.clement@cirad.fr>
##
## Octobre 2021
##
##=====================================================

## gdal library is needed to run this script
## http://www.gdal.org/

## GRASS GIS 7.x.x is also needed to run this script
## https://grass.osgeo.org/

## Read argument for download
## Set "down" to TRUE if you want to download the sources. Otherwise, the data already provided in the data_raw folder will be used.
arg <- commandArgs(trailingOnly = TRUE)
down <- TRUE
if (length(arg) > 0) {
  down <- arg[1]
}

# Libraries
library(glue)
library(here)
library(sf)
library(stars)
library(rgdal)
library(rgrass7)

## gdalwrap options
Extent <- readLines(here("output/reExtent_short.txt"))
Res <- "1000"
nodat <- "-9999"
proj.s <- "EPSG:4326"
proj.t <- "EPSG:3165"

#====== Elevation, slope aspect, roughness #======
# SRTM at 90m resolution from https://dwtkns.com/srtm/ version 4.1
## Download and unzip CGIAR-CSI 90m DEM data 
# "69_16","70_16",
tiles <- c( "69_17", "70_17")
for (i in 1:length(tiles)) {
  dst <- paste0(here("data_raw", "srtm_v1_4_90m", "srtm_"), tiles[i], ".zip")
  if (down) {
    url.tile <- paste0("https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/srtm_", tiles[i], ".zip")
    download.file(url = url.tile, destfile = dst, method = "wget", quiet = TRUE)
  }
  unzip(dst, exdir = here("data_raw", "srtm_v1_4_90m"), overwrite = TRUE)
}
## Merge srtm with c.stars
sourcefile <- list.files(here("data_raw", "srtm_v1_4_90m"), pattern = "*.tif", full.names = TRUE)
merge.tif <- c(read_stars(sourcefile[1]), read_stars(sourcefile[2]), along = 1)

write_stars(merge.tif, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
            dsn = here("data_raw","srtm_v1_4_90m","temp", "elevation_noExtent.tif"), layer = 1)

## Reproject from lat long to UTM58S (epsg: 3165) and reframe on New Caledonia
sourcefile <- here("data_raw", "srtm_v1_4_90m", "temp", "elevation_noExtent.tif")
destfile <- here("data_raw", "srtm_v1_4_90m", "temp", "elevation.tif")
system(glue("gdalwarp -overwrite -s_srs {proj.s} -t_srs {proj.t}  \\
        -r bilinear -tr 90 90 -ot Int16 -of GTiff \\
        {sourcefile} \\
        {destfile}"))
elev <- destfile
borders = st_bbox(soilgrids)
elev = st_crop(read_stars(elev), borders)
write_stars(elev, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
            dsn = destfile, layer = 1)
border <- st_read(paste(here("data_raw", "fao_gaul"), "gadm36_NCL.gpkg", sep = "/"), layer = "gadm36_NCL_0")
border <- st_transform(border,crs = st_crs(elev))
plot(elev,axes = TRUE, reset = FALSE)
plot(st_geometry(border), add = TRUE, reset = FALSE)


## Compute slope, aspect and roughness using gdaldem 
# compute slope
in_f <- here("data_raw", "srtm_v1_4_90m", "temp", "elevation.tif")
out_f <- here("data_raw", "srtm_v1_4_90m", "temp", "slope.tif")
cmd <- glue('gdaldem slope {in_f} {out_f} -co "COMPRESS=LZW" -co "PREDICTOR=2"')
system(cmd)
# compute aspect
out_f <- here("data_raw", "srtm_v1_4_90m", "temp", "aspect.tif")
cmd <- glue('gdaldem aspect {in_f} {out_f} -co "COMPRESS=LZW" -co "PREDICTOR=2"')
system(cmd)
# compute roughness
out_f <- here("data_raw", "srtm_v1_4_90m", "temp", "roughness.tif")
cmd <- glue('gdaldem roughness {in_f} {out_f} -co "COMPRESS=LZW" -co "PREDICTOR=2"')
system(cmd)
# Resolution from 90m x 90m to 1000m x 1000m using gdalwarp
# elevation
out_f <- here("data_raw", "srtm_v1_4_90m", "elevation_1km.tif")
cmd <- glue('gdalwarp -r bilinear -tr 1000 1000 \\
        -co "COMPRESS=LZW" -co "PREDICTOR=2" -overwrite {in_f} {out_f}')
system(cmd)
# aspect
in_f <- here("data_raw", "srtm_v1_4_90m", "temp", "aspect.tif")
out_f <-here("data_raw", "srtm_v1_4_90m", "aspect_1km.tif")
cmd <- glue('gdalwarp -r bilinear -tr 1000 1000 \\
        -co "COMPRESS=LZW" -co "PREDICTOR=2" -overwrite {in_f} {out_f}')
system(cmd)
# slope
in_f <- here("data_raw", "srtm_v1_4_90m", "temp", "slope.tif")
out_f <- here("data_raw", "srtm_v1_4_90m", "slope_1km.tif")
cmd <- glue('gdalwarp -r bilinear -tr 1000 1000 \\
        -co "COMPRESS=LZW" -co "PREDICTOR=2" -overwrite {in_f} {out_f}')
system(cmd)
# roughness
in_f <- here("data_raw", "srtm_v1_4_90m", "temp", "roughness.tif")
out_f <- here("data_raw", "srtm_v1_4_90m", "roughness_1km.tif")
cmd <- glue('gdalwarp -r bilinear -tr 1000 1000 -ot Int16 -of GTiff \\
        -co "COMPRESS=LZW" -co "PREDICTOR=2" -overwrite {in_f} {out_f}')
system(cmd)

#=== Solar radiation =====
# with r.sun at 90m resolution 
# Solar radiation (in Wh.m-2.day-1) was computed from altitude,
# slope and aspect using the function r.sun from the GRASS GIS software.
# We incorporated the shadowing effect of terrain to compute the solar radiation.
# Solar radiation was computed for the Julian day 79 (20th of March for regular years=equinox).
## Initialize GRASS
setwd(here("data_raw"))
Sys.setenv(LD_LIBRARY_PATH = paste("/usr/lib/grass80/lib", Sys.getenv("LD_LIBRARY_PATH"), sep = ":"))
# use a georeferenced raster
elevation <- here("data_raw", "srtm_v1_4_90m", "temp", "elevation.tif")
system(glue('grass -c {elevation} grassdata/environ'))
# connect to grass database
initGRASS(gisBase = "/usr/lib/grass80", 
          gisDbase = "grassdata", home = tempdir(), 
          location = "environ", mapset = "PERMANENT",
          override = TRUE)
## Import raster in grass
system(glue("r.in.gdal -e --o input={elevation} output=elevation"))
slope <- here("data_raw", "srtm_v1_4_90m", "temp", "slope.tif")
system(glue("r.in.gdal -e --o input={slope} output=slope"))
aspect <- here("data_raw", "srtm_v1_4_90m", "temp", "aspect.tif")
system(glue("r.in.gdal -e --o input={aspect} output=aspect"))
# Compute radiation
cmd <- glue("r.sun --overwrite --verbose elevation=elevation aspect=aspect slope=slope day=79 glob_rad=global_rad nprocs=4")
system(cmd)
# Export
system(glue("r.out.gdal -f --verbose --overwrite input=global_rad \\
  			 output={here('data_raw', 'srtm_v1_4_90m', 'temp', 'srad.tif')} type=Int16  \\
  			 createopt='COMPRESS=LZW' nodata={nodat}"))

# Resolution from 90m x 90m to 1000m x 1000m using gdalwarp
# srad
in_f <- here("data_raw", "srtm_v1_4_90m", "temp", "srad.tif")
out_f <- here("data_raw", "srtm_v1_4_90m", "srad_1km.tif")
# srad = read_stars(in_f)
# srad = st_crop(srad, borders)
# write_stars(srad, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
#            dsn = in_f, layer = 1)

cmd <- glue('gdalwarp -s_srs {proj.t} -t_srs {proj.t} \\
        -r bilinear -tr 1000 1000 -ot Int16 -of GTiff \\
        -co "COMPRESS=LZW" -co "PREDICTOR=2" -overwrite {in_f} {out_f}')
system(cmd)

##==========================================================
##
## Forest: percentage of forest in 1km2
##
##==========================================================

## References:
## (1) Forestatrisk project: https://forestatrisk.cirad.fr/rawdata.html

## Create directory
dir.create(here("data_raw","tmf_ec_jrc"))

# Download forest cover of New Caledonia in 2000 from tmf_ec_jrc (see note)
if (down) {
  download.file("https://drive.google.com/uc?export=download&id=1DBMIVEYRfOFgZFdk96xF16k8vhe15eFj",
                destfile = here("data_raw", "tmf_ec_jrc", "forest_t3.tif"), method = 'auto', mode = "wb")
}
## Reproject from South America Albers Equal Area Conic (ESRI:102033) to UTM32N (EPSG:3165) 
# Resampling from 30x30m to 1000x1000m using sum method 
## Reframe on New Caledonia and set no data values 
sourcefile <- here("data_raw", "tmf_ec_jrc", "forest_t3.tif")
destfile <- here("data_raw", "tmf_ec_jrc", "forest_n.tif")
system(glue("gdalwarp -overwrite -tap -r sum -tr 1000 1000  \\
         -srcnodata -9999 -ot UInt16 -of GTiff \\
        {sourcefile} \\
        {destfile}"))
# plot(read_stars(destfile), breaks = "equal") if you wanna see it
## If forest_n > 555 (more than 50% of the 1km cell is covered by land)
## Indeed, 555*30*30 = 499500 m2 and 556*30*30 = 500400 m2
sourcefile <- here("data_raw", "tmf_ec_jrc", "forest_n.tif")
destfile <- here("data_raw", "tmf_ec_jrc", "forest_1km.tif")
system(glue("gdal_calc.py --overwrite -A {sourcefile} --calc='(A>555)' --outfile={destfile}"))
r <- read_stars(destfile)
plot(r, reset = FALSE, axes = TRUE, col = c("white", "#7eb00c"), xlim = c(3760000, 4500000))
border <- st_read(here("data_raw", "fao_gaul", "gadm36_NCL.gpkg"), layer = "gadm36_NCL_0")
border <- st_transform(border, crs = st_crs(r))
plot(st_geometry(border), add = TRUE, reset = FALSE)
# Distance to coast
# gdal.proximity
# Distance to river 
# Distance to road
# Perturbation probability (forestatrisk)