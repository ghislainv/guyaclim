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
library(tmaptools)
library(terra)

## gdalwrap options
Res <- "1000"
nodat <- "-9999"
proj.s <- "EPSG:4326"
proj.t <- "EPSG:3165"
ISO_country_code = "NCL"
EPSG = 3165
nodat = -9999
proj.s <- "EPSG:4326"
proj.t <- "EPSG:3165"
Extent <- readLines(here("output/extent_short.txt"))

##==============================
##
## Borders of New Caledonia
##
##==============================
URL <- paste0("https://geodata.ucdavis.edu/gadm/gadm3.6/gpkg/gadm36_", ISO_country_code, "_gpkg.zip")
# Download file
# The coordinate reference system is longitude/latitude and the WGS84 datum.
download.file(URL, quiet = FALSE, 
              here("data_raw", "fao_gaul", paste0("gpkg_gadm36_", ISO_country_code, ".zip")))
# Unzip 
unzip(here("data_raw", "fao_gaul", paste0("gpkg_gadm36_", ISO_country_code, ".zip")),
      exdir=here("data_raw", "fao_gaul"), overwrite=TRUE)
# Read vector (level 0 for country borders)
border <- sf::st_read(here("data_raw", "fao_gaul", paste0("gadm36_", ISO_country_code, ".gpkg")),
                      layer=paste0("gadm36_", ISO_country_code, "_0"), quiet=FALSE)
border <- st_transform(border[1], crs = EPSG)


##==============================
##
## Soilgrids
##
##==============================

# file available only by 2° x 2° for soilgrids250
for (i in seq(163, 167, 2))
{
  for (j in seq(-23,-19, 2))
  {
    url = paste("https://maps.isric.org/mapserv?map=/map/wrb.map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage&COVERAGEID=MostProbable&FORMAT=image/tiff&SUBSET=long(",i, ".0000,", i+2, ".0000)&SUBSET=lat(", j, ".0000,", j+2, ".0000)&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326", sep = "")
    dest = here("data_raw", "soilgrids250_v2_0", paste("soilgrids_", j, "_", i, ".tif", sep = ""))
    download.file(url = url, destfile = dest, verbose = TRUE)
  }
}
# concatenate files in one .tif
sourcefile <- list.files(here("data_raw", "soilgrids250_v2_0"), pattern = ".tif", full.names = TRUE)
merge19 <- c(read_stars(sourcefile[1]), read_stars(sourcefile[2]), read_stars(sourcefile[3]), along = 1)
merge21 <- c(read_stars(sourcefile[4]), read_stars(sourcefile[5]), read_stars(sourcefile[6]), along = 1)
merge23 <- c(read_stars(sourcefile[7]), read_stars(sourcefile[8]), read_stars(sourcefile[9]), along = 1)
merge.tif <- c(merge19, merge21, merge23, along = 2)
# plot(merge.tif, breaks = "equal", main = "NCL")
write_stars(merge.tif, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
            dsn = here("data_raw","soilgrids250_v2_0", "soilgrids_NoExtent.tif"))

# resize, project and define accurate to 1km x 1km
sourcefile <- here("data_raw", "soilgrids250_v2_0", "soilgrids_NoExtent.tif")
destfile <- here("data_raw", "soilgrids250_v2_0", "soilgrids_1km.tif")
system(glue("gdalwarp -overwrite -s_srs {proj.s} -t_srs {proj.t} \\
        -r bilinear -tr 1000 1000 -ot Int16 -srcnodata {nodat} -of GTiff \\
        {sourcefile} \\
        {destfile}"))
# read_stars(destfile)
# plot(read_stars(destfile), breaks = "equal")

#Redimension others tif files
# soilgrids = read_stars(here("data_raw", "soilgrids250_v2_0", "soilgrids_1km.tif"))
# d = st_dimensions(soilgrids)
# offset = c(d[["x"]]$offset, d[["y"]]$offset)
# reExtent = paste0(c(xmin = round(offset[1]),
#                     ymin = round(offset[2]),
#                     xmax = round(offset[1] + d[["x"]]$to * d[["x"]]$delta),
#                     ymax = round(offset[2] - d[["y"]]$to * d[["y"]]$delta)), collapse = " ")
# writeLines(reExtent, here("output", "reExtent_short.txt"))

soilgrids <- read_stars(here("data_raw", "soilgrids250_v2_0", "soilgrids_1km.tif"))
soilgrids <- st_crop(soilgrids, border)
write_stars(soilgrids, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
            dsn = here("data_raw", "soilgrids250_v2_0", "soilgrids_1km.tif"), layer = 1)
# define border smaller
border = st_crop(border, st_bbox(soilgrids))


##==============================
##
##
## Tropical Moist Forest
##
## https://forobs.jrc.ec.europa.eu/TMF/download/TMF_DataUsersGuide.pdf page 12
##==============================

# file available only by 10° x 10° for moist forest
download.file("https://ies-ows.jrc.ec.europa.eu/iforce/tmf_v1/download.py?type=tile&dataset=AnnualChange_2010&lat=S20&lon=E160",
              destfile = here("data_raw", "tmf_ec_jrc", "TMF_20_160noExtent.tif"), method = 'auto', mode = "wb")

sourcefile <- here("data_raw", "tmf_ec_jrc", "TMF_20_160noExtent.tif")
destfile <- here("data_raw", "tmf_ec_jrc", "TMF_20_160_1km.tif")
system(glue("gdalwarp -overwrite -s_srs {proj.s} -t_srs {proj.t} \\
        -r average -tr 100 100 -te {Extent} -ot Int16 -of GTiff -srcnodata 0 -dstnodata {nodat} \\
        {sourcefile} {destfile}"))
forest <- st_as_stars(read_stars(destfile))
forest[forest == 2] <- 1
forest[forest == 4] <- 1
forest[forest != 1] <- 0
# forest <- st_normalize(st_crop(forest, border))
# forest <- forest %in% c(1,2,4)
write_stars(obj = forest, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
            dsn = destfile)
sourcefile <-  here("data_raw", "tmf_ec_jrc", "TMF_20_160_1km.tif")
destfile <- here("data_raw", "tmf_ec_jrc", "TMF_20_160_1kmBool.tif")
system(glue("gdalwarp -overwrite -tr 1000 1000 -r average  -ot Int16 -of GTiff -srcnodata {nodat} -dstnodata {nodat} \\
        {sourcefile} {destfile}"))
forest <- read_stars(destfile)
forest <- st_normalize(st_crop(forest, border))
# set water values to NA
write_stars(obj = forest, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
            dsn = destfile)

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

##=====================================
##
## Merge environmental variables in one .tif
##
##=====================================

forest <- read_stars(here("data_raw", "tmf_ec_jrc", "TMF_20_160_1km.tif"))
aspect <- read_stars(here("data_raw", "srtm_v1_4_90m", "aspect_1km.tif"))
elevation <- read_stars(here("data_raw", "srtm_v1_4_90m", "elevation_1km.tif"))
roughness <- read_stars(here("data_raw", "srtm_v1_4_90m", "roughness_1km.tif"))
slope <- read_stars(here("data_raw", "srtm_v1_4_90m", "slope_1km.tif"))
srad <- read_stars(here("data_raw", "srtm_v1_4_90m", "srad_1km.tif"))
soilgrids <- read_stars(here("data_raw", "soilgrids250_v2_0", "soilgrids_1km.tif"))

border <- soilgrids # st_bbox(forest)
aspect <- st_normalize(st_crop(aspect,border))
elevation <- st_normalize(st_crop(elevation, border))
roughness <- st_normalize(st_crop(roughness, border))
slope <- st_normalize(st_crop(slope, border))
srad <- st_normalize(st_crop(srad, border))
soilgrids <- st_normalize(st_crop(soilgrids, border))
forest <- st_normalize(st_crop(forest, border))

x_to <- min(st_dimensions(aspect)[['x']]$to, st_dimensions(soilgrids)[['x']]$to, st_dimensions(srad)[['x']]$to, st_dimensions(forest)[['x']]$to)
y_to <- min(st_dimensions(aspect)[['y']]$to, st_dimensions(soilgrids)[['y']]$to, st_dimensions(srad)[['y']]$to, st_dimensions(forest)[['y']]$to)

aspect <-    st_normalize(aspect   [, 1 : x_to, (1 + st_dimensions(aspect)[['y']]$to - y_to) : st_dimensions(aspect)[['y']]$to])
elevation <- st_normalize(elevation[, 1 : x_to, (1 + st_dimensions(elevation)[['y']]$to - y_to) : st_dimensions(elevation)[['y']]$to])
roughness <- st_normalize(roughness[, 1 : x_to, (1 + st_dimensions(roughness)[['y']]$to - y_to) : st_dimensions(roughness)[['y']]$to])
slope <-     st_normalize(slope    [, 1 : x_to, (1 + st_dimensions(slope)[['y']]$to - y_to) : st_dimensions(slope)[['y']]$to])
srad <-      st_normalize(srad     [, 1 : x_to, (1 + st_dimensions(srad)[['y']]$to - y_to) : st_dimensions(srad)[['y']]$to])
soilgrids <- st_normalize(soilgrids[, 1 : x_to, (1 + st_dimensions(soilgrids)[['y']]$to - y_to) : st_dimensions(soilgrids)[['y']]$to])
forest <-    st_normalize(forest   [, 1 : x_to, (1 + st_dimensions(forest)[['y']]$to - y_to) : st_dimensions(forest)[['y']]$to])

st_dimensions(aspect)[["x"]]$offset = round(st_bbox(border)$xmin)
st_dimensions(aspect)[["y"]]$offset = round(st_bbox(border)$ymax)
st_dimensions(elevation)[["x"]]$offset = round(st_bbox(border)$xmin)
st_dimensions(elevation)[["y"]]$offset = round(st_bbox(border)$ymax)
st_dimensions(roughness)[["x"]]$offset = round(st_bbox(border)$xmin)
st_dimensions(roughness)[["y"]]$offset = round(st_bbox(border)$ymax)
st_dimensions(slope)[["x"]]$offset = round(st_bbox(border)$xmin)
st_dimensions(slope)[["y"]]$offset = round(st_bbox(border)$ymax)
st_dimensions(srad)[["x"]]$offset = round(st_bbox(border)$xmin)
st_dimensions(srad)[["y"]]$offset = round(st_bbox(border)$ymax)
st_dimensions(soilgrids)[["x"]]$offset = round(st_bbox(border)$xmin)
st_dimensions(soilgrids)[["y"]]$offset = round(st_bbox(border)$ymax)
st_dimensions(forest)[["x"]]$offset = round(st_bbox(border)$xmin)
st_dimensions(forest)[["y"]]$offset = round(st_bbox(border)$ymax)

r <- c(aspect, elevation, roughness, slope, srad, soilgrids, forest, along = "band", try_hard = TRUE)
r <- split(r)
names(r) <- c("aspect", "elevation", "roughness", "slope", "srad", "soilgrids", "forest")
r <- merge(r)
writeLines(paste0(st_bbox(r), collapse = " "), here("output", "reExtent_short.txt"))
write_stars(r, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
            dsn = here("data_raw", "srtm_v1_4_90m", "environ.tif"))


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


### Tropical Moist Forest
download.file("https://ies-ows.jrc.ec.europa.eu/iforce/tmf_v1/download.py?type=tile&dataset=AnnualChange_2010&lat=S20&lon=E160",
              destfile = here("data_raw", "tmf_ec_jrc", "TMF_20_160noExtent.tif"), method = 'auto', mode = "wb")
sourcefile <- here("data_raw", "tmf_ec_jrc", "TMF_20_160noExtent.tif")
destfile <- here("data_raw", "tmf_ec_jrc", "TMF_20_160.tif")
system(glue("gdalwarp -overwrite -s_srs {proj.s} -t_srs {proj.t} \\
        -r bilinear -tr 1000 1000 -te {Extent} -ot Int16 -of GTiff -srcnodata 0 -dstnodata {nodat} \\
        {sourcefile} {destfile} \\
        "))
forest <- st_normalize(st_crop(read_stars(destfile), border))
write_stars(obj = forest, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
            dsn = destfile)
plot(forest)
