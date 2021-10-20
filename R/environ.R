##=====================================================
##
## Environmental data for Guyane
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
arg <- commandArgs(trailingOnly=TRUE)
down <- TRUE
if (length(arg)>0) {
  down <- arg[1]
}

# Libraries
library(glue)
library(here)
library(sf)
library(raster)
library(rgdal)
library(rgrass7)

## gdalwrap options
# from output/extent.txt
Extent <- "100000 230000 435000 642000"
Res <- "1000"
nodat <- "-9999"
proj.s <- "EPSG:4326"
proj.t <- "EPSG:32622"

# Download forest cover of Guyane in 2000 from tmf_ec_jrc (see note)
download.file("https://drive.google.com/uc?export=download&id=1FBL_Jy8QRi-wtG3rcsKZoJldzkwHyKn9",
              destfile=here("data_raw", "tmf_ec_jrc", "forest_t3.tif"), method = 'auto', mode="wb")

# SRTM at 90m resolution from https://dwtkns.com/srtm/ version 4.1
## Download and unzip CGIAR-CSI 90m DEM data
tiles <- c("26_11","26_12")
for (i in 1:length(tiles)) {
  dst <- paste0(here("data_raw", "srtm_v1_4_90m","srtm_"),tiles[i],".zip")
  if (down) {
    url.tile <- paste0("https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/srtm_",tiles[i],".zip")
    download.file(url=url.tile,destfile=dst,method="auto",quiet=TRUE)
  }
  unzip(dst,exdir=here("data_raw", "srtm_v1_4_90m"),overwrite=TRUE)
}
## Mosaic with gdalbuildvrt
destfile <- here("data_raw", "srtm_v1_4_90m", "temp", "elevation.vrt")
sourcefile <- here("data_raw", "srtm_v1_4_90m", "temp", "*.tif")
system(glue("gdalbuildvrt {destfile} {sourcefile}"))

## Reproject from lat long to UTM32N (epsg: 32622) (dstnodata need to be set to 32767 as we pass from Int16 (nodata=-32768) to INT2S (nodata=-32767) in R)
border <- st_read(here("data_raw","fao_gaul","FAO_GAUL_GUY.kml"))
border <- st_transform(border,crs=proj.t)
destfile <- here("data_raw","fao_gaul","FAO_GAUL_GUY_epsg32622.kml")
sourcefile <- here("data_raw", "srtm_v1_4_90m", "temp", "elevation.vrt")
destfile <- here("data_raw", "srtm_v1_4_90m", "temp", "elevation.tif")
system(glue("gdalwarp -overwrite -s_srs {proj.s} -t_srs {proj.t} -srcnodata -32768 -dstnodata -32767 \\
        -r bilinear -tr 90 90 -te {Extent} -ot Int16 -of GTiff \\
        {sourcefile} \\
        {destfile}"))
# plot(raster(here("data_raw", "srtm_v1_4_90m", "temp", "elevation.tif")))
# plot(border$geometry, add=TRUE)

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
# resolution from 90m x 90m to 1000m x 1000m using gdalwarp
# elevation
out_f <- here("data_raw", "srtm_v1_4_90m", "elevation_1km.tif")
cmd <- glue('gdalwarp -srcnodata -32767 -dstnodata -32767 -r bilinear -tr 1000 1000 -te {Extent} \\
        -co "COMPRESS=LZW" -co "PREDICTOR=2" -overwrite {in_f} {out_f}')
system(cmd)
# aspect
in_f <- here("data_raw", "srtm_v1_4_90m", "temp", "aspect.tif")
out_f <-here("data_raw", "srtm_v1_4_90m", "aspect_1km.tif")
cmd <- glue('gdalwarp -srcnodata -9999 -dstnodata -32767 -r bilinear -tr 1000 1000 -te {Extent} \\
        -co "COMPRESS=LZW" -co "PREDICTOR=2" -overwrite {in_f} {out_f}')
system(cmd)
# slope
in_f <- here("data_raw", "srtm_v1_4_90m", "temp", "slope.tif")
out_f <- here("data_raw", "srtm_v1_4_90m", "slope_1km.tif")
cmd <- glue('gdalwarp -srcnodata -9999 -dstnodata -32767 -r bilinear -tr 1000 1000 -te {Extent} \\
        -co "COMPRESS=LZW" -co "PREDICTOR=2" -overwrite {in_f} {out_f}')
system(cmd)
# roughness
in_f <- here("data_raw", "srtm_v1_4_90m", "temp", "roughness.tif")
out_f <- here("data_raw", "srtm_v1_4_90m", "roughness_1km.tif")
cmd <- glue('gdalwarp -srcnodata -9999 -dstnodata -32767 \\
        -r bilinear -tr 1000 1000 -te {Extent} -ot Int16 -of GTiff \\
        -co "COMPRESS=LZW" -co "PREDICTOR=2" -overwrite {in_f} {out_f}')
system(cmd)

# global radiation of wordclim or recompute ? 

## Reproject
# unzip
# fichier vrt mosaicing 
# gdalwarp
# rm intermediaire files 