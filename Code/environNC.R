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

## GRASS GIS 8.x.x is also needed to run this script
## https://grass.osgeo.org/

# Libraries
library(glue)
library(here)
library(sf)
library(stars)
library(rgrass7)
library(osmextract)
library(RCurl)
library(wdpar)
library(countrycode)
library(stringr)

source(here("Code", "gaulNC.R"))
## gdalwrap options
EPSG <- 3163
nodat <- -9999
proj.s <- "EPSG:4326"
proj.t <- paste("EPSG:", EPSG, sep = "")
ISO_country_code <- "NCL"
full_name <- countrycode(ISO_country_code, origin = "iso3c", destination = "country.name")
Extent <- readLines(here("output/extent_short.txt"))

##==============================
##
## Borders of New Caledonia
##
##==============================

dir.create(here("data_raw", "fao_gaul"))
URL <- paste0("https://geodata.ucdavis.edu/gadm/gadm3.6/gpkg/gadm36_", ISO_country_code, "_gpkg.zip")
# Download file
# The coordinate reference system is longitude/latitude and the WGS84 datum.
download.file(URL, quiet = FALSE, 
              here("data_raw", "fao_gaul", paste0("gpkg_gadm36_", ISO_country_code, ".zip")))
# Unzip 
unzip(here("data_raw", "fao_gaul", paste0("gpkg_gadm36_", ISO_country_code, ".zip")),
      exdir = here("data_raw", "fao_gaul"), overwrite=TRUE)
# Read vector (level 0 for country borders)
border <- sf::st_read(here("data_raw", "fao_gaul", paste0("gadm36_", ISO_country_code, ".gpkg")),
                      layer=paste0("gadm36_", ISO_country_code, "_0"), quiet=FALSE)
border <- st_transform(border[1], crs = EPSG)


##==============================
##
## Soilgrids
##
##==============================

dir.create(here("data_raw", "soilgrids250_v2_0"))
dir.create(here("data_raw", "soilgrids250_v2_0", "temp"))

# file available only by 2° x 2° for soilgrids250
for (i in seq(163, 167, 2))
{
  for (j in seq(-23,-19, 2))
  {
    url = paste("https://maps.isric.org/mapserv?map=/map/wrb.map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage&COVERAGEID=MostProbable&FORMAT=image/tiff&SUBSET=long(",i, ".0000,", i+2, ".0000)&SUBSET=lat(", j, ".0000,", j+2, ".0000)&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326", sep = "")
    dest = here("data_raw", "soilgrids250_v2_0", "temp", paste("soilgrids_", j, "_", i, ".tif", sep = ""))
    download.file(url = url, destfile = dest, verbose = TRUE)
  }
}
# concatenate files in one .tif
sourcefile <- list.files(here("data_raw", "soilgrids250_v2_0", "temp"), pattern = ".tif", full.names = TRUE)
merge19 <- c(read_stars(sourcefile[1]), read_stars(sourcefile[2]), read_stars(sourcefile[3]), along = 1)
merge21 <- c(read_stars(sourcefile[4]), read_stars(sourcefile[5]), read_stars(sourcefile[6]), along = 1)
merge23 <- c(read_stars(sourcefile[7]), read_stars(sourcefile[8]), read_stars(sourcefile[9]), along = 1)
merge.tif <- c(merge19, merge21, merge23, along = 2)
write_stars(merge.tif, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
            dsn = here("data_raw","soilgrids250_v2_0", "temp", "soilgrids_NoExtent.tif"))

# resize, project and define accurate to 1km x 1km
sourcefile <- here("data_raw", "soilgrids250_v2_0", "temp", "soilgrids_NoExtent.tif")
destfile <- here("data_raw", "soilgrids250_v2_0", "soilgrids_1km.tif")
system(glue("gdalwarp -overwrite -s_srs {proj.s} -t_srs {proj.t} \\
        -r bilinear -tr 1000 1000 -ot Int16 -srcnodata {nodat} -of GTiff \\
        {sourcefile} \\
        {destfile}"))
soilgrids <- read_stars(here("data_raw", "soilgrids250_v2_0", "soilgrids_1km.tif"))
soilgrids <- st_crop(soilgrids, border)
write_stars(soilgrids, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
            dsn = here("data_raw", "soilgrids250_v2_0", "soilgrids_1km.tif"), layer = 1)
unlink(here("data_raw", "soilgrids250_v2_0", "temp"), recursive = TRUE)
rm("merge.tif", "merge19", "merge21", "merge23", "soilgrids")

##==============================
##
##
## Tropical Moist Forest
##
## https://forobs.jrc.ec.europa.eu/TMF/download/TMF_DataUsersGuide.pdf page 12
##
##==============================

dir.create(here("data_raw", "tmf_ec_jrc"))
dir.create(here("data_raw", "tmf_ec_jrc", "temp"))
# file available only by 10° x 10° for Tropical Moist Forest
download.file("https://ies-ows.jrc.ec.europa.eu/iforce/tmf_v1/download.py?type=tile&dataset=AnnualChange_2010&lat=S20&lon=E160",
              destfile = here("data_raw", "tmf_ec_jrc", "temp", "TMF_20_160noExtent.tif"), method = 'auto', mode = "wb")

sourcefile <- here("data_raw", "tmf_ec_jrc", "temp", "TMF_20_160noExtent.tif")
destfile <- here("data_raw", "tmf_ec_jrc", "TMF_20_160_1km.tif")
system(glue("gdalwarp -overwrite -s_srs {proj.s} -t_srs {proj.t} \\
        -r average -tr 100 100 -te {Extent} -ot Int16 -of GTiff -srcnodata 0 -dstnodata {nodat} \\
        {sourcefile} {destfile}"))
forest <- st_as_stars(read_stars(destfile))
# merge differents type of forest in forest = TRUE
forest[forest == 2] <- 1
forest[forest == 4] <- 1
forest[forest != 1] <- 0
write_stars(obj = forest, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
            dsn = destfile)
sourcefile <-  here("data_raw", "tmf_ec_jrc", "TMF_20_160_1km.tif")
destfile <- here("data_raw", "tmf_ec_jrc", "TMF_20_160_1km%.tif")
system(glue("gdalwarp -overwrite -tr 1000 1000 -r average  -ot Float64 -of GTiff -srcnodata {nodat} -dstnodata {nodat} \\
        {sourcefile} {destfile}"))

forest <- read_stars(destfile)
forest <- st_normalize(st_crop(st_as_stars(round(forest *100)), border))
write_stars(obj = forest, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
            type = "Int16", dsn = destfile)
destfile <- here("data_raw", "tmf_ec_jrc", "TMF_20_160_1kmBool.tif")
system(glue("gdalwarp -overwrite -tr 1000 1000 -r average  -ot Int16 -of GTiff -srcnodata {nodat} -dstnodata {nodat} \\
        {sourcefile} {destfile}"))
unlink(here("data_raw", "tmf_ec_jrc", "temp"), recursive = TRUE)
rm("forest")

##==============================
##
## SRTM at 90m resolution from 
## Elevation, slope aspect, roughness
## 
## https://dwtkns.com/srtm/ version 4.1
##
##==============================

dir.create(here("data_raw", "srtm_v1_4_90m"))
dir.create(here("data_raw", "srtm_v1_4_90m", "temp"))
tiles_srtm <- function(extent_latlong)
{
  # Compute lat/long tiles for SRTM data from an extent.
  # This function computes lat/long tiles for SRTM data from an extent
  # in lat/long. See `<http://dwtkns.com/srtm/>`_. SRTM tiles are 5x5
  # degrees. x: -180/+180, y: +60/-60.
  # :param extent_latlong: Extent in lat/long: (xmin, ymin, xmax, ymax).
  # :return: A tuple of two strings indicating tile numbers for lat and long.

  # Tiles for SRTM data
  xmin_latlong <- floor(extent_latlong[1] / 5) * 5
  ymin_latlong <- floor(extent_latlong[3] / 5) * 5
  xmax_latlong <- ceiling(extent_latlong[2] / 5) * 5
  ymax_latlong <- ceiling(extent_latlong[4] / 5) * 5

  x <- xmin_latlong
  y <- ymin_latlong
  tileslat <- NULL
  tileslong <- NULL
  tiles <- NULL
  repeat {
    repeat {
      tileslat <- c(tileslat, x / 5 + 37)
      tileslong <- c(tileslong, -y / 5 + 12)
      y <- y + 5
      if (y > ymax_latlong) {
        break
      }
    }
    x <- x + 5
    y <- ymin_latlong
    if (x > xmax_latlong) {
      break
    }
  }
  for (i in 1:length(tileslat))
  {
    tiles <- c(tiles, paste(str_pad(tileslat[i], 2, side = "left", pad = "0"),
                            str_pad(tileslong[i], 2, side = "left", pad = "0"), sep = "_"))
    # add 0 if value between 1-9
  }
  return(tiles)
}

tiles <- tiles_srtm(as.numeric(strsplit(readLines(here("output/extent_short_latlong.txt")), " ")[[1]]))
# tiles <- c( "69_17", "70_17")

for (i in tiles) {
  if ( url.exists(paste0("https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/srtm_", i, ".zip")) )
  {
    dst <- paste0(here("data_raw", "srtm_v1_4_90m", "temp", "srtm_"), i, ".zip")
    url.tile <- paste0("https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/srtm_", i, ".zip")
    download.file(url = url.tile, destfile = dst, method = "wget", quiet = TRUE)
    unzip(dst, exdir = here("data_raw", "srtm_v1_4_90m", "temp"), overwrite = TRUE)
  }
}

# Merge and Reproject with EPSG 3163
sourcefile <- list.files(here("data_raw", "srtm_v1_4_90m", "temp"), pattern = "*.tif", full.names = TRUE)
file<-file(here("data_raw", "srtm_v1_4_90m", "temp", "sourcefilevrt.txt"))
writeLines(sourcefile, file)
close(file)

destfile <- here("data_raw", "srtm_v1_4_90m", "temp", "srtm.vrt")
system(glue('gdalbuildvrt {destfile} -input_file_list {here("data_raw", "srtm_v1_4_90m", "temp", "sourcefilevrt.txt")}'))
system(glue('gdalwarp -overwrite -t_srs {proj.t} -tap -r bilinear \\
            -co "COMPRESS=LZW" -co "PREDICTOR=2" -te {Extent} -ot Int16 -of GTiff \\
            -tr 90 90 {destfile} {here("data_raw", "srtm_v1_4_90m", "temp", "elevation.tif")}'))

## Compute slope, aspect and roughness using gdaldem 
# compute slope
in_f <- here("data_raw", "srtm_v1_4_90m", "temp", "elevation.tif")
out_f <- here("data_raw", "srtm_v1_4_90m", "temp", "slope.tif")
system(glue('gdaldem slope {in_f} {out_f} -co "COMPRESS=LZW" -co "PREDICTOR=2"'))
# compute aspect
out_f <- here("data_raw", "srtm_v1_4_90m", "temp", "aspect.tif")
system(glue('gdaldem aspect {in_f} {out_f} -co "COMPRESS=LZW" -co "PREDICTOR=2"'))
# compute roughness
out_f <- here("data_raw", "srtm_v1_4_90m", "temp", "roughness.tif")
system(glue('gdaldem roughness {in_f} {out_f} -co "COMPRESS=LZW" -co "PREDICTOR=2"'))

# Resolution from 90m x 90m to 1000m x 1000m using gdalwarp
# elevation
out_f <- here("data_raw", "srtm_v1_4_90m", "elevation_1km.tif")
system(glue('gdalwarp -r bilinear -tr 1000 1000 -ot Int16 -of GTiff \\
        -co "COMPRESS=LZW" -co "PREDICTOR=2" -overwrite {in_f} {out_f}'))
# aspect
in_f <- here("data_raw", "srtm_v1_4_90m", "temp", "aspect.tif")
out_f <-here("data_raw", "srtm_v1_4_90m", "aspect_1km.tif")
system(glue('gdalwarp -r bilinear -tr 1000 1000 -ot Int16 -of GTiff \\
        -co "COMPRESS=LZW" -co "PREDICTOR=2" -overwrite {in_f} {out_f}'))
# slope
in_f <- here("data_raw", "srtm_v1_4_90m", "temp", "slope.tif")
out_f <- here("data_raw", "srtm_v1_4_90m", "slope_1km.tif")
system(glue('gdalwarp -r bilinear -tr 1000 1000 -ot Int16 -of GTiff \\
        -co "COMPRESS=LZW" -co "PREDICTOR=2" -overwrite {in_f} {out_f}'))
# roughness
in_f <- here("data_raw", "srtm_v1_4_90m", "temp", "roughness.tif")
out_f <- here("data_raw", "srtm_v1_4_90m", "roughness_1km.tif")
system(glue('gdalwarp -r bilinear -tr 1000 1000 -ot Int16 -of GTiff \\
        -co "COMPRESS=LZW" -co "PREDICTOR=2" -overwrite {in_f} {out_f}'))
rm("merge_all")

##==============================
##
## Solar radiation
##
#### with r.sun at 90m resolution
## Solar radiation (in Wh.m-2.day-1) was computed from altitude,
## slope and aspect using the function r.sun from the GRASS GIS software.
## We incorporated the shadowing effect of terrain to compute the solar radiation.
## Solar radiation was computed for the Julian day 79 (20th of March for regular years=equinox).
##
##==============================

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
system(glue("r.sun --overwrite --verbose elevation=elevation aspect=aspect slope=slope day=79 glob_rad=global_rad"))

# Export
system(glue("r.out.gdal -f --verbose --overwrite input=global_rad \\
  			 output={here('data_raw', 'srtm_v1_4_90m', 'temp', 'srad.tif')} type=Int16  \\
  			 createopt='COMPRESS=LZW' nodata={nodat}"))

# Resolution from 90m x 90m to 1000m x 1000m using gdalwarp
# srad
in_f <- here("data_raw", "srtm_v1_4_90m", "temp", "srad.tif")
out_f <- here("data_raw", "srtm_v1_4_90m", "srad_1km.tif")
system(glue('gdalwarp  -t_srs {proj.t} \\
        -r bilinear -tr 1000 1000 -ot Int16 -of GTiff \\
        -co "COMPRESS=LZW" -co "PREDICTOR=2" -overwrite {in_f} {out_f}'))

# unlink(here("data_raw", "srtm_v1_4_90m", "temp"), recursive = TRUE)

##===========================
##
## Distance to Forest
##
##===========================

sourcefile <- here("data_raw", "tmf_ec_jrc", "TMF_20_160_1kmBool.tif")
destfile <- here("data_raw", "tmf_ec_jrc", "distForest.tif")
system(glue("gdal_proximity.py -ot Int16 -of GTiff -nodata {nodat} \\
        -values {1} -distunits GEO -use_input_nodata YES {sourcefile} {destfile}"))
# Distance max set to 100km, otherwise islands without forest get big values 
write_stars(st_crop(read_stars(destfile),border), options = c("COMPRESS=LZW", "PREDICTOR=2"), NA_value = nodat,
            dsn = destfile)
 plot(read_stars(destfile))

##===========================
##
## Distance to Sea
##
##===========================
seaBool <- read_stars(destfile) == nodat
write_stars(seaBool, options = c("COMPRESS=LZW", "PREDICTOR=2"), NA_value = nodat,
            dsn = here("data_raw", "tmf_ec_jrc", "Sea_1kmBool.tif"))
sourcefile <- here("data_raw", "tmf_ec_jrc", "Sea_1kmBool.tif")
destfile <- here("data_raw", "tmf_ec_jrc", "distSea.tif")
system(glue("gdal_proximity.py -ot Int16 -of GTiff -nodata {nodat} \\
        -values {nodat} -distunits GEO -use_input_nodata YES {sourcefile} {destfile}"))
write_stars(st_crop(read_stars(destfile), border), options = c("COMPRESS=LZW", "PREDICTOR=2"), NA_value = nodat,
            dsn = destfile)
# plot(read_stars(destfile), breaks = "equal")
rm("seaBool")

##=========================
##
## WDPA : World Database Protected Areas
## International Union for Conservation of Nature
## UNEP-WCMC (2022). Protected Area Profile for New Caledonia from the World Database of Protected Areas, May 2022.
## Available at: www.protectedplanet.net
##
##=========================

dir.create(here("data_raw", "WDPA"))
dir.create(here("data_raw", "WDPA", "temp"))
download.file(paste("https://d1gam3xoknrgr2.cloudfront.net/current/WDPA_WDOECM_May2022_Public_", ISO_country_code, ".zip", sep = ""),
              destfile = here("data_raw", "WDPA","temp", paste("WDPA_WDOECM_May2022_Public_", ISO_country_code, ".zip")), method = 'auto', mode = "wb")
WDPA <- wdpa_read(here("data_raw", "WDPA", "temp", paste("WDPA_WDOECM_May2022_Public_", ISO_country_code, ".zip")))
WDPA = st_as_stars(WDPA[3])
st_set_crs(WDPA, st_crs(border))
WDPA <- st_transform_proj(WDPA, proj.t)
WDPA <- st_combine(WDPA)
WDPA <- st_rasterize(st_as_sf(WDPA), dx = 1000, dy = 1000)
WDPA <- st_normalize(st_crop(WDPA, border))
write_stars(WDPA, options = c("COMPRESS=LZW", "PREDICTOR=2"), NA_value = nodat,
            dsn = here("data_raw", "WDPA", "WDPA_1kmBool.tif" ))
unlink(here("data_raw", "WDPA", "temp"), recursive = TRUE)
rm("WDPA")

##=========================
##
## Open Street Map : distance from cities, roads, rivers
##
##=========================

dir.create(here("data_raw", "OSM"))
dir.create(here("data_raw", "OSM", "temp"))
osm_NCL <- oe_match(full_name)
oe_download(
  file_url = osm_NCL$url, 
  file_size = osm_NCL$file_size,
  force_download = TRUE,
  max_file_size = osm_NCL$file_size + 1,
  download_directory = here("data_raw", "OSM", "temp"))

download_file <- list.files(here("data_raw", "OSM", "temp"), pattern = "osm.pbf", full.names = TRUE)
type_object <- c("lines", "points", "lines", "multipolygons", "multipolygons")
file_name <- c("roads", "place", "river", "lake", "reservoir")
osm_key <- c("highway", "place", "waterway", "natural", "natural")
osm_value <- c( 'highway=motorway or highway=trunk or highway=primary or highway=secondary or highway=primary_link or highway=secondary_link or highway=tertiary or highway=motorway_link)',
                'place=city or place=town or place=village',
                'waterway=river',
                'water=lake',
                'water=reservoir and reservoir_type!=sewage and reservoir_type!=water_storage')
destfile <- paste(substring(download_file, 1, nchar(download_file)-8), ".o5m", sep ="")
system(glue('osmconvert {download_file} -o={destfile}'))
for (i in 1:length(osm_key))
{
  osm_file <- here("data_raw", "OSM" , "temp", paste(file_name[i], ".osm", sep = ""))
  shpfile  <- here("data_raw", "OSM" , "temp", paste(file_name[i], "NoProj.shp", sep = ""))
  projshp  <- here("data_raw", "OSM" , "temp", paste(file_name[i], ".shp", sep = ""))
  file.tif <- here("data_raw", "OSM" , "temp", paste(file_name[i], ".tif", sep = ""))
  distance.tif <- here("data_raw", "OSM", paste(file_name[i], "distance", ".tif", sep = ""))
  distance_1km.tif <- here("data_raw", "OSM", paste(file_name[i], "distance_1km.tif", sep = ""))
  system(glue("osmfilter {destfile}  --keep='{osm_value[i]}' > {osm_file}"))
  system(glue("ogr2ogr -overwrite -skipfailures -f 'ESRI Shapefile' -progress \\
              -sql 'SELECT osm_id, name,{osm_key[i]}  FROM {type_object[i]} WHERE {osm_key[i]} IS NOT NULL' \\
              -lco ENCODING=UTF-8  {shpfile} {osm_file}"))
  system(glue("ogr2ogr -overwrite -s_srs EPSG:4326 -t_srs {proj.t} -f 'ESRI Shapefile' \\
              -lco ENCODING=UTF-8 {projshp} {shpfile} "))
  system(glue("gdal_rasterize  {projshp} -te {Extent} -tap -burn 1 -co 'COMPRESS=LZW' -co 'PREDICTOR=2' \\
              -ot Byte -of GTiff -a_nodata {nodat} -a_srs {proj.t} -tr 100 100 {file.tif}"))
  system(glue("gdal_proximity.py {file.tif} {distance.tif} -f -overwrite -co 'COMPRESS=LZW' -co 'PREDICTOR=2' \\
              -values 1 -ot Int16 -of GTiff -distunits GEO "))
  system(glue("gdalwarp -overwrite -r average -tr 1000 1000 -ot Int16 -srcnodata {nodat} -of GTiff \\
              -dstnodata {nodat} {distance.tif} {distance_1km.tif}"))
}
unlink(here("data_raw", "OSM", "temp"), recursive = TRUE) # delete temporary files

for (i in list.files(here("data_raw", "OSM"), pattern = "*distance_1km.tif", full.names = TRUE))
{
  write_stars(st_crop(read_stars(i), border), i)
}
file.remove(list.files(here("data_raw", "OSM"), pattern = "distance.tif", full.name = TRUE)) # delete temporary files
water <- paste("lake", "reservoir", "river", sep = "|")
watering_place <- list.files(here("data_raw","OSM"), pattern = water, full.names = TRUE)
dim_matrix <- dim(read_stars(watering_place[1])[[1]])[1]
watering_place.tif <- pmin(read_stars(watering_place[1])[[1]], 
                           read_stars(watering_place[2])[[1]],
                           read_stars(watering_place[3])[[1]])
lake <- read_stars(here("data_raw", "OSM", "lakedistance_1km.tif"))
watering_place.tif <- st_as_stars(watering_place.tif, dimension = st_dimensions(lake))
write_stars(watering_place.tif, here("data_raw", "OSM", "wateringplacedistance_1km.tif"))
file.remove(list.files(here("data_raw","OSM"), pattern = water, full.names = TRUE))

##=====================================
##
## Merge environmental variables in one .tif
##
##=====================================

system(glue('gdal_merge.py -ot Int16 -of GTiff -o {here("output", "environNC.tif")} -a_nodata {nodat} -separate \\
            -co "COMPRESS=LZW" -co "PREDICTOR=2" \\
            {here("data_raw", "tmf_ec_jrc", "TMF_20_160_1km%.tif")} {here("data_raw", "srtm_v1_4_90m", "aspect_1km.tif")} \\
            {here("data_raw", "srtm_v1_4_90m", "elevation_1km.tif")} {here("data_raw", "srtm_v1_4_90m", "roughness_1km.tif")} \\
            {here("data_raw", "srtm_v1_4_90m", "slope_1km.tif")} {here("data_raw", "srtm_v1_4_90m", "srad_1km.tif")} \\
            {here("data_raw", "soilgrids250_v2_0", "soilgrids_1km.tif")} {here("data_raw", "tmf_ec_jrc", "distForest.tif")} \\
            {here("data_raw", "tmf_ec_jrc", "distSea.tif")} {here("data_raw", "OSM", "roadsdistance_1km.tif")} \\
            {here("data_raw", "OSM", "placedistance_1km.tif")} {here("data_raw", "OSM", "wateringplacedistance_1km.tif")} \\
            {here("data_raw", "WDPA", "WDPA_1kmBool.tif")} '))
environ <- split(read_stars(here("output", "environNC.tif")))
names(environ) <- c("forest", "aspect", "elevation", "roughness", "slope", "srad", "soilgrids", 
                    "distanceForest", "distanceSea", "distanceRoad", "distancePlace", "distancewater", "WDPA")
write_stars(st_crop(merge(environ), border), dsn = here("output", "environNC.tif"))
