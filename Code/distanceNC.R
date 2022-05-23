library(stars)
library(here)
library(glue)
library(osmextract)
library(rgdal)
library(RCurl)

ISO_country_code = "NCL"
full_name = " New Caledonia"
EPSG = 3163
proj.s <- "EPSG:4326"
proj.t <- "EPSG:3163"
nodat <- -9999
Extent <- readLines(here("output/extent_short.txt"))
border <- sf::st_read(here("data_raw", "fao_gaul", paste0("gadm36_", ISO_country_code, ".gpkg")),
                      layer=paste0("gadm36_", ISO_country_code, "_0"), quiet=FALSE)
border <- st_transform(border[1], crs = EPSG)

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
plot(read_stars(destfile), breaks = "equal")

##=========================
##
## WDPA : World Database Protected Areas
## International Union for Conservation of Nature
## UNEP-WCMC (2022). Protected Area Profile for New Caledonia from the World Database of Protected Areas, May 2022.
## Available at: www.protectedplanet.net
##
##=========================
dir.create(here("data_raw", "WDPA"))
download.file(paste("https://d1gam3xoknrgr2.cloudfront.net/current/WDPA_WDOECM_May2022_Public_", ISO_country_code, ".zip", sep = ""),
              destfile = here("data_raw", "WDPA", paste("WDPA_WDOECM_May2022_Public_", ISO_country_code, ".zip")), method = 'auto', mode = "wb")
WDPA <- wdpa_read(here("data_raw", "WDPA",paste("WDPA_WDOECM_May2022_Public_", ISO_country_code, ".zip")))
WDPA = st_as_stars(WDPA[3])
st_set_crs(WDPA, st_crs(border))
WDPA <- st_transform_proj(WDPA, proj.t)
WDPA <- st_combine(WDPA)
WDPA <- st_rasterize(st_as_sf(WDPA), dx = 1000, dy = 1000)
WDPA <- st_normalize(st_crop(WDPA, border))
write_stars(WDPA, options = c("COMPRESS=LZW", "PREDICTOR=2"), NA_value = nodat,
            dsn = here("data_raw", "WDPA", "WDPA_1kmBool.tif" ))

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
osm_value <- c( 'highway=motorway or highway=trunk or highway=primary or highway=secondary or highway=primary_link or highway=secondary_link or highway=motorway_link)',
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
              {distance.tif} {distance_1km.tif}"))
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
