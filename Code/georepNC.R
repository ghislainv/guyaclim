# Libraries
library(here)
library(glue)
library(stars)

# Create directories
dir.create(here("data_raw", "georep"))
dir.create(here("data_raw", "georep", "peridotite"))
dir.create(here("data_raw", "georep", "geol"))
dir.create(here("data_raw", "georep", "peridotite", "temp"))
dir.create(here("data_raw", "georep", "geol", "temp"))

EPSG = 3165
nodat = -9999
proj.s <- "EPSG:4326"
proj.t <- paste("EPSG:", EPSG, sep = "")
ISO_country_code = "NCL"
Extent <- readLines(here("output/extent_short.txt"))

# Download Geology and Peridotites
dst <- paste0(here("data_raw", "georep", "geol", "temp/"), "geol", ".7z")
url.tile <- "https://sig-public.gouv.nc//plateforme_telechargement/Geologie_SHP_Lyr_50000.7z"
download.file(url = url.tile, destfile = dst, method = "wget", quiet = TRUE)
system('7z e  -odata_raw/georep/geol/temp/ data_raw/georep/geol/temp/geol.7z')

dst <- paste0(here("data_raw", "georep", "peridotite", "temp/"), "peridotite", ".zip")
url.tile <- "https://opendata.arcgis.com/api/v3/datasets/4daa93c2634048549d0a43c1ba7fe4aa_3/downloads/data?format=shp&spatialRefId=3163"
download.file(url = url.tile, destfile = dst, method = "auto", quiet = TRUE)
unzip(dst, exdir = here("data_raw", "georep", "peridotite", "temp"), overwrite = TRUE)

# Convert peridotites from .shp to .tif

sourcefile <- here("data_raw", "georep", "peridotite", "temp", "Massifs_de_p%C3%A9ridotites.shp")
destfile <- here("data_raw", "georep", "peridotite", "peridotite.shp")
system(glue("ogr2ogr -overwrite -f 'ESRI Shapefile' \\
              -lco ENCODING=UTF-8 {destfile} {sourcefile}  "))
sourcefile <- here("data_raw", "georep", "peridotite", "peridotite.shp")
destfile <- here("data_raw", "georep", "peridotite", "peridotite.tif")
system(glue("gdal_rasterize  {sourcefile} -te {Extent} -tap -burn 1 -co 'COMPRESS=LZW' -co 'PREDICTOR=2' \\
              -ot Byte -of GTiff -a_nodata {nodat} -a_srs {proj.t} -tr 1000 1000 {destfile}"))

# Convert geology from .shp to .tif
sourcefile <- here("data_raw", "georep", "geol", "temp", "SurfaceGeologique_50000_Vbeta.shp")
destfile <- here("data_raw", "georep", "geol", "ContourGeo.shp")
system(glue("ogr2ogr -overwrite -skipfailures -f 'ESRI Shapefile' -progress \\
              -sql 'SELECT unite, lithologie, legende  FROM SurfaceGeologique_50000_Vbeta' \\
            -lco ENCODING=UTF-8 {destfile} {sourcefile}  "))

sourcefile <- here("data_raw", "georep", "geol", "ContourGeo.shp")
destfile <- here("data_raw", "georep", "geol", "ContourGeo.tif")
system(glue("gdal_rasterize  {sourcefile} -a lithologie -co 'COMPRESS=LZW' -co 'PREDICTOR=2' \\
              -ot Byte -of GTiff -a_nodata {nodat} -tr 100 100 {destfile}"))


# not working
