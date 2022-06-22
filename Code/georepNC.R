# Libraries
library(here)
library(glue)
library(stars)
library(rgdal)
library(archive)

# Create directories
dir.create(here("data_raw", "georep"))
dir.create(here("data_raw", "georep", "peridotite"))
dir.create(here("data_raw", "georep", "geol"))
dir.create(here("data_raw", "georep", "peridotite", "temp"))
dir.create(here("data_raw", "georep", "geol", "temp"))

EPSG = 3163
nodat = -9999
proj.s <- "EPSG:4326"
proj.t <- paste0("EPSG:", EPSG)
ISO_country_code = "NCL"
Extent <- readLines(here("output/extent_short.txt"))

border <- sf::st_read(here("data_raw", "fao_gaul", paste0("gadm36_", ISO_country_code, ".gpkg")),
                      layer=paste0("gadm36_", ISO_country_code, "_0"), quiet=FALSE)
border <- st_transform(border[1], crs = EPSG)
borderExtent <- as.integer(strsplit(Extent, split = " ")[[1]])
names(borderExtent) <- c("xmin", "ymin", "xmax", "ymax")
border <- st_crop(border, st_bbox(borderExtent))

# Download Geology and Peridotites
dst <- here("data_raw", "georep", "geol", "temp", "geol.7z")
url.tile <- "https://sig-public.gouv.nc//plateforme_telechargement/Geologie_SHP_Lyr_50000.7z"
download.file(url = url.tile, destfile = dst, method = "wget", quiet = TRUE)
archive_extract(dst, dir = here("data_raw", "georep", "geol", "temp"))
# system('7z e -o{here("data_raw","georep", "geol", "temp", "geol.7z")}')

dst <- here("data_raw", "georep", "peridotite", "temp", "peridotite.zip")
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
              -ot Byte -of GTiff -a_srs {proj.t} -tr 1000 1000 {destfile}"))

# Convert geology from .shp to .tif
sourcefile <- here("data_raw", "georep", "geol", "temp", "Geologie_SHP_Lyr_50000", "SurfaceGeologique_50000_Vbeta.shp")
destfile <- here("data_raw", "georep", "geol", "temp" )
litho <- readOGR(sourcefile)
names_lithologie <- unique(litho$lithologie)
nb <- 1:length(names_lithologie)
for (i in 1:length(litho$code)){
  litho$code[i] <- nb[litho$lithologie[i] == names_lithologie][1]
}
litho$code <- as.integer(litho$code)
writeOGR(litho, destfile, layer = "geology", driver = "ESRI Shapefile")
rm("litho")

sourcefile <- here("data_raw", "georep", "geol", "temp", "geology.shp")
destfile <- here("data_raw", "georep", "geol", "ContourGeo.shp")
system(glue("ogr2ogr -overwrite -skipfailures -f 'ESRI Shapefile' -progress \\
              -sql 'SELECT unite, lithologie, legende, code  FROM geology ' \\
            -lco ENCODING=UTF-8 {destfile} {sourcefile}  "))


sourcefile <- here("data_raw", "georep", "geol", "ContourGeo.shp")
destfile <- here("data_raw", "georep", "geol", "ContourGeo.tif")
system(glue("gdal_rasterize  {sourcefile} -a code -co 'COMPRESS=LZW' -co 'PREDICTOR=2' \\
              -ot Byte -of GTiff -a_srs {proj.t} -a_nodata 0 -tr 1000 1000 {destfile}"))

##=====================================
##
## Merge Goblal environ data with NC specific data
##
##===================================== 

system(glue('gdal_merge.py -o {here("output", "environ_allNC.tif")} -of GTiff -ot Int16 -co "COMPRESS=LZW" \\
            -co "PREDICTOR=2" -separate -a_nodata {nodat} {here("output", "environNC.tif")} \\
            {here("data_raw", "georep", "geol", "ContourGeo.tif")} {here("data_raw", "georep", "peridotite", "peridotite.tif")}'))

write_stars(obj = st_crop(read_stars(here("output", "environ_allNC.tif")), border),
            options = c("COMPRESS=LZW","PREDICTOR=2"), type = "Int16",
            NA_value = nodat, dsn = here("output", "environ_allNC.tif"))
