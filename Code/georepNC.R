# Libraries
library(here)
library(glue)
library(stars)
library(rgdal)

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
dst <- paste0(here("data_raw", "georep", "geol", "temp/"), "geol", ".7z")
url.tile <- "https://sig-public.gouv.nc//plateforme_telechargement/Geologie_SHP_Lyr_50000.7z"
download.file(url = url.tile, destfile = dst, method = "wget", quiet = TRUE)
system('7z e  -odata_raw/georep/geol/temp/ data_raw/georep/geol/temp/geol.7z')

dst <- paste0(here("data_raw", "georep", "peridotite", "temp/"), "peridotite", ".zip")
url.tile <- "https://opendata.arcgis.com/api/v3/datasets/4daa93c2634048549d0a43c1ba7fe4aa_3/downloads/data?format=shp&spatialRefId=3163"
download.file(url = url.tile, destfile = dst, method = "auto", quiet = TRUE)
unzip(dst, exdir = here("data_raw", "georep", "peridotite", "temp"), overwrite = TRUE)

# Convert peridotites from .shp to .tif

sourcefile <- here("data_raw", "georep", "peridotite", "temp", "Massifs_de_peridotites.shp")
destfile <- here("data_raw", "georep", "peridotite", "peridotite.shp")
system(glue("ogr2ogr -overwrite -f 'ESRI Shapefile' \\
              -lco ENCODING=UTF-8 {destfile} {sourcefile}  "))
sourcefile <- here("data_raw", "georep", "peridotite", "peridotite.shp")
destfile <- here("data_raw", "georep", "peridotite", "peridotite.tif")
system(glue("gdal_rasterize  {sourcefile} -te {Extent} -tap -burn 1 -co 'COMPRESS=LZW' -co 'PREDICTOR=2' \\
              -ot Byte -of GTiff -a_nodata {0} -a_srs {proj.t} -tr 1000 1000 {destfile}"))

# Convert geology from .shp to .tif
sourcefile <- here("data_raw", "georep", "geol", "temp", "SurfaceGeologique_50000_Vbeta.shp")
destfile <- here("data_raw", "georep", "geol", "temp" )
litho <- readOGR(sourcefile)
names_lithologie <- unique(litho$lithologie)
nb <- 1:length(names_lithologie)
for (i in 1:length(litho$code)){
litho$code[i] <- nb[litho$lithologie[i] == names_lithologie][1]
}
litho$code <- as.integer(litho$code)
writeOGR(litho, destfile, layer = "geol", driver = "ESRI Shapefile")
rm("litho")

sourcefile <- here("data_raw", "georep", "geol", "temp", "geol.shp")
destfile <- here("data_raw", "georep", "geol", "ContourGeo.shp")
system(glue("ogr2ogr -overwrite -skipfailures -f 'ESRI Shapefile' -progress \\
              -sql 'SELECT unite, lithologie, legende, code  FROM geol ' \\
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

geol <- read_stars(here("data_raw", "georep", "geol", "ContourGeo.tif"))
peridotite <- read_stars(here("data_raw", "georep", "peridotite", "peridotite.tif"))
peridotite <- st_normalize(st_crop(peridotite, geol))
forest <- read_stars(here("data_raw", "finalfiles", "TMF_NC.tif"))
aspect <- read_stars(here("data_raw", "finalfiles", "aspectNC.tif"))
elevation <- read_stars(here("data_raw", "finalfiles", "elevationNC.tif"))
roughness <- read_stars(here("data_raw", "finalfiles", "roughnessNC.tif"))
slope <- read_stars(here("data_raw", "finalfiles", "slopeNC.tif"))
srad <- read_stars(here("data_raw", "finalfiles", "sradNC.tif"))
soilgrids <- read_stars(here("data_raw", "finalfiles", "soilgridsNC.tif"))
distanceForest <- read_stars(here("data_raw", "finalfiles", "distForestNC.tif"))
distanceSea <- read_stars(here("data_raw", "finalfiles", "distSeaNC.tif"))
distanceRoad <- read_stars(here("data_raw", "finalfiles", "distanceRoadsNC.tif"))
distancePlace <- read_stars(here("data_raw", "finalfiles", "distancePlaceNC.tif"))
distancewater <- read_stars(here("data_raw", "finalfiles", "distancewaterNC.tif"))
WDPA <- read_stars(here("data_raw", "finalfiles", "WDPA_NC.tif"))
# environ <- read_stars(here("output", "environNC.tif"))

x_to <- min(st_dimensions(peridotite)[['x']]$to, st_dimensions(geol)[['x']]$to, st_dimensions(aspect)[['x']]$to)
y_to <- min(st_dimensions(peridotite)[['y']]$to, st_dimensions(geol)[['y']]$to, st_dimensions(aspect)[['y']]$to)

geol       <- st_normalize(geol      [, 1 : x_to, (1 + st_dimensions(geol)[['y']]$to - y_to) : st_dimensions(geol)[['y']]$to])
peridotite <- st_normalize(peridotite[, 1 : x_to, (1 + st_dimensions(peridotite)[['y']]$to - y_to) : st_dimensions(peridotite)[['y']]$to])
forest         <- st_normalize(forest        [, 1 : x_to, (1 + st_dimensions(forest)[['y']]$to - y_to) : st_dimensions(forest)[['y']]$to])
aspect         <- st_normalize(aspect        [, 1 : x_to, (1 + st_dimensions(aspect)[['y']]$to - y_to) : st_dimensions(aspect)[['y']]$to])
elevation      <- st_normalize(elevation     [, 1 : x_to, (1 + st_dimensions(elevation)[['y']]$to - y_to) : st_dimensions(elevation)[['y']]$to])
roughness      <- st_normalize(roughness     [, 1 : x_to, (1 + st_dimensions(roughness)[['y']]$to - y_to) : st_dimensions(roughness)[['y']]$to])
slope          <- st_normalize(slope         [, 1 : x_to, (1 + st_dimensions(slope)[['y']]$to - y_to) : st_dimensions(slope)[['y']]$to])
srad           <- st_normalize(srad          [, 1 : x_to, (1 + st_dimensions(srad)[['y']]$to - y_to) : st_dimensions(srad)[['y']]$to])
soilgrids      <- st_normalize(soilgrids     [, 1 : x_to, (1 + st_dimensions(soilgrids)[['y']]$to - y_to) : st_dimensions(soilgrids)[['y']]$to])
distanceForest <- st_normalize(distanceForest[, 1 : x_to, (1 + st_dimensions(distanceForest)[['y']]$to - y_to) : st_dimensions(distanceForest)[['y']]$to])
distanceSea    <- st_normalize(distanceSea   [, 1 : x_to, (1 + st_dimensions(distanceSea)[['y']]$to - y_to) : st_dimensions(distanceSea)[['y']]$to])
distanceRoad   <- st_normalize(distanceRoad [, 1 : x_to, (1 + st_dimensions(distanceRoad)[['y']]$to - y_to) : st_dimensions(distanceRoad)[['y']]$to])
distancePlace  <- st_normalize(distancePlace [, 1 : x_to, (1 + st_dimensions(distancePlace)[['y']]$to - y_to) : st_dimensions(distancePlace)[['y']]$to])
distancewater  <- st_normalize(distancewater [, 1 : x_to, (1 + st_dimensions(distancewater)[['y']]$to - y_to) : st_dimensions(distancewater)[['y']]$to])
WDPA           <- st_normalize(WDPA          [, 1 : x_to, (1 + st_dimensions(WDPA)[['y']]$to - y_to) : st_dimensions(WDPA)[['y']]$to])
# environ    <- st_normalize(environ   [, 1 : x_to, (1 + st_dimensions(environ)[['y']]$to - y_to) : st_dimensions(environ)[['y']]$to])

st_dimensions(geol)[["x"]]$offset = round(st_bbox(border)$xmin)
st_dimensions(geol)[["y"]]$offset = round(st_bbox(border)$ymax)
st_dimensions(peridotite)[["x"]]$offset = round(st_bbox(border)$xmin)
st_dimensions(peridotite)[["y"]]$offset = round(st_bbox(border)$ymax)
st_dimensions(forest)[["x"]]$offset = round(st_bbox(border)$xmin)
st_dimensions(forest)[["y"]]$offset = round(st_bbox(border)$ymax)
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
st_dimensions(distanceForest)[["x"]]$offset = round(st_bbox(border)$xmin)
st_dimensions(distanceForest)[["y"]]$offset = round(st_bbox(border)$ymax)
st_dimensions(distanceSea)[["x"]]$offset = round(st_bbox(border)$xmin)
st_dimensions(distanceSea)[["y"]]$offset = round(st_bbox(border)$ymax)
st_dimensions(distanceRoad)[["x"]]$offset = round(st_bbox(border)$xmin)
st_dimensions(distanceRoad)[["y"]]$offset = round(st_bbox(border)$ymax)
st_dimensions(distancePlace)[["x"]]$offset = round(st_bbox(border)$xmin)
st_dimensions(distancePlace)[["y"]]$offset = round(st_bbox(border)$ymax)
st_dimensions(distancewater)[["x"]]$offset = round(st_bbox(border)$xmin)
st_dimensions(distancewater)[["y"]]$offset = round(st_bbox(border)$ymax)
st_dimensions(WDPA)[["x"]]$offset = round(st_bbox(border)$xmin)
st_dimensions(WDPA)[["y"]]$offset = round(st_bbox(border)$ymax)
# st_dimensions(environ)[["x"]]$offset = round(st_bbox(border)$xmin)
# st_dimensions(environ)[["y"]]$offset = round(st_bbox(border)$ymax)
r <- c(forest, aspect, elevation, roughness, slope, srad, soilgrids, distanceForest, distanceSea, 
       distanceRoad, distancePlace, distancewater, WDPA, geol, peridotite, along = "band", try_hard = TRUE)
r <- split(r)
names(r) <- c("forest", "aspect", "elevation", "roughness", "slope", "srad", "soilgrids", "distanceForest",
              "distanceSea", "distanceRoad", "distancePlace", "distancewater", "WDPA", "geology", "peridotites")
r <- merge(r)
write_stars(r, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
            dsn = here("output", "environNC.tif"))
plot(r)
dev.print(device = png, file = here("output", "environNC.png"), width = 1000)
dev.off()
