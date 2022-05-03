library(sf)
library(rgdal)
library(here)
library(stars)
library(glue)
library(tmaptools)
library(terra)

ISO_country_code = "NCL"
EPSG = 3165
nodat = -9999
proj.s <- "EPSG:4326"
proj.t <- "EPSG:3165"

## Useless in code, but still lat long of New Caledonia
# bb_ll <- c(163, -23, 169, -17)
for (i in seq(163, 167, 2))
{
  for (j in seq(-23,-19, 2))
  {
    url = paste("https://maps.isric.org/mapserv?map=/map/wrb.map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage&COVERAGEID=MostProbable&FORMAT=image/tiff&SUBSET=long(",i, ".0000,", i+2, ".0000)&SUBSET=lat(", j, ".0000,", j+2, ".0000)&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326", sep = "")
    dest = here("data_raw", "soilgrids250_v2_0", paste("soilgrids_", j, "_", i, ".tif", sep = ""))
    download.file(url = url, destfile = dest, verbose = TRUE)
  }
}
sourcefile <- list.files(here("data_raw", "soilgrids250_v2_0"), pattern = ".tif", full.names = TRUE)
merge19 <- c(read_stars(sourcefile[1]), read_stars(sourcefile[2]), read_stars(sourcefile[3]), along = 1)
merge21 <- c(read_stars(sourcefile[4]), read_stars(sourcefile[5]), read_stars(sourcefile[6]), along = 1)
merge23 <- c(read_stars(sourcefile[7]), read_stars(sourcefile[8]), read_stars(sourcefile[9]), along = 1)
merge.tif <- c(merge19, merge21, merge23, along = 2)
# plot(merge.tif, breaks = "equal", main = "NCL")
write_stars(merge.tif, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
            dsn = here("data_raw","soilgrids250_v2_0", "soilgrids_NoExtent.tif"))

sourcefile <- here("data_raw", "soilgrids250_v2_0", "soilgrids_NoExtent.tif")
destfile <- here("data_raw", "soilgrids250_v2_0", "soilgrids_1km.tif")
system(glue("gdalwarp -overwrite -s_srs {proj.s} -t_srs {proj.t} \\
        -r bilinear -tr 1000 1000 -ot Int16 -srcnodata {nodat} -of GTiff \\
        {sourcefile} \\
        {destfile}"))
read_stars(destfile)
plot(read_stars(destfile), breaks = "equal") # ferralsol == ultramafic

#Redimension others tif files
soilgrids = read_stars(here("data_raw", "soilgrids250_v2_0", "soilgrids_1km.tif"))
d = st_dimensions(soilgrids)
offset = c(d[["x"]]$offset, d[["y"]]$offset)
reExtent = paste0(c(xmin = round(offset[1]),
                    ymin = round(offset[2]),
                    xmax = round(offset[1] + d[["x"]]$to * d[["x"]]$delta),
                    ymax = round(offset[2] - d[["y"]]$to * d[["y"]]$delta)), collapse = " ")
writeLines(reExtent, here("output", "reExtent_short.txt"))

###### not working

## Merge soil with elevation, aspect, slope, solar radiation(srad) and roughness
aspect <- read_stars(here("data_raw", "srtm_v1_4_90m", "aspect_1km.tif"))
elevation <- read_stars(here("data_raw", "srtm_v1_4_90m", "elevation_1km.tif"))
roughness <- read_stars(here("data_raw", "srtm_v1_4_90m", "roughness_1km.tif"))
slope <- read_stars(here("data_raw", "srtm_v1_4_90m", "slope_1km.tif"))
srad <- read_stars(here("data_raw", "srtm_v1_4_90m", "srad_1km.tif"))
# xmin = round(max(st_bbox(aspect)$xmin, st_bbox(elevation)$xmin, st_bbox(roughness)$xmin, 
#                  st_bbox(slope)$xmin, st_bbox(srad)$xmin))
# xmax = round(min(st_bbox(aspect)$xmax, st_bbox(elevation)$xmax, st_bbox(roughness)$xmax, 
#                  st_bbox(slope)$xmax, st_bbox(srad)$xmax))
# ymin =round(max(st_bbox(aspect)$ymin, st_bbox(elevation)$ymin, st_bbox(roughness)$ymin, 
#                 st_bbox(slope)$ymin, st_bbox(srad)$ymin))
# ymax = round(min(st_bbox(aspect)$ymax, st_bbox(elevation)$ymax, st_bbox(roughness)$ymax, 
#                  st_bbox(slope)$ymax, st_bbox(srad)$ymax))
# borders = c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax)
# borders = bb(borders)
# st_crs(borders) = st_crs(aspect)
borders = st_bbox(aspect)
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
border = st_transform(border[1], crs = 3165)
border = st_crop(border, borders)
soilgrids = st_crop(soilgrids, border)
write_stars(soilgrids, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
            dsn = here("data_raw", "soilgrids250_v2_0", "soilgrids_1km.tif"), layer = 1)
soilgrids <- read_stars(here("data_raw", "soilgrids250_v2_0", "soilgrids_1km.tif"))

borders <- st_bbox(soilgrids)
aspect <- st_normalize(st_crop(aspect,borders))
elevation <- st_normalize(st_crop(elevation, borders))
roughness <- st_normalize(st_crop(roughness, borders))
slope <- st_normalize(st_crop(slope, borders))
srad <- st_normalize(st_crop(srad, borders))

st_dimensions(aspect)[["x"]]$offset = round(borders$xmin)
st_dimensions(aspect)[["y"]]$offset = round(borders$ymax)
st_dimensions(elevation)[["x"]]$offset = round(borders$xmin)
st_dimensions(elevation)[["y"]]$offset = round(borders$ymax)
st_dimensions(roughness)[["x"]]$offset = round(borders$xmin)
st_dimensions(roughness)[["y"]]$offset = round(borders$ymax)
st_dimensions(slope)[["x"]]$offset = round(borders$xmin)
st_dimensions(slope)[["y"]]$offset = round(borders$ymax)
st_dimensions(srad)[["x"]]$offset = round(borders$xmin)
st_dimensions(srad)[["y"]]$offset = round(borders$ymax)
st_dimensions(soilgrids)[["x"]]$offset = round(borders$xmin)
st_dimensions(soilgrids)[["y"]]$offset = round(borders$ymax)

x_to <- min(st_dimensions(aspect)[['x']]$to, st_dimensions(soilgrids)[['x']]$to, st_dimensions(srad)[['x']]$to)
y_to <- min(st_dimensions(aspect)[['y']]$to, st_dimensions(soilgrids)[['y']]$to, st_dimensions(srad)[['y']]$to)

aspect <- aspect[, 1:x_to, 1:y_to]
elevation <- elevation[, 1:x_to, 1:y_to]
roughness <- roughness[, 1:x_to, 1:y_to]
slope <- slope[, 1:x_to, 1:y_to]
srad <- srad[, 1:x_to, 1:y_to]

r <- c(aspect, elevation, roughness, slope, srad, soilgrids, along = "band", try_hard = TRUE)
write_stars(r, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
            dsn = here("data_raw", "srtm_v1_4_90m", "environ.tif"))

