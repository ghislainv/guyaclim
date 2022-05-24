## Futur climate data for New Caledonia
##
## Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
## Jeanne Clément <jeanne.clement@cirad.fr>
##
## May 2022
##

## gdal library is needed to run this script
## http://www.gdal.org/

# Libraries
library(glue)
library(here) 
library(gtools) # for function mixedsort to order files 
library(sf) # for spatial sf objects 
library(stars) # for spatial stars objects 
library(raster)
library(rgdal) # to use gdal  
library(insol) # for function daylength
library(rgrass7) # for function initGRASS
library(dismo) # for function bioclim

## gdalwrap options
# from output/extent.txt
Extent <- readLines(here("output/extent_short.txt"))
Res <- "1000"
nodat <- -9999
proj.s <- "EPSG:4326"
proj.t <- "EPSG:3165" 
ISO_country_code <- "NCL"
EPSG <- 3165

## Border for resizing 

border <- sf::st_read(here("data_raw", "fao_gaul", paste0("gadm36_", ISO_country_code, ".gpkg")),
                      layer=paste0("gadm36_", ISO_country_code, "_0"), quiet=FALSE)
border <- st_transform(border[1], crs = EPSG)
borderExtent <- as.integer(strsplit(Extent, split = " ")[[1]])
names(borderExtent) <- c("xmin", "ymin", "xmax", "ymax")
border <- st_crop(border, st_bbox(borderExtent))

## Initialize GRASS
setwd(here("data_raw"))
dir.create("grassdata")
Sys.setenv(LD_LIBRARY_PATH=paste("/usr/lib/grass80/lib", Sys.getenv("LD_LIBRARY_PATH"),sep=":"))
# use a georeferenced raster
system(glue('grass -c {files.tif[1]} grassdata/climate'))
# connect to grass database
initGRASS(gisBase="/usr/lib/grass80", 
          gisDbase="grassdata", home=tempdir(), 
          location="climate", mapset="PERMANENT",
          override=TRUE)

dir.create(here("data_raw","chelsa_v2_1", "futur")) ## folder for futur climatic data
 
for (model in c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")){
  dir.create(here("data_raw","chelsa_v2_1", "futur", model))
  dir.create(here("data_raw","chelsa_v2_1", "futur", model, "temp"))
  for (m in c(paste0('0',1:9),10:12)){
    
    ## Monthly minimum temperature (°C).
    download.file(paste("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/2071-2100/", model,
                        "/ssp585/tasmin/CHELSA_", tolower(model), "_r1i1p1f1_w5e5_ssp585_tasmin_", m, "_2071_2100_norm.tif", sep = ""),
                  destfile = here("data_raw", "chelsa_v2_1", "futur", model, "temp", paste0("tasmin_",m,".tif")), method = 'wget')
    ## Monthly maximum temperature (°C).
    download.file(paste("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/2071-2100/", model,
                        "/ssp585/tasmax/CHELSA_", tolower(model), "_r1i1p1f1_w5e5_ssp585_tasmax_", m, "_2071_2100_norm.tif", sep = ""),
                  destfile = here("data_raw","chelsa_v2_1", "futur", model, "temp", paste0("tasmax_",m,".tif")), method = 'wget')
    ## Monthly average temperature (°C).
    download.file(paste("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/2071-2100/", model,
                        "/ssp585/tas/CHELSA_", tolower(model), "_r1i1p1f1_w5e5_ssp585_tas_", m, "_2071_2100_norm.tif", sep = ""),
                  destfile = here("data_raw", "chelsa_v2_1", "futur", model, "temp", paste0("tas_",m,".tif")), method = 'wget')
    ## Monthly precipitation (mm ~ kg/m2).
    download.file(paste("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/2071-2100/", model,
                        "/ssp585/pr/CHELSA_", tolower(model), "_r1i1p1f1_w5e5_ssp585_pr_", m, "_2071_2100_norm.tif", sep = ""),
                  destfile = here("data_raw", "chelsa_v2_1", "futur", model, "temp", paste0("pr_",m,".tif")), method = 'wget')
  }
  for(i in 1:19){
    # 19 Bioclimatic variables 
    # See https://chelsa-climate.org/wp-admin/download-page/CHELSA_tech_specification_V2.pdf for details 
    download.file(paste("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/2071-2100/", model, 
                        "/ssp585/bio/CHELSA_bio", i, "_2071-2100_", tolower(model), "_ssp585_V.2.1.tif", sep = ""),
                  destfile=here("data_raw", "chelsa_v2_1", "futur", model, "temp", paste0("bio",i,".tif")), method = 'wget')
  }
}

for (model in c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")){
  for(var in c("tasmin", "tasmax", "tas_", "pr", "bio"))
  {
    files.tif <- list.files(here("data_raw", "chelsa_v2_1", "futur", model, "temp"), pattern = var, full.names = TRUE)
    for(i in 1:length(files.tif))
    {
      sourcefile <- files.tif[i]
      destfile <- gsub(".tif", "_1km.tif", files.tif[i])
      system(glue("gdalwarp -overwrite -s_srs {proj.s} -t_srs {proj.t} \\
        -r bilinear -tr 1000 1000 -te {Extent} -ot Int16 -of GTiff -srcnodata 0 -dstnodata {nodat} \\
        {sourcefile} {destfile}"))
      reSizeFiles <- st_crop(read_stars(destfile), border)
      write_stars(obj = reSizeFiles, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
                  dsn = destfile)
      # file.remove(sourcefile)
    }
    files.tif <- list.files(here("data_raw", "chelsa_v2_1", "futur", model, "temp"), pattern = var, full.names = TRUE)
    files.tif <- files.tif[grep("[[:digit:]]_1km", files.tif)] # remove original file but not delete it
    r <- read_stars(gtools::mixedsort(files.tif), along="band", NA_value = nodat)
    r <- split(r)
    names(r) <- c(paste(var, 1:length(files.tif), sep = ""))
    r <- merge(r)
    # Define no data values out of New Caledonia border according to elevation map
    write_stars(obj = r, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
                dsn = here("data_raw","chelsa_v2_1", "futur", model, paste0(var,"_1km.tif")))
    # file.remove(files.tif)
  }
  # Stack Tasmin, Tasmax, Tas, Pr & Bio
  files.tif <- here("data_raw", "chelsa_v2_1", "futur", model, paste0(c("tasmin","tasmax","tas_","pr", "bio"),"_1km.tif"))
  r <- c(read_stars(files.tif[1]), read_stars(files.tif[2]), read_stars(files.tif[3]), 
         read_stars(files.tif[4]), read_stars(files.tif[5]), along = "band")
  write_stars(obj = r, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
              dsn = here("output", paste0(model, "_1km.tif")))
              # file.remove(files.tif)
}


