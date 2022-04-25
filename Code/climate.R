## Climate data for French Guyana
##
## Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
## Jeanne Clément <jeanne.clement@cirad.fr>
##
## Octobre 2021
##

## gdal library is needed to run this script
## http://www.gdal.org/

# Libraries
library(glue)
library(here) 
library(gtools) # for function mixedsort to order files 
library(sf) # for spatial sf objects 
library(stars) # for spatial stars objects 
#library(raster)
library(rgdal) # to use gdal  
library(insol) # for function daylength
library(rgrass7) # for function initGRASS
library(dismo) # for function bioclim

## gdalwrap options
# from output/extent.txt
Extent <- readLines(here("output/extent_short.txt"))
Res <- 1000
nodat <- -9999
proj.s <- "EPSG:4326"
proj.t <- "EPSG:2972"
# "EPSG:32622" not legal projection

#=========== Worldclim v2.1 =================
#===== Current climate #================
## Create some directories
dir.create(here("data_raw","worldclim_v2_1")) ## folder for climatic data
dir.create(here("data_raw","worldclim_v2_1","temp")) ## Temporary folder


# Download "zip" files containing 12 GeoTiff (.tif) files, one for each month of the year (January is 1; December is 12).
# They are the average for the years 1970-2000 at 30 seconds (~1 km2) From WorldClim version 2.1.
## Monthly minimum temperature (°C).
download.file('http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_tmin.zip',
              destfile=here("data_raw","worldclim_v2_1","wc2.1_30s_tmin.zip"), method = 'wget')
## Monthly maximum temperature (°C).
download.file('http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_tmax.zip',
              destfile=here("data_raw","worldclim_v2_1","wc2.1_30s_tmax.zip"), method = 'wget')
## Average temperature (°C).
download.file('http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_tavg.zip',
              destfile=here("data_raw","worldclim_v2_1","wc2.1_30s_tavg.zip"), method = 'wget')
## Precipitation (mm).
download.file('http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_prec.zip',
              destfile=here("data_raw","worldclim_v2_1","wc2.1_30s_prec.zip"), method = 'wget')
## Wind speed (m s^-1)
download.file('http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_wind.zip',
              destfile=here("data_raw", "worldclim_v2_1", "wc2.1_30s_wind.zip"), method = 'wget')
## Water vapor pressure (kPa)
download.file('http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_vapr.zip',
              destfile=here("data_raw", "worldclim_v2_1", "wc2.1_30s_vapr.zip"), method = 'wget')

## Unzip worldclim 30s monthly Tmin, Tmax, Tavg and Prec.
## Reframe on French Guyana with a little margin.
guy <- sf::st_read(here("data_raw", "fao_gaul", "FAO_GAUL_GUY.kml"))
bb_ll <- st_bbox(guy, crs=st_crs(4326))
bb_ll_large <- c(floor(bb_ll[c("xmin", "ymin")]), ceiling(bb_ll[c("xmax", "ymax")]))
files.zip <- c("wc2.1_30s_tmax.zip", "wc2.1_30s_tmin.zip", "wc2.1_30s_tavg.zip", "wc2.1_30s_prec.zip")
for (i in 1:length(files.zip)){
  dst <- paste(here("data_raw","worldclim_v2_1"), files.zip[i], sep="/")
  unzip(dst, exdir=here("data_raw", "worldclim_v2_1", "temp"), overwrite=TRUE)
  files.tif <- grep("_crop_", list.files(here("data_raw", "worldclim_v2_1", "temp"),
                                         pattern="tif", full.names = TRUE),
                    invert=TRUE, value=TRUE)
  for(j in 1:length(files.tif)){
    sourcefile <- files.tif[j]
    destfile <- gsub("wc2.1_30s_", "wc2.1_30s_crop_", files.tif[j])
    cmd <- glue("gdal_translate -projwin {bb_ll_large['xmin']}  {bb_ll_large['ymax']} {bb_ll_large['xmax']} {bb_ll_large['ymin']} \\
                -projwin_srs {proj.s} {sourcefile} {destfile}")
    system(cmd)
    file.remove(sourcefile)
  }
  files.tif <- list.files(here("data_raw", "worldclim_v2_1", "temp"), pattern="tif", full.names = TRUE)
  r <- read_stars(files.tif, along="band")
  write_stars(obj=r, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value=nodat,
              dsn=paste(here("data_raw","worldclim_v2_1"),  
                        gsub(".zip",".tif", gsub("wc2.1_30s_","wc2.1_30s_crop_", files.zip[i])),sep="/"))
  file.remove(files.tif)
}

#====== 19 Bioclimatic variables =========
# Compute standard (19) WorldClim Bioclimatic variables from monthly Tmin, Tmax, Tavg and Prec of WorldClim version 2.1.
# Using the function r.bioclim from the GRASS GIS software.
# BIO1 = Annual Mean Temperature
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# BIO3 = Isothermality (BIO2/BIO7) (×100)
# BIO4 = Temperature Seasonality (standard deviation ×100)
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# BIO8 = Mean Temperature of Wettest Quarter
# BIO9 = Mean Temperature of Driest Quarter
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO12 = Annual Precipitation
# BIO13 = Precipitation of Wettest Month
# BIO14 = Precipitation of Driest Month
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
# BIO16 = Precipitation of Wettest Quarter
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter

## Initialize GRASS
setwd(here("data_raw"))
dir.create("grassdata")
Sys.setenv(LD_LIBRARY_PATH=paste("/usr/lib/grass78/lib", Sys.getenv("LD_LIBRARY_PATH"),sep=":"))
# use a georeferenced raster
files.tif <- list.files(here("data_raw", "worldclim_v2_1"), pattern="crop", full.names = TRUE)
system(glue('grass -c {files.tif[1]} grassdata/climate'))
# connect to grass database
initGRASS(gisBase="/usr/lib/grass78", 
          gisDbase="grassdata", home=tempdir(), 
          location="climate", mapset="PERMANENT",
          override=TRUE)
## Import raster in grass
files.tif <- list.files(here("data_raw", "worldclim_v2_1"), pattern="crop")
for(i in 1:length(files.tif)){
  system(glue("r.in.gdal --o input={here('data_raw', 'worldclim_v2_1',files.tif[i])} \\
            output={gsub('.tif', '', gsub('wc2.1_30s_crop_', '', files.tif[i]))}"))
}

# Install the r.bioclim extension
# In the terminal:
# > grass  data_raw/grassdata/climate/PERMANENT
# > g.extension r.bioclim url=https://svn.osgeo.org/grass/grass-addons/grass7/
# Check installation
# > g.extension -a
## Compute standard (19) WorldClim Bioclimatic using the function r.bioclim
# toutscale=10 : scale factor for output temperature output in (C°x10)
system("r.bioclim --verbose --overwrite tmin=`g.list type=rast pat=tmin.* map=. sep=,` \\
       tmax=`g.list type=rast pat=tmax.* map=. sep=,` \\
       tavg=`g.list type=rast pat=tavg.* map=. sep=,` \\
       prec=`g.list type=rast pat=prec.* map=. sep=,` \\
       output=wc2.1_30s_ workers=4 quartals=12 tinscale=1 toutscale=10")
# BIO3 = Isothermality (BIO2/BIO7) ...
# Mean temperature for each quarter year ...
# Traceback (most recent call last):
#   File "/home/clement/.grass7/addons/scripts/r.bioclim", line 618, in <module>
#   main()
# File "/home/clement/.grass7/addons/scripts/r.bioclim", line 296, in main
# > 296
# input = "%s,%s,%s" % (tavgl[m1], tavgl[m2], tavgl[m3]),
# TypeError: list indices must be integers or slices, not float 
# Copy and paste r.bioclim.py from https://github.com/OSGeo/grass-addons/blob/master/grass7/raster/r.bioclim/r.bioclim.py in local file grass7/addons/scripts/r.bioclim to solve error 

# Export bio1-19
system(glue("g.list --overwrite type=rast pat=wc2.1_30s_bio[0-1][0-9] map=. output={here('data_raw','worldclim_v2_1','files')}"))
files <- read.table(here('data_raw','worldclim_v2_1','files'), header = FALSE, sep = "", dec = ".")
for (i in 1:nrow(files)){
  system(glue("r.out.gdal -f --overwrite input={files[i,]} \\
  			 output={paste0(here('data_raw', 'worldclim_v2_1','temp'),'/',files[i,],'.tif')} \\
  			 type=Int16 nodata={nodat} \\
  			 createopt='compress=lzw,predictor=2'"))
}

## Reproject from lat long to UTM32N (epsg: 2972), set resolution to 1km and reframe on French Guyana
# bio1-19
files.tif <- list.files(here("data_raw", "worldclim_v2_1", "temp"), pattern="tif", full.names = TRUE)
for(i in 1:length(files.tif)){
  sourcefile <- files.tif[i]
  destfile <- gsub("wc2.1_30s_", "wc2.1_1km_", files.tif[i])
  system(glue("gdalwarp -overwrite -s_srs {proj.s} -t_srs {proj.t} \\
        -r bilinear -tr {Res} {Res} -te {Extent} -ot Int16 -of GTiff  \\
        -srcnodata {nodat} -dstnodata {nodat} \\
        {sourcefile} \\
        {destfile}"))
  file.remove(sourcefile)
}
files.tif <- list.files(here("data_raw", "worldclim_v2_1", "temp"), pattern="tif", full.names=TRUE)
(r <- read_stars(files.tif, along="band"))
write_stars(obj=r, type="Int16", options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value=nodat,
            dsn=paste(here("data_raw","worldclim_v2_1"),"wc2.1_1km_bio.tif",sep="/"))
file.remove(files.tif)
# Compute bio1-19 using dismo::biovars() function 
# tmin=stack(here("data_raw","worldclim_v2_1","wc2.1_30s_crop_tmin.tif"))
# tmax=stack(here("data_raw","worldclim_v2_1","wc2.1_30s_crop_tmax.tif"))
# tavg=stack(here("data_raw","worldclim_v2_1","wc2.1_30s_crop_tavg.tif"))
# prec=stack(here("data_raw","worldclim_v2_1","wc2.1_30s_crop_prec.tif"))
# bioclim <- dismo::biovars(tmin, tmax, prec)
# writeRaster(bioclim, type="Int16", options = c("COMPRESS=LZW","PREDICTOR=2"),
#             filename=here("data_raw","worldclim_v2_1","wc2.1_30s_crop_bio.tif"))

## Reproject from lat long to UTM32N (epsg: 2972), set resolution to 1km and reframe on French Guyana
# Tmin, Tmax, Tavg and Prec
files.tif <- list.files(here("data_raw", "worldclim_v2_1"), pattern="crop", full.names = TRUE)
for(i in 1:length(files.tif)){
  sourcefile <- files.tif[i]
  destfile <- gsub("worldclim_v2_1/wc2.1_30s_crop_", "worldclim_v2_1/temp/wc2.1_1km_", files.tif[i])
  system(glue("gdalwarp -overwrite -s_srs {proj.s} -t_srs {proj.t} -srcnodata {nodat} -dstnodata {nodat} \\
        -r bilinear -tr {Res} {Res} -te {Extent} -of GTiff  \\
        {sourcefile} \\
        {destfile}"))
  #file.remove(sourcefile)
}
files.tif <- list.files(here("data_raw", "worldclim_v2_1", "temp"), pattern="tif", full.names = TRUE)
(r <- read_stars(files.tif[4:1], along="band"))
write_stars(obj=r, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value=nodat,
            dsn=paste(here("data_raw","worldclim_v2_1"), "wc2.1_1km_clim.tif", sep="/"))
file.remove(files.tif)

#==== PET, CWD and NDM ====
## function to compute PET, CWD and NDM
## PET: potential evapotranspiration (Thornthwaite equation,1948) (or (Priestley-Taylor equation,1972)?)
## CWD: climatic water deficit
## NDM: number of dry months
#== Thornthwaite functions
thorn.indices <- function(clim){
  I <- rep(0,ncell(clim))
  for (i in 1:12) {
    Tmin <- clim[[1]][,,i]
    Tmax <- clim[[1]][,,i+12]
    Tavg <- clim[[1]][,,i+24]
    I <- I+(Tavg/5)^(1.514)
  }
  alpha <- (6.75e-7)*I^3-(7.71e-5)*I^2+(1.792e-2)*I+0.49239
  return(list(I=I,alpha=alpha))
}
thorn.f <- function(Tm,I,alpha,Jday,lat,long) {
  L <- insol::daylength(lat,long,Jday,tmz=3)[,3]
  PET <- 1.6*(L/12)*(10*Tm/I)^alpha
  return(PET)
}

pet.cwd.ndm.f <- function(clim){
  # get latitude in radians
  long_deg <- st_coordinates(st_transform(clim[,,,1],
                                          crs = 4326))[,1]
  lat_deg <- st_coordinates(st_transform(clim[,,,1],
                                         crs = 4326))[,2]
  # initialize
  cwd <- rep(0,ncell(clim)) 
  ndm <- rep(0,ncell(clim))
  pet <- rep(0,ncell(clim))
  # thorn.index
  ind <- thorn.indices(clim)
  # loop on months
  for (i in 1:12){
    cat(paste("Month: ",i,"\n",sep=""))
    evap.thorn <- clim[,,,1] # Evap Thornthwaite
    Tmin <- c(clim[[1]][,,i])
    Tmax <- c(clim[[1]][,,i+12])
    Tavg <- c(clim[[1]][,,i+24])
    Prec <- c(clim[[1]][,,i+36])
    d <- data.frame(day=(30*i)-15,Tmin,Tmax,Tavg,lat_deg,long_deg)
    d[is.na(d)] <- 0
    ## Thornthwaite
    pet.thorn <- thorn.f(Tm=d$Tavg,lat=d$lat_deg,long=d$long_deg,
                         I=ind$I,alpha=ind$alpha,Jday=d$day)*10
    pet.thorn[is.na(Tmin)] <- NA # to correct for NA values
    evap.thorn[[1]][,,1] <- pet.thorn
    if (i==1) {
      PET12.thorn <- c(evap.thorn)
    }
    if (i>1) {
      PET12.thorn <- c(PET12.thorn, evap.thorn, along='band')
    }
    pet <- pet+pet.thorn # annual PET
    pe.diff <- Prec-pet.thorn
    cwd <- cwd + pmin(pe.diff,0.0) # climatic water deficit
    dm <- rep(0,ncell(clim)) # dry month
    dm[pe.diff<0] <- 1
    ndm <- ndm+dm
  }
  # make rasters
  PET <- CWD <- NDM <- clim[,,,1]
  PET[[1]][,,1] <- pet
  CWD[[1]][,,1] <- -cwd
  NDM[[1]][,,1] <- ndm
  NDM[is.na(PET)] <- NA # to account for NA values
  return (list(PET12=PET12.thorn,PET=PET,CWD=CWD,NDM=NDM))
}

# pet.cwd.ndm
clim <- read_stars(here("data_raw", "worldclim_v2_1","wc2.1_1km_clim.tif"), along="band")
pet.cwd.ndm <- pet.cwd.ndm.f(clim)
write_stars(obj=c(pet.cwd.ndm$PET12,pet.cwd.ndm$PET,pet.cwd.ndm$CWD,pet.cwd.ndm$NDM, along="band"),
            overwrite=TRUE, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value=nodat,
            dsn=paste(here("data_raw","worldclim_v2_1","wc2.1_1km_pet_cwd_ndm.tif")))

#===== Output stack ======
clim <- read_stars(here("data_raw","worldclim_v2_1","wc2.1_1km_clim.tif"), along="band")
bioclim <- read_stars(here("data_raw","worldclim_v2_1","wc2.1_1km_bio.tif"), along="band")
pet.cwd.ndm <- read_stars(here("data_raw","worldclim_v2_1","wc2.1_1km_pet_cwd_ndm.tif"), along="band")
os <- c(clim,bioclim,pet.cwd.ndm, along="band")
os <- split(os, "band")
names(os) <- c(paste("tmin",1:12,sep=""),paste("tmax",1:12,sep=""),paste("tavg",1:12,sep=""),paste("prec",1:12,sep=""),
               paste("bio",1:19,sep=""),paste("pet",1:12,sep=""),"pet","cwd","ndm")
# T from °C to 10x°C and round prec, cwd and pet to get integers 
os_int <- round(c(os[1:36,,]*10, os[37:82,,]))

write_stars(obj=merge(os_int), type="Int16", overwrite=TRUE, options = c("COMPRESS=LZW","PREDICTOR=2"),
            NA_value=nodat, 
            dsn=paste(here("output"),"current_wc.tif", sep="/"))
##= Plot
os <- read_stars(paste(here("output"),"current_wc.tif", sep="/"), along="band")
# check if all bands have same number of NAs
# tmin, tmax, tavg, prec OK,
# bioclim less NAs 
apply(is.na(os)[[1]],3,sum)
# chack overlap with border
border <- st_read(here("data_raw","fao_gaul","FAO_GAUL_GUY.kml"))
border <- st_transform(border,crs=st_crs(os))
plot(os[,,,1], axes=TRUE, reset = FALSE)
plot(border, add=TRUE)
## PET
# Pet 1:12 ok but pet seems wrong
png(here("output","pet_wc.png"),width=600,height=600,res=72,pointsize=16)
plot(os[,,,68:80])
dev.off()

## NDM
png(here("output","ndm_wc.png"),width=600,height=600,res=72,pointsize=16)
par(mar=c(3,3,1,1))
plot(os[,,,82], main="ndm")
dev.off()

## CWD
png(here("output","cwd_wc.png"),width=600,height=600,res=72,pointsize=16)
par(mar=c(3,3,1,1))
plot(os[,,,81], main="cwd")
dev.off()

#=========== Chelsa v2.1 =================
#===== Current climate #================
# Website : https://chelsa-climate.org/
## References : 
# Scientific publication: 
# Karger, D.N., Conrad, O., Böhner, J., Kawohl, T., Kreft, H., Soria-Auza, R.W., Zimmermann, 
# N.E., Linder, H.P. & Kessler, M. (2017) Climatologies at high resolution for the earth’s 
# land surface areas. Scientific Data 4, 170122. https://doi.org/10.1038/sdata.2017.122 
# 
# Data citation: 
# Karger, D.N., Conrad, O., Böhner, J., Kawohl, T., Kreft, H., Soria-Auza, R.W., Zimmermann, 
# N.E., Linder, H.P. & Kessler, M. (2021) Climatologies at high resolution for the earth’s 
# land surface areas. EnviDat. https://doi.org/ 10.16904/envidat.228.v2.1

## Create some directories
dir.create(here("data_raw","chelsa_v2_1")) ## folder for climatic data
dir.create(here("data_raw","chelsa_v2_1","temp")) ## Temporary folder

for(m in c(paste0('0',1:9),10:12)){
  ## Monthly minimum temperature (°C).
  download.file(paste0('https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/tasmin/CHELSA_tasmin_',m,'_1981-2010_V.2.1.tif'),
                destfile=here("data_raw", "chelsa_v2_1", "temp", paste0("tasmin_",m,".tif")), method = 'wget')
  ## Monthly maximum temperature (°C).
  download.file(paste0('https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/tasmax/CHELSA_tasmax_',m,'_1981-2010_V.2.1.tif'),
                destfile=here("data_raw","chelsa_v2_1", "temp", paste0("tasmax_",m,".tif")), method = 'wget')
  ## Monthly average temperature (°C).
  download.file(paste0('https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/tas/CHELSA_tas_',m,'_1981-2010_V.2.1.tif'),
                destfile=here("data_raw", "chelsa_v2_1", "temp", paste0("tas_",m,".tif")), method = 'wget')
  ## Monthly precipitation (mm ~ kg/m2).
  download.file(paste0('https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/pr/CHELSA_pr_',m,'_1981-2010_V.2.1.tif'),
                destfile=here("data_raw", "chelsa_v2_1", "temp", paste0("pr_",m,".tif")), method = 'wget')
}

#====== 19 Bioclimatic variables =========
for(i in 1:19){
  # 19 Bioclimatic variables 
  # See https://chelsa-climate.org/wp-admin/download-page/CHELSA_tech_specification_V2.pdf for details 
  download.file(paste0('https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio',i,'_1981-2010_V.2.1.tif'),
                destfile=here("data_raw", "chelsa_v2_1", "temp", paste0("bio",i,".tif")), method = 'wget')
}

# Reproject from lat long to UTM32N (epsg: 2972), set resolution to 1km
# Reframe on French Guyana and set no data values in sea
elev <- read_stars(here("data_raw", "srtm_v1_4_90m", "elevation_1km.tif"))
for(var in c("tasmin","tasmax","tas","pr", "bio")){
  files.tif <- list.files(here("data_raw", "chelsa_v2_1", "temp"), pattern=var, full.names = TRUE)
  for(i in 1:length(files.tif)){
    sourcefile <- files.tif[i]
    destfile <- gsub(".tif", "_1km.tif", files.tif[i])
    system(glue("gdalwarp -overwrite -s_srs {proj.s} -t_srs {proj.t} \\
        -r bilinear -tr {Res} {Res} -te {Extent} -ot Int16 -of GTiff -srcnodata 0 -dstnodata {nodat} \\
        {sourcefile} \\
        {destfile}"))
  file.remove(sourcefile)
  }
  files.tif <- list.files(here("data_raw", "chelsa_v2_1", "temp"), pattern=var, full.names = TRUE)
  r <- read_stars(gtools::mixedsort(files.tif), along="band") 
  # Define no data values out of French Guyana border 
  #borders <- sf::st_read(here("data_raw", "fao_gaul", paste0("gadm36_", ISO_country_code, ".gpkg")),
  #                       layer=paste0("gadm36_", ISO_country_code, "_0"), quiet=!verbose)
  #st_transform(borders, crs=EPSG)
  # st_crop(r, borders)
  # Define no data values in sea according to elevation map 
  r[[1]][is.na(elev[[1]])] <- nodat
  write_stars(obj=r, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value=nodat,
              dsn=here("data_raw","chelsa_v2_1", paste0(var,"_1km.tif")))
 file.remove(files.tif)
}

# Stack Tasmin, Tasmax, Tas and Pr
files.tif <- here("data_raw", "chelsa_v2_1", paste0(c("tasmin","tasmax","tas","pr"),"_1km.tif"))
(r <- read_stars(files.tif, along="band"))
write_stars(obj=r, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value=nodat,
            dsn=paste(here("data_raw","chelsa_v2_1"),"clim_1km.tif", sep="/"))
file.remove(files.tif)

#==== PET, CWD and NDM ====

# pet.cwd.ndm
clim <- read_stars(here("data_raw", "chelsa_v2_1","clim_1km.tif"), along="band")
pet.cwd.ndm <- pet.cwd.ndm.f(clim)
write_stars(obj=c(pet.cwd.ndm$PET12,pet.cwd.ndm$PET,pet.cwd.ndm$CWD,pet.cwd.ndm$NDM, along="band"),
            overwrite=TRUE, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value=nodat,
            dsn=paste(here("data_raw","chelsa_v2_1","pet_cwd_ndm_1km.tif")))

#===== Output stack ======
clim <- read_stars(here("data_raw","chelsa_v2_1","clim_1km.tif"), along="band")
bioclim <- read_stars(here("data_raw","chelsa_v2_1","bio_1km.tif"), along="band")
pet.cwd.ndm <- read_stars(here("data_raw","chelsa_v2_1","pet_cwd_ndm_1km.tif"), along="band")
os <- c(clim,bioclim,pet.cwd.ndm, along="band")
os <- split(os, "band")
names(os) <- c(paste("tmin",1:12,sep=""),paste("tmax", 1:12, sep=""), paste("tavg",1:12, sep=""),paste("prec",1:12,sep=""),
               paste("bio",1:19,sep=""), paste("pet", 1:12, sep=""),"pet","cwd","ndm")
# T from °C to 10x°C and round prec, cwd and pet to get integers 
os_int <- round(c(os[1:36,,]*10, os[37:48,,], os[49:59,,]*10, os[60:82,,]))
write_stars(obj=merge(os_int), overwrite=TRUE, options = c("COMPRESS=LZW","PREDICTOR=2"),
            type="Int16", NA_value=nodat, dsn=paste(here("output"),"current_chelsa.tif", sep="/"))
# type="Int16" to reduce size => transfrom Temperatures in °C*10 or find a way to keep scale and offset from Chelsa 
##= Plot
os <- read_stars(paste(here("output"),"current_chelsa.tif", sep="/"), along="band")
# check if all bands have same number of NAs 
apply(is.na(os)[[1]],3,sum)
# check overlap with border
border <- st_read(here("data_raw","fao_gaul","FAO_GAUL_GUY.kml"))
border <- st_transform(border,crs=st_crs(os))
plot(os[,,,1],axes=TRUE, reset = FALSE)
plot(border, add=TRUE)
## PET
# Pet 1:12 monthly and annual pet 
png(here("output","pet_chelsa.png"), width=600,
    height=600, res=72, pointsize=16)
plot(os[,,,68:80])
dev.off()

## NDM
png(here("output","ndm_chelsa.png"), width=600,
    height=600, res=72, pointsize=16)
par(mar=c(3,3,1,1))
plot(os[,,,82], main="ndm")
dev.off()

## CWD
png(here("output","cwd_chelsa.png"), width=600, height=600, res=72, pointsize=16)
par(mar=c(3,3,1,1))
plot(os[,,,81], main="cwd")
dev.off()


#===== Future climate #================

