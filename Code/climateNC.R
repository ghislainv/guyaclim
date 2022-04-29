## Climate data for New Caledonia
##
## Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
## Jeanne Clément <jeanne.clement@cirad.fr>
##
## April 2022
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
Extent <- readLines(here("output/reExtent_short.txt"))
Res <- "1000"
nodat <- -9999
proj.s <- "EPSG:4326"
proj.t <- "EPSG:3165" 

## Initialize GRASS
setwd(here("data_raw"))
dir.create("grassdata")
Sys.setenv(LD_LIBRARY_PATH=paste("/usr/lib/grass80/lib", Sys.getenv("LD_LIBRARY_PATH"),sep=":"))
# use a georeferenced raster
files.tif <- list.files(here("data_raw", "worldclim_v2_1"), pattern="crop", full.names = TRUE)
system(glue('grass -c {files.tif[1]} grassdata/climate'))
# connect to grass database
initGRASS(gisBase="/usr/lib/grass80", 
          gisDbase="grassdata", home=tempdir(), 
          location="climate", mapset="PERMANENT",
          override=TRUE)

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

# Reproject from lat long to UTM58S (epsg: 3165), set resolution to 1km
# Reframe on New Caledonia and set no data values in sea
elev <- read_stars(here("data_raw", "srtm_v1_4_90m", "elevation_1km.tif"))
for(var in c("tasmin", "tasmax", "tas_", "pr", "bio"))
    { 
      files.tif <- list.files(here("data_raw", "chelsa_v2_1", "temp"), pattern=var, full.names = TRUE)
      for(i in 1:length(files.tif))
      {
        sourcefile <- files.tif[i]
        destfile <- gsub(".tif", "_1km.tif", files.tif[i])
        system(glue("gdalwarp -overwrite -s_srs {proj.s} -t_srs {proj.t} \\
        -r bilinear -tr 1000 1000 -te {Extent} -ot Int16 -of GTiff -srcnodata 0 -dstnodata {nodat} \\
        {sourcefile} \\
        {destfile}"))
        # file.remove(sourcefile)
      }
      files.tif <- list.files(here("data_raw", "chelsa_v2_1", "temp"), pattern = var, full.names = TRUE)
      files.tif <- files.tif[grep("[[:digit:]]_1km", files.tif)] # remove original file but not delete it
      r <- read_stars(gtools::mixedsort(files.tif), along="band", NA_value = nodat) 
      # Define no data values out of New Caledonia border according to elevation map 
      #r[[1]][is.na(elev[[1]])] <- nodat
      write_stars(obj=r, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value=nodat,
                  dsn=here("data_raw","chelsa_v2_1", paste0(var,"_1km.tif")))
      # file.remove(files.tif)
}
###
# file.remove(list.files(here("data_raw", "chelsa_v2_1", "temp"), pattern = "_1km", full.names = TRUE))
### plot(read_stars(here("data_raw", "chelsa_v2_1", "temp", "bio12_1km.tif")))
# Stack Tasmin, Tasmax, Tas and Pr
files.tif <- here("data_raw", "chelsa_v2_1", paste0(c("tasmin","tasmax","tas_","pr"),"_1km.tif"))
(r <- read_stars(files.tif, along = "band"))
write_stars(obj = r, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
            dsn = paste(here("data_raw","chelsa_v2_1"),"clim_1km.tif", sep = "/"))
# file.remove(files.tif)

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
write_stars(obj=merge(os), overwrite=TRUE, options = c("COMPRESS=LZW","PREDICTOR=2"),
            NA_value=nodat, dsn=paste(here("output"),"current_chelsa.tif", sep="/"))
# type="Int16" to reduce size => transfrom Temperatures in °C*10 or find a way to keep scale and offset from Chelsa 
##= Plot
os <- read_stars(paste(here("output"),"current_chelsa.tif", sep="/"), along="band")
# check if all bands have same number of NAs 
apply(is.na(os)[[1]],3,sum)
# chack overlap with border
border <- st_read(here("data_raw","fao_gaul","gadm36_NCL.gpkg"), layer = "gadm36_NCL_0")# FAO_GAUL_NCL.kml
border <- st_transform(border,crs = st_crs(os))
plot(os[,,,1],axes = TRUE, reset = FALSE)
plot(border, add=TRUE, col = NA)
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

