## Climate data for New Caledonia
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
library(sf) # for spatial sf objects 
library(stars) # for spatial stars objects 
library(insol) # for function daylength
library(stringr)

Extent <- readLines(here("output/extent_short.txt"))
nodat <- -9999
EPSG <- 3163
proj.s <- "EPSG:4326"
proj.t <- paste0("EPSG:", EPSG) 
ISO_country_code <- "NCL"

##==============================
##
## Chelsa v2.1
## Current climate
##Website : https://chelsa-climate.org/
##============================== 
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
                destfile = here("data_raw", "chelsa_v2_1", "temp", paste0("tasmin_",m,".tif")), method = 'wget')
  ## Monthly maximum temperature (°C).
  download.file(paste0('https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/tasmax/CHELSA_tasmax_',m,'_1981-2010_V.2.1.tif'),
                destfile = here("data_raw","chelsa_v2_1", "temp", paste0("tasmax_",m,".tif")), method = 'wget')
  ## Monthly average temperature (°C).
  download.file(paste0('https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/tas/CHELSA_tas_',m,'_1981-2010_V.2.1.tif'),
                destfile = here("data_raw", "chelsa_v2_1", "temp", paste0("tas_",m,".tif")), method = 'wget')
  ## Monthly precipitation (mm ~ kg/m2).
  download.file(paste0('https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/pr/CHELSA_pr_',m,'_1981-2010_V.2.1.tif'),
                destfile = here("data_raw", "chelsa_v2_1", "temp", paste0("pr_",m,".tif")), method = 'wget')
  ## Monthly cloud cover
  download.file(paste0('https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/tcc/CHELSA_tcc_',m,'_1981-2010_V.2.1.tif'),
                destfile = here("data_raw", "chelsa_v2_1", "temp", paste0("tcc_mean_",m,".tif")), method = 'wget')
  ## Monthly Pet_ penman
  # https://www.fao.org/3/x0490e/x0490e06.htm
  download.file(paste0('https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/pet/CHELSA_pet_penman_',m,'_1981-2010_V.2.1.tif'),
                destfile = here("data_raw", "chelsa_v2_1", "temp", paste0("pet_penman",m,".tif")), method = 'wget')
}

##==============================
##
## 19 Bioclimatic variables
##
##==============================

for(i in 1:19){
  # 19 Bioclimatic variables
  # See https://chelsa-climate.org/wp-admin/download-page/CHELSA_tech_specification_V2.pdf for details
  download.file(paste0('https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio',i,'_1981-2010_V.2.1.tif'),
                destfile=here("data_raw", "chelsa_v2_1", "temp", paste0("bio", str_pad(i, 2, pad = "0"), ".tif")), method = 'wget')
}

# Reproject from lat long to UTM58S (epsg: 3163), set resolution to 1km
# Reframe on New Caledonia and set no data values in sea


border <- sf::st_read(here("data_raw", "fao_gaul", paste0("gadm36_", ISO_country_code, ".gpkg")),
                      layer=paste0("gadm36_", ISO_country_code, "_0"), quiet=FALSE)
border <- st_transform(border[1], crs = EPSG)
borderExtent <- as.integer(strsplit(Extent, split = " ")[[1]])
names(borderExtent) <- c("xmin", "ymin", "xmax", "ymax")
border <- st_crop(border, st_bbox(borderExtent))

for(var in c("tasmin", "tasmax", "tas_", "pr", "bio", "tcc", "pet_penman"))
{ 
  files.tif <- list.files(here("data_raw", "chelsa_v2_1", "temp"), pattern = var, full.names = TRUE)
  for(i in 1:length(files.tif))
  {
    sourcefile <- files.tif[i]
    destfile <- gsub(".tif", "_1km.tif", files.tif[i])
    system(glue("gdalwarp -overwrite -s_srs {proj.s} -t_srs {proj.t} \\
        -r bilinear -tr 1000 1000 -te {Extent} -ot Int16 -of GTiff -srcnodata 0 -dstnodata {nodat} \\
        {sourcefile} {destfile}"))
    reSizeFiles <- st_crop(read_stars(destfile), border)
    if (var %in% c("tasmin", "tasmax", "tas_") | (var == "bio" & i <= 11))
    {
      # stock °C as integer to reduce size
      # °C * 10 to keep information
      reSizeFiles <- round(reSizeFiles * 10)
    }
    write_stars(obj = reSizeFiles, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
                type = "Int16", dsn = destfile)
    # file.remove(sourcefile)
  }
  files.tif <- list.files(here("data_raw", "chelsa_v2_1", "temp"), pattern = var, full.names = TRUE)
  files.tif <- files.tif[grep("[[:digit:]]_1km", files.tif)] # remove original file but not delete it
  r <- read_stars(sort(files.tif), along="band", NA_value = nodat) 
  r <- split(r)
  names(r) <- c(paste(var, 1:length(names(r)), sep = ""))
  r <- merge(r)
  # Define no data values out of New Caledonia border according to elevation map 
  write_stars(obj = r, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
              type = "Int16", dsn = here("data_raw","chelsa_v2_1", paste0(var,"_1km.tif")))
  # file.remove(files.tif)
}
###
# file.remove(list.files(here("data_raw", "chelsa_v2_1", "temp"), pattern = "_1km", full.names = TRUE))
###

# Stack Tasmin, Tasmax, Tas, Pr, Tcc, Pet Penman & bio
files.tif <- here("data_raw", "chelsa_v2_1", paste0(c("tasmin","tasmax","tas_","pr", "tcc", "pet_penman", "bio"),"_1km.tif"))
r <- c(read_stars(files.tif[1]), read_stars(files.tif[2]), read_stars(files.tif[3]), read_stars(files.tif[4]),
       read_stars(files.tif[5]), read_stars(files.tif[6]), read_stars(files.tif[7]), along = "band")
write_stars(obj = r, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
            type = "Int16", dsn = here("data_raw","chelsa_v2_1", "clim_1km.tif"))
# file.remove(files.tif)
rm(r)

##==============================
##
## CWD & NDM
## 
## CWD: climatic water deficit
## NDM: number of dry months
##==============================

pr_file <- here("data_raw", "chelsa_v2_1","pr_1km.tif")
pet_penman_file <- here("data_raw", "chelsa_v2_1","pet_penman_1km.tif")

for (i in 1:12)
{
  # CWD = PET_PENMAN - Pr
  # CWD is a positive values for a lack of water 
   system(glue('gdal_calc.py -A {pet_penman_file} --A_band={i} -B {pr_file} --B_band={i} --quiet --type=Int16 \\
               --creation-option="COMPRESS=LZW" --creation-option="PREDICTOR=2"  --calc="A-B" --NoDataValue={nodat} \\
               --outfile={here("data_raw", "chelsa_v2_1", paste0("cwd", i, "_1km.tif"))} --overwrite'))
}

for (i in 1:12)
{
  # Number of Dry Month ie sum(CWD > 0)
  cwd_file <- here("data_raw", "chelsa_v2_1", paste0("cwd", i, "_1km.tif"))
  ndm_file <- here("data_raw", "chelsa_v2_1", "temp", paste0("ndm", i, "_1km.tif"))
  system(glue('gdal_calc.py -A {cwd_file} --A_band={1} --quiet --type=Int16 \\
              --creation-option="COMPRESS=LZW" --creation-option="PREDICTOR=2" \\
              --outfile={ndm_file} --calc="A>0" --overwrite --NoDataValue={nodat}'))
}
ndm_files <- list.files(here("data_raw", "chelsa_v2_1", "temp"), pattern = "ndm", full.names = TRUE)
system(glue('gdal_calc.py -A {ndm_files[1]} -B {ndm_files[2]} -C {ndm_files[3]} -D {ndm_files[4]}  -E {ndm_files[5]} \\
            -F {ndm_files[6]} -G {ndm_files[7]} -H {ndm_files[8]} -I {ndm_files[9]} -J {ndm_files[10]} -K {ndm_files[11]} \\
            -L {ndm_files[12]} --quiet --type=Int16 --creation-option="COMPRESS=LZW" --creation-option="PREDICTOR=2" \\
            --outfile={here("data_raw", "chelsa_v2_1", "ndm_1km.tif")} --NoDataValue={nodat} \\
            --calc="A+B+C+D+E+F+G+H+I+J+K+L" --overwrite'))

cwd_files <- list.files(here("data_raw", "chelsa_v2_1"), pattern = "cwd", full.names = TRUE)
system(glue('gdal_calc.py -A {cwd_files[1]} -B {cwd_files[2]} -C {cwd_files[3]} -D {cwd_files[4]} -E {cwd_files[5]} \\
            -F {cwd_files[6]} -G {cwd_files[7]} -H {cwd_files[8]} -I {cwd_files[9]} -J {cwd_files[10]} -K {cwd_files[11]} \\
            -L {cwd_files[12]} --quiet --type=Int16 --creation-option="COMPRESS=LZW" --creation-option="PREDICTOR=2" \\
            --outfile={here("data_raw", "chelsa_v2_1", "cwd_1km.tif")} --NoDataValue={nodat} \\
            --calc="numpy.maximum(A,0)+numpy.maximum(B,0)+numpy.maximum(C,0)+numpy.maximum(D,0)+numpy.maximum(E,0)+numpy.maximum(F,0) \\
            +numpy.maximum(G,0)+numpy.maximum(H,0)+numpy.maximum(I,0)+numpy.maximum(J,0)+numpy.maximum(K,0)+numpy.maximum(L,0)" --overwrite'))

system(glue('gdal_merge.py -o {here("output", "current_chelsaNC.tif")} -of GTiff -ot Int16 -co "COMPRESS=LZW" \\
            -co "PREDICTOR=2" -separate -a_nodata {nodat} {here("data_raw", "chelsa_v2_1", "clim_1km.tif")} \\
            {here("data_raw", "chelsa_v2_1", "cwd_1km.tif")} {here("data_raw", "chelsa_v2_1", "ndm_1km.tif")}'))
