library(here)
library(sp)
library(jSDM)
library(stars)
library(glue)
library(rgdal)
library(readr)
library("rnaturalearth") # plotting maps
library("rnaturalearthdata")
library("rnaturalearthhires")
library(ggplot2)
library(viridis)
library(rgrass7)
library(terra)

EPSG <- 3163
nodat <- -9999
ISO_country_code <- "NCL"
proj.t <- paste0("EPSG:", EPSG) 
PA <- read.csv2(here("data_raw", "NCpippn", "Presence_Absence.csv"), sep = ",")
PA$X <- NULL
Ab <- read.csv2(here("data_raw", "NCpippn", "Abondance.csv"), sep = ",")
Ab$X <- NULL

# Merge climate and environ in one tif file

# system(glue('gdal_merge.py -o {here("output", "dataNC.tif")} -of GTiff -ot Int16 -co "COMPRESS=LZW" \\
#             -co "PREDICTOR=2" -separate -a_nodata {nodat} {here("output", "environ_allNC.tif")} \\
#             {here("output", "current_chelsaNC.tif")} '))

write_stars(c(read_stars(here("output", "environ_allNC.tif")), read_stars(here("output", "current_chelsaNC.tif")), 
              along = "band"), dsn = here("output", "dataNC.tif"), 
            options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat)

system(glue('gdal_translate -b 15 -b 88 -b 91 -b 102 -b 107 \\
            {here("output","dataNC.tif")} {here("output", "jSDM_data_model_miss_pr.tif")}'))

system(glue('gdal_calc.py -A {here("output", "dataNC.tif")} --A_band={52} \\
                          -B {here("output", "dataNC.tif")} --B_band={53} \\
                          -C {here("output", "dataNC.tif")} --C_band={54} \\
                          -D {here("output", "dataNC.tif")} --D_band={55} \\
                          -E {here("output", "dataNC.tif")} --E_band={56} \\
                          -F {here("output", "dataNC.tif")} --F_band={57} \\
                          -G {here("output", "dataNC.tif")} --G_band={58} \\
                          -H {here("output", "dataNC.tif")} --H_band={59} \\
                          -I {here("output", "dataNC.tif")} --I_band={60} \\
                          -J {here("output", "dataNC.tif")} --J_band={61} \\
                          -K {here("output", "dataNC.tif")} --K_band={62} \\
                          -L {here("output", "dataNC.tif")} --L_band={63} \\
              --creation-option="COMPRESS=LZW" --creation-option="PREDICTOR=2" --quiet --type=Int16 \\
              --outfile={here("output", "prec_annual.tif")} --calc="A+B+C+D+E+F+G+H+I+J+K+L" --overwrite --NoDataValue={nodat}'))

system(glue('gdal_merge.py -o {here("output", "jSDM_data_model_original.tif")} -of GTiff -ot Int16 -co "COMPRESS=LZW" \\
            -co "PREDICTOR=2" -separate -a_nodata {nodat} \\
            {here("output", "jSDM_data_model_miss_pr.tif")} {here("output", "prec_annual.tif")}'))
system(glue('gdal_calc.py -A {here("output", "jSDM_data_model_original.tif")} --A_band=2 \\
                          -B {here("output", "jSDM_data_model_original.tif")} --B_band=3 \\
                          -C {here("output", "jSDM_data_model_original.tif")} --C_band=4 \\
                          -D {here("output", "jSDM_data_model_original.tif")} --D_band=5 \\
                          -E {here("output", "jSDM_data_model_original.tif")} --E_band=6 \\
                          --calc="A^2" --calc="B^2" --calc="C^2" --calc="D^2" --calc="E^2" \\
             --creation-option="COMPRESS=LZW" --creation-option="PREDICTOR=2" --quiet --type=Int16 \\
             --outfile={here("output", "jSDM_data_model_square.tif")}  --overwrite --NoDataValue={nodat}'))

# system(glue('gdal_merge.py -o {here("output", "jSDM_data_model_combine.tif")} -of GTiff -ot Int16 -co "COMPRESS=LZW" \\
#             -co "PREDICTOR=2" -separate -a_nodata {nodat} \\
#             {here("output", "jSDM_data_model_original.tif")} {here("output", "jSDM_data_model_square.tif")}'))

write_stars(c(read_stars(here("output", "jSDM_data_model_original.tif")), read_stars(here("output", "jSDM_data_model_square.tif")), 
              along = "band"), dsn = here("output", "jSDM_data_model_combine.tif"), 
            options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat)

unlink(here("output", paste0("jSDM_data_model_", c("miss_pr.tif", "original.tif", "square.tif"))))

data <- split(read_stars(here("output", "jSDM_data_model_combine.tif")))
names(data) <- c("ultramafic", "bio1", "bio4", "bio15", "cwd", "prec", "bio1^2", "bio4^2", "bio15^2", "cwd^2", "prec^2")
dataF <- merge(data)
dataF <- split(st_apply(dataF, "attributes", scale))
dataF <- c(data[1,,], dataF[2:11,,])

write_stars(merge(dataF), dsn = here("output", "jSDM_data_final.tif"), 
            options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat)

latlong <- read.csv2(here("data_raw", "NCpippn", "coord_site.csv"), sep = ",")
latlong$X <- NULL
coord <- matrix(as.numeric(unlist(latlong)), nrow = nrow(latlong))[,2:3]
data_site <- matrix(0, nrow = dim(coord)[1], ncol = 11)

for (i in 1:dim(coord)[1]) 
{
  data_site[i,] <- as.integer(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "jSDM_data_model_combine.tif")} \\
                                             -wgs84 {coord[i,1]} {coord[i,2]}'), intern = TRUE))
}
data_site <- data.frame(data_site)
colnames(data_site) <- c("ultramafic", "bio1", "bio4", "bio15", "cwd", "prec", "bio1^2", "bio4^2", "bio15^2", "cwd^2", "prec^2")

# bio1, bio4, bio15, cwd, sum(pr), peridotites           
var_jSDM <- data.matrix(data_site)
var_jSDM[,1] <- as.numeric(var_jSDM[,1] == 1)
write.csv(var_jSDM, here("data_raw", "NCpippn", "var_site.csv"))

# Take main species for losing less time first 300 for Abundance

# nb_species <- colSums(PA)
# names(nb_species) <- NULL
# nb_min_species <- sort(nb_species, decreasing = TRUE)[501]
# PA <- PA[,colSums(PA) >= nb_min_species]
# PA$X <- NULL

nb_species <- colSums(Ab)
# names(nb_species) <- NULL
nb_min_species <- sort(nb_species, decreasing = TRUE)[301]
Ab <- Ab[,colSums(Ab) >= nb_min_species]
Ab$X<- NULL
rm("current", "environ")

##=================
##
## jSDM Poisson Log
##
##=================

Sys.time()
jSDM_pois_log <- jSDM_poisson_log(
  burnin = 5000,
  mcmc = 10000,
  thin = 10,
  count_data = data.matrix(Ab),
  site_data = var_jSDM,
  site_formula = ~. ,
  n_latent = 2,
  site_effect = "random",
  V_lambda = 1,
  V_beta = 1,
  shape = 0.1,
  rate = 0.1,
  seed = 1234,
  verbose = 1
)
Sys.time()
save(jSDM_pois_log, file = here("output", "jSDM_pois_log.RData"))
load(here("output", "jSDM_pois_log.RData"))

## beta_j of the five first species
top_species <- which(colSums(Ab) >= sort(colSums(Ab), decreasing = TRUE)[5])
colSums(PA)[top_species] 
names(top_species) <- NULL
par(mfrow=c(3,2))
np <- nrow(jSDM_pois_log$model_spec$beta_start)
for (i in top_species) {
  for (p in 1:np) {
    coda::traceplot(coda::as.mcmc(jSDM_pois_log$mcmc.sp[[i]][,p]))
    coda::densplot(coda::as.mcmc(jSDM_pois_log$mcmc.sp[[i]][,p]), 
                   main = paste(colnames(jSDM_pois_log$mcmc.sp[[i]])[p], 
                                ", species : ",  names(top_species[top_species == i])))
  }
}


## lambda_j of the first five species
n_latent <- jSDM_pois_log$model_spec$n_latent
par(mfrow=c(2,2))
for (j in top_species) {
  for (l in 1:n_latent) {
    coda::traceplot(jSDM_pois_log$mcmc.sp[[j]][,np+l])
    coda::densplot(jSDM_pois_log$mcmc.sp[[j]][,np+l], 
                   main = paste(colnames(jSDM_pois_log$mcmc.sp[[j]])
                                [np+l], ", species : ",names(top_species[top_species == j])))
  }
}

## Latent variables W_i for the first two sites

par(mfrow=c(2,2))
for (l in 1:n_latent) {
  for (i in 1:2) {
    coda::traceplot(jSDM_pois_log$mcmc.latent[[paste0("lv_",l)]][,i],
                    main = paste0("Latent variable W_", l, ", site ", i))
    coda::densplot(jSDM_pois_log$mcmc.latent[[paste0("lv_",l)]][,i],
                   main = paste0("Latent variable W_", l, ", site ", i))
  }
}


## alpha_i of the first two sites
plot(coda::as.mcmc(jSDM_pois_log$mcmc.alpha[,1:2]))

## V_alpha
par(mfrow=c(2,2))
coda::traceplot(jSDM_pois_log$mcmc.V_alpha)
coda::densplot(jSDM_pois_log$mcmc.V_alpha)
## Deviance
coda::traceplot(jSDM_pois_log$mcmc.Deviance)
coda::densplot(jSDM_pois_log$mcmc.Deviance)


## probit_theta
par (mfrow=c(2,1))
hist(jSDM_pois_log$log_theta_latent, main = "Predicted log theta", xlab ="predicted log theta")
hist(jSDM_pois_log$theta_latent, main = "Predicted theta", xlab ="predicted theta", breaks = 20)

# do not run
## plot_residual_cor(jSDM_pois_log, tl.cex=0.5)
##=================
##
## jSDM Binomial Probit
## about 1h10 to run
##=================

Sys.time()
jSDM_binom_pro <- jSDM_binomial_probit(
  burnin = 5000,
  mcmc = 10000,
  thin = 10,
  presence_data = data.matrix(PA),
  site_formula = ~.,
  site_data = var_jSDM,
  n_latent = 2,
  site_effect = "random",
  V_lambda = 1,
  V_beta = 1,
  shape = 0.1,
  rate = 0.1,
  seed = 1234,
  verbose = 1
)
Sys.time()
save(jSDM_binom_pro, file = here("output", "jSDM_binom_pro.RData"))
load(here("output", "jSDM_binom_pro.RData"))

top_species <- which(colSums(PA) >= sort(colSums(PA), decreasing = TRUE)[5])
np <- nrow(jSDM_binom_pro$model_spec$beta_start)

## beta_j of the top five species
par(mfrow=c(3,2))
for (j in top_species) {
  for (p in 1:np) {
    coda::traceplot(coda::as.mcmc(jSDM_binom_pro$mcmc.sp[[j]][,p]))
    coda::densplot(coda::as.mcmc(jSDM_binom_pro$mcmc.sp[[j]][,p]), 
                   main = paste(colnames(jSDM_binom_pro$mcmc.sp[[j]])[p],
                                ", species : ", names(top_species[top_species == j])))
  }
}

## lambda_j of the top five species
n_latent <- jSDM_binom_pro$model_spec$n_latent
par(mfrow=c(2,2))
for (j in top_species) {
  for (l in 1:n_latent) {
    coda::traceplot(coda::as.mcmc(jSDM_binom_pro$mcmc.sp[[j]][,np+l]))
    coda::densplot(coda::as.mcmc(jSDM_binom_pro$mcmc.sp[[j]][,np+l]), 
                   main = paste(colnames(jSDM_binom_pro$mcmc.sp[[j]])
                                [np+l],", species : ", names(top_species[top_species == j])))
  }
}

## Latent variables W_i for the first two sites
par(mfrow=c(2,2))
for (l in 1:n_latent) {
  for (i in 1:2) {
    coda::traceplot(jSDM_binom_pro$mcmc.latent[[paste0("lv_",l)]][,i],
                    main = paste0("Latent variable W_", l, ", site ", i))
    coda::densplot(jSDM_binom_pro$mcmc.latent[[paste0("lv_",l)]][,i],
                   main = paste0("Latent variable W_", l, ", site ", i))
  }
}

## alpha_i of the first two sites
plot(coda::as.mcmc(jSDM_binom_pro$mcmc.alpha[,1:2]))

## V_alpha
par(mfrow=c(2,2))
coda::traceplot(jSDM_binom_pro$mcmc.V_alpha)
coda::densplot(jSDM_binom_pro$mcmc.V_alpha)
## Deviance
coda::traceplot(jSDM_binom_pro$mcmc.Deviance)
coda::densplot(jSDM_binom_pro$mcmc.Deviance)

## probit_theta
par (mfrow=c(2,1))
hist(jSDM_binom_pro$probit_theta_latent,
     main = "Predicted probit theta", xlab ="predicted probit theta")
hist(jSDM_binom_pro$theta_latent,
     main = "Predicted theta", xlab ="predicted theta")

##===================
##
## Plotting occurences for each species
##
##===================

species_to_plot <- "2321" # Calophyllum Caledonium Vieil.
NC_PIPPN <- read.csv2(here("data_raw", "NCpippn", "NC_PIPPN_2022.csv"), sep = ",")
name_species <- NC_PIPPN$taxaname[NC_PIPPN$id_taxon_ref == species_to_plot][1]
species_to_plot <- paste0("X", species_to_plot)
coord <- read.csv2(here("data_raw", "NCpippn", "coord_site.csv"), sep = ",")
coord$longitude <- as.numeric(coord$longitude)
coord$latitude <- as.numeric(coord$latitude)
coord$X <- NULL
coord[, species_to_plot] <- PA[, species_to_plot]

NC <- ne_countries(scale = 10, returnclass = "sf")
pres <- data.frame(coord[coord[, species_to_plot] == 1, 2:3])
abs <- data.frame(coord[coord[, species_to_plot] == 0, 2:3])

ggplot(data = NC) +
  geom_sf() +
  geom_point(data = coord[,c("longitude", "latitude", species_to_plot)], 
             aes(x = longitude, y = latitude, shape = as.logical(coord[,species_to_plot])), size = 1) +
  ggtitle(paste0('Observed occurences of ', name_species)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_shape_manual(labels = c("Absence", "Presence"), values = c(1, 16), name = "Species state") +
  theme(legend.position = c(0.2, 0.2)) +
  coord_sf(xlim = c(163, 169 ), ylim = c(-22.84805679, -20), expand = FALSE)
# Reel extent in lat/long
# xlim = c(155.86889648, 172.09008789 ), ylim = c(-22.84805679, -17.39916611)

cuts <- c("[0,0.1]", "[0.1,0.2]", "[0.2,0.3]", "[0.3,0.4]", "[0.4,0.5]", "[0.5,0.6]", "[0.6,0.7]", "[0.7,0.8]", "[0.8,0.9]", "[0.9,1]")
data_theta <- data.frame(coord[,c("longitude", "latitude")])
data_theta$theta <- jSDM_binom_pro$theta_latent[, species_to_plot]

ggplot(data = NC) +
  geom_sf() +
  geom_point(data = data_theta, 
             aes(x = longitude, y = latitude, colour = factor(floor(theta * 10) / 10)), size = 1) +
  ggtitle(paste0("Estimated probabilities of presence for ", name_species)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_colour_manual(labels = cuts, values = rev(rocket(10)), name = "Prediction") +
  theme(legend.position = c(0.2, 0.3)) +
  coord_sf(xlim = c(163, 169 ), ylim = c(-22.84805679, -20), expand = FALSE)

##===============
##
##Plotting Species Richness
##
##===============

# sites richness 
species_richness <- data.frame(richness = rowSums(PA))
species_richness$latitude <- coord$latitude
species_richness$longitude <- coord$longitude
species_richness$richness[species_richness$richness >= 28] <- 35 # fix max number of species, can be change
species_richness$estimate <- rowSums(jSDM_binom_pro$theta_latent)
species_richness$estimate[species_richness$estimate >= 28] <- 35 # fix max number of species, can be change
cuts <- c("[0,5]","[5, 10]", "[10, 15]", "[15,20]", "[20,25]", "[25,30]", "[30,35]", paste0("[35,",dim(PA)[2],"]"))
col <- rev(rocket(length(cuts)))
ggplot(data = NC) +
  geom_sf() +
  geom_point(data = species_richness, 
             aes(x = longitude, y = latitude, colour = factor(floor(richness / 4) * 4)), size = 1) +
  ggtitle("Observed species richness") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_colour_manual(labels = cuts, values = rev(rocket(length(cuts))), name = "Number of species") +
  theme(legend.position = c(0.2, 0.3)) +
  coord_sf(xlim = c(163, 169 ), ylim = c(-22.84805679, -20), expand = FALSE)

# Estimated richness

ggplot(data = NC) +
  geom_sf() +
  geom_point(data = species_richness, 
             aes(x = longitude, y = latitude, colour = factor(floor(estimate / 4) * 4)), size = 1) +
  ggtitle("Estimated species richness") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_colour_manual(labels = cuts, values = rev(rocket(length(cuts))), name = "Number of species") +
  theme(legend.position = c(0.2, 0.3)) +
  coord_sf(xlim = c(163, 169 ), ylim = c(-22.84805679, -20), expand = FALSE)
par(mfrow = c(1,1))
plot(rowSums(PA), rowSums(jSDM_binom_pro$theta_latent), main = "Species richness for area with less than 50 differents species",
     xlab = "observed", ylab = "estimated", 
     cex.main = 1.4, cex.lab = 1.3, xlim = c(0, 50), ylim = c(0, 50))
abline(a = 0, b = 1, col = 'red')

##====================
##
## Spatial Interpolation
##
##====================

## Initialize GRASS
setwd(here("output"))
Sys.setenv(LD_LIBRARY_PATH = paste("/usr/lib/grass80/lib", Sys.getenv("LD_LIBRARY_PATH"), sep = ":"))
# use a georeferenced raster
system(glue('grass -c {here("output", "jSDM_data_final.tif")} grassdata/plot'))
# connect to grass database
initGRASS(gisBase = "/usr/lib/grass80", 
          gisDbase = "grassdata", home = tempdir(), 
          location = "plot", mapset = "PERMANENT",
          override = TRUE)


# Import sites parameters
params_sites <- read.csv(here("data_raw", "NCpippn", "var_site.csv"), sep = ",")
coord <- read.csv(here("data_raw", "NCpippn", "coord_site.csv"), sep = ",")
params_sites$latitude <- coord$latitude
params_sites$longitude <- coord$longitude
longlat <- SpatialPoints(params_sites[, c("longitude", "latitude")])
proj4string(longlat) <- CRS("+proj=longlat +ellps=GRS80 +units=m")
# lat-long to UTM58S projection 
xy = spTransform(longlat, CRS("EPSG:3163"))

# alpha 
#alpha_sp <- SpatialPointsDataFrame(xy, data.frame(alpha = jSDM_binom_pro$mcmc.alpha[1000,]))
#alpha_sp <- st_as_sf(x=data.frame(alpha = colMeans(jSDM_binom_pro$mcmc.alpha),
#                                      x=xy@coords[,1], y = xy@coords[,2]), coords=c("x","y"), crs=st_crs(3163))
alpha_sp <- terra::vect(x = data.frame(alpha = colMeans(jSDM_binom_pro$mcmc.alpha),
                                       x = xy@coords[,1], y = xy@coords[,2]), geom = c("x", "y"),
                        crs = crs("EPSG:3163")@projargs)
# writeOGR(obj = alpha_sp, dsn = here("output"), driver = "ESRI Shapefile", layer = "alpha")
rgrass7::write_VECT(alpha_sp, "alpha")


# alpha <- here("output", "alpha.shp")
# system(glue('grass -c {alpha} grassdata/plot'))

# Import Madagascar shape file
clim <- merge(split(read_stars(here("output", "jSDM_data_final.tif")))[2:11,,])

NC <- rast(clim[[1]], crs = "EPSG:3163")
rgrass7::write_RAST(NC, "NC")

# Re-sample with RST
# Note: use maskmap=Madagascar to save computation time
# for punctual data use function v.surf.rst
system("v.surf.rst --overwrite --verbose -t tension=3 input=alpha zcolumn=alpha \\
       smooth=0.0 elevation=alpha_rst")
# Export
system(glue('r.out.gdal --overwrite input=alpha_rst \\
             output={here("output", "alphas.tif")} type=Float32 \\
             createopt="compress=lzw,predictor=2"'))

# Representation
test_in <- read_stars(here("output", "alphas.tif"))
plot(test_in, main="Site effect alpha interpolated by RST")
test_xy <- rep(0, dim(coord)[1])
for (i in 1:dim(coord)[1]) 
{
  test_xy[i] <- as.integer(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "alphas.tif")} \\
                                             -wgs84 {coord[i,"longitude"]} {coord[i,"latitude"]}'), intern = TRUE))
}
plot(test_xy, colMeans(jSDM_binom_pro$mcmc.alpha), xlab = "alpha interpolated by RST", ylab = "alpha estimated by JSDM",
     main = "Random site effect")
abline(a = 0, b = 1, col = 'red')
