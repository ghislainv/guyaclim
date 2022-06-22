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
library(readr)
library(stringr)
library(ade4)
library(parallel)
library(doParallel)

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

system(glue('gdal_translate -b 15 -b 88 -b 91 -b 99 -b 102 -b 107 \\
            {here("output","dataNC.tif")} {here("output", "jSDM_data_model_original.tif")}'))

jSDM_data_model_square <- read_stars(here("output", "jSDM_data_model_original.tif"))
jSDM_data_model_original <- read_stars(here("output", "jSDM_data_model_original.tif"))
for (i in 2:6) {
  jSDM_data_model_square[[1]][,,i] <- jSDM_data_model_square[[1]][,,i] ^ 2  
}

dataF <- split(c(jSDM_data_model_original, jSDM_data_model_square[,,,2:6], along = "band"))
data <- split(read_stars(here("output", "jSDM_data_model_original.tif")))
names(dataF) <- c("ultramafic", "bio1", "bio4", "bio12^2", "bio15", "cwd", "bio1^2", "bio4^2", "bio12^2", "bio15^2", "cwd^2")
dataF <- merge(dataF)
dataF <- split(st_apply(dataF, "attributes", scale))
dataF <- c(data[1,,], dataF[2:11,,])

write_stars(merge(dataF), dsn = here("output", "jSDM_data_final.tif"), 
            options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat)
unlink(here("output", "jSDM_data_model_original.tif"))
rm(jSDM_data_model_square, jSDM_data_model_original)

latlong <- read.csv2(here("data_raw", "NCpippn", "coord_site.csv"), sep = ",")
latlong$X <- NULL
coord <- matrix(as.numeric(unlist(latlong)), nrow = nrow(latlong))[,2:3]
data_site <- matrix(0, nrow = dim(coord)[1], ncol = 11)

for (i in 1:dim(coord)[1]) 
{
  data_site[i,] <- as.integer(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "jSDM_data_final.tif")} \\
                                             -wgs84 {coord[i,1]} {coord[i,2]}'), intern = TRUE))
}
data_site <- data.frame(data_site)
colnames(data_site) <- c("ultramafic", "bio1", "bio4", "bio12", "bio15", "cwd", "bio1^2", "bio4^2", "bio12^2", "bio15^2", "cwd^2")
data_site$ultramafic <- as.numeric(data_site$ultramafic == 1)
write.csv(data_site, here("output", "data_site.csv"))

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


##=================
##
## jSDM Poisson Log
## about 11h10 to run
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
## Spatial Interpolation for alpha and latents variables
##
##====================

border <- sf::st_read(here("data_raw", "fao_gaul", paste0("gadm36_", ISO_country_code, ".gpkg")),
                      layer=paste0("gadm36_", ISO_country_code, "_0"), quiet=FALSE)
border <- st_transform(border[1], crs = EPSG)

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

#====
# alpha 
#====

alpha_sp <- terra::vect(x = data.frame(alpha = colMeans(jSDM_binom_pro$mcmc.alpha),
                                       x = xy@coords[,1], y = xy@coords[,2]), geom = c("x", "y"),
                        crs = xy@proj4string@projargs)
# crs("EPSG:3163")@projargs
rgrass7::write_VECT(alpha_sp, "alpha")

# # Import New Caledonia shape file
# clim <- read_stars(here("output", "jSDM_data_final.tif"))
# NC <- rast(clim[[1]], crs = "EPSG:3163")
# rgrass7::write_RAST(NC, "NC")

# Re-sample with RST
# for punctual data use function v.surf.rst
system("v.surf.rst --overwrite --verbose -t tension=3 input=alpha zcolumn=alpha \\
       smooth=0.0 elevation=alpha_rst")

# Export
system(glue('r.out.gdal --overwrite input=alpha_rst \\
             output={here("output", "alphas.tif")} type=Float32 \\
             createopt="compress=lzw,predictor=2"'))

# Representation
alpha_in <- read_stars(here("output", "alphas.tif"))
plot(st_crop(alpha_in, border), main = "Site effect alpha interpolated by RST")
alpha_xy <- rep(0, dim(coord)[1])
for (i in 1:dim(coord)[1]) 
{
  alpha_xy[i] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "alphas.tif")} \\
                                             -wgs84 {coord[i,"longitude"]} {coord[i,"latitude"]}'), intern = TRUE))
}
plot(alpha_xy, colMeans(jSDM_binom_pro$mcmc.alpha), xlab = "alpha interpolated by RST", ylab = "alpha estimated by JSDM",
     main = "Random site effect")
abline(a = 0, b = 1, col = 'red')

#====
# W1 
#====

W1_sp <- terra::vect(x = data.frame(W1 = colMeans(jSDM_binom_pro$mcmc.latent$lv_1),
                                    x = xy@coords[,1], y = xy@coords[,2]), geom = c("x", "y"),
                     crs = xy@proj4string@projargs)
# crs("EPSG:3163")@projargs
rgrass7::write_VECT(W1_sp, "W1")

# Re-sample with RST
# Note: use maskmap=Madagascar to save computation time
# for punctual data use function v.surf.rst
system("v.surf.rst --overwrite --verbose -t tension=3 input=W1 zcolumn=W1 \\
       smooth=0.0 elevation=W1_rst ")
# Export
system(glue('r.out.gdal --overwrite input=W1_rst \\
             output={here("output", "lv_W1.tif")} type=Float32 \\
             createopt="compress=lzw,predictor=2"'))
# Representation
W1_in <- read_stars(here("output", "lv_W1.tif"))
plot(st_crop(W1_in, border), main = "Latent variable W1 interpolated by RST")
W1_xy <- rep(0, dim(coord)[1])
for (i in 1:dim(coord)[1]) 
{
  W1_xy[i] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "lv_W1.tif")} \\
                                             -wgs84 {coord[i,"longitude"]} {coord[i,"latitude"]}'), intern = TRUE))
}
plot(W1_xy, colMeans(jSDM_binom_pro$mcmc.latent$lv_1), xlab = "W1 interpolated by RST", ylab = "W1 estimated by JSDM", 
     main = "First latent axe")
abline(a = 0, b = 1, col = 'red')

#====
# W2 
#====

W2_sp <- terra::vect(x = data.frame(W2 = colMeans(jSDM_binom_pro$mcmc.latent$lv_2),
                                    x = xy@coords[,1], y = xy@coords[,2]), geom = c("x", "y"),
                     crs = xy@proj4string@projargs)
# crs("EPSG:3163")@projargs
rgrass7::write_VECT(W2_sp, "W2")

# Re-sample with RST
# Note: use maskmap=Madagascar to save computation time
# for punctual data use function v.surf.rst
system("v.surf.rst --overwrite --verbose -t tension=3 input=W2 zcolumn=W2 \\
       smooth=0.0 elevation=W2_rst ")
# Export
system(glue('r.out.gdal --overwrite input=W2_rst \\
             output={here("output", "lv_W2.tif")} type=Float32 \\
             createopt="compress=lzw,predictor=2"'))
# Representation
W2_in <- read_stars(here("output", "lv_W2.tif"))
plot(st_crop(W2_in, border), main = "Latent variable W2 interpolated by RST")
W2_xy <- rep(0, dim(coord)[1])
for (i in 1:dim(coord)[1]) 
{
  W2_xy[i] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "lv_W2.tif")} \\
                                             -wgs84 {coord[i,"longitude"]} {coord[i,"latitude"]}'), intern = TRUE))
}
plot(W2_xy, colMeans(jSDM_binom_pro$mcmc.latent$lv_2), xlab = "W2 interpolated by RST", ylab = "W2 estimated by JSDM", 
     main = "First latent axe")
abline(a = 0, b = 1, col = 'red')

##====================
##
## Compute probabilities for each species
##
##====================

rst_alpha <- read_stars(here("output", "alphas.tif"))
rst_W1 <- read_stars(here("output", "lv_W1.tif"))
rst_W2 <- read_stars(here("output", "lv_W2.tif"))
load(here("output", "jSDM_binom_pro.RData"))
data_site <- read.csv2(here("output", "data_site.csv"), sep = ",")[, -1]

# Species parameters 
### fixed species effect lambdas and betas 

# params_species <- data.frame(jSDM_binom_pro$mcmc.sp)

# Current climatic variables 
X <- data.frame(intercep = rep(1, nrow = data_site), data_site)
scaled_clim_var <- read_stars(here("output", "jSDM_data_final.tif"))
scaled_clim_var[[1]][1,,] = as.numeric(!is.na(scaled_clim_var[[1]][1,,]))
scaled_clim_var <- split(scaled_clim_var)
names(scaled_clim_var) <- c("ultramafic", "bio1", "bio4", "bio15", "cwd", "prec", "bio1^2", "bio4^2", "bio15^2", "cwd^2", "prec^2")

np <- st_dimensions(merge(scaled_clim_var))$attributes$to
nsp <- length(colnames(jSDM_binom_pro$theta_latent))
npart <- 30
species <- colnames(jSDM_binom_pro$model_spec$presence_data)
n_latent <- length(jSDM_binom_pro$mcmc.latent)
lambdas <- matrix(0, nsp, n_latent)
betas <- matrix(0, nsp, np + 1)
for (j in 1:nsp){
  for (l in 1:n_latent){
    lambdas[j, l] <- mean(jSDM_binom_pro$mcmc.sp[[j]][, np + l])
  }
  for (p in 1:(np + 1)){
    betas[j, p] <- mean(jSDM_binom_pro$mcmc.sp[[j]][, p])
  }
}
colnames(betas) <- colnames(jSDM_binom_pro$mcmc.sp[[1]])[1:(np + 1)]
params_species <- data.frame(species = species,
                             betas, lambda_1 = lambdas[, 1], lambda_2 = lambdas[, 2])
write.csv(params_species, file = here("output", "params_species.csv"), row.names = F)
params_species <- read.csv2(here("output", "params_species.csv"), sep = ",", dec = ".")

# Function to compute probit_theta at New Caledonia scale 
predfun <- function(scaled_clim_var, params_species, rst_alpha, rst_W1, rst_W2, species.range){
  # Get lambda values for each species
  lambda_1 <- as.matrix(params_species[, "lambda_1"])
  lambda_2 <- as.matrix(params_species[, "lambda_2"])
  beta <- as.matrix(params_species[,2:13])
  
  ## Xbeta_1
  np <- st_dimensions(merge(scaled_clim_var))$attributes$to
  Xbeta_1 <- st_as_stars(matrix(beta[1,1][[1]], ncol = ncol(rst_alpha), nrow = nrow(rst_alpha)))

  st_crs(Xbeta_1) <- st_crs(rst_alpha)
  st_dimensions(Xbeta_1)$X1$delta <- 1000
  st_dimensions(Xbeta_1)$X2$delta <- -1000
  st_set_bbox(Xbeta_1, st_bbox(rst_alpha))
  st_dimensions(Xbeta_1)[["X1"]]$offset <- st_dimensions(rst_alpha)[['x']]$offset
  st_dimensions(Xbeta_1)[["X2"]]$offset <- st_dimensions(rst_alpha)[['y']]$offset
  for (p in 1:(np - 1)) {
    Xbeta_1 <- Xbeta_1 + scaled_clim_var[[p]] * beta[1, p + 1] 
  }
  
  ## Wlambda_1
  Wlambda_1 <- rst_W1*lambda_1[1] + rst_W2*lambda_2[1]
  ## probit_theta_1
  probit_theta_1 <- Xbeta_1 + Wlambda_1 + rst_alpha
  probit_theta <- probit_theta_1
  remove(list = c("probit_theta_1","Wlambda_1"))
  
  ## Other species
  for (j in (species.range[1] + 1):species.range[2]) {
    
    ## Xbeta_j
    Xbeta_j <- Xbeta_1
    Xbeta_j[[1]] <- rep(beta[j,1][[1]], ncell(Xbeta_j))
    for (p in 1:(np - 1)) {
      Xbeta_j <- Xbeta_j + scaled_clim_var[[p]] * beta[j, p + 1] 
    }
    
    ## Wlambda_j
    Wlambda_j <- rst_W1 * lambda_1[j] + rst_W2 * lambda_2[j] 
    
    ## probit_theta_j
    probit_theta_j <- Xbeta_j + Wlambda_j + rst_alpha
    probit_theta <- c(probit_theta, probit_theta_j)
    remove(list=c("probit_theta_j", "Xbeta_j", "Wlambda_j"))
  }
  names(probit_theta) <- params_species$species[species.range[1]:species.range[2]]
  return(probit_theta)
}


Sys.time() # 6 min to run 
# Compute theta in n parts because it's too large

dir.create(here("output", "theta"))
first.species <- seq(1, nsp, by = floor(nsp / npart) + 1)
for (n in 1:npart){
  probit_theta <- predfun(scaled_clim_var, params_species, rst_alpha, rst_W1, rst_W2,
                          species.range = c(first.species[n], min(nsp, first.species[n] + floor(nsp / npart))))
  write_stars(merge(probit_theta), options = c("COMPRESS=LZW", "PREDICTOR=2"),
              dsn = here("output", "theta", paste0("RST_probit_theta_", str_pad(n, width = 2, pad = "0"), ".tif")))
  
  # Compute probabilities of presence theta
  theta <- merge(probit_theta)
  remove(probit_theta)
  for (j in 1:st_dimensions(theta)$attributes$to) {
    theta[[1]][,,j] <- pnorm(theta[[1]][,,j])
  }
  write_stars(theta, options = c("COMPRESS=LZW", "PREDICTOR=2"), 
              dsn = here("output", "theta", paste0("RST_theta_", str_pad(n, width = 2, pad = "0"), ".tif")))
  remove(theta)
}
Sys.time()

##====================
##
## Predictive maps of presence probabilities
##
##====================

species_to_plot <- 1644
# params_species <- data.frame(jSDM_binom_pro$mcmc.sp[species_to_plot])
# nsp <- nrow(params_species)
theta <- read_stars(here("output", "theta", "RST_theta_01.tif"))
latlong <- read.csv2(here("data_raw", "NCpippn", "coord_site.csv"), sep = ",")
latlong$X <- NULL
coord <- matrix(as.numeric(unlist(latlong)), nrow = nrow(latlong))[,2:3]
NC_PIPPN <- read.csv2(here("data_raw", "NCpippn", "data_clear.csv"), sep = ",")
name_species <- NC_PIPPN$nom_taxon[NC_PIPPN$id_taxon_ref == species_to_plot][1]

# Observed presence absence
id_pres <- which(PA[, paste0("X", species_to_plot)] == 1)
obs_pres <- data.frame(lat = coord[id_pres, 2], long = coord[id_pres, 1])
obs_pres <- st_as_sf(obs_pres, coords = c("long", "lat"), crs = 3163)

# Plotting for one specie
ggplot() + 
  geom_stars(data = theta[,,,1]) +
  ggtitle(paste('Interpolated current probabilities of presence \n for', name_species)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = c(0.2, 0.3)) +
  scale_fill_gradientn(colours = rev(rocket(5)), na.value = "transparent") +
  geom_sf(data = obs_pres, shape = 3, color = "white", size = 2)
  
##=============
##
## Predictive maps of species richness
##
##=============

list_theta <- list.files(here("output", "theta"), pattern = "RST_theta_", full.names = TRUE)
theta_sum <- sum(rast(list_theta[1]))
for (i in list_theta[2:npart])
{
 theta_sum <- theta_sum + sum(rast(i))
}
terra::writeRaster(theta_sum, here("output", "theta_sum.tif"))
theta_sum <- read_stars(here("output", "theta_sum.tif"))

ggplot() + 
  geom_stars(data = theta_sum) +
  ggtitle("Estimated current species richness") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = c(0.2, 0.3)) +
  scale_fill_gradientn(colours = rev(rocket(5)), na.value = "transparent") +
  theme_bw()
  
forest <- split(read_stars(here("output", "environNC.tif")))[1,,]
forest[forest < 50] <- NA
theta_forest <- theta_sum
theta_forest[is.na(forest)] <- NA

ggplot() + 
  geom_stars(data = theta_forest) +
  ggtitle("Estimated current species richness on forest area") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = c(0.2, 0.3)) +
  scale_fill_gradientn(colours = rev(rocket(5)), na.value = "transparent") 

##===============
##
## Predictive maps of species turnover
## diversity beta
##
##===============

# Sample 10000 cells of presence probabilities raster
npart <- 30
theta <- terra::rast(here("output", "theta", "RST_theta_01.tif"))
samp_theta <- spatSample(theta, size = 10000, cells = TRUE, na.rm = TRUE, method = "random")
samp_cells <- samp_theta[, "cell"]
samp_theta <- samp_theta[, -1]
for(n in 2:npart){
  theta <- terra::rast(here("output", "theta", paste0("RST_theta_",str_pad(n, 2, pad = "0") , ".tif")))
  samp_theta <- cbind(samp_theta, theta[samp_cells])
}
params_species <- read.csv(file = here("output", "params_species.csv"))
colnames(samp_theta) <- params_species$species

# Data frame to compute PCA on pixels dissimilarity 
theta_df <- data.frame(samp_theta)
pca_theta <- ade4::dudi.pca(theta_df, center = TRUE, scale = TRUE, nf = 3, scannf = FALSE)
save(pca_theta, samp_cells, file = here("output", "PCA_theta.RData"))
str(pca_theta)
fviz_eig(pca_theta)

# Find species with more effective in PCA
nb_species_effect <- 5
pca_powerfull_species <- matrix(0, ncol = 3, nrow = nb_species_effect)
colnames(pca_powerfull_species) <- c("Axis 1 (R)", "Axis 2 (G)", "Axis 3 (B)")
for (i in 1:3) {
  pca_powerfull_species[, i] <- names(sort(get_pca_var(pca_theta)$contrib[, i], decreasing = TRUE)[1:nb_species_effect])
  
}
for (l in 1:3) {
  for (k in 1:nb_species_effect) {
    pca_powerfull_species[k, l] <- NC_PIPPN$nom_taxon[NC_PIPPN$id_taxon_ref == substring(pca_powerfull_species[k, l], 2)][1]
  }
}

# Perform projections of PCA results on all raster's rows in parallel
# Coordinates on the 3 axes retained in the PCA
coords <- theta[[c(1, 2, 3)]]
values(coords) <- 0
names(coords) <- colnames(pca_theta$li)
coords <- st_as_stars(coords)
Sys.time() # 5 min to run

for(k in 1:nrow(theta)){
  theta_k <- NULL
  cat(k, "/", nrow(theta),"\n")
  # theta <- terra::rast(here("output", "theta", "RST_theta_01.tif"))
  ## Make a cluster for parallel computation
  # detect the number of CPU cores on the current host
  # ncores <- parallel::detectCores()
  # clust <- makeCluster(ncores - 2)
  # registerDoParallel(clust)
  # theta_k <- foreach(n = 1:npart, .combine = "cbind", .packages = c("here", "stringr")) %dopar%{
  for (n in 1:npart) {
    theta <- terra::rast(here("output","theta", paste0("RST_theta_", str_pad(n, 2, pad = "0"), ".tif")))
    theta_k <- cbind(theta_k, terra::values(theta, row = k, nrow = 1))
  }
  # return(theta_k)
  # }
  # stopCluster(clust)
  colnames(theta_k) <- params_species$species
  theta_k <- data.frame(theta_k)
  coords[[1]][,k,] <- as.matrix(ade4::suprow(pca_theta,theta_k)$lisup)
}
Sys.time()
write_stars(coords, dsn = here("output", "coords.tif"), options = c("COMPRESS=LZW", "PREDICTOR=2"))

# Change the coordinate scale for [0.255]. 
## Min reduced to 0

coords_RGB <- coords
for (l in 1:3) {
  Min <- min(read_stars(here("output", "coords.tif"))[,,,l][[1]], na.rm = TRUE)
  coords_RGB[[1]][,,l] <- coords[[1]][,,l]- Min
}

Min2 <- min(coords_RGB[[1]], na.rm = TRUE)
Min2
## Max at 255
remove(coords)

for (l in 1:3) {
  Max <- max(coords_RGB[,,,l][[1]], na.rm = TRUE)
  Max
  coords_RGB[[1]][,,l] <- round((coords_RGB[[1]][,,l] / Max)*255)
}
Max2 <- max(coords_RGB[[1]], na.rm = TRUE)
Max2

# Coloration RGB
write_stars(coords_RGB,dsn = here("output", "species_turnover.tif"), options = c("COMPRESS=LZW", "PREDICTOR=2"))

##================
##
## Current species turn over restricted to forest cover
##
##================

forest <- split(read_stars(here("output", "environNC.tif")))$forest
# current species turnover
# change resolution of species richness raster from 1000x1000m to 30x30m
# in_f <- here("output", "species_turnover.tif")
# out_f <- here("output", "species_turnover_30m.tif")
# system(glue('gdalwarp -tr 30 30 -te {xmin(forest)} {ymin(forest)} {xmax(forest)} {ymax(forest)} -r near  \\
#             -co "COMPRESS=LZW" -co "PREDICTOR=2" -overwrite {in_f} {out_f}'))

# remove species richness values where there was no forest in 2000s
species_turnover <- read_stars(here("output", "species_turnover.tif"))
species_turnover_forest <- species_turnover
forest[forest < 50] <- NA
for ( l in 1:3)
{
  species_turnover_forest[[1]][,,l][is.na(forest)] <- NA
}
write_stars(species_turnover_forest, dsn = here("output", "species_turnover_forest.tif"),
            options = c("COMPRESS=LZW", "PREDICTOR=2"), overwrite=TRUE)

# Representation of species turn over restricted to forest cover in 2000's

species_turnover <- terra::rast(here("output", "species_turnover.tif"))
species_turnover_forest <- terra::rast(here("output", "species_turnover_forest.tif"))
# par(mfrow=c(1,2), mar=c(1,1,2,1))
plotRGB(species_turnover)
title(main = "Current species turn over restricted to forest cover",
      cex.main = 1.5, line = -2)

