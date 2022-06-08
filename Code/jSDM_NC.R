library(here)
library(jSDM)
library(stars)
library(glue)

EPSG <- 3163
nodat <- -9999
proj.t <- paste0("EPSG:", EPSG) 
PA <- read.csv2(here("data_raw", "NCpippn", "Presence_Absence.csv"), sep = ",")
Ab <- read.csv2(here("data_raw", "NCpippn", "Abondance.csv"), sep = ",")

# Merge climate and environ in one tif file

latlong <- read.csv2(here("data_raw", "NCpippn", "coord_site.csv"), sep = ",")
latlong$X <- NULL
current <- split(read_stars(here("output", "current_chelsaNC.tif")))
names(current) <-  c(paste("tmin",1:12,sep=""),paste("tmax", 1:12, sep=""), paste("tas",1:12, sep=""), paste("prec",1:12,sep=""),
                     paste("tcc", 1:12, sep=""), paste("pet", 1:12, sep=""), paste("bio",1:19,sep=""), "cwd", "ndm")

current <- merge(current)
environ <- split(read_stars(here("output", "environ_allNC.tif")))
names(environ) <- c("forest", "aspect", "elevation", "roughness", "slope", "srad", "soilgrids", "distanceForest",
                    "distanceSea", "distanceRoad", "distancePlace", "distancewater", "WDPA", "geology", "ultramafic")
environ <- merge(environ)
coord <- matrix(as.numeric(unlist(latlong)), nrow = nrow(latlong))[,2:3]
current_site <- matrix(0, nrow = dim(coord)[1], ncol = st_dimensions(current)$attribute$to)
environ_site <- matrix(0, nrow = dim(coord)[1], ncol = st_dimensions(environ)$attribute$to)

for (i in 1:dim(coord)[1]) 
{
  current_site[i,] <- as.integer(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "current_chelsaNC.tif")} \\
                                             -wgs84 {coord[i,1]} {coord[i,2]}'), intern = TRUE))
  environ_site[i,] <- as.integer(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "environ_allNC.tif")} \\
                                             -wgs84 {coord[i,1]} {coord[i,2]}'), intern = TRUE))
}
current_site <- data.frame(current_site)
environ_site <- data.frame(environ_site)
colnames(current_site) <- c(paste("tmin",1:12,sep=""),paste("tmax", 1:12, sep=""), paste("tas",1:12, sep=""), paste("prec",1:12,sep=""),
                            paste("tcc", 1:12, sep=""), paste("pet", 1:12, sep=""), paste("bio",1:19,sep=""), "cwd", "ndm")
colnames(environ_site) <- c("forest", "aspect", "elevation", "roughness", "slope", "srad", "soilgrids", "distanceForest",
                            "distanceSea", "distanceRoad", "distancePlace", "distancewater", "WDPA", "geology", "ultramafic")

# bio1, bio4, bio15, cwd, sum(pr), peridotites           
var_jSDM <- matrix(c(current_site$bio1, current_site$bio1^2, current_site$bio4, current_site$bio4^2, current_site$bio15,
                     current_site$bio15^2, current_site$cwd, current_site$cwd^2, rowSums(current_site[, paste0("prec", 1:12)]),
                     rowSums(current_site[, paste0("prec", 1:12)])^2, as.numeric(environ_site$ultramafic == 1)), nrow = dim(coord)[1])
write_stars(var_jSDM, options = c("COMPRESS=LZW","PREDICTOR=2"), NA_value = nodat,
            type = "Int16", dsn = here("output", "var_jSDM.tif"))
# Take main species for losing less time first 300
nb_species <- colSums(PA)
names(nb_species) <- NULL
nb_min_species <- sort(nb_species, decreasing = TRUE)[501]
PA <- PA[,colSums(PA) >= nb_min_species]
PA$X <- NULL

nb_species <- colSums(Ab)
names(nb_species) <- NULL
nb_min_species <- sort(nb_species, decreasing = TRUE)[501]
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
  seed = 1234,
  verbose = 1
)
Sys.time()

## beta_j of the five first species
top_species <- which(colSums(PA) >= sort(colSums(PA), decreasing = TRUE)[5])
colSums(PA)[top_species]
par(mfrow=c(3,2))
np <- nrow(jSDM_pois_log$model_spec$beta_start)
for (i in top_species) {
  for (p in 1:np) {
    coda::traceplot(coda::as.mcmc(jSDM_pois_log$mcmc.sp[[i]][,p]))
    coda::densplot(coda::as.mcmc(jSDM_pois_log$mcmc.sp[[i]][,p]), 
                   main = paste(colnames(jSDM_pois_log$mcmc.sp[[i]])[p], 
                                ", species : ", i))
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
                                [np+l], ", species : ",j))
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
##
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
  seed = 1234,
  verbose = 1
)
Sys.time()



top_species <- which(colSums(PA) >= sort(colSums(PA), decreasing = TRUE)[5])
np <- nrow(jSDM_binom_pro$model_spec$beta_start)

## beta_j of the top five species
par(mfrow=c(3,2))
for (j in top_species) {
  for (p in 1:np) {
    coda::traceplot(coda::as.mcmc(jSDM_binom_pro$mcmc.sp[[j]][,p]))
    coda::densplot(coda::as.mcmc(jSDM_binom_pro$mcmc.sp[[j]][,p]), 
                   main = paste(colnames(jSDM_binom_pro$mcmc.sp[[j]])[p],
                                ", species : ",j))
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
                                [np+l],", species : ",j))
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

