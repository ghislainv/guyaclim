library(here)
library(stars)
library(sf)
library(sp)
library(ggplot2)
library(glue)
library(terra)
library("rnaturalearth") # plotting maps
library("rnaturalearthdata")
library("rnaturalearthhires")
library(rgrass7)
library(jSDM)
library(RColorBrewer)

dir.create(here("output", "CV"))
data <- read.csv2(here("data_raw", "NCpippn", "coord_site.csv"), sep = ",")
data$X <- NULL
border <- split(read_stars(here("output", "environNC.tif")))[1,,]
border[[1]][1100:1706,1:400] <- NA
border[[1]][1200:1706,1:570] <- NA
border[[1]][!is.na(border[[1]])] <- 1 
border <- st_as_sf(border)

# Create groups for Cross Validation Leave One Out
dist_site <- dist(data[,2:3], method = "canberra")
plot(hclust(dist_site), labels = FALSE)
nb_class <- 20
PCA_site_group <- cutree(hclust(dist_site), k = nb_class)
sort(table(PCA_site_group))
data$group <- PCA_site_group

# create color for each group  
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[1:nb_class]

# plot groups 
NC <- ne_countries(scale = 10, returnclass = "sf")
ggplot(data = NC) +
  geom_sf() +
  geom_point(data = data, aes(x = as.numeric(longitude), y = as.numeric(latitude),color = factor(group)), size = 1) +
  scale_color_manual(values = col_vector) +
  ggtitle("Groups of sites for CV") +
  xlab("Longitude") + 
  ylab("Latitude") +
  theme_bw() +
  coord_sf(xlim = c(163, 169 ), ylim = c(-22.84805679, -19.5), expand = FALSE) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

# import files for jSDM and CV

set.seed(1234)
EPSG <- 3163
nodat <- -9999
ISO_country_code <- "NCL"
proj.t <- paste0("EPSG:", EPSG) 
PA <- read.csv2(here("data_raw", "NCpippn", "Presence_Absence.csv"), sep = ",")
PA$X <- NULL
data_site <- read.csv2(here("output", "data_site.csv"), sep =",")
data_site$X <- NULL
data_site <- cbind(data_site, data)

sens <- rep(0, nb_class)
spe <- sens
TSS <- sens
for (i in 1:nb_class) {
  var_jSDM <- data.matrix(data_site[data$group != i,1:12])
  var_jSDM[,1] <- as.numeric(var_jSDM[,1] == 1)
  for (j in 2:12) {
    var_jSDM[, j] <- scale(var_jSDM[, j])
  }
  PA_group <- PA[data$group != i, ]
  
  
  Sys.time()
  jSDM_binom_pro <- jSDM_binomial_probit(
    burnin = 5000,
    mcmc = 10000,
    thin = 10,
    presence_data = data.matrix(PA_group),
    site_formula = ~.,
    site_data = var_jSDM,
    n_latent = 2,
    site_effect = "random",
    V_lambda = 1,
    V_beta = c(rep(1, 12), 0.001),
    mu_beta = c(rep(0, 12), 0.25),
    shape = 0.1,
    rate = 0.1,
    seed = 1234,
    verbose = 1
  )
  Sys.time()
  save(jSDM_binom_pro, file = here("output", "CV", paste0("jSDM_binom_pro_group_", i, ".RData")))
  load(here("output", "CV", paste0("jSDM_binom_pro_group_", i, ".RData")))
  
  ##====================
  ##
  ## Spatial Interpolation for alpha and latents variables
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
  params_sites <- data.frame(var_jSDM)
  coord <- data_site[data_site$group != i, 13:15]
  params_sites$latitude <- as.numeric(coord$latitude)
  params_sites$longitude <- as.numeric(coord$longitude)
  longlat <- SpatialPoints(params_sites[, c("longitude", "latitude")])
  proj4string(longlat) <- CRS("+proj=longlat +ellps=GRS80 +units=m")
  # lat-long to UTM58S projection 
  xy = spTransform(longlat, CRS("EPSG:3163"))
  
  #====
  # alpha 
  #====
  
  alpha_sp <- terra::vect(x = data.frame(alpha = colMeans(jSDM_binom_pro$mcmc.alpha[data$group != i,]),
                                         x = xy@coords[,1], y = xy@coords[,2]), geom = c("x", "y"),
                          crs = xy@proj4string@projargs)
  rgrass7::write_VECT(alpha_sp, "alpha")
  
  # Re-sample with RST
  # for punctual data use function v.surf.rst
  system("v.surf.rst --overwrite --verbose -t tension=3 input=alpha zcolumn=alpha \\
       smooth=0.0 elevation=alpha_rst")
  
  # Export
  system(glue('r.out.gdal --overwrite input=alpha_rst \\
             output={here("output", "CV", paste0("alphas_", i, ".tif"))} type=Float32 \\
             createopt="compress=lzw,predictor=2"'))
  
  # Representation
  alpha_in <- read_stars(here("output", "CV", paste0("alphas_", i, ".tif")))
  plot(st_crop(alpha_in, border), main = "Site effect alpha interpolated by RST")
  alpha_xy <- rep(0, dim(data_site)[1])
  for (j in 1:dim(data_site)[1]) 
  {
    alpha_xy[j] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "CV", paste0("alphas_", i, ".tif"))} \\
                                             -wgs84 {data_site[j,"longitude"]} {data_site[j,"latitude"]}'), intern = TRUE))
  }
  plot(alpha_xy[data_site$group != i], colMeans(jSDM_binom_pro$mcmc.alpha), xlab = "alpha interpolated by RST", ylab = "alpha estimated by JSDM",
       main = "Random site effect")
  abline(a = 0, b = 1, col = 'red')
  
  #====
  # W1 
  #====
  
  W1_sp <- terra::vect(x = data.frame(W1 = colMeans(jSDM_binom_pro$mcmc.latent$lv_1),
                                      x = xy@coords[,1], y = xy@coords[,2]), geom = c("x", "y"),
                       crs = xy@proj4string@projargs)
  rgrass7::write_VECT(W1_sp, "W1")
  
  # Re-sample with RST
  system("v.surf.rst --overwrite --verbose -t tension=3 input=W1 zcolumn=W1 \\
       smooth=0.0 elevation=W1_rst ")
  # Export
  system(glue('r.out.gdal --overwrite input=W1_rst \\
             output={here("output", "CV", paste0("lv_W1_", i, ".tif"))} type=Float32 \\
             createopt="compress=lzw,predictor=2"'))
  # Representation
  W1_in <- read_stars(here("output", "CV", paste0("lv_W1_", i, ".tif")))
  plot(st_crop(W1_in, border), main = "Latent variable W1 interpolated by RST")
  W1_xy <- rep(0, dim(data_site)[1])
  for (j in 1:dim(data_site)[1]) 
  {
    W1_xy[j] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "CV", paste0("lv_W1_", i, ".tif"))} \\
                                             -wgs84 {data_site[j,"longitude"]} {data_site[j,"latitude"]}'), intern = TRUE))
  }
  plot(W1_xy[data_site$group != i], colMeans(jSDM_binom_pro$mcmc.latent$lv_1), xlab = "W1 interpolated by RST", ylab = "W1 estimated by JSDM", 
       main = "First latent axe")
  abline(a = 0, b = 1, col = 'red')
  
  #====
  # W2 
  #====
  
  W2_sp <- terra::vect(x = data.frame(W2 = colMeans(jSDM_binom_pro$mcmc.latent$lv_2),
                                      x = xy@coords[,1], y = xy@coords[,2]), geom = c("x", "y"),
                       crs = xy@proj4string@projargs)
  rgrass7::write_VECT(W2_sp, "W2")
  
  # Re-sample with RST
  system("v.surf.rst --overwrite --verbose -t tension=3 input=W2 zcolumn=W2 \\
       smooth=0.0 elevation=W2_rst ")
  # Export
  system(glue('r.out.gdal --overwrite input=W2_rst \\
             output={here("output", "CV", paste0("lv_W2_", i, ".tif"))} type=Float32 \\
             createopt="compress=lzw,predictor=2"'))
  # Representation
  W2_in <- read_stars(here("output", "CV", paste0("lv_W2_", i, ".tif")))
  plot(st_crop(W2_in, border), main = "Latent variable W2 interpolated by RST")
  W2_xy <- rep(0, dim(data_site)[1])
  for (j in 1:dim(data_site)[1]) 
  {
    W2_xy[j] <- as.numeric(system(glue('gdallocationinfo -l_srs {proj.t} -valonly {here("output", "CV", paste0("lv_W2_", i, ".tif"))} \\
                                             -wgs84 {data_site[j,"longitude"]} {data_site[j,"latitude"]}'), intern = TRUE))
  }
  plot(W2_xy[data_site$group != i], colMeans(jSDM_binom_pro$mcmc.latent$lv_2), xlab = "W2 interpolated by RST", ylab = "W2 estimated by JSDM", 
       main = "Second latent axe")
  abline(a = 0, b = 1, col = 'red')
  
  # Lambda & Beta
  n_latent <- 2
  np <- 11
  nsp <- length(colnames(jSDM_binom_pro$theta_latent))
  lambdas <- matrix(0, nsp, n_latent)
  betas <- matrix(0, nsp, np + 2)
  for (j in 1:nsp){
    for (l in 1:n_latent){
      lambdas[j, l] <- mean(jSDM_binom_pro$mcmc.sp[[j]][, np + 2 + l])
    }
    for (p in 1:(np + 2)){
      betas[j, p] <- mean(jSDM_binom_pro$mcmc.sp[[j]][, p])
    }
  }
  colnames(betas) <- colnames(jSDM_binom_pro$mcmc.sp[[1]])[1:(np + 2)]
  
  var_pred <- data_site[data_site$group == i, ]
  result <- matrix(0, ncol = nsp, nrow = dim(var_pred)[1])
  for (k in 1:dim(var_pred)[1]) {
    result[k,] <- pnorm(betas[, 2:13] %*% as.numeric(var_pred[k, 1:12]) + alpha_xy[k] + betas[k, 1] + W1_xy[k] * lambdas[,1] + W2_xy[k] * lambdas[,2])
  }
  
  Sensitivity_CV <- function(PA, theta){
    n_sites <- nrow(theta)
    score <- rep(0, n_sites)
    for(i in 1:n_sites){
      # Sensitivity 
      obs_sp <- which(PA[i,] > 0)
      nobs_sp <- length(obs_sp)
      pred_sp <- which(theta[i,] >= sort(theta[i,], decreasing = TRUE)[nobs_sp])
      score[i] <- sum(pred_sp %in% obs_sp) / ifelse(nobs_sp != 0, nobs_sp, 1)
    }
    return(score)
  }
  Specificity_CV <- function(PA, theta){
    n_sites <- nrow(theta)
    score <- rep(0, n_sites)
    for(i in 1:n_sites){
      # Specificity 
      abs_sp <- which(PA[i,] == 0)
      nabs_sp <- length(abs_sp)
      pred_abs_sp <- which(theta[i,] <= sort(theta[i,])[nabs_sp])
      score[i] <- sum(pred_abs_sp %in% abs_sp) / ifelse(nabs_sp !=0, nabs_sp, 1)
    }
    return(score)
  }
  sens[i] = sum(Sensitivity_CV(PA, result)) / nrow(result)
  spe[i] = sum(Specificity_CV(PA, result)) / nrow(result)
  TSS[i] = sens[i] + spe[i] - 1
  Sys.time()
}
Sys.time()

##===================
##
## Analisys of each jSDM model 
##
##===================

for (i in 1:nb_class){
  readline()
  load(here("output", "CV", paste0("jSDM_binom_pro_group_", i, ".RData")))
  
  top_species <- which(colSums(PA) >= sort(colSums(PA), decreasing = TRUE)[5])
  np <- nrow(jSDM_binom_pro$model_spec$beta_start)
  
  # ## beta_j of the top five species
  # par(mfrow = c(3, 2))
  # for (j in top_species) {
  #   for (p in 1:np) {
  #     coda::traceplot(coda::as.mcmc(jSDM_binom_pro$mcmc.sp[[j]][, p]))
  #     coda::densplot(coda::as.mcmc(jSDM_binom_pro$mcmc.sp[[j]][, p]), 
  #                    main = paste(colnames(jSDM_binom_pro$mcmc.sp[[j]])[p],
  #                                 ", species : ", names(top_species[top_species == j])))
  #   }
  # }
  
  ## lambda_j of the top five species
  n_latent <- jSDM_binom_pro$model_spec$n_latent
  par(mfrow = c(2, 2))
  for (j in top_species) {
    for (l in 1:n_latent) {
      coda::traceplot(coda::as.mcmc(jSDM_binom_pro$mcmc.sp[[j]][, np + l]))
      coda::densplot(coda::as.mcmc(jSDM_binom_pro$mcmc.sp[[j]][, np + l]), 
                     main = paste(colnames(jSDM_binom_pro$mcmc.sp[[j]])
                                  [np + l],", species : ", names(top_species[top_species == j])))
    }
  }
  
  ## Latent variables W_i for the first two sites
  par(mfrow = c(2, 2))
  for (l in 1:n_latent) {
    for (i in 1:2) {
      coda::traceplot(jSDM_binom_pro$mcmc.latent[[paste0("lv_", l)]][, i],
                      main = paste0("Latent variable W_", l, ", site ", i))
      coda::densplot(jSDM_binom_pro$mcmc.latent[[paste0("lv_", l)]][, i],
                     main = paste0("Latent variable W_", l, ", site ", i))
    }
  }
  
  ## alpha_i of the first two sites
  plot(coda::as.mcmc(jSDM_binom_pro$mcmc.alpha[, 1:2]))
  
  ## V_alpha
  par(mfrow = c(2, 2))
  coda::traceplot(jSDM_binom_pro$mcmc.V_alpha)
  coda::densplot(jSDM_binom_pro$mcmc.V_alpha)
  ## Deviance
  coda::traceplot(jSDM_binom_pro$mcmc.Deviance)
  coda::densplot(jSDM_binom_pro$mcmc.Deviance)
  
  ## probit_theta
  par (mfrow = c(2, 1))
  hist(jSDM_binom_pro$probit_theta_latent,
       main = "Predicted probit theta", xlab = "predicted probit theta")
  hist(jSDM_binom_pro$theta_latent,
       main = "Predicted theta", xlab = "predicted theta")
  
}
