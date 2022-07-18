library(stats)
library(here)
library(stars)
library(terra)
library(stringr)
library(EMCluster)
library(factoextra)
library(terrainr)
library(grid)

set.seed(1234)
load(here("output", "PCA_theta.RData"))
dir.create(here("output", "classification"))
##================
##
## Dendogram
## determine number of class
##================
dist_col <- dist(pca_theta$co)
plot(hclust(dist_col), labels = FALSE) # 5 groupes
nb_class <- 5
PCA_group <- cutree(hclust(dist_col), k = nb_class)
data_clear <- read.csv2(here("data_raw", "NCpippn", "data_clear.csv"), sep =",")
for (j in 1:nb_class) {
  group <- names(PCA_group)[PCA_group == j]
  for (k in 1:length(group)) {
    group[k] <- data_clear$nom_taxon[data_clear$id_taxon_ref == substring(group[k], 2)][1]
    writeLines(sort(group), here("output", "classification", paste0("group_PCA_", j, ".txt")))
  }
}
print(table(PCA_group))

##================
##
## Expectation Maximisation
##
##================
nb_class <- 5
init_EM <- rand.EM(pca_theta$co, nclass = nb_class, min.n = 10)
EM_class <- assign.class(pca_theta$co, emcluster(pca_theta$co, init_EM))
EM_group <- EM_class$class
names(EM_group) <- rownames(pca_theta$co)
data_clear <- read.csv2(here("data_raw", "NCpippn", "data_clear.csv"), sep =",")
for (j in 1:nb_class) {
  group <- names(EM_group)[EM_group == j]
  for (k in 1:length(group)) {
    group[k] <- data_clear$nom_taxon[data_clear$id_taxon_ref == substring(group[k], 2)][1]
    writeLines(sort(group), here("output", "classification", paste0("group_EM_", j, ".txt")))
  }
}
print(table(EM_group))

##================
##
## K-means
##
##================
nb_class <- 5
KM_class <- kmeans(pca_theta$co, centers = nb_class, iter.max = 20, nstart = 1000)
KM_group <- KM_class$cluster
data_clear <- read.csv2(here("data_raw", "NCpippn", "data_clear.csv"), sep =",")
for (j in 1:nb_class) {
  group <- names(KM_group)[KM_group == j]
  for (k in 1:length(group)) {
    group[k] <- data_clear$nom_taxon[data_clear$id_taxon_ref == substring(group[k], 2)][1]
    writeLines(sort(group), here("output", "classification", paste0("group_KM_", j, ".txt")))
  }
}
print(KM_class$size)

##================
##
## Comparaison between EM, KM and PCA
##
##================

## EM -> KM
similarity_EM_KM <- matrix(0, ncol = nb_class, nrow = 4)
colnames(similarity_EM_KM) <- paste0("EM_", 1:nb_class)
rownames(similarity_EM_KM) <- c("number of species", "number of species in common", " percentage of EM group", "number of KM group")
for (j in 1:nb_class) {
  EM <- readLines(here("output", "classification", paste0("group_EM_", j, ".txt")))
  nb_common <- rep(0, nb_class)
  for (k in 1:nb_class){
    KM <- readLines(here("output", "classification", paste0("group_KM_", k, ".txt")))
    nb_common[k] <- as.numeric(max(length(intersect(EM, KM)), 0))
  }
  similarity_EM_KM[1, j] <- length(EM)
  similarity_EM_KM[2, j] <- as.integer(max(nb_common, na.rm = TRUE))
  similarity_EM_KM[3, j] <- as.integer(similarity_EM_KM[2, j] / length(EM) * 100) 
  similarity_EM_KM[4, j] <- as.integer(which.max(nb_common))
}
print(similarity_EM_KM)

## KM -> EM
similarity_KM_EM <- matrix(0, ncol = nb_class, nrow = 4)
colnames(similarity_KM_EM) <- paste0("KM_", 1:nb_class)
rownames(similarity_KM_EM) <- c("number of species", "number of species in common", " percentage of KM group", "number of EM group")
for (j in 1:nb_class) {
  KM <- readLines(here("output", "classification", paste0("group_KM_", j, ".txt")))
  nb_common <- rep(0, nb_class)
  for (k in 1:nb_class){
    EM <- readLines(here("output", "classification", paste0("group_EM_", k, ".txt")))
    nb_common[k] <- as.numeric(max(length(intersect(EM, KM)), 0))
  }
  similarity_KM_EM[1, j] <- length(KM)
  similarity_KM_EM[2, j] <- as.integer(max(nb_common, na.rm = TRUE))
  similarity_KM_EM[3, j] <- as.integer(similarity_KM_EM[2, j] / length(KM) * 100) 
  similarity_KM_EM[4, j] <- as.integer(which.max(nb_common))
}
print(similarity_KM_EM)

## PCA -> EM
similarity_PCA_EM <- matrix(0, ncol = nb_class, nrow = 4)
colnames(similarity_PCA_EM) <- paste0("PCA_", 1:nb_class)
rownames(similarity_PCA_EM) <- c("number of species", "number of species in common", " percentage of PCA group", "number of EM group")
for (j in 1:nb_class) {
  PCA <- readLines(here("output", "classification", paste0("group_PCA_", j, ".txt")))
  nb_common <- rep(0, nb_class)
  for (k in 1:nb_class){
    EM <- readLines(here("output", "classification", paste0("group_EM_", k, ".txt")))
    nb_common[k] <- as.numeric(max(length(intersect(EM, PCA)), 0))
  }
  similarity_PCA_EM[1, j] <- length(PCA)
  similarity_PCA_EM[2, j] <- as.integer(max(nb_common, na.rm = TRUE))
  similarity_PCA_EM[3, j] <- as.integer(similarity_PCA_EM[2, j] / length(PCA) * 100) 
  similarity_PCA_EM[4, j] <- as.integer(which.max(nb_common))
}
print(similarity_PCA_EM)

## EM -> PCA
similarity_EM_PCA <- matrix(0, ncol = nb_class, nrow = 4)
colnames(similarity_EM_PCA) <- paste0("EM_", 1:nb_class)
rownames(similarity_EM_PCA) <- c("number of species", "number of species in common", " percentage of EM group", "number of PCA group")
for (j in 1:nb_class) {
  EM <- readLines(here("output", "classification", paste0("group_EM_", j, ".txt")))
  nb_common <- rep(0, nb_class)
  for (k in 1:nb_class){
    PCA <- readLines(here("output", "classification", paste0("group_PCA_", k, ".txt")))
    nb_common[k] <- as.numeric(max(length(intersect(EM, PCA)), 0))
  }
  similarity_EM_PCA[1, j] <- length(EM)
  similarity_EM_PCA[2, j] <- as.integer(max(nb_common, na.rm = TRUE))
  similarity_EM_PCA[3, j] <- as.integer(similarity_EM_PCA[2, j] / length(EM) * 100) 
  similarity_EM_PCA[4, j] <- as.integer(which.max(nb_common))
}
print(similarity_EM_PCA)

## PCA -> KM
similarity_PCA_KM <- matrix(0, ncol = nb_class, nrow = 4)
colnames(similarity_PCA_KM) <- paste0("PCA_", 1:nb_class)
rownames(similarity_PCA_KM) <- c("number of species", "number of species in common", " percentage of PCA group", "number of KM group")
for (j in 1:nb_class) {
  PCA <- readLines(here("output", "classification", paste0("group_PCA_", j, ".txt")))
  nb_common <- rep(0, nb_class)
  for (k in 1:nb_class){
    KM <- readLines(here("output", "classification", paste0("group_KM_", k, ".txt")))
    nb_common[k] <- as.numeric(max(length(intersect(KM, PCA)), 0))
  }
  similarity_PCA_KM[1, j] <- length(PCA)
  similarity_PCA_KM[2, j] <- as.integer(max(nb_common, na.rm = TRUE))
  similarity_PCA_KM[3, j] <- as.integer(similarity_PCA_KM[2, j] / length(PCA) * 100) 
  similarity_PCA_KM[4, j] <- as.integer(which.max(nb_common))
}
print(similarity_PCA_KM)

## KM -> PCA
similarity_KM_PCA <- matrix(0, ncol = nb_class, nrow = 4)
colnames(similarity_KM_PCA) <- paste0("KM_", 1:nb_class)
rownames(similarity_KM_PCA) <- c("number of species", "number of species in common", " percentage of KM group", "number of PCA group")
for (j in 1:nb_class) {
  KM <- readLines(here("output", "classification", paste0("group_KM_", j, ".txt")))
  nb_common <- rep(0, nb_class)
  for (k in 1:nb_class){
    PCA <- readLines(here("output", "classification", paste0("group_PCA_", k, ".txt")))
    nb_common[k] <- as.numeric(max(length(intersect(KM, PCA)), 0))
  }
  similarity_KM_PCA[1, j] <- length(KM)
  similarity_KM_PCA[2, j] <- as.integer(max(nb_common, na.rm = TRUE))
  similarity_KM_PCA[3, j] <- as.integer(similarity_KM_PCA[2, j] / length(KM) * 100) 
  similarity_KM_PCA[4, j] <- as.integer(which.max(nb_common))
}
print(similarity_KM_PCA)

##================
##
## get species who are always in same group
##
##================
print(similarity_EM_KM)
print(similarity_KM_EM)
print(similarity_PCA_KM)
print(similarity_KM_PCA)
print(similarity_PCA_EM)
print(similarity_EM_PCA)

nb_class <- 5
final_groups <- list(NULL)
nb_species_by_group <- rep(0, nb_class)
for (i in 1:5) {
  pca <- readLines(here("output", "classification", paste0("group_PCA_", i, ".txt")))
  em_group <- similarity_PCA_EM[4, i]
  km_group <- similarity_PCA_KM[4, i]
  em <- readLines(here("output", "classification", paste0("group_EM_", em_group, ".txt")))
  km <- readLines(here("output", "classification", paste0("group_KM_", km_group, ".txt")))
  common_species <- intersect(intersect(pca, em), km)
  nb_species_by_group[i] <- length(common_species)
  final_groups[[i]] <- common_species
  writeLines(sort(final_groups[[i]]), here("output", "classification", paste0("final_group_", i, ".txt")))
}
nb_species_by_group
sum(nb_species_by_group)


# methode bourine qui donne 28 groupes
# nb_class <- 5
# final_groups <- list(NULL)
# nb_species_by_group <- array(rep(0,125), c(5, 5, 5))
# for (i in 1:5) {
#   pca <- readLines(here("output", "classification", paste0("group_PCA_", i, ".txt")))
#   for (j in 1:5) {
#     em <- readLines(here("output", "classification", paste0("group_EM_", j, ".txt")))
#     for (k in 1:5) {
#       km <- readLines(here("output", "classification", paste0("group_KM_", k, ".txt")))
#       common_species <- intersect(intersect(pca, em), km)
#       nb_species_by_group[i, j, k] <- length(common_species)
#       # final_groups[[i]] <- common_species
#     }
#   }
# }

##================
##
## Color of each group 
##
##================

npart <- 30
nb_pixel_color <- 600
col_group <- matrix(0, nrow = nb_class, ncol = 3)
hex_color <- rep(0, nb_class)
for (i in 1:nb_class) {
  group_i <- readLines(here("output", "classification", paste0("final_group_", i, ".txt")))
  for (j in 1:length(group_i)){
    group_i[j] <- paste0("X", data_clear$id_taxon_ref[data_clear$nom_taxon == group_i[j]][1])
  }
  theta_group <- terra::rast(here("output", "theta", "RST_theta_forest_01.tif"))
  sum_theta_group <- theta_group[["X1644"]] - theta_group[["X1644"]] # init rast empty with right crs, extent,...
  for (k in 1:npart) {
    theta_group <- terra::rast(here("output", "theta", paste0("RST_theta_forest_", str_pad(k, width = 2, pad = "0"), ".tif")))
    common <- intersect(names(theta_group), group_i)
    if (length(common) != 0)
    {
      # stack probabilities of each species of each group
      sum_theta_group <- sum_theta_group + sum(theta_group[[common]])
    }
  }
  terra::writeRaster(sum_theta_group, here("output", "classification", paste0("theta_group_", i, ".tif")), overwrite = TRUE)
  sum_theta_group <- st_as_stars(sum_theta_group)
  # value_min <- sort(sum_theta_group[[1]], decreasing = TRUE)[nb_pixel_color]
  nb_in_group <- length(group_i)
  sum_theta_group[[1]][sum_theta_group[[1]] < 0.6 * nb_in_group] <- NA
  sum_theta_group_sf <- st_as_sf(sum_theta_group)
  write_stars(st_crop(read_stars(here("output", "species_turnover.tif")), sum_theta_group_sf),
              dsn = here("output", "classification", "species_group.tif"))
  
  plotRGB(terra::rast(here("output", "classification", "species_group.tif")), main = paste0("Main area of species from group ", i))
  RGB_map <- read_stars(here("output", "species_turnover.tif"))
  R <- RGB_map[[1]][,,1]
  G <- RGB_map[[1]][,,2]
  B <- RGB_map[[1]][,,3]
  # get values for each band RGB
  R <- R[!is.na(sum_theta_group[[1]])]
  G <- G[!is.na(sum_theta_group[[1]])]
  B <- B[!is.na(sum_theta_group[[1]])]
  
  # Find nearest color of mean color
  mindist <- which.min(sum(colMeans(cbind(R, G, B))^2) - rowSums(cbind(R, G, B)^2))

  R_near <- R[mindist] 
  G_near <- G[mindist]
  B_near <- B[mindist] 
  col_group[i, ] <- c(R_near, G_near, B_near)
  hex_color[i] <- rgb(red = col_group[i,1], green = col_group[i,2], blue = col_group[i,3], maxColorValue = 255)
}

ggplot() +
  geom_spatial_rgb(data = terra::rast(here("output", "species_turnover.tif")), 
                   aes(x = x, y = y, r = red, g = green, b = blue)) +
  geom_label(aes(x = 200000, y = 260000, label = "Group 1", size = 5), fill = hex_color[1]) +
  geom_label(aes(x = 200000, y = 245000, label = "Group 2", size = 5), fill = hex_color[2]) +
  geom_label(aes(x = 200000, y = 230000, label = "Group 3", size = 5), fill = hex_color[3]) +
  geom_label(aes(x = 200000, y = 215000, label = "Group 4", size = 5), fill = hex_color[4]) +
  geom_label(aes(x = 200000, y = 200000, label = "Group 5", size = 5), fill = hex_color[5]) +
  theme_bw() +
  theme(legend.position = "none")

# sum(!is.na(read_stars(here("output", "classification", "theta_group_1.tif"))[[1]])) = 6548

##================
##
## CAH on probabilities for each species and each pixels
##
##================

list_theta <- list.files(here("output", "theta"), pattern = "RST_theta_forest", full.names = TRUE)
load(here("output", "cell.RData"))
PA <- read.csv2(here("data_raw", "NCpippn", "Presence_Absence.csv"), sep = ",")
PA$X <- NULL
prob_all <- matrix(0, ncol = dim(PA)[2], nrow = dim(cell)[1])
for (j in 1:length(list_theta)) {
  theta <- read_stars(list_theta[j])[[1]]
  for (i in 1:dim(cell)[1]) {
    prob_all[i, (1 + dim(theta)[3] * (j - 1)): min(dim(theta)[3] * j, dim(PA)[2])] <- theta[cell[j, 1], cell[j, 2], ]
  }
}