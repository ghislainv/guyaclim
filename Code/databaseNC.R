library(here)

dir.create(here("data_raw", "NCpippn"))
# url <- "insert url"
# download.file(url, destfile, method = "auto")
destfile <- here("data_raw", "NCpippn", "NC_PIPPN_2022.csv")

dbs <-read.csv2(destfile, sep = "," )
data <- dbs[,c("id_individu", "id_locality", "id_taxon_ref", "plot", "dbh", "longitude",
               "latitude", "nom_taxon", "iscurrent", "isliving")]
# rm(dbs)
data <- data[data$isliving,]
data <- data[data$iscurrent, c("id_individu", "id_locality", "id_taxon_ref",
                               "plot", "dbh", "longitude", "latitude", "nom_taxon")]
format <- unique(data$plot)
aires <- round(c(20*20, pi*11.28^2, 100*100, pi*10^2, 20*120, 30*100,
                 70*70, 50*50, 40*60, 10*10, 10*20, 50*100, 20*50))
for ( i in 1:length(data$plot)){
  data$aire[i] <- aires[data$plot[i] == format]
}
data$dbh <- as.numeric(data$dbh)
data <- data[!(is.na(data$id_taxon_ref) * (data$nom_taxon == "")),]
data <- data[!(data$dbh < 10),]
data <- data[!is.na(data$latitude),]
aires <- sort(unique(aires))
nb_aire <- 1:length(aires)
nb_indi_aire <- 1:length(aires)

for (i in 1:length(aires))
{
  nb_indi_aire[i] <- sum(data$aire == aires[i])  
}

latlong <- data[!duplicated(data$id_locality), c("id_locality", "longitude", "latitude")]
rownames(latlong) <- NULL
write.csv(latlong, here("data_raw", "NCpippn", "coord_site.csv"))
# rm("data", "dbs")

##=====================
##
## Abundance & Presence Absence
##
##=====================

locality <- unique(data$id_locality)
taxon <- unique(data$id_taxon_ref)
Ab <- matrix(0, nrow = length(locality), ncol = length(taxon))
for (i in 1:length(locality))
{
  for (j in 1:length(taxon)) 
  {
    Ab[i,j] <- Ab[i,j] + sum(data$id_taxon_ref == taxon[j] & data$id_locality == locality[i])
  }
  print(i)
}
Ab <- data.frame(Ab)
colnames(Ab) <- taxon
rownames(Ab) <- locality

PA <- Ab
PA <- 1 * (PA != 0)
write.csv(Ab, file = here("data_raw", "NCpippn", "Abondance.csv"))
write.csv(PA, file = here("data_raw", "NCpippn", "Presence_Absence.csv"))

##=====================
##
## Ajust to same Area
##
##=====================
# 
# # roundAb <- Ab
# area_locality <- data$aire[!duplicated(data$id_locality)]
# 
# # Proportion to m² then round 
# # roundAb <- data.frame(diag(min(aires) / area_locality) %*% as.matrix(roundAb))
# # roundAb <- round(roundAb)
# # colnames(roundAb) <- taxon
# # rownames(roundAb) <- locality
# # sum(Ab != 0) = 3009
# # sum(roundAb != 0) = 11377
# # loss 73% of individus
# 
# # Mean individus by area
# Rsum <- rowSums(roundAb)
# mhu <- (nb_indi_aire / aires * min(aires)) / table(area_locality)
# sigma <- rep(0, length(aires))
# for (i in 1:length(aires)) {
#   sigma[i] <- sum((Rsum[area_locality == aires[i]] - mhu[i])^2) / table(area_locality)[i]
# }
# AreaTotal <- matrix(c(aires, table(area_locality), nb_indi_aire , mhu, sigma), ncol = 5, byrow = FALSE)
# colnames(AreaTotal) <- c("area", "number of places", "number of individus", "mean of individus by 10*10",  "Sd")
# 
# plot(AreaTotal[,1], AreaTotal[,4] * AreaTotal[,1] / 100, xlab = "Aire in m²", ylab = "Mean number of individus")
