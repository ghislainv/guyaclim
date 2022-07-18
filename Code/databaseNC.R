library(here)

dir.create(here("data_raw", "NCpippn"))
data <-read.csv2(here("data_raw", "NCpippn", "NC_PIPPN_latest.csv"), sep = "," )
data$X <- NULL
table(round(as.numeric(data$AREA)))
latlong <- data[!duplicated(data$plot_name), c("plot_name", "longitude", "latitude")]
write.csv(data, here("data_raw", "NCpippn", "data_clear.csv"))
write.csv(latlong, here("data_raw", "NCpippn", "coord_site.csv"))

##=====================
##
## Abundance & Presence Absence
##
##=====================

locality <- unique(data$plot_name)
taxon <- unique(data$taxaname)
Ab <- matrix(0, nrow = length(locality), ncol = length(taxon))
for (i in 1:length(locality))
{
  for (j in 1:length(taxon)) 
  {
    Ab[i,j] <- Ab[i,j] + sum(data$taxaname == taxon[j] & data$plot_name == locality[i])
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