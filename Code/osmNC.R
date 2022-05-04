library(osmdata)
get_overpass_url()
new_url <- "https://overpass.openstreetmap.ie/api/interpreter"
set_overpass_url(new_url)
bb <- getbb('Nouvelle Calédonie')
q <- opq(bbox = bb) %>%
  add_osm_feature(key = 'natural', value = 'water')