library(here)
# Download forest cover of Guyane in 2000 from tmf_ec_jrc
download.file('https://drive.google.com/uc?export=download&id=1FBL_Jy8QRi-wtG3rcsKZoJldzkwHyKn9',
              destfile=here("data_raw","tmf_ec_jrc","forest_t3.tif"), method = 'wget', mode="wb")
