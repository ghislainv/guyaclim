# Libraries
library(here)

# Download forest cover of Guyane in 2000 from tmf_ec_jrc (see note)
download.file("https://drive.google.com/uc?export=download&id=1FBL_Jy8QRi-wtG3rcsKZoJldzkwHyKn9",
              destfile=here("data_raw", "tmf_ec_jrc", "forest_t3.tif"), method = 'auto', mode="wb")

# SRTM at 90m resolution from https://dwtkns.com/srtm/ version 4.1
download.file("https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/srtm_26_11.zip",
              destfile=here("data_raw", "srtm_v1_4_90m", "srtm_26_11.zip"), method = 'auto')
download.file("https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/srtm_26_12.zip",
              destfile=here("data_raw", "srtm_v1_4_90m", "srtm_26_12.zip"), method = 'auto')

# Reproject
#