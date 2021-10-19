library(here)

# Standard (19) WorldClim Bioclimatic variables for WorldClim version 2.1. 
# They are the average for the years 1970-2000 at 30 seconds (~1 km2).
# Download a "zip" file containing 19 GeoTiff (.tif) files, one for each variable : 
# BIO1 = Annual Mean Temperature
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# BIO3 = Isothermality (BIO2/BIO7) (×100)
# BIO4 = Temperature Seasonality (standard deviation ×100)
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# BIO8 = Mean Temperature of Wettest Quarter
# BIO9 = Mean Temperature of Driest Quarter
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO12 = Annual Precipitation
# BIO13 = Precipitation of Wettest Month
# BIO14 = Precipitation of Driest Month
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
# BIO16 = Precipitation of Wettest Quarter
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter 
# BIO19 = Precipitation of Coldest Quarter
download.file('http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_bio.zip',
              destfile=here("data_raw","worldclim_v2_1","wc2.1_30s_bio.zip"), method = 'wget')

# Download "zip" files containing 12 GeoTiff (.tif) files, one for each month of the year (January is 1; December is 12).
# From WorldClim version 2.1 at 30 seconds (~1 km2).
## Monthly minimum temperature (°C).
download.file('http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_tmin.zip',
              destfile=here("data_raw","worldclim_v2_1","wc2.1_30s_tmin.zip"), method = 'wget')
## Monthly maximum temperature (°C).
download.file('http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_tmax.zip',
              destfile=here("data_raw","worldclim_v2_1","wc2.1_30s_tmax.zip"), method = 'wget')
## Average temperature (°C).
download.file('http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_tavg.zip',
              destfile=here("data_raw","worldclim_v2_1","wc2.1_30s_tavg.zip"), method = 'wget')
## Precipitation (mm).
download.file('http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_prec.zip',
              destfile=here("data_raw","worldclim_v2_1","wc2.1_30s_prec.zip"), method = 'wget')
## Solar radiation (kJ m^-2 day^-1).
download.file('http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_srad.zip',
              destfile=here("data_raw","worldclim_v2_1","wc2.1_30s_srad.zip"), method = 'wget')
## Wind speed (m s^-1)
download.file('http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_wind.zip',
              destfile=here("data_raw","worldclim_v2_1","wc2.1_30s_wind.zip"), method = 'wget')
## Water vapor pressure (kPa)
download.file('http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_vapr.zip',
              destfile=here("data_raw","worldclim_v2_1","wc2.1_30s_vapr.zip"), method = 'wget')