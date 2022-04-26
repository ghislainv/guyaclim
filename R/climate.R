## ==============================================================================
## author          :Ghislain Vieilledent, Jeanne Clément
## email           :ghislain.vieilledent@cirad.fr, jeanne.clement16@laposte.net
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

#' @name climate
#' @aliases climate 
#' @title Download and resample the climatic data for a given country or area. 
#' @description The \code{climate} Download and resample the climatic data for a given country or area.
#' @param ISO_country_code Specify the ISO 3166 code of the country whose climate data is required.
#' @param EPSG wanted projection for output climate data
#' @param area_borders shapefile giving the borders of the area whose climate data is required.
#' @param extent if specified ISO_country_code and area_borders are ignored, give extent corresponding to \code{EPSG} 
#' @param source Worldclim / Chelsa
#' @param future A Boolean which determines whether or not future climate data should be downloaded and resampled
#' @param GCM global climate model (GCM) data from CMIP5 (IPPC Fifth Assessment) 
#' @param RCP Future conditions are for two CO2 emission scenarios (RCP 4.5 and 8.5) 
#' @param period time periods 2050s and 2080s 
#' @param output_file if null return stars object 
#' @param verbose A Boolean which determines whether or not the download and resampling progress is printed on the screen. Default is \code{TRUE}.
#' @return A file in format .tif and type Int16 including : 
#'  \item{tmin.1-12}{.}
#'  \item{tmax.1-12}{.}
#'  \item{tavg.1-12}{.}
#'  \item{prec.1-12}{.}
#'  \item{bio.1-19}{}
#'  \item{pet.1-12}{.} 
#'  \item{pet}{}
#'  \item{ndm}{}
#'  \item{cwd}{}
#'
#' @details 
#' 
#' 
#' @examples 
#' #==============================================
#' # climate()
#' # Example on French Guiana
#' #==============================================
#' 
#' #=================
#' #== Load libraries
#' library(guyaclim)
#' 
#' @references
#' worldclim 
#'
#' Chelsa
#' 
#' @author 
#' Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
#'  
#' Jeanne Clément <jeanne.clement16@laposte.net>
#' 
#' @seealso \code{\link{get_extent}}, \code{\link{environ}}
#' @keywords climatic data 
#' @export

climate <- function() {
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
                destfile=here("data_raw","worldclim_v2_1","wc2.1_30s_bio.zip"), method = 'auto')
  
  # Download "zip" files containing 12 GeoTiff (.tif) files, one for each month of the year (January is 1; December is 12).
  # From WorldClim version 2.1 at 30 seconds (~1 km2).
  ## Monthly minimum temperature (°C).
  download.file('http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_tmin.zip',
                destfile=here("data_raw","worldclim_v2_1","wc2.1_30s_tmin.zip"), method = 'auto')
  ## Monthly maximum temperature (°C).
  download.file('http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_tmax.zip',
                destfile=here("data_raw","worldclim_v2_1","wc2.1_30s_tmax.zip"), method = 'auto')
  ## Average temperature (°C).
  download.file('http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_tavg.zip',
                destfile=here("data_raw","worldclim_v2_1","wc2.1_30s_tavg.zip"), method = 'auto')
  ## Precipitation (mm).
  download.file('http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_prec.zip',
                destfile=here("data_raw","worldclim_v2_1","wc2.1_30s_prec.zip"), method = 'auto')
  ## Solar radiation (kJ m^-2 day^-1).
  download.file('http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_srad.zip',
                destfile=here("data_raw","worldclim_v2_1","wc2.1_30s_srad.zip"), method = 'auto')
  ## Wind speed (m s^-1)
  download.file('http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_wind.zip',
                destfile=here("data_raw","worldclim_v2_1","wc2.1_30s_wind.zip"), method = 'auto')
  ## Water vapor pressure (kPa)
  download.file('http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_vapr.zip',
                destfile=here("data_raw","worldclim_v2_1","wc2.1_30s_vapr.zip"), method = 'auto')
}