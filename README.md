
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `guyaclim` R Package <img src="man/figures/logo.svg" align="right" alt="" width="120" />

[![R-CMD-check](https://github.com/ghislainv/guyaclim/workflows/R-CMD-check/badge.svg)](https://github.com/ghislainv/guyaclim/actions)

`guyaclim` provides a set of climatic and environmental data for French
Guiana that can be used for ecological analyses. Data come from a
variety of sources (WorldClim, Chelsa, and Metéo France for the climate;
SRTM, BRGM, and various research projects for the environment). Climatic
and environmental data are available as multiband raster files at 1km
resolution in the geographic projection UTM 22N (epsg: 32622) and cover
all the territory of French Guiana (extent of the rasters is
xmin=100000, ymin=230000, xmax=435000, ymax=642000).

## System requirements

Make sure GDAL and GRASS GIS are installed on your system.

## Installation

You can install **guyaclim** from
[GitHub](https://github.com/ghislainv/guyaclim) with:

``` r
devtools::install_github("ghislainv/guyaclim")
```

## Contributing

The `guyaclim` R package is Open Source and released under the [GNU GPL
version 3](https://www.gnu.org/licenses/gpl-3.0.en.html) license.
Anybody who is interested can contribute to the package development
following our [Community guidelines](articles/Contributing.html). Every
contributor must agree to follow the project’s [Code of
Conduct](articles/Code_of_conduct.html).
