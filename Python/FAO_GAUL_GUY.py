#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ecology.ghislainv.fr
# python_version  :>=3
# license         :GPLv3
# ==============================================================================


# Third party imports
import ee
ee.Initialize()

countries = ee.FeatureCollection("FAO/GAUL/2015/level0")
guyane = countries.filter(ee.Filter.eq("ADM0_CODE", 86))

# Export the FeatureCollection to a KML file.
task = ee.batch.Export.table.toDrive(
  collection=guyane,
  folder="Metradica",
  fileNamePrefix="FAO_GAUL_GUY",
  description="FAO_GAUL_GUY_py",
  fileFormat="KML")
task.start()

URL = guyane.getDownloadURL(filetype="KML", filename="FAO_GAUL_GUY")
# print(URL)
# https://earthengine.googleapis.com/v1alpha/projects/earthengine-legacy/tables/b41df7168964496b4375fb8f4b8d3631-3f6b51d518de231ca4fac96a1eca37e6:getFeatures
