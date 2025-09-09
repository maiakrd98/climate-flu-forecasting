# This is me trying (to figure out Climate R)
library(AOI)
library(terra)
library(climateR)
library(tidyterra)
library(ggplot2)
library(tidyr)
library(sf)

cat <- catalog
AOI = aoi_get(state = "CO")
plot(AOI$geometry)

# this does work!! (for NC)
max_temps <- getGridMET(AOI, varname = "tmmx", startDate = "2012-10-01", endDate = "2012-10-30",dryrun = FALSE)
max_temps_by_date <- terra::extract(max_temps$daily_maximum_temperature,AOI)
###
#rmax

max_rhum <- getGridMET(AOI, varname = "rmax", startDate = "2012-10-01", endDate = "2012-10-30",dryrun = FALSE)
max_rhum_by_date <- terra::extract(max_rhum$daily_maximum_relative_humidity,AOI)
min_rhum <- getGridMET(AOI, varname = "rmin", startDate = "2012-10-01", endDate = "2012-10-30",dryrun = FALSE)
min_rhum_by_date <- terra::extract(min_rhum$daily_minimum_relative_humidity,AOI)

climater_dap(id, args, verbose, dryrun, print.arg = FALSE)

climater_filter(
  id = NULL,
  asset = NULL,
  AOI = NULL,
  startDate = NULL,
  endDate = NULL,
  varname = NULL,
  model = NULL,
  scenario = NULL,
  ensemble = NULL
)

get_data(dap)

# not working yet 
dap(
  AOI = AOI,
  startDate = "2014-10-29",
  endDate = "2015-10-29",
  varname = rhum,
  ID = esrlIcoads1ge
)

dap <- dap_crop(URL = 'https://coastwatch.pfeg.noaa.gov/erddap/griddap/esrlIcoads1ge', catalog = cat, AOI = AOI, 
                startDate = "2014-10-29", endDate = "2015-10-29", start = NULL, 
                end = NULL, varname = rhum, verbose = TRUE)
x = dap_get(dap)

max_rel_hum <- getMACA(
  AOI,
  varname = "rhsmax",
  timeRes = "month",
  model = "CCSM4",
  scenario = "rcp45",
  startDate = "2014-10-29",
  endDate = "2015-10-29",
  verbose = TRUE,
  ID = NULL,
  dryrun = FALSE
)