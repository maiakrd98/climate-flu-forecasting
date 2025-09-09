### Use weather and population data to predict triangle for the 2022-2023 season and use ARIMA to forecast deviations ###

# ==============================================================================
# Import libraries
# ==============================================================================

library(tidyverse) 
library(purrr)
library(plotly)
library(forecast)
#library(segmented)
library(ggplot2)
library(tidyr)
#library(urca)
#library(signal)
#library(ggdist) 
#library(viridis)
#library(Metrics)
#library(pairsD3)
library(AOI)
library(terra)
library(climateR)
library(tidyterra)
library(tidycensus)

# ==============================================================================
# Import data
# ==============================================================================

# percent visits due to ILI
ILI <- read.csv("ILINetByState.csv", header=TRUE, skip = 1) 
# percent of flu tests positive post 2015
PercentPos1 <- read.csv("PercentPosByState(post2016).csv", header=TRUE, skip = 1) 
# percent of flu tests positive pre 2015
PercentPos2 <- read.csv("PercentPosByState(pre2016).csv", header=TRUE, skip = 1) 
# 2010 population data
PopData2010 <- read.csv("2010PopData.csv", header=TRUE)

states <- unique(PercentPos1$REGION)
S = length(states)-2 # don't include Virgin Islands or New York City

StartYear = 2010
EndYear = 2019
N = EndYear-StartYear+1

# state areas in alphabetical order in mi^2
state_areas <- append(state.area, 3424) # add Puerto Rico
state_areas <- append(state_areas, 61, after=8) # add DC
# state latitudes in alphabetical order
latitudes <- state.center$y
latitudes <- append(latitudes, 18.2208) # add Puerto Rico
latitudes <- append(latitudes, 38.9072, after=8) # add DC

# import data from 2010-2019
flu_data <- read.csv("flu_data(3).csv")

# ==============================================================================
# Import climate and population data for 2022-2023 season
# ==============================================================================

# set years, states, and latitudes
# [(S)*(n-StartYear-2)+s] minus 2 since we're skipping two covid years
for (n in 2022:2022) {
  for (s in 1:(S)) {
    flu_data <- flu_data %>% add_row(Year = n, State = states[s], latitudes[s])
  }
}

# get temperature data by year and state
for (n in 2022:2022) {
  print(n)
  for (s in 1:(S)) {
    print(states[s])
    if (s==9) {
      AOI = aoi_get(fip=11)
    }
    else if (s==52) {
      AOI = aoi_get(country=states[s])
    }
    else {
      AOI = aoi_get(state=states[s])
    }
    # skip Alaska, Hawaii, and Puerto Rico for now (b/c GridMET data is non contiguous US phobic)
    if (s==2 | s==12 | s==52) {
      flu_data$MeanMaxTemp[(S)*(n-StartYear-2)+s] <- NA
      flu_data$OctMaxTemp[(S)*(n-StartYear-2)+s] <- NA
      flu_data$NovMaxTemp[(S)*(n-StartYear-2)+s] <- NA
      flu_data$DecMaxTemp[(S)*(n-StartYear-2)+s] <- NA
      flu_data$JanMaxTemp[(S)*(n-StartYear-2)+s] <- NA
      flu_data$FebMaxTemp[(S)*(n-StartYear-2)+s] <- NA
      flu_data$MarMaxTemp[(S)*(n-StartYear-2)+s] <- NA
      flu_data$AprMaxTemp[(S)*(n-StartYear-2)+s] <- NA
      flu_data$MayMaxTemp[(S)*(n-StartYear-2)+s] <- NA
      flu_data$JunMaxTemp[(S)*(n-StartYear-2)+s] <- NA
      flu_data$JulMaxTemp[(S)*(n-StartYear-2)+s] <- NA
      flu_data$AugMaxTemp[(S)*(n-StartYear-2)+s] <- NA
      flu_data$SepMaxTemp[(S)*(n-StartYear-2)+s] <- NA
    }
    else {
      max_temps <- getGridMET(AOI, varname = "tmmx", startDate = paste(as.character(n),sep="","-10-01"), endDate = paste(as.character(n+1),sep="","-09-30"),dryrun = FALSE)
      max_temps_by_date <- terra::extract(max_temps$daily_maximum_temperature,AOI)
      max_temps_by_date <- max_temps_by_date[,2:(366+(n %% 4 == 3))]
      #assign(paste("max_temps",as.character(n),states[s],sep="_"), max_temps_by_date) #get separate data frames for each year/state so that you can go back and subset them by month
      max_temps_oct <- max_temps_by_date[,1:31]
      max_temps_nov <- max_temps_by_date[,32:61]
      max_temps_dec <- max_temps_by_date[,62:92]
      max_temps_jan <- max_temps_by_date[,93:123]
      max_temps_feb <- max_temps_by_date[,124:(151+(n %% 4 == 3))]
      max_temps_mar <- max_temps_by_date[,(152+(n %% 4 == 3)):(182+(n %% 4 == 3))]
      max_temps_apr <- max_temps_by_date[,(183+(n %% 4 == 3)):(212+(n %% 4 == 3))]
      max_temps_may <- max_temps_by_date[,(213+(n %% 4 == 3)):(243+(n %% 4 == 3))]
      max_temps_jun <- max_temps_by_date[,(244+(n %% 4 == 3)):(273+(n %% 4 == 3))]
      max_temps_jul <- max_temps_by_date[,(274+(n %% 4 == 3)):(304+(n %% 4 == 3))]
      max_temps_aug <- max_temps_by_date[,(305+(n %% 4 == 3)):(335+(n %% 4 == 3))]
      max_temps_sep <- max_temps_by_date[,(336+(n %% 4 == 3)):(365+(n %% 4 == 3))]
      flu_data$MeanMaxTemp[(S)*(n-StartYear-2)+s] <- mean(unlist(max_temps_by_date),na.rm=TRUE)-273.15
      flu_data$OctMaxTemp[(S)*(n-StartYear-2)+s] <- mean(unlist(max_temps_oct),na.rm=TRUE)-273.15
      flu_data$NovMaxTemp[(S)*(n-StartYear-2)+s] <- mean(unlist(max_temps_nov),na.rm=TRUE)-273.15
      flu_data$DecMaxTemp[(S)*(n-StartYear-2)+s] <- mean(unlist(max_temps_dec),na.rm=TRUE)-273.15
      flu_data$JanMaxTemp[(S)*(n-StartYear-2)+s] <- mean(unlist(max_temps_jan),na.rm=TRUE)-273.15
      flu_data$FebMaxTemp[(S)*(n-StartYear-2)+s] <- mean(unlist(max_temps_feb),na.rm=TRUE)-273.15
      flu_data$MarMaxTemp[(S)*(n-StartYear-2)+s] <- mean(unlist(max_temps_mar),na.rm=TRUE)-273.15
      flu_data$AprMaxTemp[(S)*(n-StartYear-2)+s] <- mean(unlist(max_temps_apr),na.rm=TRUE)-273.15
      flu_data$MayMaxTemp[(S)*(n-StartYear-2)+s] <- mean(unlist(max_temps_may),na.rm=TRUE)-273.15
      flu_data$JunMaxTemp[(S)*(n-StartYear-2)+s] <- mean(unlist(max_temps_jun),na.rm=TRUE)-273.15
      flu_data$JulMaxTemp[(S)*(n-StartYear-2)+s] <- mean(unlist(max_temps_jul),na.rm=TRUE)-273.15
      flu_data$AugMaxTemp[(S)*(n-StartYear-2)+s] <- mean(unlist(max_temps_aug),na.rm=TRUE)-273.15
      flu_data$SepMaxTemp[(S)*(n-StartYear-2)+s] <- mean(unlist(max_temps_sep),na.rm=TRUE)-273.15
    }
  }
  #write.csv(data,"flu_data.csv",row.names = FALSE)
}

# get max relative humidity data by year and state
for (n in 2022:2022) {
  print(n)
  for (s in 1:(S)) {
    print(states[s])
    if (s==9) {
      AOI = aoi_get(fip=11)
    }
    else if (s==52) {
      AOI = aoi_get(country=states[s])
    }
    else {
      AOI = aoi_get(state=states[s])
    }
    # skip Alaska, Hawaii, and Puerto Rico for now (b/c GridMET data is non contiguous US phobic)
    if (s==2 | s==12 | s==52) {
      flu_data$MeanMaxRelHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$OctMaxRelHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$NovMaxRelHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$DecMaxRelHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$JanMaxRelHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$FebMaxRelHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$MarMaxRelHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$AprMaxRelHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$MayMaxRelHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$JunMaxRelHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$JulMaxRelHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$AugMaxRelHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$SepMaxRelHum[(S)*(n-StartYear-2)+s] <- NA
    }
    else {
      max_rhum <- getGridMET(AOI, varname = "rmax", startDate = paste(as.character(n),sep="","-10-01"), endDate = paste(as.character(n+1),sep="","-09-30"), dryrun = FALSE)
      max_rhum_by_date <- terra::extract(max_rhum$daily_maximum_relative_humidity,AOI)
      max_rhum_by_date <- max_rhum_by_date[,2:(366+(n %% 4 == 3))]
      #assign(paste("max_rhum",as.character(n),states[s],sep="_"), max_rhum_by_date)
      max_rhum_oct <- max_rhum_by_date[,1:31]
      max_rhum_nov <- max_rhum_by_date[,32:61]
      max_rhum_dec <- max_rhum_by_date[,62:92]
      max_rhum_jan <- max_rhum_by_date[,93:123]
      max_rhum_feb <- max_rhum_by_date[,124:(151+(n %% 4 == 3))]
      max_rhum_mar <- max_rhum_by_date[,(152+(n %% 4 == 3)):(182+(n %% 4 == 3))]
      max_rhum_apr <- max_rhum_by_date[,(183+(n %% 4 == 3)):(212+(n %% 4 == 3))]
      max_rhum_may <- max_rhum_by_date[,(213+(n %% 4 == 3)):(243+(n %% 4 == 3))]
      max_rhum_jun <- max_rhum_by_date[,(244+(n %% 4 == 3)):(273+(n %% 4 == 3))]
      max_rhum_jul <- max_rhum_by_date[,(274+(n %% 4 == 3)):(304+(n %% 4 == 3))]
      max_rhum_aug <- max_rhum_by_date[,(305+(n %% 4 == 3)):(335+(n %% 4 == 3))]
      max_rhum_sep <- max_rhum_by_date[,(336+(n %% 4 == 3)):(365+(n %% 4 == 3))]
      flu_data$MeanMaxRelHum[(S)*(n-StartYear-2)+s] <- mean(unlist(max_rhum_by_date),na.rm=TRUE)
      flu_data$OctMaxRelHum[(S)*(n-StartYear-2)+s] <- mean(unlist(max_rhum_oct),na.rm=TRUE)
      flu_data$NovMaxRelHum[(S)*(n-StartYear-2)+s] <- mean(unlist(max_rhum_nov),na.rm=TRUE)
      flu_data$DecMaxRelHum[(S)*(n-StartYear-2)+s] <- mean(unlist(max_rhum_dec),na.rm=TRUE)
      flu_data$JanMaxRelHum[(S)*(n-StartYear-2)+s] <- mean(unlist(max_rhum_jan),na.rm=TRUE)
      flu_data$FebMaxRelHum[(S)*(n-StartYear-2)+s] <- mean(unlist(max_rhum_feb),na.rm=TRUE)
      flu_data$MarMaxRelHum[(S)*(n-StartYear-2)+s] <- mean(unlist(max_rhum_mar),na.rm=TRUE)
      flu_data$AprMaxRelHum[(S)*(n-StartYear-2)+s] <- mean(unlist(max_rhum_apr),na.rm=TRUE)
      flu_data$MayMaxRelHum[(S)*(n-StartYear-2)+s] <- mean(unlist(max_rhum_may),na.rm=TRUE)
      flu_data$JunMaxRelHum[(S)*(n-StartYear-2)+s] <- mean(unlist(max_rhum_jun),na.rm=TRUE)
      flu_data$JulMaxRelHum[(S)*(n-StartYear-2)+s] <- mean(unlist(max_rhum_jul),na.rm=TRUE)
      flu_data$AugMaxRelHum[(S)*(n-StartYear-2)+s] <- mean(unlist(max_rhum_aug),na.rm=TRUE)
      flu_data$SepMaxRelHum[(S)*(n-StartYear-2)+s] <- mean(unlist(max_rhum_sep),na.rm=TRUE)
    }
    #write.csv(data,"flu_data.csv",row.names = FALSE)
  }
}


# get min relative humidity data by year and state
for (n in 2022:2022) {
  print(n)
  for (s in 1:(S)) {
    print(states[s])
    if (s==9) {
      AOI = aoi_get(fip=11)
    }
    else if (s==52) {
      AOI = aoi_get(country=states[s])
    }
    else {
      AOI = aoi_get(state=states[s])
    }
    # skip Alaska, Hawaii, and Puerto Rico for now (b/c GridMET data is non contiguous US phobic)
    if (s==2 | s==12 | s==52) {
      flu_data$MeanMinRelHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$OctMinRelHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$NovMinRelHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$DecMinRelHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$JanMinRelHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$FebMinRelHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$MarMinRelHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$AprMinRelHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$MayMinRelHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$JunMinRelHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$JulMinRelHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$AugMinRelHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$SepMinRelHum[(S)*(n-StartYear-2)+s] <- NA
    }
    else {
      min_rhum <- getGridMET(AOI, varname = "rmin", startDate = paste(as.character(n),sep="","-10-01"), endDate = paste(as.character(n+1),sep="","-09-30"), dryrun = FALSE)
      min_rhum_by_date <- terra::extract(min_rhum$daily_minimum_relative_humidity,AOI)
      min_rhum_by_date <- min_rhum_by_date[,2:(366+(n %% 4 == 3))]
      #assign(paste("min_rhum",as.character(n),states[s],sep="_"), min_rhum_by_date)
      min_rhum_oct <- min_rhum_by_date[,1:31]
      min_rhum_nov <- min_rhum_by_date[,32:61]
      min_rhum_dec <- min_rhum_by_date[,62:92]
      min_rhum_jan <- min_rhum_by_date[,93:123]
      min_rhum_feb <- min_rhum_by_date[,124:(151+(n %% 4 == 3))]
      min_rhum_mar <- min_rhum_by_date[,(152+(n %% 4 == 3)):(182+(n %% 4 == 3))]
      min_rhum_apr <- min_rhum_by_date[,(183+(n %% 4 == 3)):(212+(n %% 4 == 3))]
      min_rhum_may <- min_rhum_by_date[,(213+(n %% 4 == 3)):(243+(n %% 4 == 3))]
      min_rhum_jun <- min_rhum_by_date[,(244+(n %% 4 == 3)):(273+(n %% 4 == 3))]
      min_rhum_jul <- min_rhum_by_date[,(274+(n %% 4 == 3)):(304+(n %% 4 == 3))]
      min_rhum_aug <- min_rhum_by_date[,(305+(n %% 4 == 3)):(335+(n %% 4 == 3))]
      min_rhum_sep <- min_rhum_by_date[,(336+(n %% 4 == 3)):(365+(n %% 4 == 3))]
      flu_data$MeanMinRelHum[(S)*(n-StartYear-2)+s] <- mean(unlist(min_rhum_by_date),na.rm=TRUE)
      flu_data$OctMinRelHum[(S)*(n-StartYear-2)+s] <- mean(unlist(min_rhum_oct),na.rm=TRUE)
      flu_data$NovMinRelHum[(S)*(n-StartYear-2)+s] <- mean(unlist(min_rhum_nov),na.rm=TRUE)
      flu_data$DecMinRelHum[(S)*(n-StartYear-2)+s] <- mean(unlist(min_rhum_dec),na.rm=TRUE)
      flu_data$JanMinRelHum[(S)*(n-StartYear-2)+s] <- mean(unlist(min_rhum_jan),na.rm=TRUE)
      flu_data$FebMinRelHum[(S)*(n-StartYear-2)+s] <- mean(unlist(min_rhum_feb),na.rm=TRUE)
      flu_data$MarMinRelHum[(S)*(n-StartYear-2)+s] <- mean(unlist(min_rhum_mar),na.rm=TRUE)
      flu_data$AprMinRelHum[(S)*(n-StartYear-2)+s] <- mean(unlist(min_rhum_apr),na.rm=TRUE)
      flu_data$MayMinRelHum[(S)*(n-StartYear-2)+s] <- mean(unlist(min_rhum_may),na.rm=TRUE)
      flu_data$JunMinRelHum[(S)*(n-StartYear-2)+s] <- mean(unlist(min_rhum_jun),na.rm=TRUE)
      flu_data$JulMinRelHum[(S)*(n-StartYear-2)+s] <- mean(unlist(min_rhum_jul),na.rm=TRUE)
      flu_data$AugMinRelHum[(S)*(n-StartYear-2)+s] <- mean(unlist(min_rhum_aug),na.rm=TRUE)
      flu_data$SepMinRelHum[(S)*(n-StartYear-2)+s] <- mean(unlist(min_rhum_sep),na.rm=TRUE)
    }
    #write.csv(data,"flu_data.csv",row.names = FALSE)
  }
}


# get mean absolute humidity data by year and state
for (n in 2022:2022) {
  print(n)
  for (s in 1:(S)) {
    print(states[s])
    if (s==9) {
      AOI = aoi_get(fip=11)
    }
    else if (s==52) {
      AOI = aoi_get(country=states[s])
    }
    else {
      AOI = aoi_get(state=states[s])
    }
    # skip Alaska, Hawaii, and Puerto Rico for now (b/c GridMET data is non contiguous US phobic)
    if (s==2 | s==12 | s==52) {
      flu_data$MeanMeanAbsHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$OctMeanAbsHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$NovMeanAbsHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$DecMeanAbsHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$JanMeanAbsHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$FebMeanAbsHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$MarMeanAbsHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$AprMeanAbsHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$MayMeanAbsHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$JunMeanAbsHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$JulMeanAbsHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$AugMeanAbsHum[(S)*(n-StartYear-2)+s] <- NA
      flu_data$SepMeanAbsHum[(S)*(n-StartYear-2)+s] <- NA
    }
    else {
      mean_ahum <- getGridMET(AOI, varname = "sph", startDate = paste(as.character(n),sep="","-10-01"), endDate = paste(as.character(n+1),sep="","-09-30"), dryrun = FALSE)
      mean_ahum_by_date <- terra::extract(mean_ahum$daily_mean_specific_humidity,AOI)
      mean_ahum_by_date <- mean_ahum_by_date[,2:(366+(n %% 4 == 3))]
      #assign(paste("max_rhum",as.character(n),states[s],sep="_"), max_rhum_by_date)
      mean_ahum_oct <- mean_ahum_by_date[,1:31]
      mean_ahum_nov <- mean_ahum_by_date[,32:61]
      mean_ahum_dec <- mean_ahum_by_date[,62:92]
      mean_ahum_jan <- mean_ahum_by_date[,93:123]
      mean_ahum_feb <- mean_ahum_by_date[,124:(151+(n %% 4 == 3))]
      mean_ahum_mar <- mean_ahum_by_date[,(152+(n %% 4 == 3)):(182+(n %% 4 == 3))]
      mean_ahum_apr <- mean_ahum_by_date[,(183+(n %% 4 == 3)):(212+(n %% 4 == 3))]
      mean_ahum_may <- mean_ahum_by_date[,(213+(n %% 4 == 3)):(243+(n %% 4 == 3))]
      mean_ahum_jun <- mean_ahum_by_date[,(244+(n %% 4 == 3)):(273+(n %% 4 == 3))]
      mean_ahum_jul <- mean_ahum_by_date[,(274+(n %% 4 == 3)):(304+(n %% 4 == 3))]
      mean_ahum_aug <- mean_ahum_by_date[,(305+(n %% 4 == 3)):(335+(n %% 4 == 3))]
      mean_ahum_sep <- mean_ahum_by_date[,(336+(n %% 4 == 3)):(365+(n %% 4 == 3))]
      flu_data$MeanMeanAbsHum[(S)*(n-StartYear-2)+s] <- mean(unlist(mean_ahum_by_date),na.rm=TRUE)
      flu_data$OctMeanAbsHum[(S)*(n-StartYear-2)+s] <- mean(unlist(mean_ahum_oct),na.rm=TRUE)
      flu_data$NovMeanAbsHum[(S)*(n-StartYear-2)+s] <- mean(unlist(mean_ahum_nov),na.rm=TRUE)
      flu_data$DecMeanAbsHum[(S)*(n-StartYear-2)+s] <- mean(unlist(mean_ahum_dec),na.rm=TRUE)
      flu_data$JanMeanAbsHum[(S)*(n-StartYear-2)+s] <- mean(unlist(mean_ahum_jan),na.rm=TRUE)
      flu_data$FebMeanAbsHum[(S)*(n-StartYear-2)+s] <- mean(unlist(mean_ahum_feb),na.rm=TRUE)
      flu_data$MarMeanAbsHum[(S)*(n-StartYear-2)+s] <- mean(unlist(mean_ahum_mar),na.rm=TRUE)
      flu_data$AprMeanAbsHum[(S)*(n-StartYear-2)+s] <- mean(unlist(mean_ahum_apr),na.rm=TRUE)
      flu_data$MayMeanAbsHum[(S)*(n-StartYear-2)+s] <- mean(unlist(mean_ahum_may),na.rm=TRUE)
      flu_data$JunMeanAbsHum[(S)*(n-StartYear-2)+s] <- mean(unlist(mean_ahum_jun),na.rm=TRUE)
      flu_data$JulMeanAbsHum[(S)*(n-StartYear-2)+s] <- mean(unlist(mean_ahum_jul),na.rm=TRUE)
      flu_data$AugMeanAbsHum[(S)*(n-StartYear-2)+s] <- mean(unlist(mean_ahum_aug),na.rm=TRUE)
      flu_data$SepMeanAbsHum[(S)*(n-StartYear-2)+s] <- mean(unlist(mean_ahum_sep),na.rm=TRUE)
    }
    #write.csv(data,"flu_data.csv",row.names = FALSE)
  }
}

# get census data by year and state
for (n in 2022:2022) {
  for (s in 1:(S)) {
    acs_data_n <- get_acs(geography = "state", variables = c("S0101_C01_001","S0101_C01_022"), 
                          state = states[s], year=n, survey = 'acs1')
    total_pop <- acs_data_n$estimate[1]
    under_18 <- acs_data_n$estimate[2]
    flu_data$Under18[(S)*(n-StartYear-2)+s] = under_18/total_pop
    flu_data$PopDensity[(S)*(n-StartYear-2)+s] = total_pop/state_areas[s]
  }
}
