# Get climate and population data for each year and state
# This includes max temp, max relative humidity, min relative humidity, mean absolute humidity,
# population density, percent of population under 18, and latitude 

# to recover data after crash: data <- read.csv("../data/flu_data.csv")

library(segmented)
library(tidycensus)
library(tidyverse)
library(AOI)
library(terra)
library(climateR)
library(tidyterra)
library(ggplot2)
library(tidyr)
library(sf)

#census_api_key("YOUR API KEY GOES HERE")

# state areas in alphabetical order in mi^2
state_areas <- append(state.area, 3424) # add Puerto Rico
state_areas <- append(state_areas, 61, after=8) # add DC
# state latitudes in alphabetical order
latitudes <- state.center$y
latitudes <- append(latitudes, 18.2208) # add Puerto Rico
latitudes <- append(latitudes, 38.9072, after=8) # add DC

states <- unique(PercentPos1$REGION)
S = length(states)-2 # don't include Virgin Islands or New York City


# Note: we will call the 2014-2015 season the 2014 season, etc. 
# That is, the n season lasts from week 40 in year n to week 39 in year n+1

# start and end years
StartYear = 2010
EndYear = 2019
N = EndYear-StartYear+1
# data frame to hold peak weeks, end weeks, slopes, census, climate, etc.
data <- data.frame(Year = rep(0, N*S),
                   State = rep(states[1:S], N),
                   PopDensity = rep(0, N*S),
                   Under18 = rep(0, N*S),
                   Lat = rep(latitudes, N),
                   MeanMaxTemp = rep(0, N*S), # in Celsius
                   OctMaxTemp = rep(0, N*S), # mean max temp for Oct
                   NovMaxTemp = rep(0, N*S), # mean max temp for Nov
                   DecMaxTemp = rep(0, N*S), # mean max temp for Dec
                   JanMaxTemp = rep(0, N*S), # mean max temp for Jan
                   FebMaxTemp = rep(0, N*S), # mean max temp for Feb
                   MarMaxTemp = rep(0, N*S), # mean max temp for Mar
                   AprMaxTemp = rep(0, N*S), # mean max temp for Apr
                   MayMaxTemp = rep(0, N*S), # mean max temp for May
                   JunMaxTemp = rep(0, N*S), # mean max temp for Jun
                   JulMaxTemp = rep(0, N*S), # mean max temp for Jul
                   AugMaxTemp = rep(0, N*S), # mean max temp for Aug
                   SepMaxTemp = rep(0, N*S), # mean max temp for Sep
                   MeanMaxRelHum = rep(0, N*S), # Max relative humidity
                   OctMaxRelHum = rep(0, N*S), 
                   NovMaxRelHum = rep(0, N*S), 
                   DecMaxRelHum = rep(0, N*S), 
                   JanMaxRelHum = rep(0, N*S), 
                   FebMaxRelHum = rep(0, N*S), 
                   MarMaxRelHum = rep(0, N*S), 
                   AprMaxRelHum = rep(0, N*S), 
                   MayMaxRelHum = rep(0, N*S),
                   JunMaxRelHum = rep(0, N*S), 
                   JulMaxRelHum = rep(0, N*S),
                   AugMaxRelHum = rep(0, N*S),
                   SepMaxRelHum = rep(0, N*S),
                   MeanMinRelHum = rep(0, N*S), # Min relative humidity
                   OctMinRelHum = rep(0, N*S), 
                   NovMinRelHum = rep(0, N*S), 
                   DecMinRelHum = rep(0, N*S), 
                   JanMinRelHum = rep(0, N*S), 
                   FebMinRelHum = rep(0, N*S), 
                   MarMinRelHum = rep(0, N*S), 
                   AprMinRelHum = rep(0, N*S), 
                   MayMinRelHum = rep(0, N*S),
                   JunMinRelHum = rep(0, N*S), 
                   JulMinRelHum = rep(0, N*S),
                   AugMinRelHum = rep(0, N*S),
                   SepMinRelHum = rep(0, N*S)
)

# set years
for (n in StartYear:EndYear) {
  for (s in 1:(S)) {
    data$Year[(S)*(n-StartYear)+s] = n
  }
}
# save to csv in case of crash
write.csv(data,"../data/flu_data.csv",row.names = FALSE)

# get temperature data by year and state
for (n in StartYear:EndYear) {
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
      data$MeanMaxTemp[(S)*(n-StartYear)+s] <- NA
      data$OctMaxTemp[(S)*(n-StartYear)+s] <- NA
      data$NovMaxTemp[(S)*(n-StartYear)+s] <- NA
      data$DecMaxTemp[(S)*(n-StartYear)+s] <- NA
      data$JanMaxTemp[(S)*(n-StartYear)+s] <- NA
      data$FebMaxTemp[(S)*(n-StartYear)+s] <- NA
      data$MarMaxTemp[(S)*(n-StartYear)+s] <- NA
      data$AprMaxTemp[(S)*(n-StartYear)+s] <- NA
      data$MayMaxTemp[(S)*(n-StartYear)+s] <- NA
      data$JunMaxTemp[(S)*(n-StartYear)+s] <- NA
      data$JulMaxTemp[(S)*(n-StartYear)+s] <- NA
      data$AugMaxTemp[(S)*(n-StartYear)+s] <- NA
      data$SepMaxTemp[(S)*(n-StartYear)+s] <- NA
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
      data$MeanMaxTemp[(S)*(n-StartYear)+s] <- mean(unlist(max_temps_by_date),na.rm=TRUE)-273.15
      data$OctMaxTemp[(S)*(n-StartYear)+s] <- mean(unlist(max_temps_oct),na.rm=TRUE)-273.15
      data$NovMaxTemp[(S)*(n-StartYear)+s] <- mean(unlist(max_temps_nov),na.rm=TRUE)-273.15
      data$DecMaxTemp[(S)*(n-StartYear)+s] <- mean(unlist(max_temps_dec),na.rm=TRUE)-273.15
      data$JanMaxTemp[(S)*(n-StartYear)+s] <- mean(unlist(max_temps_jan),na.rm=TRUE)-273.15
      data$FebMaxTemp[(S)*(n-StartYear)+s] <- mean(unlist(max_temps_feb),na.rm=TRUE)-273.15
      data$MarMaxTemp[(S)*(n-StartYear)+s] <- mean(unlist(max_temps_mar),na.rm=TRUE)-273.15
      data$AprMaxTemp[(S)*(n-StartYear)+s] <- mean(unlist(max_temps_apr),na.rm=TRUE)-273.15
      data$MayMaxTemp[(S)*(n-StartYear)+s] <- mean(unlist(max_temps_may),na.rm=TRUE)-273.15
      data$JunMaxTemp[(S)*(n-StartYear)+s] <- mean(unlist(max_temps_jun),na.rm=TRUE)-273.15
      data$JulMaxTemp[(S)*(n-StartYear)+s] <- mean(unlist(max_temps_jul),na.rm=TRUE)-273.15
      data$AugMaxTemp[(S)*(n-StartYear)+s] <- mean(unlist(max_temps_aug),na.rm=TRUE)-273.15
      data$SepMaxTemp[(S)*(n-StartYear)+s] <- mean(unlist(max_temps_sep),na.rm=TRUE)-273.15
    }
  }
  #write.csv(data,"../data/flu_data.csv",row.names = FALSE)
}

# save to csv in case of crash
#write.csv(data,"../data/flu_data.csv",row.names = FALSE)

# get max relative humidity data by year and state
for (n in StartYear:EndYear) {
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
      data$MeanMaxRelHum[(S)*(n-StartYear)+s] <- NA
      data$OctMaxRelHum[(S)*(n-StartYear)+s] <- NA
      data$NovMaxRelHum[(S)*(n-StartYear)+s] <- NA
      data$DecMaxRelHum[(S)*(n-StartYear)+s] <- NA
      data$JanMaxRelHum[(S)*(n-StartYear)+s] <- NA
      data$FebMaxRelHum[(S)*(n-StartYear)+s] <- NA
      data$MarMaxRelHum[(S)*(n-StartYear)+s] <- NA
      data$AprMaxRelHum[(S)*(n-StartYear)+s] <- NA
      data$MayMaxRelHum[(S)*(n-StartYear)+s] <- NA
      data$JunMaxRelHum[(S)*(n-StartYear)+s] <- NA
      data$JulMaxRelHum[(S)*(n-StartYear)+s] <- NA
      data$AugMaxRelHum[(S)*(n-StartYear)+s] <- NA
      data$SepMaxRelHum[(S)*(n-StartYear)+s] <- NA
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
      data$MeanMaxRelHum[(S)*(n-StartYear)+s] <- mean(unlist(max_rhum_by_date),na.rm=TRUE)
      data$OctMaxRelHum[(S)*(n-StartYear)+s] <- mean(unlist(max_rhum_oct),na.rm=TRUE)
      data$NovMaxRelHum[(S)*(n-StartYear)+s] <- mean(unlist(max_rhum_nov),na.rm=TRUE)
      data$DecMaxRelHum[(S)*(n-StartYear)+s] <- mean(unlist(max_rhum_dec),na.rm=TRUE)
      data$JanMaxRelHum[(S)*(n-StartYear)+s] <- mean(unlist(max_rhum_jan),na.rm=TRUE)
      data$FebMaxRelHum[(S)*(n-StartYear)+s] <- mean(unlist(max_rhum_feb),na.rm=TRUE)
      data$MarMaxRelHum[(S)*(n-StartYear)+s] <- mean(unlist(max_rhum_mar),na.rm=TRUE)
      data$AprMaxRelHum[(S)*(n-StartYear)+s] <- mean(unlist(max_rhum_apr),na.rm=TRUE)
      data$MayMaxRelHum[(S)*(n-StartYear)+s] <- mean(unlist(max_rhum_may),na.rm=TRUE)
      data$JunMaxRelHum[(S)*(n-StartYear)+s] <- mean(unlist(max_rhum_jun),na.rm=TRUE)
      data$JulMaxRelHum[(S)*(n-StartYear)+s] <- mean(unlist(max_rhum_jul),na.rm=TRUE)
      data$AugMaxRelHum[(S)*(n-StartYear)+s] <- mean(unlist(max_rhum_aug),na.rm=TRUE)
      data$SepMaxRelHum[(S)*(n-StartYear)+s] <- mean(unlist(max_rhum_sep),na.rm=TRUE)
    }
    #write.csv(data,"../data/flu_data.csv",row.names = FALSE)
  }
}

# get min relative humidity data by year and state
for (n in StartYear:EndYear) {
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
      data$MeanMinRelHum[(S)*(n-StartYear)+s] <- NA
      data$OctMinRelHum[(S)*(n-StartYear)+s] <- NA
      data$NovMinRelHum[(S)*(n-StartYear)+s] <- NA
      data$DecMinRelHum[(S)*(n-StartYear)+s] <- NA
      data$JanMinRelHum[(S)*(n-StartYear)+s] <- NA
      data$FebMinRelHum[(S)*(n-StartYear)+s] <- NA
      data$MarMinRelHum[(S)*(n-StartYear)+s] <- NA
      data$AprMinRelHum[(S)*(n-StartYear)+s] <- NA
      data$MayMinRelHum[(S)*(n-StartYear)+s] <- NA
      data$JunMinRelHum[(S)*(n-StartYear)+s] <- NA
      data$JulMinRelHum[(S)*(n-StartYear)+s] <- NA
      data$AugMinRelHum[(S)*(n-StartYear)+s] <- NA
      data$SepMinRelHum[(S)*(n-StartYear)+s] <- NA
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
      data$MeanMinRelHum[(S)*(n-StartYear)+s] <- mean(unlist(min_rhum_by_date),na.rm=TRUE)
      data$OctMinRelHum[(S)*(n-StartYear)+s] <- mean(unlist(min_rhum_oct),na.rm=TRUE)
      data$NovMinRelHum[(S)*(n-StartYear)+s] <- mean(unlist(min_rhum_nov),na.rm=TRUE)
      data$DecMinRelHum[(S)*(n-StartYear)+s] <- mean(unlist(min_rhum_dec),na.rm=TRUE)
      data$JanMinRelHum[(S)*(n-StartYear)+s] <- mean(unlist(min_rhum_jan),na.rm=TRUE)
      data$FebMinRelHum[(S)*(n-StartYear)+s] <- mean(unlist(min_rhum_feb),na.rm=TRUE)
      data$MarMinRelHum[(S)*(n-StartYear)+s] <- mean(unlist(min_rhum_mar),na.rm=TRUE)
      data$AprMinRelHum[(S)*(n-StartYear)+s] <- mean(unlist(min_rhum_apr),na.rm=TRUE)
      data$MayMinRelHum[(S)*(n-StartYear)+s] <- mean(unlist(min_rhum_may),na.rm=TRUE)
      data$JunMinRelHum[(S)*(n-StartYear)+s] <- mean(unlist(min_rhum_jun),na.rm=TRUE)
      data$JulMinRelHum[(S)*(n-StartYear)+s] <- mean(unlist(min_rhum_jul),na.rm=TRUE)
      data$AugMinRelHum[(S)*(n-StartYear)+s] <- mean(unlist(min_rhum_aug),na.rm=TRUE)
      data$SepMinRelHum[(S)*(n-StartYear)+s] <- mean(unlist(min_rhum_sep),na.rm=TRUE)
    }
    #write.csv(data,"../data/flu_data.csv",row.names = FALSE)
  }
}

# get mean absolute humidity data by year and state
for (n in StartYear:EndYear) {
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
      data$MeanMeanAbsHum[(S)*(n-StartYear)+s] <- NA
      data$OctMeanAbsHum[(S)*(n-StartYear)+s] <- NA
      data$NovMeanAbsHum[(S)*(n-StartYear)+s] <- NA
      data$DecMeanAbsHum[(S)*(n-StartYear)+s] <- NA
      data$JanMeanAbsHum[(S)*(n-StartYear)+s] <- NA
      data$FebMeanAbsHum[(S)*(n-StartYear)+s] <- NA
      data$MarMeanAbsHum[(S)*(n-StartYear)+s] <- NA
      data$AprMeanAbsHum[(S)*(n-StartYear)+s] <- NA
      data$MayMeanAbsHum[(S)*(n-StartYear)+s] <- NA
      data$JunMeanAbsHum[(S)*(n-StartYear)+s] <- NA
      data$JulMeanAbsHum[(S)*(n-StartYear)+s] <- NA
      data$AugMeanAbsHum[(S)*(n-StartYear)+s] <- NA
      data$SepMeanAbsHum[(S)*(n-StartYear)+s] <- NA
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
      data$MeanMeanAbsHum[(S)*(n-StartYear)+s] <- mean(unlist(mean_ahum_by_date),na.rm=TRUE)
      data$OctMeanAbsHum[(S)*(n-StartYear)+s] <- mean(unlist(mean_ahum_oct),na.rm=TRUE)
      data$NovMeanAbsHum[(S)*(n-StartYear)+s] <- mean(unlist(mean_ahum_nov),na.rm=TRUE)
      data$DecMeanAbsHum[(S)*(n-StartYear)+s] <- mean(unlist(mean_ahum_dec),na.rm=TRUE)
      data$JanMeanAbsHum[(S)*(n-StartYear)+s] <- mean(unlist(mean_ahum_jan),na.rm=TRUE)
      data$FebMeanAbsHum[(S)*(n-StartYear)+s] <- mean(unlist(mean_ahum_feb),na.rm=TRUE)
      data$MarMeanAbsHum[(S)*(n-StartYear)+s] <- mean(unlist(mean_ahum_mar),na.rm=TRUE)
      data$AprMeanAbsHum[(S)*(n-StartYear)+s] <- mean(unlist(mean_ahum_apr),na.rm=TRUE)
      data$MayMeanAbsHum[(S)*(n-StartYear)+s] <- mean(unlist(mean_ahum_may),na.rm=TRUE)
      data$JunMeanAbsHum[(S)*(n-StartYear)+s] <- mean(unlist(mean_ahum_jun),na.rm=TRUE)
      data$JulMeanAbsHum[(S)*(n-StartYear)+s] <- mean(unlist(mean_ahum_jul),na.rm=TRUE)
      data$AugMeanAbsHum[(S)*(n-StartYear)+s] <- mean(unlist(mean_ahum_aug),na.rm=TRUE)
      data$SepMeanAbsHum[(S)*(n-StartYear)+s] <- mean(unlist(mean_ahum_sep),na.rm=TRUE)
    }
    #write.csv(data,"../data/flu_data.csv",row.names = FALSE)
  }
}

# get census data by year and state
#for (n in StartYear:EndYear) {
for (n in StartYear:EndYear) {
  for (s in 1:S) {
    if (n<2017) {
      acs_data_n <- get_acs(geography = "state", variables = c("S0101_C01_001","S0101_C01_025"), 
                            state = states[s], year=n, survey = 'acs1')
      total_pop <- acs_data_n$estimate[1]
      percent_18_and_over_18 <- acs_data_n$estimate[2]
      data$Under18[(S)*(n-StartYear)+s] = (100-percent_18_and_over_18)/100
      data$PopDensity[(S)*(n-StartYear)+s] = total_pop/state_areas[s]
    }
    else {
      acs_data_n <- get_acs(geography = "state", variables = c("S0101_C01_001","S0101_C01_022"), 
                            state = states[s], year=n, survey = 'acs1')
      total_pop <- acs_data_n$estimate[1]
      under_18 <- acs_data_n$estimate[2]
      data$Under18[(S)*(n-StartYear)+s] = under_18/total_pop
      data$PopDensity[(S)*(n-StartYear)+s] = total_pop/state_areas[s]
    }
  }
}
#write.csv(data,"../data/flu_data.csv",row.names = FALSE)