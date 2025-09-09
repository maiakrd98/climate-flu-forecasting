# Fitting a piecewise linear function to flu incidence and looking for correlations
# of parameters with population data

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

# Note: we will call the 2014-2015 season the 2014 season, etc. 
# That is, the n season lasts from week 40 in year n to week 39 in year n+1

# start and end years
StartYear = 2010
EndYear = 2019
N = EndYear-StartYear+1
# data frame to hold peak weeks, end weeks, slopes, census, climate, etc.
data <- data.frame(Year = rep(0, N*S),
                   State = rep(states[1:S], N),
                   Break1 = rep(0, N*S),
                   Break2 = rep(0, N*S),
                   Slope1 = rep(0, N*S),
                   Slope2 = rep(0, N*S),
                   Slope3 = rep(0, N*S),
                   Val1 = rep(0, N*S),
                   Val2 = rep(0, N*S),
                   Val3 = rep(0, N*S),
                   PeakWeek = rep(0, N*S),
                   UpSlope = rep(0, N*S),
                   DownSlope = rep(0, N*S),
                   PeakVal = rep(0, N*S),
                   PopDensity = rep(0, N*S),
                   Under18 = rep(0, N*S))

# get census data by year and state
for (n in StartYear:EndYear) {
  for (s in 1:(S)) {
    acs_data_n <- get_acs(geography = "state", variables = c("B01001_002","B01001_003","B01001_004","B01001_005",
                                                             "B01001_006","B01001_007","B01001_026","B01001_027",
                                                             "B01001_028","B01001_029","B01001_030","B01001_031"), 
                          state = states[s], year=n)
    total_pop <-acs_data_n$estimate[1]+acs_data_n$estimate[7]
    under_18 <- sum(acs_data_n$estimate[2:5])+sum(acs_data_n$estimate[8:11])
    data$Under18[(S)*(n-StartYear)+s] = under_18/total_pop
    data$PopDensity[(S)*(n-StartYear)+s] = total_pop/state_areas[s]
  }
}

# fit and plot data by year and state
for (n in StartYear:EndYear) {
  for (s in 1:(S)) {
    data$Year[(S)*(n-StartYear)+s] = n
    ILIn <- ILI[((ILI$YEAR==n & ILI$WEEK>=40) | (ILI$YEAR==n+1 & ILI$WEEK<40)) & ILI$REGION==states[s],]
    if (n<2015) {
      PercentPosn <- PercentPos2[((PercentPos2$YEAR==n & PercentPos2$WEEK>=40) | (PercentPos2$YEAR==n+1 & PercentPos2$WEEK<40)) & PercentPos2$REGION==states[s],]
    }
    else{
      PercentPosn <- PercentPos1[((PercentPos1$YEAR==n & PercentPos1$WEEK>=40) | (PercentPos1$YEAR==n+1 & PercentPos1$WEEK<40)) & PercentPos1$REGION==states[s],]
    }
    # skip if more than 5 missing data points?!?!? or more that 40 weeks of zeros??
    if(sum(PercentPosn$PERCENT.POSITIVE=="X")>5 | sum(PercentPosn$PERCENT.POSITIVE==0)>40) {
      data$Break1[(S)*(n-StartYear)+s]=NA
      data$Break2[(S)*(n-StartYear)+s]=NA
      data$Slope1[(S)*(n-StartYear)+s]=NA
      data$Slope2[(S)*(n-StartYear)+s]=NA
      data$Slope3[(S)*(n-StartYear)+s]=NA
      data$Val1[(S)*(n-StartYear)+s]=NA
      data$Val2[(S)*(n-StartYear)+s]=NA
      data$Val3[(S)*(n-StartYear)+s]=NA
      data$PeakWeek[(S)*(n-StartYear)+s]=NA
      data$UpSlope[(S)*(n-StartYear)+s]=NA
      data$DownSlope[(S)*(n-StartYear)+s]=NA
      data$PeakVal[(S)*(n-StartYear)+s]=NA
    }
    
    else{
      IncidenceN <- as.numeric(ILIn$X.UNWEIGHTED.ILI)*as.numeric(PercentPosn$PERCENT.POSITIVE)
      LogIncidenceN <- log(IncidenceN+0.1*(IncidenceN==0))
      week <- 1:nrow(ILIn)
      
      # delete -Inf's and NA's (bad if too many?!?!?)
      NAIndicies <- c()
      InfIndicies <- c()
      for (w in week){
        if (is.na(LogIncidenceN[w])) {NAIndicies <- c(NAIndicies,w)}
        else if (LogIncidenceN[w]==-Inf) {InfIndicies <- c(InfIndicies,w)}
      }
      
      if (!is.null(NAIndicies) & !is.null(InfIndicies)) {
        week <- week[-c(NAIndicies,InfIndicies)]
        LogIncidenceN <- LogIncidenceN[-c(NAIndicies,InfIndicies)]
      }
      else if (!is.null(NAIndicies)) {
        week <- week[-NAIndicies]
        LogIncidenceN <- LogIncidenceN[-NAIndicies]
      }
      else if (!is.null(InfIndicies)) {
        week <- week[-InfIndicies]
        LogIncidenceN <- LogIncidenceN[-InfIndicies]
      }
      
      # first fit regular linear regression
      lmfitN <- lm(LogIncidenceN~week)
      
      #fit piecewise regression model to original model
      segmented.fitN <- segmented(lmfitN, seg.Z = ~week, npsi=2, control=seg.control(quant=TRUE))
      
      #keep track of breakpoints (peak/end week), slopes
      sum <- summary.segmented(segmented.fitN)
      data$Break1[(S)*(n-StartYear)+s] <- sum$psi[1,2]
      data$Break2[(S)*(n-StartYear)+s] <- sum$psi[2,2]
      data$Slope1[(S)*(n-StartYear)+s] <- sum$coefficients[2,1]
      data$Slope2[(S)*(n-StartYear)+s] <- sum$coefficients[2,1]+sum$coefficients[3,1]
      data$Slope3[(S)*(n-StartYear)+s] <- sum$coefficients[2,1]+sum$coefficients[3,1]+sum$coefficients[4,1]
      data$Val1[(S)*(n-StartYear)+s] <- sum$coefficients[1,1]
      
      #calculate peak/end value
      data$Val2[(S)*(n-StartYear)+s] <- data$Val1[(S)*(n-StartYear)+s]+data$Slope1[(S)*(n-StartYear)+s]*data$Break1[(S)*(n-StartYear)+s]
      data$Val3[(S)*(n-StartYear)+s] <- data$Val2[(S)*(n-StartYear)+s]+data$Slope2[(S)*(n-StartYear)+s]*(data$Break2[(S)*(n-StartYear)+s]-data$Break1[(S)*(n-StartYear)+s])
      data$PeakVal[(S)*(n-StartYear)+s] <- max(data$Val2[(S)*(n-StartYear)+s],data$Val3[(S)*(n-StartYear)+s])
      if (data$PeakVal[(S)*(n-StartYear)+s]==data$Val2[(S)*(n-StartYear)+s]) {
        data$UpSlope[(S)*(n-StartYear)+s] <- data$Slope1[(S)*(n-StartYear)+s]
        data$DownSlope[(S)*(n-StartYear)+s] <- data$Slope2[(S)*(n-StartYear)+s]
        data$PeakWeek[(S)*(n-StartYear)+s] <- data$Break1[(S)*(n-StartYear)+s]
      }
      if (data$PeakVal[(S)*(n-StartYear)+s]==data$Val3[(S)*(n-StartYear)+s]) {
        data$UpSlope[(S)*(n-StartYear)+s] <- data$Slope2[(S)*(n-StartYear)+s]
        data$DownSlope[(S)*(n-StartYear)+s] <- data$Slope3[(S)*(n-StartYear)+s]
        data$PeakWeek[(S)*(n-StartYear)+s] <- data$Break2[(S)*(n-StartYear)+s]
      }
      
      #plot original data
      #plot(week, LogIncidenceN, main = sprintf("Flu incidence in %s over the %i-%i season", states[s], n,n+1))
      
      #add segmented regression model
      #plot(segmented.fitN, add=T)
    }
  }
}

# check for correlations
# also need to test significance and do a linear model or something...
cor(data[, c('PeakWeek','UpSlope','DownSlope','PeakVal','PopDensity','Under18')], use = "complete.obs")

library(Hmisc)
res <- rcorr(as.matrix(data[, c('PeakWeek','UpSlope','DownSlope','PeakVal','PopDensity','Under18')])) # rcorr() accepts matrices only

# display p-values (rounded to 3 decimals)
round(res$P, 4)

# scatterplots
install.packages("pairsD3")
library(pairsD3)

colors = rep("black",S)
colors[7] = "red"
pairs(data[, c('PeakWeek','UpSlope','DownSlope','PeakVal','PopDensity','Under18')],col=hcl.colors(52, "Temps"))
pairs(data[, c('PeakWeek','UpSlope','DownSlope','PeakVal','PopDensity','Under18')],col=colors)
shinypairs(data[, c('PeakWeek','UpSlope','DownSlope','PeakVal','PopDensity','Under18','State','Year')])

# compare variation b/w and w/i years/states
sd_by_state = data.frame('PeakWeek' = rep(0,S+2),
                         'UpSlope' = rep(0,S+2),
                         'DownSlope' = rep(0,S+2),
                         'PeakVal' = rep(0,S+2),
                         'State' = c(states[1:S],'Average','Overall'))
for (s in 1:S) {
  sd_by_state$PeakWeek[s]=sd(data[data$State == states[s],'PeakWeek'],na.rm = TRUE)
  sd_by_state$UpSlope[s]=sd(data[data$State == states[s],'UpSlope'],na.rm = TRUE)
  sd_by_state$DownSlope[s]=sd(data[data$State == states[s],'DownSlope'],na.rm = TRUE)
  sd_by_state$PeakVal[s]=sd(data[data$State == states[s],'PeakVal'],na.rm = TRUE)
}
# average sd within states
sd_by_state$PeakWeek[S+1]=mean(sd_by_state$PeakWeek[1:S],na.rm = TRUE)
sd_by_state$UpSlope[S+1]=mean(sd_by_state$UpSlope[1:S],na.rm = TRUE)
sd_by_state$DownSlope[S+1]=mean(sd_by_state$DownSlope[1:S],na.rm = TRUE)
sd_by_state$PeakVal[S+1]=mean(sd_by_state$PeakVal[1:S],na.rm = TRUE)
# overall sd
sd_by_state$PeakWeek[S+2]=sd(data[,'PeakWeek'],na.rm = TRUE)
sd_by_state$UpSlope[S+2]=sd(data[,'UpSlope'],na.rm = TRUE)
sd_by_state$DownSlope[S+2]=sd(data[,'DownSlope'],na.rm = TRUE)
sd_by_state$PeakVal[S+2]=sd(data[,'PeakVal'],na.rm = TRUE)

sd_by_year = data.frame('PeakWeek' = rep(0,N+2),
                        'UpSlope' = rep(0,N+2),
                        'DownSlope' = rep(0,N+2),
                        'PeakVal' = rep(0,N+2),
                        'Year' = c(StartYear:EndYear,'Average','Overall'))
for (n in StartYear:EndYear) {
  sd_by_year$PeakWeek[n-StartYear+1]=sd(data[data$Year == n,'PeakWeek'],na.rm = TRUE)
  sd_by_year$UpSlope[n-StartYear+1]=sd(data[data$Year == n,'UpSlope'],na.rm = TRUE)
  sd_by_year$DownSlope[n-StartYear+1]=sd(data[data$Year == n,'DownSlope'],na.rm = TRUE)
  sd_by_year$PeakVal[n-StartYear+1]=sd(data[data$Year == n,'PeakVal'],na.rm = TRUE)
}
# average sd within years
sd_by_year$PeakWeek[N+1]=mean(sd_by_year$PeakWeek[1:N],na.rm = TRUE)
sd_by_year$UpSlope[N+1]=mean(sd_by_year$UpSlope[1:N],na.rm = TRUE)
sd_by_year$DownSlope[N+1]=mean(sd_by_year$DownSlope[1:N],na.rm = TRUE)
sd_by_year$PeakVal[N+1]=mean(sd_by_year$PeakVal[1:N],na.rm = TRUE)
# overall sd
sd_by_year$PeakWeek[N+2]=sd(data[,'PeakWeek'],na.rm = TRUE)
sd_by_year$UpSlope[N+2]=sd(data[,'UpSlope'],na.rm = TRUE)
sd_by_year$DownSlope[N+2]=sd(data[,'DownSlope'],na.rm = TRUE)
sd_by_year$PeakVal[N+2]=sd(data[,'PeakVal'],na.rm = TRUE)