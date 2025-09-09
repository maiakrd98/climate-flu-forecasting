# Fitting a piecewise linear function to flu incidence by state and extracting breakpoints
# and slopes

library(segmented)

# percent visits due to ILI
ILI <- read.csv("ILINetByState.csv", header=TRUE, skip = 1) 
# percent of flu tests positive post 2015
PercentPos1 <- read.csv("PercentPosByState(post2016).csv", header=TRUE, skip = 1) 
# percent of flu tests positive pre 2015
PercentPos2 <- read.csv("PercentPosByState(pre2016).csv", header=TRUE, skip = 1) 
# 2010 population data
PopData2010 <- read.csv("2010PopData.csv", header=TRUE)

states <- unique(PercentPos1$REGION)
S = length(states)

# Note: we will call the 2014-2015 season the 2014 season, etc. 
# That is, the n season lasts from week 40 in year n to week 39 in year n+1

# start and end years
StartYear = 2010
EndYear = 2019
N = EndYear-StartYear+1
# data frame to hold peak weeks, end weeks, slopes, etc.
data <- data.frame(Year = rep(0, N*S),
                   State = rep(states, N),
                   PeakWeek = rep(0, N*S),
                   EndWeek = rep(0, N*S),
                   UpSlope = rep(0, N*S),
                   DownSlope = rep(0, N*S),
                   FlatSlope = rep(0, N*S),
                   StartVal = rep(0, N*S),
                   PeakVal = rep(0, N*S),
                   EndVal = rep(0, N*S),
                   PopDensity = rep(0, N*S),
                   Under18 = rep(0, N*S))

names(PopData2010)[names(PopData2010) == 'District.of.Columbia..Total..Estimate'] <- 'District of Columbia..Total..Estimate'
names(PopData2010)[names(PopData2010) == 'New.Hampshire..Total..Estimate'] <- 'New Hampshire..Total..Estimate'
names(PopData2010)[names(PopData2010) == 'New.Jersey..Total..Estimate'] <- 'New Jersey..Total..Estimate'
names(PopData2010)[names(PopData2010) == 'New.Mexico..Total..Estimate'] <- 'New Mexico..Total..Estimate'
names(PopData2010)[names(PopData2010) == 'New.York..Total..Estimate'] <- 'New York..Total..Estimate'
names(PopData2010)[names(PopData2010) == 'North.Carolina..Total..Estimate'] <- 'North Carolina..Total..Estimate'
names(PopData2010)[names(PopData2010) == 'North.Dakota..Total..Estimate'] <- 'North Dakota..Total..Estimate'
names(PopData2010)[names(PopData2010) == 'Rhode.Island..Total..Estimate'] <- 'Rhode Island..Total..Estimate'
names(PopData2010)[names(PopData2010) == 'South.Carolina..Total..Estimate'] <- 'South Carolina..Total..Estimate'
names(PopData2010)[names(PopData2010) == 'South.Dakota..Total..Estimate'] <- 'South Dakota..Total..Estimate'
names(PopData2010)[names(PopData2010) == 'West.Virginia..Total..Estimate'] <- 'West Virginia..Total..Estimate'
names(PopData2010)[names(PopData2010) == 'Puerto.Rico..Total..Estimate'] <- 'Puerto Rico..Total..Estimate'
names(PopData2010)[names(PopData2010) == 'United.States.Virgin.Islands..Total..Estimate'] <- 'Virgin Islands..Total..Estimate'

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
      data$PeakWeek[(S)*(n-StartYear)+s]=NA
      data$EndWeek[(S)*(n-StartYear)+s]=NA
      data$UpSlope[(S)*(n-StartYear)+s]=NA
      data$DownSlope[(S)*(n-StartYear)+s]=NA
      data$FlatSlope[(S)*(n-StartYear)+s]=NA
      data$StartVal[(S)*(n-StartYear)+s]=NA
      data$PeakVal[(S)*(n-StartYear)+s]=NA
      data$EndVal[(S)*(n-StartYear)+s]=NA
    }
    
    else{
      IncidenceN <- as.numeric(ILIn$X.UNWEIGHTED.ILI)*as.numeric(PercentPosn$PERCENT.POSITIVE)
      LogIncidenceN <- log(IncidenceN+0.01)
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
      data$PeakWeek[(S)*(n-StartYear)+s] <- sum$psi[1,2]
      data$EndWeek[(S)*(n-StartYear)+s] <- sum$psi[2,2]
      data$UpSlope[(S)*(n-StartYear)+s] <- sum$coefficients[2,1]
      data$DownSlope[(S)*(n-StartYear)+s] <- sum$coefficients[2,1]+sum$coefficients[3,1]
      data$FlatSlope[(S)*(n-StartYear)+s] <- sum$coefficients[2,1]+sum$coefficients[3,1]+sum$coefficients[4,1]
      data$StartVal[(S)*(n-StartYear)+s] <- sum$coefficients[1,1]
      
      #calculate peak/end value
      data$PeakVal[(S)*(n-StartYear)+s] <- data$StartVal[(S)*(n-StartYear)+s]+data$UpSlope[(S)*(n-StartYear)+s]*data$PeakWeek[(S)*(n-StartYear)+s]
      data$EndVal[(S)*(n-StartYear)+s] <- data$PeakVal[(S)*(n-StartYear)+s]+data$DownSlope[(S)*(n-StartYear)+s]*(data$EndWeek[(S)*(n-StartYear)+s]-data$PeakWeek[(S)*(n-StartYear)+s])
      
      #plot original data
      plot(week, LogIncidenceN, main = sprintf("Flu incidence in %s over the %i-%i season", states[s], n,n+1))
      
      #add segmented regression model
      plot(segmented.fitN, add=T)
    }
  }
}