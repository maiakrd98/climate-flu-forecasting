# Fitting a piecewise linear function to flu incidence

library(segmented)

# percent visits due to ILI
ILI <- read.csv("ILINetNat.csv", header=TRUE, skip = 1) 
# percent of flu tests positive post 2015
PercentPos1 <- read.csv("PercentPosNat(2016+).csv", header=TRUE, skip = 1) 
# percent of flu tests positive pre 2015
PercentPos2 <- read.csv("PercentPosNat(pre2016).csv", header=TRUE, skip = 1) 

# Note: we will call the 2014-2015 season the 2014 season, etc. 
# That is, the n season lasts from week 40 in year n to week 39 in year n+1

# start and end years
StartYear = 2010
EndYear = 2023
N = EndYear-StartYear
# vector to hold peak weeks, end weeks, and slopes
PeakWeek <- rep(0, N) 
EndWeek <- rep(0, N)
UpSlope <- rep(0, N)
DownSlope <- rep(0, N)
FlatSlope <- rep(0, N)
StartVal <- rep(0, N)
PeakVal <- rep(0, N)
EndVal <- rep(0, N)

# fit and plot data by year
for (n in StartYear:EndYear) {
  ILIn <- ILI[(ILI$YEAR==n & ILI$WEEK>=40) | (ILI$YEAR==n+1) & ILI$WEEK<40,]
  if (n<2015) {
    PercentPosn <- PercentPos2[(PercentPos2$YEAR==n & PercentPos2$WEEK>=40) | (PercentPos2$YEAR==n+1) & PercentPos2$WEEK<40,]
  }
  else{
    PercentPosn <- PercentPos1[(PercentPos1$YEAR==n & PercentPos1$WEEK>=40) | (PercentPos1$YEAR==n+1) & PercentPos1$WEEK<40,]
  }
  IncidenceN <- ILIn$X.UNWEIGHTED.ILI*PercentPosn$PERCENT.POSITIVE
  LogIncidenceN <- log(IncidenceN)
  week <- 1:nrow(ILIn)
  
  # first fit regular linear regression
  lmfitN <- lm(LogIncidenceN~week)
  
  #fit piecewise regression model to original model, estimating a breakpoint at week 20
  segmented.fitN <- segmented(lmfitN, seg.Z = ~week, npsi=2)
  
  #view summary of segmented model
  summary(segmented.fitN)
  
  #keep track of breakpoints (peak/end week), slopes
  sum <- summary.segmented(segmented.fitN)
  PeakWeek[n-StartYear+1] <- sum$psi[1,2]
  EndWeek[n-StartYear+1] <- sum$psi[2,2]
  UpSlope[n-StartYear+1] <- sum$coefficients[2,1]
  DownSlope[n-StartYear+1] <- sum$coefficients[2,1]+sum$coefficients[3,1]
  FlatSlope[n-StartYear+1] <- sum$coefficients[2,1]+sum$coefficients[3,1]+sum$coefficients[4,1]
  StartVal[n-StartYear+1] <- sum$coefficients[1,1]
  
  #calculate peak/end value
  PeakVal[n-StartYear+1] <- StartVal[n-StartYear+1]+UpSlope[n-StartYear+1]*PeakWeek[n-StartYear+1]
  EndVal[n-StartYear+1] <- PeakVal[n-StartYear+1]+DownSlope[n-StartYear+1]*(EndWeek[n-StartYear+1]-PeakWeek[n-StartYear+1])
  
  #plot original data
  plot(week, LogIncidenceN, main = sprintf("Flu incidence over the %i-%i season", n,n+1))
  
  #add segmented regression model
  plot(segmented.fitN, add=T)
  
  #plot lines from extracted parameters to check
  est <- rep(0, rev(week)[1])
  for (w in week) {
    if (w<PeakWeek[n-StartYear+1]){
      est[w]=StartVal[n-StartYear+1]+UpSlope[n-StartYear+1]*w
    }
    if (PeakWeek[n-StartYear+1]<w && w<EndWeek[n-StartYear+1]){
      est[w]=PeakVal[n-StartYear+1]+DownSlope[n-StartYear+1]*(w-PeakWeek[n-StartYear+1])
    }
    if (w>EndWeek[n-StartYear+1]){
      est[w]=EndVal[n-StartYear+1]+FlatSlope[n-StartYear+1]*(w-EndWeek[n-StartYear+1])
    }
  }
  lines(week,est,col='blue')
  }