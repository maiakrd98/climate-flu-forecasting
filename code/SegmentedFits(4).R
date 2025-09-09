# Testing other ways to fit piecewise line (trapezoid)
library(segmented)
library(tidyverse)
library(ggplot2)
library(tidyr)

#LogILIplus_test <- read.csv("log_ILIplus(2).csv")

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

# create csv with log(ILI+) for each week for each year and state
LogILIplus <- data.frame(Year = rep(0, N*S),
                         State = rep(states[1:S], N))
for (n in StartYear:EndYear) {
  for (s in 1:S) {
    LogILIplus$Year[S*(n-StartYear)+s] = n
    ILIn <- ILI[((ILI$YEAR==n & ILI$WEEK>=40) | (ILI$YEAR==n+1 & ILI$WEEK<40)) & ILI$REGION==states[s],]
    if (n<2015) {
      PercentPosn <- PercentPos2[((PercentPos2$YEAR==n & PercentPos2$WEEK>=40) | (PercentPos2$YEAR==n+1 & PercentPos2$WEEK<40)) & PercentPos2$REGION==states[s],]
    }
    else{
      PercentPosn <- PercentPos1[((PercentPos1$YEAR==n & PercentPos1$WEEK>=40) | (PercentPos1$YEAR==n+1 & PercentPos1$WEEK<40)) & PercentPos1$REGION==states[s],]
    }
    if(sum(PercentPosn$PERCENT.POSITIVE=="X")>5 | sum(PercentPosn$PERCENT.POSITIVE==0)>40) {
    }
    else {
      IncidenceN <- as.numeric(ILIn$X.UNWEIGHTED.ILI)*as.numeric(PercentPosn$PERCENT.POSITIVE)
      # make zeros equal to the minimum nonzero value
      IncidenceN <- IncidenceN + min(IncidenceN[IncidenceN!=0], na.rm = TRUE)*(IncidenceN==0)
      LogIncidenceN <- log(IncidenceN)
      #LogIncidenceN <- log(IncidenceN+0.1*(IncidenceN==0))
      
      LogILIplus[S*(n-StartYear)+s,3:54] = LogIncidenceN
      #if(length(LogIncidenceN)!=52){print(n)}
    }
  }
}
# (4) since this is the third method of dealing with zeros/fitting piecewise line
# probably not necessary since we're not changing the zeros
write.csv(LogILIplus,"log_ILIplus(4).csv",row.names = FALSE)
LogILIplus <- read.csv("log_ILIplus(2).csv")


data <- read.csv("flu_data(2).csv")
data$NumBreaks <- rep(NA, N*S)
data$Break1 <- rep(NA, N*S)
data$Break2 <- rep(NA, N*S)
data$Break3 <- rep(NA, N*S)
data$Slope1 <- rep(NA, N*S)
data$Slope2 <- rep(NA, N*S)
data$Slope3 <- rep(NA, N*S)
data$Slope4 <- rep(NA, N*S)
data$Val1 <- rep(NA, N*S)
data$Val2 <- rep(NA, N*S)
data$Val3 <- rep(NA, N*S)
data$Val4 <- rep(NA, N*S)
data$PeakWeek <- rep(NA, N*S)
data$UpSlope <- rep(NA, N*S)
data$DownSlope <- rep(NA, N*S)
data$PeakVal <- rep(NA, N*S)

for (n in StartYear:EndYear) {
  for (s in 1:S) {
    # data$Year[S*(n-StartYear)+s] = n
    ILIn <- ILI[((ILI$YEAR==n & ILI$WEEK>=40) | (ILI$YEAR==n+1 & ILI$WEEK<40)) & ILI$REGION==states[s],]
    if (n<2015) {
      PercentPosn <- PercentPos2[((PercentPos2$YEAR==n & PercentPos2$WEEK>=40) | (PercentPos2$YEAR==n+1 & PercentPos2$WEEK<40)) & PercentPos2$REGION==states[s],]
    } else{
      PercentPosn <- PercentPos1[((PercentPos1$YEAR==n & PercentPos1$WEEK>=40) | (PercentPos1$YEAR==n+1 & PercentPos1$WEEK<40)) & PercentPos1$REGION==states[s],]
    }
    # skip if more than 5 missing data points?!?!? or more that 40 weeks of zeros??
    if(sum(PercentPosn$PERCENT.POSITIVE=="X")>5 | sum(PercentPosn$PERCENT.POSITIVE==0)>40) {
      # do nothing (i.e. leave NAs)
    } else {
      IncidenceN <- as.numeric(ILIn$X.UNWEIGHTED.ILI)*as.numeric(PercentPosn$PERCENT.POSITIVE)
      # make zeros equal to the minimum nonzero value
      IncidenceN <- IncidenceN + min(IncidenceN[IncidenceN!=0], na.rm=TRUE)*(IncidenceN==0)
      LogIncidenceN <- log(IncidenceN)
      #LogIncidenceN <- log(IncidenceN+0.1*(IncidenceN==0))
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
      } else if (!is.null(NAIndicies)) {
        week <- week[-NAIndicies]
        LogIncidenceN <- LogIncidenceN[-NAIndicies]
      } else if (!is.null(InfIndicies)) {
        week <- week[-InfIndicies]
        LogIncidenceN <- LogIncidenceN[-InfIndicies]
      }
      
      # first fit regular linear regression
      lmfitN <- lm(LogIncidenceN~week)
      
      # fit piecewise regression model to original model
      segmented.fitN <- segreg(LogIncidenceN ~ seg(week,npsi=3, est=c(1,0,1,1)),control=seg.control(fc=0.9))
        #selgmented(lmfitN, seg.Z = ~week, Kmax = 3, stop.if=7, type="bic")
      slopes <- slope(segmented.fitN)$week[,1]
      #slopes <- segmented.fitN$coefficients
      intercept <- segmented.fitN$coefficients[1]
      
      #plot first guess at fit
      png(paste("~/FluForecasting/FluForecasting/plots_trapezoid_test/",states[s],n,".png",sep = ""))
      plot(week, LogIncidenceN, main = sprintf("Flu incidence in %s over the %i-%i season (?)", states[s], n,n+1))
      plot(segmented.fitN, add=T)
      dev.off()
      
      # # if there are less than 3 breakpoints, leave the fit
      # if (dim(segmented.fitN$psi)[1]>=3) {
      #   if (slopes[2]<0 | slopes[3]>0) { # if 2nd slope is neg or 3rd slope is pos
      #     # go back to (at most) two breakpoints
      #     segmented.fitN <- selgmented(lmfitN, seg.Z = ~week, Kmax = 2, stop.if=7, type="bic")
      #   } else {
      #     break1val <- intercept+slopes[1]*segmented.fitN$psi[1,2]
      #     peakval <- break1val+slopes[2]*(segmented.fitN$psi[2,2]-segmented.fitN$psi[1,2])
      #     break3val <- peakval+slopes[3]*(segmented.fitN$psi[3,2]-segmented.fitN$psi[2,2])
      #     if ((peakval-break1val)/(peakval-break3val)>2 | (peakval-break1val)/(peakval-break3val)<0.5) {
      #       # go back to (at most) two breakpoints
      #       segmented.fitN <- selgmented(lmfitN, seg.Z = ~week, Kmax = 2, stop.if=7, type="bic")
      #     }
      #   }
      # }
      # 
      # #keep track of breakpoints (peak/end week), slopes
      # slopes <- slope(segmented.fitN)$week[,1]
      # sum <- summary.segmented(segmented.fitN)
      # intercept <- segmented.fitN$coefficients[1]
      # data$NumBreaks[S*(n-StartYear)+s] <- dim(segmented.fitN$psi)[1]
      # data$Break1[S*(n-StartYear)+s] <- sum$psi[1,2]
      # data$Slope1[S*(n-StartYear)+s] <- slopes[1]
      # data$Slope2[S*(n-StartYear)+s] <- slopes[2]
      # data$Val1[S*(n-StartYear)+s] <- intercept
      # data$Val2[S*(n-StartYear)+s] <- data$Val1[S*(n-StartYear)+s]+data$Slope1[S*(n-StartYear)+s]*data$Break1[S*(n-StartYear)+s]
      # if (dim(segmented.fitN$psi)[1]>1) {
      #   data$Break2[S*(n-StartYear)+s] <- sum$psi[2,2]
      #   data$Slope3[S*(n-StartYear)+s] <- slopes[3]
      #   data$Val3[S*(n-StartYear)+s] <- data$Val2[S*(n-StartYear)+s]+data$Slope2[S*(n-StartYear)+s]*(data$Break2[S*(n-StartYear)+s]-data$Break1[S*(n-StartYear)+s])
      #   if (dim(segmented.fitN$psi)[1]>2) {
      #     data$Break3[S*(n-StartYear)+s] <- sum$psi[3,2]
      #     data$Slope4[S*(n-StartYear)+s] <- slopes[4]
      #     data$Val4[(S)*(n-StartYear)+s] <- data$Val3[(S)*(n-StartYear)+s]+data$Slope3[(S)*(n-StartYear)+s]*(data$Break3[(S)*(n-StartYear)+s]-data$Break2[(S)*(n-StartYear)+s])
      #   }
      # }
      # 
      # 
      # #calculate peak/end value
      # 
      # data$PeakVal[S*(n-StartYear)+s] <- max(data$Val2[S*(n-StartYear)+s],data$Val3[S*(n-StartYear)+s],na.rm=TRUE)
      # if (data$PeakVal[S*(n-StartYear)+s]==data$Val2[S*(n-StartYear)+s]) {
      #   data$UpSlope[S*(n-StartYear)+s] <- data$Slope1[S*(n-StartYear)+s]
      #   data$DownSlope[S*(n-StartYear)+s] <- data$Slope2[S*(n-StartYear)+s]
      #   data$PeakWeek[S*(n-StartYear)+s] <- data$Break1[S*(n-StartYear)+s]
      # } else if (data$PeakVal[S*(n-StartYear)+s]==data$Val3[S*(n-StartYear)+s]) {
      #   data$UpSlope[S*(n-StartYear)+s] <- data$Slope2[S*(n-StartYear)+s]
      #   data$DownSlope[S*(n-StartYear)+s] <- data$Slope3[S*(n-StartYear)+s]
      #   data$PeakWeek[S*(n-StartYear)+s] <- data$Break2[S*(n-StartYear)+s]
      # }
      # 
      # png(paste("~/FluForecasting/FluForecasting/plots_segmented_test/",states[s],n,".png",sep = ""))
      # #plot original data
      # plot(week, LogIncidenceN, main = sprintf("Flu incidence in %s over the %i-%i season", states[s], n,n+1))
      # #add segmented regression model
      # plot(segmented.fitN, add=T)
      # dev.off()
    }  
  }
}

data1 <- data[,c(1:2,70,3:4,71,5:7,72,8:10,73,11:69)]
write.csv(data1,"flu_data(3).csv",row.names = FALSE)
