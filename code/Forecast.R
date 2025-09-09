### Using ARIMA to forecast ###

# ==============================================================================
# Import libraries
# ==============================================================================

library(tidyverse) 
library(purrr)
library(plotly)
library(forecast)
library(segmented)
library(ggplot2)
library(tidyr)
library(urca)
library(signal)
library(ggdist) 
library(viridis)
library(Metrics)
library(pairsD3)

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


# ==============================================================================
# Save fits and residuals in a tibble
# ==============================================================================
# ?????
fit_data <- tibble(Year = rep(NA, N*S),
                   State = rep(NA, N*S),
                   Fit = rep(NA, N*S),
                   Residuals = rep(NA, N*S))

for (n in StartYear:EndYear) {
  for (s in 1:S) {
    fit_data$Year[S*(n-StartYear)+s] = n
    fit_data$State[S*(n-StartYear)+s] = states[s]
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
      segmented.fitN <- selgmented(lmfitN, seg.Z = ~week, Kmax = 3, stop.if=7, type="bic")
      slopes <- slope(segmented.fitN)$week[,1]
      #slopes <- segmented.fitN$coefficients
      intercept <- segmented.fitN$coefficients[1]
      
      #plot first guess at fit
      plot(week, LogIncidenceN, main = sprintf("Flu incidence in %s over the %i-%i season (?)", states[s], n,n+1))
      plot(segmented.fitN, add=T)
      
      # if there are less than 3 breakpoints, leave the fit
      if (dim(segmented.fitN$psi)[1]>=3) {
        if (slopes[2]<0 | slopes[3]>0) { # if 2nd slope is neg or 3rd slope is pos
          # go back to (at most) two breakpoints
          segmented.fitN <- selgmented(lmfitN, seg.Z = ~week, Kmax = 2, stop.if=7, type="bic")
        } else {
          break1val <- intercept+slopes[1]*segmented.fitN$psi[1,2]
          peakval <- break1val+slopes[2]*(segmented.fitN$psi[2,2]-segmented.fitN$psi[1,2])
          break3val <- peakval+slopes[3]*(segmented.fitN$psi[3,2]-segmented.fitN$psi[2,2])
          if ((peakval-break1val)/(peakval-break3val)>2 | (peakval-break1val)/(peakval-break3val)<0.5) {
            # go back to (at most) two breakpoints
            segmented.fitN <- selgmented(lmfitN, seg.Z = ~week, Kmax = 2, stop.if=7, type="bic")
          }
        }
      }
      
      fit_data$Fit[S*(n-StartYear)+s] <- list(segmented.fitN)
      fit_data$Residuals[S*(n-StartYear)+s] <- list(LogIncidenceN-predict(fit_data$Fit[S*(n-StartYear)+s][[1]]))
    }  
  }
}

# ==============================================================================
# Let's try for CO 2014 (using best fit triangle, not predicted triangle)
# ==============================================================================
n = 2014
s = 6
ILIn <- ILI[((ILI$YEAR==n & ILI$WEEK>=40) | (ILI$YEAR==n+1 & ILI$WEEK<40)) & ILI$REGION==states[s],]
if (n<2015) {
  PercentPosn <- PercentPos2[((PercentPos2$YEAR==n & PercentPos2$WEEK>=40) | (PercentPos2$YEAR==n+1 & PercentPos2$WEEK<40)) & PercentPos2$REGION==states[s],]
} else{
  PercentPosn <- PercentPos1[((PercentPos1$YEAR==n & PercentPos1$WEEK>=40) | (PercentPos1$YEAR==n+1 & PercentPos1$WEEK<40)) & PercentPos1$REGION==states[s],]
}

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
segmented.fitN <- selgmented(lmfitN, seg.Z = ~week, Kmax = 3, stop.if=7, type="bic")
slopes <- slope(segmented.fitN)$week[,1]
#slopes <- segmented.fitN$coefficients
intercept <- segmented.fitN$coefficients[1]

#plot first guess at fit
plot(week, LogIncidenceN, main = sprintf("Flu incidence in %s over the %i-%i season (?)", states[s], n,n+1))
plot(segmented.fitN, add=T)

# if there are less than 3 breakpoints, leave the fit
if (dim(segmented.fitN$psi)[1]>=3) {
  if (slopes[2]<0 | slopes[3]>0) { # if 2nd slope is neg or 3rd slope is pos
    # go back to (at most) two breakpoints
    segmented.fitN <- selgmented(lmfitN, seg.Z = ~week, Kmax = 2, stop.if=7, type="bic")
  } else {
    break1val <- intercept+slopes[1]*segmented.fitN$psi[1,2]
    peakval <- break1val+slopes[2]*(segmented.fitN$psi[2,2]-segmented.fitN$psi[1,2])
    break3val <- peakval+slopes[3]*(segmented.fitN$psi[3,2]-segmented.fitN$psi[2,2])
    if ((peakval-break1val)/(peakval-break3val)>2 | (peakval-break1val)/(peakval-break3val)<0.5) {
      # go back to (at most) two breakpoints
      segmented.fitN <- selgmented(lmfitN, seg.Z = ~week, Kmax = 2, stop.if=7, type="bic")
    }
  }
}

# plot final fit
plot(week, LogIncidenceN, main = sprintf("Flu incidence in %s over the %i-%i season", states[s], n,n+1))
plot(segmented.fitN, add=T)
# plot residuals
plot(week, LogIncidenceN-predict(segmented.fitN), main = sprintf("Residuals in %s over the %i-%i season", states[s], n,n+1))
# plot ACF, PACF
ggtsdisplay(LogIncidenceN-predict(segmented.fitN))

# fit ARIMA(1,0,0) to all the data
arima100 <- Arima(LogIncidenceN-predict(segmented.fitN), order=c(1,0,0),method="ML",include.mean = FALSE)
checkresiduals(arima100) # what are the residuals exactly??
autoplot(arima100) # come back to this
plot(forecast(arima100))
ggplot() + theme_classic() +
  geom_line(aes(x=week, y=LogIncidenceN-predict(segmented.fitN)-arima100$residuals),color='green4', linetype='dashed') +
  geom_line(aes(x=week, y=LogIncidenceN-predict(segmented.fitN)),color='slateblue')

# fit ARIMA(4,0,2) to all the data
arima402 <- Arima(LogIncidenceN-predict(segmented.fitN), order=c(4,0,2),method="ML",include.mean = FALSE)
checkresiduals(arima402)
autoplot(arima402)
plot(forecast(arima402))
ggplot() + theme_classic() +
  geom_line(aes(x=week, y=LogIncidenceN-predict(segmented.fitN)-arima402$residuals),color='green4', linetype='dashed') +
  geom_line(aes(x=week, y=LogIncidenceN-predict(segmented.fitN)),color='slateblue')


# fit ARIMA(2,0,1) to all the data
arima201 <- Arima(LogIncidenceN-predict(segmented.fitN), order=c(2,0,1),method="ML",include.mean = FALSE)
checkresiduals(arima201)
autoplot(arima201)
plot(forecast(arima201))
ggplot() + theme_classic() +
  geom_line(aes(x=week, y=LogIncidenceN-predict(segmented.fitN)-arima201$residuals),color='green4', linetype='dashed') +
  geom_line(aes(x=week, y=LogIncidenceN-predict(segmented.fitN)),color='slateblue')


# fit ARIMA(1,0,0) to previous data at each time step
forecast_arima100 <- data.frame(week = rep(NA, length(week)-2),
                                mean = rep(NA, length(week)-2),
                                low95 = rep(NA, length(week)-2),
                                high95 = rep(NA, length(week)-2))
for (i in 3:length(week)) {
  forecast_arima100$week[i-2] <- i
  arima_temp <- Arima((LogIncidenceN-predict(segmented.fitN))[1:(i-1)], order=c(1,0,0),method="ML",include.mean = FALSE)
  forecast_arima100$mean[i-2] <- forecast(arima_temp)$mean[1]
  forecast_arima100$low95[i-2] <- forecast(arima_temp)$lower[11]
  forecast_arima100$high95[i-2] <- forecast(arima_temp)$upper[11]
}

ggplot() + theme_classic() +
  geom_line(aes(x=week, y=c(0,0,forecast_arima100$mean)),color='green4',linetype='dashed') +
  geom_ribbon(aes(x=week, ymin=c(0,0,forecast_arima100$low95), ymax=c(0,0,forecast_arima100$high95)),alpha=0.2, fill='green4') +
  geom_line(aes(x=week, y=LogIncidenceN-predict(segmented.fitN)),color='slateblue')


# fit ARIMA(4,0,2) to previous data at each time step
forecast_arima402 <- data.frame(week = rep(NA, length(week)-2),
                                mean = rep(NA, length(week)-2),
                                low95 = rep(NA, length(week)-2),
                                high95 = rep(NA, length(week)-2))
for (i in 3:length(week)) {
  forecast_arima402$week[i-2] <- i
  arima_temp <- Arima((LogIncidenceN-predict(segmented.fitN))[1:(i-1)], order=c(4,0,2),method="ML",include.mean = FALSE)
  forecast_arima402$mean[i-2] <- forecast(arima_temp)$mean[1]
  forecast_arima402$low95[i-2] <- forecast(arima_temp)$lower[11]
  forecast_arima402$high95[i-2] <- forecast(arima_temp)$upper[11]
}

ggplot() + theme_classic() +
  geom_line(aes(x=week, y=c(0,0,forecast_arima402$mean)),color='green4',linetype='dashed') +
  geom_ribbon(aes(x=week, ymin=c(0,0,forecast_arima402$low95), ymax=c(0,0,forecast_arima402$high95)),alpha=0.2,fill='green4') +
  geom_line(aes(x=week, y=LogIncidenceN-predict(segmented.fitN)),color='slateblue')


# fit ARIMA(2,0,1) to previous data at each time step
forecast_arima201 <- data.frame(week = rep(NA, length(week)-2),
                                mean = rep(NA, length(week)-2),
                                low95 = rep(NA, length(week)-2),
                                high95 = rep(NA, length(week)-2))
for (i in 3:length(week)) {
  forecast_arima201$week[i-2] <- i
  arima_temp <- Arima((LogIncidenceN-predict(segmented.fitN))[1:(i-1)], order=c(2,0,1),method="ML",include.mean = FALSE)
  forecast_arima201$mean[i-2] <- forecast(arima_temp)$mean[1]
  forecast_arima201$low95[i-2] <- forecast(arima_temp)$lower[11]
  forecast_arima201$high95[i-2] <- forecast(arima_temp)$upper[11]
}

ggplot() + theme_classic() +
  geom_line(aes(x=week, y=c(0,0,forecast_arima201$mean)),color='green4',linetype='dashed') +
  geom_ribbon(aes(x=week, ymin=c(0,0,forecast_arima201$low95), ymax=c(0,0,forecast_arima201$high95)),alpha=0.2,fill='green4') +
  geom_line(aes(x=week, y=LogIncidenceN-predict(segmented.fitN)),color='slateblue')


# ==============================================================================
# Compare coefficients across years and states
# ==============================================================================

arima_coeffs <- data.frame(Year = rep(NA, N*S),
                           State = rep(NA, N*S),
                           arima100_ar1 = rep(NA, N*S),
                           arima201_ar1 = rep(NA, N*S),
                           arima201_ar2 = rep(NA, N*S),
                           arima201_ma1 = rep(NA, N*S))

for (n in StartYear:EndYear) {
  for (s in 1:S) {
    arima_coeffs$Year[S*(n-StartYear)+s] = n
    arima_coeffs$State[S*(n-StartYear)+s] = states[s]
    
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
      segmented.fitN <- selgmented(lmfitN, seg.Z = ~week, Kmax = 3, stop.if=7, type="bic")
      slopes <- slope(segmented.fitN)$week[,1]
      #slopes <- segmented.fitN$coefficients
      intercept <- segmented.fitN$coefficients[1]
      
      #plot first guess at fit
      plot(week, LogIncidenceN, main = sprintf("Flu incidence in %s over the %i-%i season (?)", states[s], n,n+1))
      plot(segmented.fitN, add=T)
      
      # if there are less than 3 breakpoints, leave the fit
      if (dim(segmented.fitN$psi)[1]>=3) {
        if (slopes[2]<0 | slopes[3]>0) { # if 2nd slope is neg or 3rd slope is pos
          # go back to (at most) two breakpoints
          segmented.fitN <- selgmented(lmfitN, seg.Z = ~week, Kmax = 2, stop.if=7, type="bic")
        } else {
          break1val <- intercept+slopes[1]*segmented.fitN$psi[1,2]
          peakval <- break1val+slopes[2]*(segmented.fitN$psi[2,2]-segmented.fitN$psi[1,2])
          break3val <- peakval+slopes[3]*(segmented.fitN$psi[3,2]-segmented.fitN$psi[2,2])
          if ((peakval-break1val)/(peakval-break3val)>2 | (peakval-break1val)/(peakval-break3val)<0.5) {
            # go back to (at most) two breakpoints
            segmented.fitN <- selgmented(lmfitN, seg.Z = ~week, Kmax = 2, stop.if=7, type="bic")
          }
        }
      }
      arima100 <- Arima((LogIncidenceN-predict(segmented.fitN)), order=c(1,0,0),method="ML",include.mean = FALSE)
      arima201 <- Arima((LogIncidenceN-predict(segmented.fitN)), order=c(2,0,1),method="ML",include.mean = FALSE)
      arima_coeffs$arima100_ar1[S*(n-StartYear)+s] = arima100$coef['ar1']
      arima_coeffs$arima201_ar1[S*(n-StartYear)+s] = arima201$coef['ar1']
      arima_coeffs$arima201_ar2[S*(n-StartYear)+s] = arima201$coef['ar2']
      arima_coeffs$arima201_ma1[S*(n-StartYear)+s] = arima201$coef['ma1']
    }
  }
}

# histograms of coefficients
arima_coeffs_long <- pivot_longer(data = arima_coeffs, cols = c('arima100_ar1','arima201_ar1','arima201_ar2','arima201_ma1'))
ggplot(data = arima_coeffs_long, aes(value, name)) + theme_classic() +
  ggdist::stat_halfeye(slab_fill='turquoise3', slab_alpha=0.5) #+
  #geom_jitter(size = 1, alpha = .15, height = 0.06)
ggplot(data = arima_coeffs_long, aes(value, name)) + theme_classic() +
  ggdist::stat_halfeye(slab_fill='snow4', slab_alpha=0.5) +
  geom_jitter(aes(color = State), size = 1, alpha = .4, height = 0.1)
ggplot(data = arima_coeffs_long[arima_coeffs_long$State=='Colorado',], aes(value, name)) + theme_classic() +
  ggdist::stat_halfeye(slab_fill='turquoise3', slab_alpha=0.5) +
  geom_jitter(size = 2, alpha = .4, height = 0.1)


# scatterplots of coefficients
pairs(arima_coeffs[3:6])
ggplot(data = arima_coeffs) + theme_classic() + geom_point(aes(x=arima201_ar1, y=arima201_ar2, color=arima201_ma1))+scale_color_viridis()
ggplot(data = arima_coeffs) + theme_classic() + geom_point(aes(x=arima201_ar1, y=arima201_ma1, color=arima201_ma1))+scale_color_viridis()
ggplot(data = arima_coeffs) + theme_classic() + geom_point(aes(x=arima201_ar2, y=arima201_ma1, color=arima201_ma1))+scale_color_viridis()
ggplot(data = arima_coeffs) + theme_classic() + geom_point(aes(x=arima100_ar1, y=arima201_ar1, color=arima201_ma1))+scale_color_viridis()
ggplot(data = arima_coeffs) + theme_classic() + geom_point(aes(x=arima100_ar1, y=arima201_ar2, color=arima201_ma1))+scale_color_viridis()
ggplot(data = arima_coeffs) + theme_classic() + geom_point(aes(x=arima100_ar1, y=arima201_ma1, color=arima201_ma1))+scale_color_viridis()


# ==============================================================================
# Calculate ARIMA model based on 2010-2017 and test it on 2018-2019 for CO/CA
# ==============================================================================

# CO first
s=6
CO_residuals <- NULL
for (n in StartYear:(EndYear-2)) {
  print(n)
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
    segmented.fitN <- selgmented(lmfitN, seg.Z = ~week, Kmax = 3, stop.if=7, type="bic")
    slopes <- slope(segmented.fitN)$week[,1]
    #slopes <- segmented.fitN$coefficients
    intercept <- segmented.fitN$coefficients[1]
    
    #plot first guess at fit
    plot(week, LogIncidenceN, main = sprintf("Flu incidence in %s over the %i-%i season (?)", states[s], n,n+1))
    plot(segmented.fitN, add=T)
    
    # if there are less than 3 breakpoints, leave the fit
    if (dim(segmented.fitN$psi)[1]>=3) {
      if (slopes[2]<0 | slopes[3]>0) { # if 2nd slope is neg or 3rd slope is pos
        # go back to (at most) two breakpoints
        segmented.fitN <- selgmented(lmfitN, seg.Z = ~week, Kmax = 2, stop.if=7, type="bic")
      } else {
        break1val <- intercept+slopes[1]*segmented.fitN$psi[1,2]
        peakval <- break1val+slopes[2]*(segmented.fitN$psi[2,2]-segmented.fitN$psi[1,2])
        break3val <- peakval+slopes[3]*(segmented.fitN$psi[3,2]-segmented.fitN$psi[2,2])
        if ((peakval-break1val)/(peakval-break3val)>2 | (peakval-break1val)/(peakval-break3val)<0.5) {
          # go back to (at most) two breakpoints
          segmented.fitN <- selgmented(lmfitN, seg.Z = ~week, Kmax = 2, stop.if=7, type="bic")
        }
      }
    }
    CO_residuals <- c(CO_residuals,(LogIncidenceN-predict(segmented.fitN)))
  }
}

ggtsdisplay(CO_residuals)

CO_arima_auto <- auto.arima(CO_residuals, stepwise = FALSE, approximation = FALSE)
autoplot(CO_arima_auto)
checkresiduals(CO_arima_auto)

CO_arima100 <- Arima(CO_residuals, order = c(1,0,0), method = 'ML', include.mean = FALSE)
autoplot(CO_arima100)
checkresiduals(CO_arima100)

CO_arima201 <- Arima(CO_residuals, order = c(2,0,1), method = 'ML', include.mean = FALSE)
autoplot(CO_arima201)
checkresiduals(CO_arima201)

CO_arima402 <- Arima(CO_residuals, order = c(4,0,2), method = 'ML', include.mean = FALSE)
autoplot(CO_arima402)
checkresiduals(CO_arima402)

CO_arima606 <- Arima(CO_residuals, order = c(6,0,6), method = 'ML', include.mean = FALSE)
autoplot(CO_arima606)
checkresiduals(CO_arima606)

CO_arima707 <- Arima(CO_residuals, order = c(7,0,7), method = 'ML', include.mean = FALSE)
autoplot(CO_arima707)
checkresiduals(CO_arima707)

plot(forecast(CO_residuals, model=CO_arima606))

# testing on 2018-2019
CO_test_residuals <- NULL

# get actual residuals
for (n in 2018:2019) {
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
    segmented.fitN <- selgmented(lmfitN, seg.Z = ~week, Kmax = 3, stop.if=7, type="bic")
    slopes <- slope(segmented.fitN)$week[,1]
    #slopes <- segmented.fitN$coefficients
    intercept <- segmented.fitN$coefficients[1]
    
    # if there are less than 3 breakpoints, leave the fit
    if (dim(segmented.fitN$psi)[1]>=3) {
      if (slopes[2]<0 | slopes[3]>0) { # if 2nd slope is neg or 3rd slope is pos
        # go back to (at most) two breakpoints
        segmented.fitN <- selgmented(lmfitN, seg.Z = ~week, Kmax = 2, stop.if=7, type="bic")
      } else {
        break1val <- intercept+slopes[1]*segmented.fitN$psi[1,2]
        peakval <- break1val+slopes[2]*(segmented.fitN$psi[2,2]-segmented.fitN$psi[1,2])
        break3val <- peakval+slopes[3]*(segmented.fitN$psi[3,2]-segmented.fitN$psi[2,2])
        if ((peakval-break1val)/(peakval-break3val)>2 | (peakval-break1val)/(peakval-break3val)<0.5) {
          # go back to (at most) two breakpoints
          segmented.fitN <- selgmented(lmfitN, seg.Z = ~week, Kmax = 2, stop.if=7, type="bic")
        }
      }
    }
    CO_test_residuals <- c(CO_test_residuals,(LogIncidenceN-predict(segmented.fitN)))
  }
}

# predict 1 week ahead
CO_test_1wk_ahead_auto <- NULL
CO_test_1wk_ahead_100 <- NULL
CO_test_1wk_ahead_201 <- NULL
CO_test_1wk_ahead_402 <- NULL
CO_test_1wk_ahead_606 <- NULL
CO_test_1wk_ahead_707 <- NULL

#training set
for (k in 1:length(CO_residuals)) {
  CO_test_1wk_ahead_auto[k] <- forecast(CO_residuals[1:(k-1)], model=CO_arima_auto)$mean[1]
  CO_test_1wk_ahead_100[k] <- forecast(CO_residuals[1:(k-1)], model=CO_arima100)$mean[1]
  CO_test_1wk_ahead_201[k] <- forecast(CO_residuals[1:(k-1)], model=CO_arima201)$mean[1]
  CO_test_1wk_ahead_402[k] <- forecast(CO_residuals[1:(k-1)], model=CO_arima402)$mean[1]
  CO_test_1wk_ahead_606[k] <- forecast(CO_residuals[1:(k-1)], model=CO_arima606)$mean[1]
  CO_test_1wk_ahead_707[k] <- forecast(CO_residuals[1:(k-1)], model=CO_arima707)$mean[1]
}
#testing set
for (k in 1:length(CO_test_residuals)) {
  CO_test_1wk_ahead_auto[k+length(CO_residuals)] <- forecast(c(CO_residuals,CO_test_residuals[1:(k-1)]), model=CO_arima_auto)$mean[1]
  CO_test_1wk_ahead_100[k+length(CO_residuals)] <- forecast(c(CO_residuals,CO_test_residuals[1:(k-1)]), model=CO_arima100)$mean[1]
  CO_test_1wk_ahead_201[k+length(CO_residuals)] <- forecast(c(CO_residuals,CO_test_residuals[1:(k-1)]), model=CO_arima201)$mean[1]
  CO_test_1wk_ahead_402[k+length(CO_residuals)] <- forecast(c(CO_residuals,CO_test_residuals[1:(k-1)]), model=CO_arima402)$mean[1]
  CO_test_1wk_ahead_606[k+length(CO_residuals)] <- forecast(c(CO_residuals,CO_test_residuals[1:(k-1)]), model=CO_arima606)$mean[1]
  CO_test_1wk_ahead_707[k+length(CO_residuals)] <- forecast(c(CO_residuals,CO_test_residuals[1:(k-1)]), model=CO_arima707)$mean[1]
}


ggplot()+theme_classic()+
  geom_point(aes(x=1:(length(CO_residuals)+length(CO_test_residuals)),y=CO_test_1wk_ahead_auto), color = '#1E3231')+
  geom_point(aes(x=1:(length(CO_residuals)+length(CO_test_residuals)),y=CO_test_1wk_ahead_100), color = '#485665')+ #88D18A
  geom_point(aes(x=1:(length(CO_residuals)+length(CO_test_residuals)),y=CO_test_1wk_ahead_201), color = '#8E7C93')+
  geom_point(aes(x=1:(length(CO_residuals)+length(CO_test_residuals)),y=CO_test_1wk_ahead_402), color = '#AF91AA')+
  geom_point(aes(x=1:(length(CO_residuals)+length(CO_test_residuals)),y=CO_test_1wk_ahead_606), color = '#D0A5C0')+
  geom_point(aes(x=1:(length(CO_residuals)+length(CO_test_residuals)),y=CO_test_1wk_ahead_707), color = '#F6C0D0')+
  geom_line(aes(x=1:(length(CO_residuals)+length(CO_test_residuals)),y=c(CO_residuals,CO_test_residuals)))+
  geom_vline(xintercept=length(CO_residuals), linetype='dashed', col = 'black')

# calculate RMSEs for testing set
rmse(CO_test_residuals,CO_test_1wk_ahead_auto[(length(CO_test_1wk_ahead_auto)-length(CO_test_residuals)+1):length(CO_test_1wk_ahead_auto)])
rmse(CO_test_residuals,CO_test_1wk_ahead_100[(length(CO_test_1wk_ahead_100)-length(CO_test_residuals)+1):length(CO_test_1wk_ahead_auto)])
rmse(CO_test_residuals,CO_test_1wk_ahead_201[(length(CO_test_1wk_ahead_201)-length(CO_test_residuals)+1):length(CO_test_1wk_ahead_auto)])
rmse(CO_test_residuals,CO_test_1wk_ahead_402[(length(CO_test_1wk_ahead_402)-length(CO_test_residuals)+1):length(CO_test_1wk_ahead_auto)])
rmse(CO_test_residuals,CO_test_1wk_ahead_606[(length(CO_test_1wk_ahead_606)-length(CO_test_residuals)+1):length(CO_test_1wk_ahead_auto)])
rmse(CO_test_residuals,CO_test_1wk_ahead_707[(length(CO_test_1wk_ahead_707)-length(CO_test_residuals)+1):length(CO_test_1wk_ahead_auto)])



# CA 
s=5
CA_residuals <- NULL
for (n in StartYear:(EndYear-2)) {
  print(n)
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
    segmented.fitN <- selgmented(lmfitN, seg.Z = ~week, Kmax = 3, stop.if=7, type="bic")
    slopes <- slope(segmented.fitN)$week[,1]
    #slopes <- segmented.fitN$coefficients
    intercept <- segmented.fitN$coefficients[1]
    
    #plot first guess at fit
    plot(week, LogIncidenceN, main = sprintf("Flu incidence in %s over the %i-%i season (?)", states[s], n,n+1))
    plot(segmented.fitN, add=T)
    
    # if there are less than 3 breakpoints, leave the fit
    if (dim(segmented.fitN$psi)[1]>=3) {
      if (slopes[2]<0 | slopes[3]>0) { # if 2nd slope is neg or 3rd slope is pos
        # go back to (at most) two breakpoints
        segmented.fitN <- selgmented(lmfitN, seg.Z = ~week, Kmax = 2, stop.if=7, type="bic")
      } else {
        break1val <- intercept+slopes[1]*segmented.fitN$psi[1,2]
        peakval <- break1val+slopes[2]*(segmented.fitN$psi[2,2]-segmented.fitN$psi[1,2])
        break3val <- peakval+slopes[3]*(segmented.fitN$psi[3,2]-segmented.fitN$psi[2,2])
        if ((peakval-break1val)/(peakval-break3val)>2 | (peakval-break1val)/(peakval-break3val)<0.5) {
          # go back to (at most) two breakpoints
          segmented.fitN <- selgmented(lmfitN, seg.Z = ~week, Kmax = 2, stop.if=7, type="bic")
        }
      }
    }
    CA_residuals <- c(CA_residuals,(LogIncidenceN-predict(segmented.fitN)))
  }
}

ggtsdisplay(CA_residuals)

CA_arima_auto <- auto.arima(CA_residuals, stepwise = FALSE, approximation = FALSE)
checkresiduals(CA_arima_auto)

CA_arima001 <- Arima(CA_residuals, order = c(0,0,1), method = 'ML', include.mean = FALSE)
checkresiduals(CA_arima001)
CA_arima001$aic

CA_arima002 <- Arima(CA_residuals, order = c(0,0,2), method = 'ML', include.mean = FALSE)
checkresiduals(CA_arima002)
CA_arima002$aic

CA_arima003 <- Arima(CA_residuals, order = c(0,0,3), method = 'ML', include.mean = FALSE)
checkresiduals(CA_arima003)
CA_arima003$aic

CA_arima004 <- Arima(CA_residuals, order = c(0,0,4), method = 'ML', include.mean = FALSE)
checkresiduals(CA_arima004)
CA_arima004$aic

CA_arima100 <- Arima(CA_residuals, order = c(1,0,0), method = 'ML', include.mean = FALSE)
autoplot(CA_arima100)
checkresiduals(CA_arima100)

CA_arima201 <- Arima(CA_residuals, order = c(2,0,1), method = 'ML', include.mean = FALSE)
autoplot(CA_arima201)
checkresiduals(CA_arima201)

CA_arima402 <- Arima(CA_residuals, order = c(4,0,2), method = 'ML', include.mean = FALSE)
autoplot(CA_arima402)
checkresiduals(CA_arima402)

CA_arima606 <- Arima(CA_residuals, order = c(6,0,6), method = 'ML', include.mean = FALSE)
autoplot(CA_arima606)
checkresiduals(CA_arima606)

CA_arima707 <- Arima(CA_residuals, order = c(7,0,7), method = 'ML', include.mean = FALSE)
autoplot(CA_arima707)
checkresiduals(CA_arima707)

# ===================================================================================
# Calculate ARIMA models based on 2010-2017 and test them on 2018-2019 for each state
# ===================================================================================

# we want to keep track of: rmse, aic?, for each state/model, plots(for each state)
# how many models should we include? all of them? (5x5)

arima_testing <- data.frame(Season = rep(NA, 25*S*2),
                            State = rep(NA, 25*S*2),
                            ModelAR = rep(NA, 25*S*2),
                            ModelMA = rep(NA, 25*S*2),
                            AIC = rep(NA, 25*S*2),
                            TestingRMSE = rep(NA, 25*S*2))

auto_arimas_by_state <- data.frame(Season = rep(NA, S*2),
                                   State = rep(NA, S*2),
                                   AutoARIMA_AR = rep(NA, S*2),
                                   AutoARIMA_MA = rep(NA, S*2),
                                   LowestAIC_AR = rep(NA, S*2),
                                   LowestAIC_MA = rep(NA, S*2),
                                   LowestRMSE_AR = rep(NA, S*2),
                                   LowestRMSE_MA = rep(NA, S*2))

# first do the full year (not just flu season)
for (s in 1:S) {
  state_residuals <- unlist(fit_data$Residuals[fit_data$State==states[s] & is.na(fit_data$Residuals)==FALSE & fit_data$Year<2018])
  if (!is.null(state_residuals)) {
    auto_arima <- auto.arima(state_residuals, stepwise = FALSE, approximation = FALSE) 
    auto_arimas_by_state$AutoARIMA_AR[s] <- auto_arima$arma[1]
    auto_arimas_by_state$AutoARIMA_MA[s] <- auto_arima$arma[2]
  }
  for (k in 0:4) {
    for (j in 0:4) {
      arima_testing$Season[25*(s-1)+5*k+j+1] = FALSE
      arima_testing$State[25*(s-1)+5*k+j+1] = states[s]
      arima_testing$ModelAR[25*(s-1)+5*k+j+1] = k
      arima_testing$ModelMA[25*(s-1)+5*k+j+1] = j
      if (!is.null(state_residuals)) {
        arima <- Arima(state_residuals, order = c(k,0,j), method = 'ML', include.mean = FALSE)
        arima_testing$AIC[25*(s-1)+5*k+j+1] <- arima$aic
        
        # rmse on testing data
        testing_residuals <- unlist(fit_data$Residuals[fit_data$State==states[s] & is.na(fit_data$Residuals)==FALSE & fit_data$Year>=2018])
        predicted_testing_residuals <- NULL
        for (i in 1:length(testing_residuals)) {
          predicted_testing_residuals[i] <- forecast(c(state_residuals, testing_residuals[1:i]), model=arima)$mean[1]
        }
        arima_testing$TestingRMSE[25*(s-1)+5*k+j+1] <- rmse(actual = testing_residuals, predicted = predicted_testing_residuals)
      }
    }
  }
}

for (s in 1:S) {
  auto_arimas_by_state$Season[s] <- FALSE
  auto_arimas_by_state$State[s] <- states[s]
  if (!is.na(arima_testing$AIC[arima_testing$Season==FALSE&arima_testing$State==states[s]][1])) {
    minAICindex <- which.min(arima_testing$AIC[arima_testing$Season==FALSE&arima_testing$State==states[s]])
    auto_arimas_by_state$LowestAIC_AR[s] <- arima_testing$ModelAR[arima_testing$Season==FALSE&arima_testing$State==states[s]][minAICindex]
    auto_arimas_by_state$LowestAIC_MA[s] <- arima_testing$ModelMA[arima_testing$Season==FALSE&arima_testing$State==states[s]][minAICindex]
  }
  if (!is.na(arima_testing$TestingRMSE[arima_testing$Season==FALSE&arima_testing$State==states[s]][1])) {
    minRMSEindex <- which.min(arima_testing$TestingRMSE[arima_testing$Season==FALSE&arima_testing$State==states[s]])
    auto_arimas_by_state$LowestRMSE_AR[s] <- arima_testing$ModelAR[arima_testing$Season==FALSE&arima_testing$State==states[s]][minRMSEindex]
    auto_arimas_by_state$LowestRMSE_MA[s] <- arima_testing$ModelMA[arima_testing$Season==FALSE&arima_testing$State==states[s]][minRMSEindex]
  }
}

# then do just flu seasons glued together

# ????
for (s in 1:S) {
  state_residuals <- NULL
  #calculate weeks of flu season (start and end of triangle)
  for (n in StartYear:2017) {
    if(!is.na(flu_data$PeakWeek[flu_data$State==states[s]&flu_data$Year==n])) {
      peak_week <- flu_data$PeakWeek[flu_data$State==states[s]&flu_data$Year==n]
      if (peak_week==flu_data$Break1[flu_data$State==states[s]&flu_data$Year==n]) { 
        season_start = 1
        if (is.na(flu_data$Break2[flu_data$State==states[s]&flu_data$Year==n])) { season_end = 52 }
        else { season_end = floor(flu_data$Break2[flu_data$State==states[s]&flu_data$Year==n]) }
      }
      else if (peak_week==flu_data$Break2[flu_data$State==states[s]&flu_data$Year==n]) { 
        season_start = ceiling(flu_data$Break1[flu_data$State==states[s]&flu_data$Year==n])
        if (is.na(flu_data$Break3[flu_data$State==states[s]&flu_data$Year==n])) { season_end = 52 }
        else { season_end = floor(flu_data$Break3[flu_data$State==states[s]&flu_data$Year==n]) }
      }
      state_residuals_by_year <- unlist(fit_data$Residuals[fit_data$State==states[s] & is.na(fit_data$Residuals)==FALSE & fit_data$Year==n])[season_start:season_end]
      state_residuals <- c(state_residuals, state_residuals_by_year)
      ###
      print(n)
      print(season_start)
      print(season_end)
      ###
    }
  }
  
  if (!is.null(state_residuals)) {
    auto_arima <- auto.arima(state_residuals, stepwise = FALSE, approximation = FALSE) 
    auto_arimas_by_state$AutoARIMA_AR[S+s] <- auto_arima$arma[1]
    auto_arimas_by_state$AutoARIMA_MA[S+s] <- auto_arima$arma[2]
  }
  for (k in 0:4) {
    for (j in 0:4) {
      arima_testing$Season[1300+25*(s-1)+5*k+j+1] = TRUE
      arima_testing$State[1300+25*(s-1)+5*k+j+1] = states[s]
      arima_testing$ModelAR[1300+25*(s-1)+5*k+j+1] = k
      arima_testing$ModelMA[1300+25*(s-1)+5*k+j+1] = j
      if (!is.null(state_residuals)) {
        arima <- Arima(state_residuals, order = c(k,0,j), method = 'ML', include.mean = FALSE)
        arima_testing$AIC[1300+25*(s-1)+5*k+j+1] <- arima$aic
        
        # rmse on testing data
        #???
        testing_state_residuals <- NULL
        for (n in 2018:EndYear) {
          if (!is.na(flu_data$PeakWeek[flu_data$State==states[s]&flu_data$Year==n])) {
            peak_week <- flu_data$PeakWeek[flu_data$State==states[s]&flu_data$Year==n]
            if (peak_week==flu_data$Break1[flu_data$State==states[s]&flu_data$Year==n]) { 
              season_start = 1
              if (is.na(flu_data$Break2[flu_data$State==states[s]&flu_data$Year==n])) { season_end = 52 } #length(???)
              else { season_end = floor(flu_data$Break2[flu_data$State==states[s]&flu_data$Year==n]) }
            }
            else if (peak_week==flu_data$Break2[flu_data$State==states[s]&flu_data$Year==n]) { 
              season_start = ceiling(flu_data$Break1[flu_data$State==states[s]&flu_data$Year==n])
              if (is.na(flu_data$Break3[flu_data$State==states[s]&flu_data$Year==n])) { season_end = 52 }
              else { season_end = floor(flu_data$Break3[flu_data$State==states[s]&flu_data$Year==n]) }
            }
            testing_state_residuals_by_year <- unlist(fit_data$Residuals[fit_data$State==states[s] & is.na(fit_data$Residuals)==FALSE & fit_data$Year==n])[season_start:season_end]
            testing_state_residuals <- c(testing_state_residuals, testing_state_residuals_by_year)
            ###
            print(n)
            print(season_start)
            print(season_end)
            ###
          }
        }
        
        predicted_testing_residuals <- NULL
        for (i in 1:length(testing_state_residuals)) {
          predicted_testing_residuals[i] <- forecast(c(state_residuals, testing_state_residuals[1:i]), model=arima)$mean[1]
        }
        arima_testing$TestingRMSE[1300+25*(s-1)+5*k+j+1] <- rmse(actual = testing_state_residuals, predicted = predicted_testing_residuals)
      }
    }
  }
}

for (s in 1:S) {
  auto_arimas_by_state$Season[S+s] <- TRUE
  auto_arimas_by_state$State[S+s] <- states[s]
  if (!is.na(arima_testing$AIC[arima_testing$Season==TRUE&arima_testing$State==states[s]][1])) {
    minAICindex <- which.min(arima_testing$AIC[arima_testing$Season==TRUE&arima_testing$State==states[s]])
    auto_arimas_by_state$LowestAIC_AR[S+s] <- arima_testing$ModelAR[arima_testing$Season==TRUE&arima_testing$State==states[s]][minAICindex]
    auto_arimas_by_state$LowestAIC_MA[S+s] <- arima_testing$ModelMA[arima_testing$Season==TRUE&arima_testing$State==states[s]][minAICindex]
  }
  if (!is.na(arima_testing$TestingRMSE[arima_testing$Season==TRUE&arima_testing$State==states[s]][1])) {
    minRMSEindex <- which.min(arima_testing$TestingRMSE[arima_testing$Season==TRUE&arima_testing$State==states[s]])
    auto_arimas_by_state$LowestRMSE_AR[S+s] <- arima_testing$ModelAR[arima_testing$Season==TRUE&arima_testing$State==states[s]][minRMSEindex]
    auto_arimas_by_state$LowestRMSE_MA[S+s] <- arima_testing$ModelMA[arima_testing$Season==TRUE&arima_testing$State==states[s]][minRMSEindex]
  }
}

# compare min RMSEs between whole year and just flu season
whole_year_rmses <- NULL
flu_season_rmses <- NULL
for (s in 1:S) {
  whole_year_rmses <- c(whole_year_rmses, min(arima_testing$TestingRMSE[arima_testing$Season==FALSE&arima_testing$State==states[s]]))
  flu_season_rmses <- c(flu_season_rmses, min(arima_testing$TestingRMSE[arima_testing$Season==TRUE&arima_testing$State==states[s]]))
}

mean(whole_year_rmses, na.rm = TRUE)
mean(flu_season_rmses, na.rm = TRUE)
plot(whole_year_rmses, flu_season_rmses)
abline(a=0,b=1)

# count number of (0,0,0), etc models chosen
dim(auto_arimas_by_state[(auto_arimas_by_state$AutoARIMA_AR==0 & auto_arimas_by_state$AutoARIMA_MA==0 & !is.na(auto_arimas_by_state$AutoARIMA_AR) & auto_arimas_by_state$Season==FALSE),])[1]

# compare AIC and RMSE
plot(arima_testing$AIC, arima_testing$TestingRMSE)
shinypairs(arima_testing)

# ===========================================================================================
# Calculate ARIMA model based on 2010-2017 (flu season only) and test it on 2018-2019 for CO
# ===========================================================================================

flu_data <- read.csv("flu_data(3).csv")
s=6
CO_flu_season_residuals <- NULL
for (n in StartYear:(EndYear-2)) {
  print(n)
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
    segmented.fitN <- selgmented(lmfitN, seg.Z = ~week, Kmax = 3, stop.if=7, type="bic")
    slopes <- slope(segmented.fitN)$week[,1]
    #slopes <- segmented.fitN$coefficients
    intercept <- segmented.fitN$coefficients[1]
    
    #plot first guess at fit
    #plot(week, LogIncidenceN, main = sprintf("Flu incidence in %s over the %i-%i season (?)", states[s], n,n+1))
    #plot(segmented.fitN, add=T)
    
    # if there are less than 3 breakpoints, leave the fit
    if (dim(segmented.fitN$psi)[1]>=3) {
      if (slopes[2]<0 | slopes[3]>0) { # if 2nd slope is neg or 3rd slope is pos
        # go back to (at most) two breakpoints
        segmented.fitN <- selgmented(lmfitN, seg.Z = ~week, Kmax = 2, stop.if=7, type="bic")
      } else {
        break1val <- intercept+slopes[1]*segmented.fitN$psi[1,2]
        peakval <- break1val+slopes[2]*(segmented.fitN$psi[2,2]-segmented.fitN$psi[1,2])
        break3val <- peakval+slopes[3]*(segmented.fitN$psi[3,2]-segmented.fitN$psi[2,2])
        if ((peakval-break1val)/(peakval-break3val)>2 | (peakval-break1val)/(peakval-break3val)<0.5) {
          # go back to (at most) two breakpoints
          segmented.fitN <- selgmented(lmfitN, seg.Z = ~week, Kmax = 2, stop.if=7, type="bic")
        }
      }
    }
    #plot final fit
    plot(week, LogIncidenceN, main = sprintf("Flu incidence in %s over the %i-%i season", states[s], n,n+1))
    plot(segmented.fitN, add=T)
    
    #calculate weeks of flu season (start and end of triangle)
    peak_week <- flu_data$PeakWeek[flu_data$State==states[s]&flu_data$Year==n]
    if (peak_week==flu_data$Break1[flu_data$State==states[s]&flu_data$Year==n]) { 
      season_start = 1
      season_end = floor(flu_data$Break2[flu_data$State==states[s]&flu_data$Year==n])
    }
    if (peak_week==flu_data$Break2[flu_data$State==states[s]&flu_data$Year==n]) { 
      season_start = ceiling(flu_data$Break1[flu_data$State==states[s]&flu_data$Year==n])
      season_end = floor(flu_data$Break3[flu_data$State==states[s]&flu_data$Year==n])
    }
    print(n)
    print(season_start)
    print(season_end)
    CO_flu_season_residuals <- c(CO_flu_season_residuals,(LogIncidenceN-predict(segmented.fitN))[season_start:season_end])
  }
}

ggtsdisplay(CO_flu_season_residuals)

CO_season_arima_auto <- auto.arima(CO_flu_season_residuals, stepwise = FALSE, approximation = FALSE)
checkresiduals(CO_season_arima_auto)

CO_season_arima100 <- Arima(CO_flu_season_residuals, order = c(1,0,0), method = 'ML', include.mean = FALSE)
checkresiduals(CO_season_arima100)

CO_season_arima201 <- Arima(CO_flu_season_residuals, order = c(2,0,1), method = 'ML', include.mean = FALSE)
checkresiduals(CO_season_arima201)

CO_season_arima402 <- Arima(CO_flu_season_residuals, order = c(4,0,2), method = 'ML', include.mean = FALSE)
checkresiduals(CO_season_arima402)

CO_season_arima606 <- Arima(CO_flu_season_residuals, order = c(6,0,6), method = 'ML', include.mean = FALSE)
checkresiduals(CO_season_arima606)

CO_season_arima707 <- Arima(CO_flu_season_residuals, order = c(7,0,7), method = 'ML', include.mean = FALSE)
checkresiduals(CO_season_arima707)

plot(forecast(CO_residuals, model=CO_arima606))

# testing on 2018-2019
CO_flu_season_test_residuals <- NULL

# get actual residuals
for (n in 2018:2019) {
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
    segmented.fitN <- selgmented(lmfitN, seg.Z = ~week, Kmax = 3, stop.if=7, type="bic")
    slopes <- slope(segmented.fitN)$week[,1]
    #slopes <- segmented.fitN$coefficients
    intercept <- segmented.fitN$coefficients[1]
    
    # if there are less than 3 breakpoints, leave the fit
    if (dim(segmented.fitN$psi)[1]>=3) {
      if (slopes[2]<0 | slopes[3]>0) { # if 2nd slope is neg or 3rd slope is pos
        # go back to (at most) two breakpoints
        segmented.fitN <- selgmented(lmfitN, seg.Z = ~week, Kmax = 2, stop.if=7, type="bic")
      } else {
        break1val <- intercept+slopes[1]*segmented.fitN$psi[1,2]
        peakval <- break1val+slopes[2]*(segmented.fitN$psi[2,2]-segmented.fitN$psi[1,2])
        break3val <- peakval+slopes[3]*(segmented.fitN$psi[3,2]-segmented.fitN$psi[2,2])
        if ((peakval-break1val)/(peakval-break3val)>2 | (peakval-break1val)/(peakval-break3val)<0.5) {
          # go back to (at most) two breakpoints
          segmented.fitN <- selgmented(lmfitN, seg.Z = ~week, Kmax = 2, stop.if=7, type="bic")
        }
      }
    }
    #calculate weeks of flu season (start and end of triangle)
    peak_week <- flu_data$PeakWeek[flu_data$State==states[s]&flu_data$Year==n]
    if (peak_week==flu_data$Break1[flu_data$State==states[s]&flu_data$Year==n]) { 
      season_start = 1
      season_end = floor(flu_data$Break2[flu_data$State==states[s]&flu_data$Year==n])
    }
    if (peak_week==flu_data$Break2[flu_data$State==states[s]&flu_data$Year==n]) { 
      season_start = ceiling(flu_data$Break1[flu_data$State==states[s]&flu_data$Year==n])
      season_end = floor(flu_data$Break3[flu_data$State==states[s]&flu_data$Year==n])
    }
    CO_flu_season_test_residuals <- c(CO_flu_season_test_residuals,(LogIncidenceN-predict(segmented.fitN))[season_start:season_end])
  }
}

# predict 1 week ahead
CO_flu_season_test_1wk_ahead_auto <- NULL
CO_flu_season_test_1wk_ahead_100 <- NULL
CO_flu_season_test_1wk_ahead_201 <- NULL
CO_flu_season_test_1wk_ahead_402 <- NULL
CO_flu_season_test_1wk_ahead_606 <- NULL
CO_flu_season_test_1wk_ahead_707 <- NULL

#training set
for (k in 1:length(CO_flu_season_residuals)) {
  CO_flu_season_test_1wk_ahead_auto[k] <- forecast(CO_flu_season_residuals[1:(k-1)], model=CO_season_arima_auto)$mean[1]
  CO_flu_season_test_1wk_ahead_100[k] <- forecast(CO_flu_season_residuals[1:(k-1)], model=CO_season_arima100)$mean[1]
  CO_flu_season_test_1wk_ahead_201[k] <- forecast(CO_flu_season_residuals[1:(k-1)], model=CO_season_arima201)$mean[1]
  CO_flu_season_test_1wk_ahead_402[k] <- forecast(CO_flu_season_residuals[1:(k-1)], model=CO_season_arima402)$mean[1]
  CO_flu_season_test_1wk_ahead_606[k] <- forecast(CO_flu_season_residuals[1:(k-1)], model=CO_season_arima606)$mean[1]
  CO_flu_season_test_1wk_ahead_707[k] <- forecast(CO_flu_season_residuals[1:(k-1)], model=CO_season_arima707)$mean[1]
}
#testing set
for (k in 1:length(CO_flu_season_test_residuals)) {
  CO_flu_season_test_1wk_ahead_auto[k+length(CO_flu_season_residuals)] <- forecast(c(CO_flu_season_residuals,CO_flu_season_test_residuals[1:(k-1)]), model=CO_season_arima_auto)$mean[1]
  CO_flu_season_test_1wk_ahead_100[k+length(CO_flu_season_residuals)] <- forecast(c(CO_flu_season_residuals,CO_flu_season_test_residuals[1:(k-1)]), model=CO_season_arima100)$mean[1]
  CO_flu_season_test_1wk_ahead_201[k+length(CO_flu_season_residuals)] <- forecast(c(CO_flu_season_residuals,CO_flu_season_test_residuals[1:(k-1)]), model=CO_season_arima201)$mean[1]
  CO_flu_season_test_1wk_ahead_402[k+length(CO_flu_season_residuals)] <- forecast(c(CO_flu_season_residuals,CO_flu_season_test_residuals[1:(k-1)]), model=CO_season_arima402)$mean[1]
  CO_flu_season_test_1wk_ahead_606[k+length(CO_flu_season_residuals)] <- forecast(c(CO_flu_season_residuals,CO_flu_season_test_residuals[1:(k-1)]), model=CO_season_arima606)$mean[1]
  CO_flu_season_test_1wk_ahead_707[k+length(CO_flu_season_residuals)] <- forecast(c(CO_flu_season_residuals,CO_flu_season_test_residuals[1:(k-1)]), model=CO_season_arima707)$mean[1]
}


ggplot()+theme_classic()+
  geom_point(aes(x=1:(length(CO_flu_season_residuals)+length(CO_flu_season_test_residuals)),y=CO_flu_season_test_1wk_ahead_auto), color = '#1E3231')+
  geom_point(aes(x=1:(length(CO_flu_season_residuals)+length(CO_flu_season_test_residuals)),y=CO_flu_season_test_1wk_ahead_100), color = '#485665')+ #88D18A
  geom_point(aes(x=1:(length(CO_flu_season_residuals)+length(CO_flu_season_test_residuals)),y=CO_flu_season_test_1wk_ahead_201), color = '#8E7C93')+
  geom_point(aes(x=1:(length(CO_flu_season_residuals)+length(CO_flu_season_test_residuals)),y=CO_flu_season_test_1wk_ahead_402), color = '#AF91AA')+
  geom_point(aes(x=1:(length(CO_flu_season_residuals)+length(CO_flu_season_test_residuals)),y=CO_flu_season_test_1wk_ahead_606), color = '#D0A5C0')+
  geom_point(aes(x=1:(length(CO_flu_season_residuals)+length(CO_flu_season_test_residuals)),y=CO_flu_season_test_1wk_ahead_707), color = '#F6C0D0')+
  geom_line(aes(x=1:(length(CO_flu_season_residuals)+length(CO_flu_season_test_residuals)),y=c(CO_flu_season_residuals,CO_flu_season_test_residuals)))+
  geom_vline(xintercept=length(CO_flu_season_residuals), linetype='dashed', col = 'black')

# calculate RMSEs for testing set
rmse(CO_flu_season_test_residuals,CO_flu_season_test_1wk_ahead_auto[(length(CO_flu_season_test_1wk_ahead_auto)-length(CO_flu_season_test_residuals)+1):length(CO_flu_season_test_1wk_ahead_auto)])
rmse(CO_flu_season_test_residuals,CO_flu_season_test_1wk_ahead_100[(length(CO_flu_season_test_1wk_ahead_100)-length(CO_flu_season_test_residuals)+1):length(CO_flu_season_test_1wk_ahead_auto)])
rmse(CO_flu_season_test_residuals,CO_flu_season_test_1wk_ahead_201[(length(CO_flu_season_test_1wk_ahead_201)-length(CO_flu_season_test_residuals)+1):length(CO_flu_season_test_1wk_ahead_auto)])
rmse(CO_flu_season_test_residuals,CO_flu_season_test_1wk_ahead_402[(length(CO_flu_season_test_1wk_ahead_402)-length(CO_flu_season_test_residuals)+1):length(CO_flu_season_test_1wk_ahead_auto)])
rmse(CO_flu_season_test_residuals,CO_flu_season_test_1wk_ahead_606[(length(CO_flu_season_test_1wk_ahead_606)-length(CO_flu_season_test_residuals)+1):length(CO_flu_season_test_1wk_ahead_auto)])
rmse(CO_flu_season_test_residuals,CO_flu_season_test_1wk_ahead_707[(length(CO_flu_season_test_1wk_ahead_707)-length(CO_flu_season_test_residuals)+1):length(CO_flu_season_test_1wk_ahead_auto)])



