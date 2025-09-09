### Arima Model! ###
### Deciding which model to use ###

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

x <- rnorm(52)+1:52
arima(x, order=c(0,0,0))
auto.arima(x)

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

arima_test <- data.frame(Year = rep(NA, N*S),
                         State = rep(NA, N*S),
                         Ljung_Box_pval = rep(NA, N*S),
                         AR = rep(NA, N*S),
                         MA = rep(NA, N*S),
                         seasonalAR = rep(NA, N*S),
                         seasonalMA = rep(NA, N*S),
                         period = rep(NA, N*S),
                         non_seasonal_differences = rep(NA, N*S),
                         seasonal_differences = rep(NA, N*S),
                         min_aic_AR = rep(NA, N*S),
                         min_aic_MA = rep(NA, N*S),
                         '(0,0,0)AIC' = rep(NA, N*S),
                         '(1,0,0)AIC' = rep(NA, N*S),
                         '(2,0,0)AIC' = rep(NA, N*S),
                         '(3,0,0)AIC' = rep(NA, N*S),
                         '(4,0,0)AIC' = rep(NA, N*S),
                         '(5,0,0)AIC' = rep(NA, N*S),
                         '(0,0,1)AIC' = rep(NA, N*S),
                         '(1,0,1)AIC' = rep(NA, N*S),
                         '(2,0,1)AIC' = rep(NA, N*S),
                         '(3,0,1)AIC' = rep(NA, N*S),
                         '(4,0,1)AIC' = rep(NA, N*S),
                         '(5,0,1)AIC' = rep(NA, N*S),
                         '(0,0,2)AIC' = rep(NA, N*S),
                         '(1,0,2)AIC' = rep(NA, N*S),
                         '(2,0,2)AIC' = rep(NA, N*S),
                         '(3,0,2)AIC' = rep(NA, N*S),
                         '(4,0,2)AIC' = rep(NA, N*S),
                         '(5,0,2)AIC' = rep(NA, N*S),
                         '(0,0,3)AIC' = rep(NA, N*S),
                         '(1,0,3)AIC' = rep(NA, N*S),
                         '(2,0,3)AIC' = rep(NA, N*S),
                         '(3,0,3)AIC' = rep(NA, N*S),
                         '(4,0,3)AIC' = rep(NA, N*S),
                         '(5,0,3)AIC' = rep(NA, N*S),
                         '(0,0,4)AIC' = rep(NA, N*S),
                         '(1,0,4)AIC' = rep(NA, N*S),
                         '(2,0,4)AIC' = rep(NA, N*S),
                         '(3,0,4)AIC' = rep(NA, N*S),
                         '(4,0,4)AIC' = rep(NA, N*S),
                         '(5,0,4)AIC' = rep(NA, N*S),
                         '(0,0,5)AIC' = rep(NA, N*S),
                         '(1,0,5)AIC' = rep(NA, N*S),
                         '(2,0,5)AIC' = rep(NA, N*S),
                         '(3,0,5)AIC' = rep(NA, N*S),
                         '(4,0,5)AIC' = rep(NA, N*S),
                         '(5,0,5)AIC' = rep(NA, N*S))

for (n in StartYear:EndYear) {
  for (s in 1:S) {
    arima_test$Year[S*(n-StartYear)+s] = n
    arima_test$State[S*(n-StartYear)+s] = states[s]
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
      
      # conduct box test to check for autocorrelation and fit auto arima
      boxtest <- Box.test(LogIncidenceN-predict(segmented.fitN), lag=16, type="Ljung-Box")
      arima <- auto.arima(LogIncidenceN-predict(segmented.fitN))
      
      arima_test$Ljung_Box_pval[S*(n-StartYear)+s] = boxtest$p.value
      arima_test$AR[S*(n-StartYear)+s] = arima$arma[1]
      arima_test$MA[S*(n-StartYear)+s] = arima$arma[2]
      arima_test$seasonalAR[S*(n-StartYear)+s] = arima$arma[3]
      arima_test$seasonalMA[S*(n-StartYear)+s] = arima$arma[4]
      arima_test$period[S*(n-StartYear)+s] = arima$arma[5]
      arima_test$non_seasonal_differences[S*(n-StartYear)+s] = arima$arma[6]
      arima_test$seasonal_differences[S*(n-StartYear)+s] = arima$arma[7]
      arima_test$X.0.0.0.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(0,0,0),method="ML",include.mean = FALSE)$aic
      arima_test$X.1.0.0.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(1,0,0),method="ML",include.mean = FALSE)$aic
      arima_test$X.2.0.0.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(2,0,0),method="ML",include.mean = FALSE)$aic
      arima_test$X.3.0.0.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(3,0,0),method="ML",include.mean = FALSE)$aic
      arima_test$X.4.0.0.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(4,0,0),method="ML",include.mean = FALSE)$aic
      arima_test$X.5.0.0.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(5,0,0),method="ML",include.mean = FALSE)$aic
      arima_test$X.0.0.1.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(0,0,1),method="ML",include.mean = FALSE)$aic
      arima_test$X.1.0.1.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(1,0,1),method="ML",include.mean = FALSE)$aic
      arima_test$X.2.0.1.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(2,0,1),method="ML",include.mean = FALSE)$aic
      arima_test$X.3.0.1.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(3,0,1),method="ML",include.mean = FALSE)$aic
      arima_test$X.4.0.1.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(4,0,1),method="ML",include.mean = FALSE)$aic
      arima_test$X.5.0.1.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(5,0,1),method="ML",include.mean = FALSE)$aic
      arima_test$X.0.0.2.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(0,0,2),method="ML",include.mean = FALSE)$aic
      arima_test$X.1.0.2.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(1,0,2),method="ML",include.mean = FALSE)$aic
      arima_test$X.2.0.2.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(2,0,2),method="ML",include.mean = FALSE)$aic
      arima_test$X.3.0.2.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(3,0,2),method="ML",include.mean = FALSE)$aic
      arima_test$X.4.0.2.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(4,0,2),method="ML",include.mean = FALSE)$aic
      arima_test$X.5.0.2.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(5,0,2),method="ML",include.mean = FALSE)$aic
      arima_test$X.0.0.3.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(0,0,3),method="ML",include.mean = FALSE)$aic
      arima_test$X.1.0.3.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(1,0,3),method="ML",include.mean = FALSE)$aic
      arima_test$X.2.0.3.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(2,0,3),method="ML",include.mean = FALSE)$aic
      arima_test$X.3.0.3.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(3,0,3),method="ML",include.mean = FALSE)$aic
      arima_test$X.4.0.3.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(4,0,3),method="ML",include.mean = FALSE)$aic
      arima_test$X.5.0.3.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(5,0,3),method="ML",include.mean = FALSE)$aic
      arima_test$X.0.0.4.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(0,0,4),method="ML",include.mean = FALSE)$aic
      arima_test$X.1.0.4.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(1,0,4),method="ML",include.mean = FALSE)$aic
      arima_test$X.2.0.4.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(2,0,4),method="ML",include.mean = FALSE)$aic
      arima_test$X.3.0.4.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(3,0,4),method="ML",include.mean = FALSE)$aic
      arima_test$X.4.0.4.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(4,0,4),method="ML",include.mean = FALSE)$aic
      arima_test$X.5.0.4.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(5,0,4),method="ML",include.mean = FALSE)$aic
      arima_test$X.0.0.5.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(0,0,5),method="ML",include.mean = FALSE)$aic
      arima_test$X.1.0.5.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(1,0,5),method="ML",include.mean = FALSE)$aic
      arima_test$X.2.0.5.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(2,0,5),method="ML",include.mean = FALSE)$aic
      arima_test$X.3.0.5.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(3,0,5),method="ML",include.mean = FALSE)$aic
      arima_test$X.4.0.5.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(4,0,5),method="ML",include.mean = FALSE)$aic
      arima_test$X.5.0.5.AIC[S*(n-StartYear)+s] = arima(LogIncidenceN-predict(segmented.fitN),order=c(5,0,5),method="ML",include.mean = FALSE)$aic
      x = which.min(arima_test[(S*(n-StartYear)+s),13:48])
      arima_test$min_aic_AR[S*(n-StartYear)+s] = (x-1)%%6
      arima_test$min_aic_MA[S*(n-StartYear)+s] = floor((x-1)/6)
    }  
  }
}

for (j in 0:5) {
  for (k in 0:5) {
    print(c(k,j,length(arima_test$min_aic_AR[arima_test$min_aic_AR==k&arima_test$min_aic_MA==j&is.na(arima_test$min_aic_AR)==FALSE])))
  }
}
for (j in 0:5) {
  for (k in 0:5) {
    print(c(k,j,mean(arima_test[[paste0("X.",k,".0.",j,".AIC")]],na.rm = TRUE)))
  }
}

# note this is not defined in this script yet...
acf(LogIncidenceN-predict(segmented.fitN))
boxtest <- Box.test(LogIncidenceN-predict(segmented.fitN), lag=16, type="Ljung-Box")
auto.arima(LogIncidenceN)
arima <- auto.arima(LogIncidenceN-predict(segmented.fitN))

# testing what ACF actually computes
x = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
xbar = mean(x)
len = length(x)
k=2

rk_den = 0
for (i in 1:len) {
  rk_den = rk_den + (x[i]-xbar)^2
}

rk_num = 0 
for (i in (k+1):len) {
  rk_num = rk_num + (x[i]-xbar)*(x[i-k]-xbar)
}

rk=rk_num/rk_den

acf(x)

rk_den_alt_1 = 0
for (i in 1:(len-k)) {
  rk_den_alt_1 = rk_den_alt_1 + (x[i]-xbar)^2
}
rk_den_alt_2 = 0
for (i in (k+1):len) {
  rk_den_alt_2 = rk_den_alt_2 + (x[i]-xbar)^2
}
rk_alt = (rk_num/(sqrt(rk_den_alt_1*rk_den_alt_2)))*(len-k-1)/(len-1) #????


### difference equations vs AR models ###

# 1st order
y1 = 10
phi = -1
c = 14
n = 200

y = rep(NA,n)
y_AR = rep(NA,n)

y[1] = y1
y_AR[1] = y1

for (i in 2:n) {
  y[i] = c + phi*y[i-1]
  y_AR[i] = c + phi*y_AR[i-1]+rnorm(1,mean=0,sd=0.1)
}

plot(y)
plot(y_AR)

# 2nd order
y1 = 1
y2 = -10
phi1 = 0.1
phi2 = 0.7
c = 3
n = 100

y = rep(NA,n)
y_AR = rep(NA,n)

y[1] = y1
y_AR[1] = y1
y[2] = y2
y_AR[2] = y2

for (i in 3:n) {
  y[i] = c + phi1*y[i-1] + phi2*y[i-2]
  y_AR[i] = c + phi1*y_AR[i-1] + phi2*y_AR[i-2] + rnorm(1,mean=0,sd=0.1)
}

plot(y)
plot(y_AR)
