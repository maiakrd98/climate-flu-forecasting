### Stan for flu forecasting ###
library('rstan')
library(tidyverse) 
library(purrr)

flu_data <- read.csv("flu_data.csv")
states = unique(flu_data$State)
years = unique(flu_data$Year)
S = length(states) # number of states
N = length(years)

# first, do a model like the school examples that doesn't include covariates (separated by state)
peak_week_aves <- data.frame(state = unique(flu_data$State),
                             peakwk = rep(0,S),
                             sd = rep(0,S))
for (s in 1:S) {
  peak_week_aves$peakwk[s]=mean(flu_data[(flu_data$State==states[s]),"PeakWeek"],na.rm = T)
  peak_week_aves$sd[s]=sd(flu_data[(flu_data$State==states[s]),"PeakWeek"],na.rm = T)
}
peak_week_aves <- peak_week_aves[is.na(peak_week_aves$sd)==F,]

J <- length(peak_week_aves$state)
y<-peak_week_aves$peakwk
sigma <- peak_week_aves$sd
peakwk_fit <- stan(file="peakwk1.stan", data=c("J","y","sigma"), iter=1000, chains=4) 
print(peakwk_fit) 
plot(peakwk_fit)

# still a model like the school examples that doesn't include covariates (separated by year)
peak_week_yraves <- data.frame(year = unique(flu_data$Year),
                             peakwk = rep(0,N),
                             sd = rep(0,N))
for (n in 1:N) {
  peak_week_yraves$peakwk[n]=mean(flu_data[(flu_data$Year==years[n]),"PeakWeek"],na.rm = T)
  peak_week_yraves$sd[n]=sd(flu_data[(flu_data$Year==years[n]),"PeakWeek"],na.rm = T)
}
peak_week_yraves <- peak_week_yraves[is.na(peak_week_yraves$sd)==F,]

J <- length(peak_week_yraves$year)
y<-peak_week_yraves$peakwk
sigma <- peak_week_yraves$sd
peakwk_fit <- stan(file="peakwk1(2).stan", data=c("J","y","sigma"), iter=1000, chains=4) 
print(peakwk_fit) 
plot(peakwk_fit)

# takes in a vector x and outputs a vector whose mean and sd have been normalized
normalize <- function(x) {
  norm_x <- (x-mean(x))/sd(x)
  return(norm_x)
}

##### Peak week linear model #####
# now, try to add in the linear model with weather data...
# this is a linear model, but not hierarchical 
peak_week_data <- flu_data[c("Year","State","PeakWeek","PopDensity","Under18","Lat","MeanMaxTemp","MeanMaxRelHum","MeanMinRelHum","MeanMeanAbsHum")]
peak_week_data <- peak_week_data[is.na(peak_week_data$PeakWeek)==F,]
peak_week_data <- peak_week_data[is.na(peak_week_data$MeanMinRelHum)==F,]
N = length(peak_week_data$State)
M = 7 # number of covariates
y <- normalize(peak_week_data$PeakWeek)
x <- array(dim=c(N,M))
x[,1] <- normalize(peak_week_data$PopDensity)
x[,2] <- normalize(peak_week_data$Under18)
x[,3] <- normalize(peak_week_data$Lat)
x[,4] <- normalize(peak_week_data$MeanMaxTemp)
x[,5] <- normalize(peak_week_data$MeanMaxRelHum)
x[,6] <- normalize(peak_week_data$MeanMinRelHum)
x[,7] <- normalize(peak_week_data$MeanMeanAbsHum)
peakwk_fit_2 <- stan(file="peakwk2(2).stan", data=c("N","M","x","y"), iter=1000, chains=4)
print(peakwk_fit_2) 
plot(peakwk_fit_2)



##### Up slope linear model #####
up_slope_data <- flu_data[c("Year","State","UpSlope","PopDensity","Under18","Lat","MeanMaxTemp","MeanMaxRelHum","MeanMinRelHum","MeanMeanAbsHum")]
up_slope_data <- up_slope_data[is.na(up_slope_data$UpSlope)==F,]
up_slope_data <- up_slope_data[is.na(up_slope_data$MeanMinRelHum)==F,]
N = length(up_slope_data$State)
M = 7 # number of covariates
y <- normalize(up_slope_data$UpSlope)
x <- array(dim=c(N,M))
x[,1] <- normalize(up_slope_data$PopDensity)
x[,2] <- normalize(up_slope_data$Under18)
x[,3] <- normalize(up_slope_data$Lat)
x[,4] <- normalize(up_slope_data$MeanMaxTemp)
x[,5] <- normalize(up_slope_data$MeanMaxRelHum)
x[,6] <- normalize(up_slope_data$MeanMinRelHum)
x[,7] <- normalize(up_slope_data$MeanMeanAbsHum)
upslope_fit <- stan(file="peakwk2(2).stan", data=c("N","M","x","y"), iter=1000, chains=4)
print(upslope_fit) 
plot(upslope_fit)



##### Down slope linear model #####
down_slope_data <- flu_data[c("Year","State","DownSlope","PopDensity","Under18","Lat","MeanMaxTemp","MeanMaxRelHum","MeanMinRelHum","MeanMeanAbsHum")]
down_slope_data <- down_slope_data[is.na(down_slope_data$DownSlope)==F,]
down_slope_data <- down_slope_data[is.na(down_slope_data$MeanMinRelHum)==F,]
N = length(down_slope_data$State)
M = 7 # number of covariates
y <- normalize(down_slope_data$DownSlope)
x <- array(dim=c(N,M))
x[,1] <- normalize(down_slope_data$PopDensity)
x[,2] <- normalize(down_slope_data$Under18)
x[,3] <- normalize(down_slope_data$Lat)
x[,4] <- normalize(down_slope_data$MeanMaxTemp)
x[,5] <- normalize(down_slope_data$MeanMaxRelHum)
x[,6] <- normalize(down_slope_data$MeanMinRelHum)
x[,7] <- normalize(down_slope_data$MeanMeanAbsHum)
downslope_fit <- stan(file="peakwk2(2).stan", data=c("N","M","x","y"), iter=1000, chains=4)
print(downslope_fit) 
plot(downslope_fit)






##### Peak value linear model #####
peak_val_data <- flu_data[c("Year","State","PeakVal","PopDensity","Under18","Lat","MeanMaxTemp","MeanMaxRelHum","MeanMinRelHum","MeanMeanAbsHum")]
peak_val_data <- peak_val_data[is.na(peak_val_data$PeakVal)==F,]
peak_val_data <- peak_val_data[is.na(peak_val_data$MeanMinRelHum)==F,]
N = length(peak_val_data$State)
M = 7 # number of covariates
y <- normalize(peak_val_data$PeakVal)
x <- array(dim=c(N,M))
x[,1] <- normalize(peak_val_data$PopDensity)
x[,2] <- normalize(peak_val_data$Under18)
x[,3] <- normalize(peak_val_data$Lat)
x[,4] <- normalize(peak_val_data$MeanMaxTemp)
x[,5] <- normalize(peak_val_data$MeanMaxRelHum)
x[,6] <- normalize(peak_val_data$MeanMinRelHum)
x[,7] <- normalize(peak_val_data$MeanMeanAbsHum)
peakval_fit <- stan(file="peakwk2(2).stan", data=c("N","M","x","y"), iter=1000, chains=4)
print(peakval_fit) 
plot(peakval_fit)





####### now try to make it a hierarchical linear model #######
peak_week_data <- flu_data[c("Year","State","PeakWeek","PopDensity","Under18","Lat","MeanMaxTemp","MeanMaxRelHum","MeanMinRelHum")]
peak_week_data <- peak_week_data[is.na(peak_week_data$PeakWeek)==F,]
peak_week_data <- peak_week_data[is.na(peak_week_data$MeanMinRelHum)==F,]
N = length(peak_week_data$State)
y <- peak_week_data$PeakWeek
x1 <- peak_week_data$PopDensity
x2 <- peak_week_data$Under18
x3 <- peak_week_data$Lat
x4 <- peak_week_data$MeanMaxTemp
x5 <- peak_week_data$MeanMaxRelHum
x6 <- peak_week_data$MeanMinRelHum
peakwk_fit_2 <- stan(file="peakwk2.stan", data=c("N","x1","x2","x3","x4","x5","x6","y"), iter=1000, chains=4)
print(peakwk_fit_2) 
plot(peakwk_fit_2)