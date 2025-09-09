# Testing different models (all multivariate)
#model_loos <- read.csv("model_loos.csv")

# ==============================================================================
# Import libraries
# ==============================================================================

library(tidyverse) 
library(purrr)
library(rstan)
library(plotly)
library(reshape2)
library(mvtnorm)
library(loo)

# ==============================================================================
# Import and prepare flu data
# ==============================================================================
flu_data <- read.csv("flu_data(3).csv")
LogILIplus <- read.csv("log_ILIplus(2).csv")

states = unique(flu_data$State)
years = unique(flu_data$Year)
S = length(states) # number of states
N = length(years)

# remove NAs
flu_data_1 <- flu_data[is.na(flu_data$PeakWeek)==F,]
flu_data_1 <- flu_data_1[is.na(flu_data_1$MeanMaxTemp)==F,]

K = length(flu_data_1$Year) # number of observations

# ==============================================================================
# Normalize variables 
# ==============================================================================

# takes in a vector x and outputs a vector whose mean and sd have been normalized
normalize <- function(x) {
  norm_x <- (x-mean(x))/sd(x)
  return(norm_x)
}

# keep track of means and standard deviations of flu data
means_and_sds <- data.frame(Variable = c("PeakWeek","UpSlope", "DownSlope", "PeakVal", 
                                         "PopDensity", "Under18", "Lat", "MeanMaxTemp",
                                         "MeanMaxRelHum", "MeanMinRelHum", "MeanMeanAbsHum"),
                            Mean = rep(0),
                            SD = rep(0))
for (i in 1:11) {
  means_and_sds$Mean[i] <- mean(flu_data_1[[means_and_sds[i,1]]])
  means_and_sds$SD[i] <- sd(flu_data_1[[means_and_sds[i,1]]])
}

# normalize flu data
flu_data_norm <- flu_data_1
for (i in 3:dim(flu_data_1)[2]) {
  flu_data_norm[,i] <- normalize(flu_data_1[,i])
}

# the things we are predicting (there's definitely a real name for that...)
y <- array(dim=c(K,4))
y[,1] <- flu_data_norm$PeakWeek
y[,2] <- flu_data_norm$UpSlope
y[,3] <- flu_data_norm$DownSlope
y[,4] <- flu_data_norm$PeakVal


# ==============================================================================
# Test models
# ==============================================================================

# all covariates
x <- array(dim=c(7,K))
x[1,] <- flu_data_norm$PopDensity
x[2,] <- flu_data_norm$Under18
x[3,] <- flu_data_norm$MeanMaxTemp
x[4,] <- flu_data_norm$MeanMaxRelHum
x[5,] <- flu_data_norm$MeanMinRelHum
x[6,] <- flu_data_norm$MeanMeanAbsHum
x[7,] <- flu_data_norm$Lat

p <- expand.grid(1:0,1:0,1:0,1:0,1:0,1:0,1:0)

model_loos <- data.frame(PopDensity = as.logical(p[,1]),
                         Under18 = as.logical(p[,2]),
                         MeanMaxTemp = as.logical(p[,3]),
                         MeanMaxRelHum = as.logical(p[,4]),
                         MeanMinRelHum = as.logical(p[,5]),
                         MeanMeanAbsHum = as.logical(p[,6]),
                         Lat = as.logical(p[,7]),
                         looic = rep(0),
                         p_loo = rep(0))

for (i in 1:127) {
  M = sum(p[i,])
  x_temp <- x[as.logical(p[i,]),]
  
  if (M==1) { x_temp <- matrix(x_temp,nrow=1,ncol=K) }
  
  mvn_model_temp <- stan(file="linear_flu_mvn.stan",data=list(K=K, M=M, x=x_temp, y=y),
                         iter=3000, chains=4)
  
  log_lik_temp <- extract_log_lik(mvn_model_temp)
  # calculate information criteria
  loo_mvn_temp <- loo(log_lik_temp)
  model_loos$looic[i] <- loo_mvn_temp$estimates[3,1]
  model_loos$p_loo[i] <- loo_mvn_temp$estimates[2,1]
  
  # save csv in case of crash
  write.csv(model_loos,"model_loos_new.csv",row.names = FALSE)
}

# apparently, the best model is PopDensity, Under18, MeanMaxTemp, and Lat?????

# ==============================================================================
# Check top models
# ==============================================================================

M = 4
x_test <- array(dim=c(4,K))
x_test[1,] <- flu_data_norm$PopDensity
x_test[2,] <- flu_data_norm$Under18
x_test[3,] <- flu_data_norm$MeanMaxTemp
x_test[4,] <- flu_data_norm$Lat

mvn_model_test <- stan(file="linear_flu_mvn.stan",data=list(K=K, M=M, x=x_test, y=y),
                       iter=3000, chains=4)
log_lik_test <- extract_log_lik(mvn_model_test)
# calculate information criteria
loo_mvn_test <- loo(log_lik_test)

mvn_model_test

# peak week
peakwk_pars <- c("alpha[1]")
for (i in 1:M) {
  peakwk_pars <- c(peakwk_pars, paste0("beta[1,",i,"]"))
}

peakwk_parameter_plot <- plot(mvn_model_test, show_density = TRUE, fill_color = "lightsteelblue1", pars=peakwk_pars) + 
  ggtitle("Posterior distributions for peak week parameters") +
  theme(plot.title = element_text(face = 'bold', size =18, hjust = 0.6),
        axis.text=element_text(size=16))

# up slope
upslope_pars <- c("alpha[2]")
for (i in 1:M) {
  upslope_pars <- c(upslope_pars, paste0("beta[2,",i,"]"))
}

upslope_parameter_plot <- plot(mvn_model_test, show_density = TRUE, fill_color = "lightsteelblue1", pars=upslope_pars) + 
  ggtitle("Posterior distributions for up slope parameters") +
  theme(plot.title = element_text(face = 'bold', size =18, hjust = 0.6),
        axis.text=element_text(size=16))

# down slope
downslope_pars <- c("alpha[3]")
for (i in 1:M) {
  downslope_pars <- c(downslope_pars, paste0("beta[3,",i,"]"))
}

downslope_parameter_plot <- plot(mvn_model_test, show_density = TRUE, fill_color = "lightsteelblue1", pars=downslope_pars) + 
  ggtitle("Posterior distributions for down slope parameters") +
  theme(plot.title = element_text(face = 'bold', size =18, hjust = 0.6),
        axis.text=element_text(size=16))

# peak value
pkval_pars <- c("alpha[4]")
for (i in 1:M) {
  pkval_pars <- c(pkval_pars, paste0("beta[4,",i,"]"))
}

pkval_parameter_plot <- plot(mvn_model_test, show_density = TRUE, fill_color = "lightsteelblue1", pars=pkval_pars) + 
  ggtitle("Posterior distributions for peak value parameters") +
  theme(plot.title = element_text(face = 'bold', size =18, hjust = 0.6),
        axis.text=element_text(size=16))

# covariance matrix
cov_parameter_plot <- plot(mvn_model_test, show_density = TRUE, fill_color = "lightsteelblue1", pars="sigma") + 
  ggtitle("Posterior distributions for covariance matrix") +
  theme(plot.title = element_text(face = 'bold', size =18, hjust = 0.6),
        axis.text=element_text(size=16))

# ==============================================================================
# Check month by month temp??
# ==============================================================================

# all covariates
z <- array(dim=c(15,K))
z[1,] <- flu_data_norm$PopDensity
z[2,] <- flu_data_norm$Under18
z[3,] <- flu_data_norm$Lat
z[4,] <- flu_data_norm$OctMaxTemp
z[5,] <- flu_data_norm$NovMaxTemp
z[6,] <- flu_data_norm$DecMaxTemp
z[7,] <- flu_data_norm$JanMaxTemp
z[8,] <- flu_data_norm$FebMaxTemp
z[9,] <- flu_data_norm$MarMaxTemp
z[10,] <- flu_data_norm$AprMaxTemp
z[11,] <- flu_data_norm$MayMaxTemp
z[12,] <- flu_data_norm$JunMaxTemp
z[13,] <- flu_data_norm$JulMaxTemp
z[14,] <- flu_data_norm$AugMaxTemp
z[15,] <- flu_data_norm$SepMaxTemp

p <- expand.grid(1:0,1:0,1:0,1:0,1:0,1:0,1:0,1:0,1:0,1:0,1:0,1:0,1:0,1:0,1:0)

months_model_loos <- data.frame(PopDensity = as.logical(p[,1]),
                                Under18 = as.logical(p[,2]),
                                Lat = as.logical(p[,3]),
                                OctMaxTemp = as.logical(p[,4]),
                                NovMaxTemp = as.logical(p[,5]),
                                DecMaxTemp = as.logical(p[,6]),
                                JanMaxTemp = as.logical(p[,7]),
                                FebMaxTemp = as.logical(p[,8]),
                                MarMaxTemp = as.logical(p[,9]),
                                AprMaxTemp = as.logical(p[,10]),
                                MayMaxTemp = as.logical(p[,11]),
                                JunMaxTemp = as.logical(p[,12]),
                                JulMaxTemp = as.logical(p[,13]),
                                AugMaxTemp = as.logical(p[,14]),
                                SepMaxTemp = as.logical(p[,15]),
                                looic = rep(0),
                                p_loo = rep(0))

for (i in 2:32767) {
  M = sum(p[i,])
  z_temp <- z[as.logical(p[i,]),]
  
  if (M==1) { z_temp <- matrix(z_temp,nrow=1,ncol=K) }
  
  mvn_months_model_temp <- stan(file="linear_flu_mvn.stan",data=list(K=K, M=M, x=z_temp, y=y),
                                iter=3000, chains=4)
  
  log_lik_months <- extract_log_lik(mvn_months_model_temp)
  # calculate information criteria
  loo_mvn_months <- loo(log_lik_months)
  months_model_loos$looic[i] <- loo_mvn_months$estimates[3,1]
  months_model_loos$p_loo[i] <- loo_mvn_months$estimates[2,1]
  
  # save csv in case of crash
  write.csv(months_model_loos,"months_model_loos.csv",row.names = FALSE)
}

# ==============================================================================
# Test duplicate Lat column
# ==============================================================================

# all covariates
x <- array(dim=c(8,K))
x[1,] <- flu_data_norm$PopDensity
x[2,] <- flu_data_norm$Under18
x[3,] <- flu_data_norm$MeanMaxTemp
x[4,] <- flu_data_norm$MeanMaxRelHum
x[5,] <- flu_data_norm$MeanMinRelHum
x[6,] <- flu_data_norm$MeanMeanAbsHum
x[7,] <- flu_data_norm$Lat
x[8,] <- flu_data_norm$Lat

p <- expand.grid(1:0,1:0,1:0,1:0,1:0,1:0,1:0,1:0)

model_loos_1 <- data.frame(PopDensity = as.logical(p[,1]),
                         Under18 = as.logical(p[,2]),
                         MeanMaxTemp = as.logical(p[,3]),
                         MeanMaxRelHum = as.logical(p[,4]),
                         MeanMinRelHum = as.logical(p[,5]),
                         MeanMeanAbsHum = as.logical(p[,6]),
                         Lat = as.logical(p[,7]),
                         Lat2 = as.logical(p[,8]),
                         looic = rep(0),
                         p_loo = rep(0))

for (i in 1) {
  M = sum(p[i,])
  x_temp <- x[as.logical(p[i,]),]
  
  if (M==1) { x_temp <- matrix(x_temp,nrow=1,ncol=K) }
  
  mvn_model_temp <- stan(file="linear_flu_mvn.stan",data=list(K=K, M=M, x=x_temp, y=y),
                         iter=3000, chains=4)
  
  log_lik_temp <- extract_log_lik(mvn_model_temp)
  # calculate information criteria
  loo_mvn_temp <- loo(log_lik_temp)
  model_loos_1$looic[i] <- loo_mvn_temp$estimates[3,1]
  model_loos_1$p_loo[i] <- loo_mvn_temp$estimates[2,1]
  
  # save csv in case of crash
  write.csv(model_loos_1,"model_loos_1.csv",row.names = FALSE)
}

# what about just a regular (non Bayesian linear model)
flu_data_norm$Lat2 <- flu_data_norm$Lat

mlm <- lm(cbind(PeakWeek, UpSlope, DownSlope, PeakVal) ~ PopDensity + Under18 + MeanMaxTemp + MeanMaxRelHum + MeanMinRelHum + MeanMeanAbsHum + Lat + Lat2, data = flu_data_norm)
summary(mlm)
coef(mlm)

# ==============================================================================
# Test linear transform of Lat column
# ==============================================================================

# all covariates
x <- array(dim=c(8,K))
x[1,] <- flu_data_norm$PopDensity
x[2,] <- flu_data_norm$Under18
x[3,] <- flu_data_norm$MeanMaxTemp
x[4,] <- flu_data_norm$MeanMaxRelHum
x[5,] <- flu_data_norm$MeanMinRelHum
x[6,] <- flu_data_norm$MeanMeanAbsHum
x[7,] <- flu_data_norm$Lat
x[8,] <- -2*flu_data_norm$Lat+1

p <- expand.grid(1:0,1:0,1:0,1:0,1:0,1:0,1:0,1:0)

model_loos_2 <- data.frame(PopDensity = as.logical(p[,1]),
                           Under18 = as.logical(p[,2]),
                           MeanMaxTemp = as.logical(p[,3]),
                           MeanMaxRelHum = as.logical(p[,4]),
                           MeanMinRelHum = as.logical(p[,5]),
                           MeanMeanAbsHum = as.logical(p[,6]),
                           Lat = as.logical(p[,7]),
                           Lat2 = as.logical(p[,8]),
                           looic = rep(0),
                           p_loo = rep(0))

for (i in 1) {
  M = sum(p[i,])
  x_temp <- x[as.logical(p[i,]),]
  
  if (M==1) { x_temp <- matrix(x_temp,nrow=1,ncol=K) }
  
  mvn_model_temp <- stan(file="linear_flu_mvn.stan",data=list(K=K, M=M, x=x_temp, y=y),
                         iter=3000, chains=4)
  
  log_lik_temp <- extract_log_lik(mvn_model_temp)
  # calculate information criteria
  loo_mvn_temp <- loo(log_lik_temp)
  model_loos_2$looic[i] <- loo_mvn_temp$estimates[3,1]
  model_loos_2$p_loo[i] <- loo_mvn_temp$estimates[2,1]
  
  # save csv in case of crash
  write.csv(model_loos_2,"model_loos_2.csv",row.names = FALSE)
}

peakwk_pars <- c("alpha[1]")
for (i in 1:M) {
  peakwk_pars <- c(peakwk_pars, paste0("beta[1,",i,"]"))
}

peakwk_parameter_plot <- plot(mvn_model_temp, show_density = TRUE, fill_color = "lightsteelblue1", pars=peakwk_pars) + 
  ggtitle("Posterior distributions for peak week parameters") +
  theme(plot.title = element_text(face = 'bold', size =18, hjust = 0.6),
        axis.text=element_text(size=16))

traceplot(mvn_model_temp, pars=peakwk_pars, inc_warmup = TRUE)

# what about just a regular (non Bayesian linear model)
flu_data_norm$Lat2 <- -2*flu_data_norm$Lat+1

mlm <- lm(cbind(PeakWeek, UpSlope, DownSlope, PeakVal) ~ PopDensity + Under18 + MeanMaxTemp + MeanMaxRelHum + MeanMinRelHum + MeanMeanAbsHum + Lat2 + Lat, data = flu_data_norm)
summary(mlm)
coef(mlm)


# ==============================================================================
# Test linear transform of Lat column plus random noise
# ==============================================================================

flu_data_norm$Lat2 <- -0.94*flu_data_norm$Lat+rnorm(K,0,0.34)

# all covariates
x <- array(dim=c(8,K))
x[1,] <- flu_data_norm$PopDensity
x[2,] <- flu_data_norm$Under18
x[3,] <- flu_data_norm$MeanMaxTemp
x[4,] <- flu_data_norm$MeanMaxRelHum
x[5,] <- flu_data_norm$MeanMinRelHum
x[6,] <- flu_data_norm$MeanMeanAbsHum
x[7,] <- flu_data_norm$Lat
x[8,] <- flu_data_norm$Lat2

p <- expand.grid(1:0,1:0,1:0,1:0,1:0,1:0,1:0,1:0)

model_loos_3 <- data.frame(PopDensity = as.logical(p[,1]),
                           Under18 = as.logical(p[,2]),
                           MeanMaxTemp = as.logical(p[,3]),
                           MeanMaxRelHum = as.logical(p[,4]),
                           MeanMinRelHum = as.logical(p[,5]),
                           MeanMeanAbsHum = as.logical(p[,6]),
                           Lat = as.logical(p[,7]),
                           Lat2 = as.logical(p[,8]),
                           looic = rep(0),
                           p_loo = rep(0))

for (i in 2:255) {
  M = sum(p[i,])
  x_temp <- x[as.logical(p[i,]),]
  
  if (M==1) { x_temp <- matrix(x_temp,nrow=1,ncol=K) }
  
  mvn_model_temp <- stan(file="linear_flu_mvn.stan",data=list(K=K, M=M, x=x_temp, y=y),
                         iter=3000, chains=4)
  
  log_lik_temp <- extract_log_lik(mvn_model_temp)
  # calculate information criteria
  loo_mvn_temp <- loo(log_lik_temp)
  model_loos_3$looic[i] <- loo_mvn_temp$estimates[3,1]
  model_loos_3$p_loo[i] <- loo_mvn_temp$estimates[2,1]
  
  # save csv in case of crash
  write.csv(model_loos_3,"model_loos_3.csv",row.names = FALSE)
}

peakwk_pars <- c("alpha[1]")
for (i in 1:M) {
  peakwk_pars <- c(peakwk_pars, paste0("beta[1,",i,"]"))
}

peakwk_parameter_plot <- plot(mvn_model_temp, show_density = TRUE, fill_color = "lightsteelblue1", pars=peakwk_pars) + 
  ggtitle("Posterior distributions for peak week parameters") +
  theme(plot.title = element_text(face = 'bold', size =18, hjust = 0.6),
        axis.text=element_text(size=16))

traceplot(mvn_model_temp, pars=peakwk_pars)

# what about just a regular (non Bayesian) linear model
mlm <- lm(cbind(PeakWeek, UpSlope, DownSlope, PeakVal) ~ PopDensity + Under18 + MeanMaxTemp + MeanMaxRelHum + MeanMinRelHum + MeanMeanAbsHum + Lat, data = flu_data_norm)
summary(mlm)
coef(mlm)

# just the runs with Lat2 and without MeanMaxTemp
model_loos_lat2 <- model_loos_3[(model_loos_3$Lat2==TRUE & model_loos_3$MeanMaxTemp==FALSE),]

# just the runs with MeanMaxTemp and without Lat2
model_loos_temp <- model_loos_3[(model_loos_3$Lat2==FALSE & model_loos_3$MeanMaxTemp==TRUE),]

# histogram of differences between IC with temp vs lat2
hist(model_loos_lat2$looic-model_loos_temp$looic)

# models with temp rather than "lat2" generally do better, as expected


# ==============================================================================
# Test linear transform of Lat column plus random noise, but subtracting Lat line
# from Lat2 and MeanMaxTemp
# ==============================================================================

flu_data_norm$Lat2 <- -0.94*flu_data_norm$Lat+rnorm(K,0,0.34)

# all covariates
x <- array(dim=c(8,K))
x[1,] <- flu_data_norm$PopDensity
x[2,] <- flu_data_norm$Under18
x[3,] <- flu_data_norm$MeanMaxTemp+0.94*flu_data_norm$Lat
x[4,] <- flu_data_norm$MeanMaxRelHum
x[5,] <- flu_data_norm$MeanMinRelHum
x[6,] <- flu_data_norm$MeanMeanAbsHum
x[7,] <- flu_data_norm$Lat
x[8,] <- flu_data_norm$Lat2+0.94*flu_data_norm$Lat

p <- expand.grid(1:0,1:0,1:0,1:0,1:0,1:0,1:0,1:0)

model_loos_4 <- data.frame(PopDensity = as.logical(p[,1]),
                           Under18 = as.logical(p[,2]),
                           MeanMaxTemp = as.logical(p[,3]),
                           MeanMaxRelHum = as.logical(p[,4]),
                           MeanMinRelHum = as.logical(p[,5]),
                           MeanMeanAbsHum = as.logical(p[,6]),
                           Lat = as.logical(p[,7]),
                           Lat2 = as.logical(p[,8]),
                           looic = rep(0),
                           p_loo = rep(0))

for (i in 2:255) {
  M = sum(p[i,])
  x_temp <- x[as.logical(p[i,]),]
  
  if (M==1) { x_temp <- matrix(x_temp,nrow=1,ncol=K) }
  
  mvn_model_temp <- stan(file="linear_flu_mvn.stan",data=list(K=K, M=M, x=x_temp, y=y),
                         iter=3000, chains=4)
  
  log_lik_temp <- extract_log_lik(mvn_model_temp)
  # calculate information criteria
  loo_mvn_temp <- loo(log_lik_temp)
  model_loos_4$looic[i] <- loo_mvn_temp$estimates[3,1]
  model_loos_4$p_loo[i] <- loo_mvn_temp$estimates[2,1]
  
  # save csv in case of crash
  write.csv(model_loos_4,"model_loos_4.csv",row.names = FALSE)
}

peakwk_pars <- c("alpha[1]")
for (i in 1:M) {
  peakwk_pars <- c(peakwk_pars, paste0("beta[1,",i,"]"))
}

peakwk_parameter_plot <- plot(mvn_model_temp, show_density = TRUE, fill_color = "lightsteelblue1", pars=peakwk_pars) + 
  ggtitle("Posterior distributions for peak week parameters") +
  theme(plot.title = element_text(face = 'bold', size =18, hjust = 0.6),
        axis.text=element_text(size=16))

traceplot(mvn_model_temp, pars=peakwk_pars)

# what about just a regular (non Bayesian) linear model
mlm <- lm(cbind(flu_data_norm$PeakWeek, flu_data_norm$UpSlope, flu_data_norm$DownSlope, flu_data_norm$PeakVal) 
          ~ x_temp[1,]+x_temp[2,]+x_temp[3,]+x_temp[4,]+x_temp[5,]+x_temp[6,]+x_temp[7,])
summary(mlm)
coef(mlm)

# just the runs with Lat2 and without MeanMaxTemp
model_loos_lat2_resid <- model_loos_4[(model_loos_4$Lat2==TRUE & model_loos_4$MeanMaxTemp==FALSE),]

# just the runs with MeanMaxTemp and without Lat2
model_loos_temp_resid <- model_loos_4[(model_loos_4$Lat2==FALSE & model_loos_4$MeanMaxTemp==TRUE),]

# histogram of differences between IC with temp vs lat2
hist(model_loos_lat2_resid$looic-model_loos_temp_resid$looic)

# ==============================================================================
# Try ignoring data points that are causing high pareto k values
# ==============================================================================

# regular temp
x <- array(dim=c(8,(K-2)))
x[1,] <- flu_data_norm$PopDensity[c(1:29,31:156,158:324)]
x[2,] <- flu_data_norm$Under18[c(1:29,31:156,158:324)]
x[3,] <- flu_data_norm$MeanMaxTemp[c(1:29,31:156,158:324)]
x[4,] <- flu_data_norm$MeanMaxRelHum[c(1:29,31:156,158:324)]
x[5,] <- flu_data_norm$MeanMinRelHum[c(1:29,31:156,158:324)]
x[6,] <- flu_data_norm$MeanMeanAbsHum[c(1:29,31:156,158:324)]
x[7,] <- flu_data_norm$Lat[c(1:29,31:156,158:324)]
x[8,] <- flu_data_norm$Lat2[c(1:29,31:156,158:324)]

# temp residuals
x <- array(dim=c(8,(K-2)))
x[1,] <- flu_data_norm$PopDensity[c(1:29,31:156,158:324)]
x[2,] <- flu_data_norm$Under18[c(1:29,31:156,158:324)]
x[3,] <- flu_data_norm$MeanMaxTemp[c(1:29,31:156,158:324)]+0.94*flu_data_norm$Lat[c(1:29,31:156,158:324)]
x[4,] <- flu_data_norm$MeanMaxRelHum[c(1:29,31:156,158:324)]
x[5,] <- flu_data_norm$MeanMinRelHum[c(1:29,31:156,158:324)]
x[6,] <- flu_data_norm$MeanMeanAbsHum[c(1:29,31:156,158:324)]
x[7,] <- flu_data_norm$Lat[c(1:29,31:156,158:324)]
x[8,] <- flu_data_norm$Lat2[c(1:29,31:156,158:324)]+0.94*flu_data_norm$Lat[c(1:29,31:156,158:324)]



y <- array(dim=c((K-2),4))
y[,1] <- flu_data_norm$PeakWeek[c(1:29,31:156,158:324)]
y[,2] <- flu_data_norm$UpSlope[c(1:29,31:156,158:324)]
y[,3] <- flu_data_norm$DownSlope[c(1:29,31:156,158:324)]
y[,4] <- flu_data_norm$PeakVal[c(1:29,31:156,158:324)]


# ==============================================================================
# Test linear transform of Lat column plus random noise, but subtracting Lat line
# from Lat2 and MeanMaxTemp, without the two data points with high k values
# ==============================================================================

flu_data_norm$Lat2 <- -0.94*flu_data_norm$Lat+rnorm(K,0,0.34)

x <- array(dim=c(8,(K-2)))
x[1,] <- flu_data_norm$PopDensity[c(1:29,31:156,158:324)]
x[2,] <- flu_data_norm$Under18[c(1:29,31:156,158:324)]
x[3,] <- flu_data_norm$MeanMaxTemp[c(1:29,31:156,158:324)]+0.94*flu_data_norm$Lat[c(1:29,31:156,158:324)]
x[4,] <- flu_data_norm$MeanMaxRelHum[c(1:29,31:156,158:324)]
x[5,] <- flu_data_norm$MeanMinRelHum[c(1:29,31:156,158:324)]
x[6,] <- flu_data_norm$MeanMeanAbsHum[c(1:29,31:156,158:324)]
x[7,] <- flu_data_norm$Lat[c(1:29,31:156,158:324)]
x[8,] <- flu_data_norm$Lat2[c(1:29,31:156,158:324)]+0.94*flu_data_norm$Lat[c(1:29,31:156,158:324)]

y <- array(dim=c((K-2),4))
y[,1] <- flu_data_norm$PeakWeek[c(1:29,31:156,158:324)]
y[,2] <- flu_data_norm$UpSlope[c(1:29,31:156,158:324)]
y[,3] <- flu_data_norm$DownSlope[c(1:29,31:156,158:324)]
y[,4] <- flu_data_norm$PeakVal[c(1:29,31:156,158:324)]

p <- expand.grid(1:0,1:0,1:0,1:0,1:0,1:0,1:0,1:0)

model_loos_5 <- data.frame(PopDensity = as.logical(p[,1]),
                           Under18 = as.logical(p[,2]),
                           MeanMaxTemp = as.logical(p[,3]),
                           MeanMaxRelHum = as.logical(p[,4]),
                           MeanMinRelHum = as.logical(p[,5]),
                           MeanMeanAbsHum = as.logical(p[,6]),
                           Lat = as.logical(p[,7]),
                           Lat2 = as.logical(p[,8]),
                           looic = rep(0),
                           p_loo = rep(0))

for (i in 80) {
  M = sum(p[i,])
  x_temp <- x[as.logical(p[i,]),]
  
  if (M==1) { x_temp <- matrix(x_temp,nrow=1,ncol=K-2) }
  
  mvn_model_temp <- stan(file="linear_flu_mvn.stan",data=list(K=K-2, M=M, x=x_temp, y=y),
                         iter=3000, chains=4)
  
  log_lik_temp <- extract_log_lik(mvn_model_temp)
  # calculate information criteria
  loo_mvn_temp <- loo(log_lik_temp)
  model_loos_5$looic[i] <- loo_mvn_temp$estimates[3,1]
  model_loos_5$p_loo[i] <- loo_mvn_temp$estimates[2,1]
  
  # save csv in case of crash
  write.csv(model_loos_5,"model_loos_5.csv",row.names = FALSE)
}

peakwk_pars <- c("alpha[1]")
for (i in 1:M) {
  peakwk_pars <- c(peakwk_pars, paste0("beta[1,",i,"]"))
}

peakwk_parameter_plot <- plot(mvn_model_temp, show_density = TRUE, fill_color = "lightsteelblue1", pars=peakwk_pars) + 
  ggtitle("Posterior distributions for peak week parameters") +
  theme(plot.title = element_text(face = 'bold', size =18, hjust = 0.6),
        axis.text=element_text(size=16))

traceplot(mvn_model_temp, pars=peakwk_pars)

# what about just a regular (non Bayesian) linear model
mlm <- lm(cbind(flu_data_norm$PeakWeek, flu_data_norm$UpSlope, flu_data_norm$DownSlope, flu_data_norm$PeakVal) 
          ~ x_temp[1,]+x_temp[2,]+x_temp[3,]+x_temp[4,]+x_temp[5,]+x_temp[6,]+x_temp[7,])
summary(mlm)
coef(mlm)

# just the runs with Lat2 and without MeanMaxTemp
model_loos_lat2_resid_5 <- model_loos_5[(model_loos_5$Lat2==TRUE & model_loos_5$MeanMaxTemp==FALSE),]

# just the runs with MeanMaxTemp and without Lat2
model_loos_temp_resid_5 <- model_loos_5[(model_loos_5$Lat2==FALSE & model_loos_5$MeanMaxTemp==TRUE),]

# histogram of differences between IC with temp vs lat2
hist(model_loos_lat2_resid_5$looic-model_loos_temp_resid_5$looic)

