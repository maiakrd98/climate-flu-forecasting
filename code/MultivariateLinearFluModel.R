# Multivariate

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
flu_data <- read.csv("flu_data(2).csv")
LogILIplus <- read.csv("log_ILIplus(2).csv")

states = unique(flu_data$State)
years = unique(flu_data$Year)
S = length(states) # number of states
N = length(years)
M = 6#7 # number of covariates

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

# covariates
x <- array(dim=c(M,K))
x[1,] <- flu_data_norm$PopDensity
x[2,] <- flu_data_norm$Under18
x[3,] <- flu_data_norm$MeanMaxTemp
x[4,] <- flu_data_norm$MeanMaxRelHum
x[5,] <- flu_data_norm$MeanMinRelHum
x[6,] <- flu_data_norm$MeanMeanAbsHum

# the things we are predicting (there's definitely a real name for that...)
y <- array(dim=c(K,4))
y[,1] <- flu_data_norm$PeakWeek
y[,2] <- flu_data_norm$UpSlope
y[,3] <- flu_data_norm$DownSlope
y[,4] <- flu_data_norm$PeakVal

# ==============================================================================
# Run multivariate linear model in stan 
# ==============================================================================
mvn_model <- stan(file="linear_flu_mvn.stan",data=list(K=K, M=M, x=x, y=y), 
                     iter=3000, chains=4)

# ==============================================================================
# Check convergence
# ==============================================================================

# look at Rhat and effective sample size
mvn_model

# plot chains
traceplot(mvn_model, pars = c("alpha", "beta", "sigma"), inc_warmup = F)

# ==============================================================================
# Cross validation with loo package
# ==============================================================================

# extract log likelihood values to use in loo package
log_lik <- extract_log_lik(mvn_model)
# calculate information criteria
waic_mvn <- waic(log_lik)
loo_mvn <- loo(log_lik)

pareto_k_table(x)
plot(loo_mvn)

# ==============================================================================
# Plot posterior parameter distributions
# ==============================================================================

# peak week
peakwk_pars <- c("alpha[1]")
for (i in 1:M) {
  peakwk_pars <- c(peakwk_pars, paste0("beta[1,",i,"]"))
}

peakwk_parameter_plot <- plot(mvn_model, show_density = TRUE, fill_color = "lightsteelblue1", pars=peakwk_pars) + 
  ggtitle("Posterior distributions for peak week parameters") +
  theme(plot.title = element_text(face = 'bold', size =18, hjust = 0.6),
        axis.text=element_text(size=16))

# up slope
upslope_pars <- c("alpha[2]")
for (i in 1:M) {
  upslope_pars <- c(upslope_pars, paste0("beta[2,",i,"]"))
}

upslope_parameter_plot <- plot(mvn_model, show_density = TRUE, fill_color = "lightsteelblue1", pars=upslope_pars) + 
  ggtitle("Posterior distributions for up slope parameters") +
  theme(plot.title = element_text(face = 'bold', size =18, hjust = 0.6),
        axis.text=element_text(size=16))

# down slope
downslope_pars <- c("alpha[3]")
for (i in 1:M) {
  downslope_pars <- c(downslope_pars, paste0("beta[3,",i,"]"))
}

downslope_parameter_plot <- plot(mvn_model, show_density = TRUE, fill_color = "lightsteelblue1", pars=downslope_pars) + 
  ggtitle("Posterior distributions for down slope parameters") +
  theme(plot.title = element_text(face = 'bold', size =18, hjust = 0.6),
        axis.text=element_text(size=16))

# peak value
pkval_pars <- c("alpha[4]")
for (i in 1:M) {
  pkval_pars <- c(pkval_pars, paste0("beta[4,",i,"]"))
}

pkval_parameter_plot <- plot(mvn_model, show_density = TRUE, fill_color = "lightsteelblue1", pars=pkval_pars) + 
  ggtitle("Posterior distributions for peak value parameters") +
  theme(plot.title = element_text(face = 'bold', size =18, hjust = 0.6),
        axis.text=element_text(size=16))

# covariance matrix
cov_parameter_plot <- plot(mvn_model, show_density = TRUE, fill_color = "lightsteelblue1", pars="sigma") + 
  ggtitle("Posterior distributions for covariance matrix") +
  theme(plot.title = element_text(face = 'bold', size =18, hjust = 0.6),
        axis.text=element_text(size=16))


# ==============================================================================
# Predict peak week, up slope, etc for each year and state
# ==============================================================================
postdraws <- extract(mvn_model)

# these are normalized!
predicted_parameters <- tibble(Year = flu_data_norm$Year,
                               State = flu_data_norm$State,
                               PeakWeek = rep(0),
                               UpSlope = rep(0),
                               DownSlope = rep(0),
                               PeakVal = rep(0))

# predict normalized peak week, etc.
for (i in 1:(dim(flu_data_norm)[1])) {
  if (predicted_parameters$Year[i]!=flu_data_norm$Year[i]) {
    print("Years are not equal!")
  }
  else if (predicted_parameters$State[i]!=flu_data_norm$State[i]) {
    print("States are not equal!")
  }
  else {
    predicted_parameters$PeakWeek[i] = mean(postdraws$alpha[,1])+
      mean(postdraws$beta[,1,1])*flu_data_norm$PopDensity[i]+
      mean(postdraws$beta[,1,2])*flu_data_norm$Under18[i]+
      mean(postdraws$beta[,1,3])*flu_data_norm$MeanMaxTemp[i]+
      mean(postdraws$beta[,1,4])*flu_data_norm$MeanMaxRelHum[i]+
      mean(postdraws$beta[,1,5])*flu_data_norm$MeanMinRelHum[i]+
      mean(postdraws$beta[,1,6])*flu_data_norm$MeanMeanAbsHum[i]
    predicted_parameters$UpSlope[i] = mean(postdraws$alpha[,2])+
      mean(postdraws$beta[,2,1])*flu_data_norm$PopDensity[i]+
      mean(postdraws$beta[,2,2])*flu_data_norm$Under18[i]+
      mean(postdraws$beta[,2,3])*flu_data_norm$MeanMaxTemp[i]+
      mean(postdraws$beta[,2,4])*flu_data_norm$MeanMaxRelHum[i]+
      mean(postdraws$beta[,2,5])*flu_data_norm$MeanMinRelHum[i]+
      mean(postdraws$beta[,2,6])*flu_data_norm$MeanMeanAbsHum[i]
    predicted_parameters$DownSlope[i] = mean(postdraws$alpha[,3])+
      mean(postdraws$beta[,3,1])*flu_data_norm$PopDensity[i]+
      mean(postdraws$beta[,3,2])*flu_data_norm$Under18[i]+
      mean(postdraws$beta[,3,3])*flu_data_norm$MeanMaxTemp[i]+
      mean(postdraws$beta[,3,4])*flu_data_norm$MeanMaxRelHum[i]+
      mean(postdraws$beta[,3,5])*flu_data_norm$MeanMinRelHum[i]+
      mean(postdraws$beta[,3,6])*flu_data_norm$MeanMeanAbsHum[i]
    predicted_parameters$PeakVal[i] = mean(postdraws$alpha[,4])+
      mean(postdraws$beta[,4,1])*flu_data_norm$PopDensity[i]+
      mean(postdraws$beta[,4,2])*flu_data_norm$Under18[i]+
      mean(postdraws$beta[,4,3])*flu_data_norm$MeanMaxTemp[i]+
      mean(postdraws$beta[,4,4])*flu_data_norm$MeanMaxRelHum[i]+
      mean(postdraws$beta[,4,5])*flu_data_norm$MeanMinRelHum[i]+
      mean(postdraws$beta[,4,6])*flu_data_norm$MeanMeanAbsHum[i]
  }
}

# un-normalize
for (i in 1:(dim(flu_data_norm)[1])) {
  predicted_parameters$PeakWeek[i] = (predicted_parameters$PeakWeek[i]*means_and_sds[means_and_sds$Variable == 'PeakWeek','SD']) + means_and_sds[means_and_sds$Variable == 'PeakWeek','Mean']
  predicted_parameters$UpSlope[i] = (predicted_parameters$UpSlope[i]*means_and_sds[means_and_sds$Variable == 'UpSlope','SD']) + means_and_sds[means_and_sds$Variable == 'UpSlope','Mean']
  predicted_parameters$DownSlope[i] = (predicted_parameters$DownSlope[i]*means_and_sds[means_and_sds$Variable == 'DownSlope','SD']) + means_and_sds[means_and_sds$Variable == 'DownSlope','Mean']
  predicted_parameters$PeakVal[i] = (predicted_parameters$PeakVal[i]*means_and_sds[means_and_sds$Variable == 'PeakVal','SD']) + means_and_sds[means_and_sds$Variable == 'PeakVal','Mean']
}

# ==============================================================================
# Predict peak week, up slope, etc for each year and state, for first 100 draws
# Also add error drawn from multi_normal(0, sigma)
# ==============================================================================
postdraws <- extract(mvn_model)

# these are normalized!
predicted_parameters_100 <- data.frame(Year = rep(flu_data_norm$Year,4),
                                       State = rep(flu_data_norm$State,4),
                                       Parameter = c(rep("PeakWeek",K),rep("UpSlope",K),rep("DownSlope",K),rep("PeakVal",K)))

# predict normalized peak week, etc.
for (i in 1:K) {
  for (d in 1:100) { # d is what draw we are on
    # generate random draws from multivariate normal distribution
    mvn_draw <- rmvnorm(1, sigma = postdraws$sigma[d,1:4,1:4])
    
    # peak week
    predicted_parameters_100[i,d+3] = postdraws$alpha[d,1]+
      postdraws$beta[d,1,1]*flu_data_norm$PopDensity[i]+
      postdraws$beta[d,1,2]*flu_data_norm$Under18[i]+
      postdraws$beta[d,1,3]*flu_data_norm$MeanMaxTemp[i]+
      postdraws$beta[d,1,4]*flu_data_norm$MeanMaxRelHum[i]+
      postdraws$beta[d,1,5]*flu_data_norm$MeanMinRelHum[i]+
      postdraws$beta[d,1,6]*flu_data_norm$MeanMeanAbsHum[i]+
      mvn_draw[1]
    
    # up slope
    predicted_parameters_100[(i+K),d+3] = postdraws$alpha[d,2]+
      postdraws$beta[d,2,1]*flu_data_norm$PopDensity[i]+
      postdraws$beta[d,2,2]*flu_data_norm$Under18[i]+
      postdraws$beta[d,2,3]*flu_data_norm$MeanMaxTemp[i]+
      postdraws$beta[d,2,4]*flu_data_norm$MeanMaxRelHum[i]+
      postdraws$beta[d,2,5]*flu_data_norm$MeanMinRelHum[i]+
      postdraws$beta[d,2,6]*flu_data_norm$MeanMeanAbsHum[i]+
      mvn_draw[2]
    
    # down slope
    predicted_parameters_100[(i+2*K),d+3] = postdraws$alpha[d,3]+      
      postdraws$beta[d,3,1]*flu_data_norm$PopDensity[i]+
      postdraws$beta[d,3,2]*flu_data_norm$Under18[i]+
      postdraws$beta[d,3,3]*flu_data_norm$MeanMaxTemp[i]+
      postdraws$beta[d,3,4]*flu_data_norm$MeanMaxRelHum[i]+
      postdraws$beta[d,3,5]*flu_data_norm$MeanMinRelHum[i]+
      postdraws$beta[d,3,6]*flu_data_norm$MeanMeanAbsHum[i]+
      mvn_draw[3]
    
    # peak value
    predicted_parameters_100[(i+3*K),d+3] = postdraws$alpha[d,4]+
      postdraws$beta[d,4,1]*flu_data_norm$PopDensity[i]+
      postdraws$beta[d,4,2]*flu_data_norm$Under18[i]+
      postdraws$beta[d,4,3]*flu_data_norm$MeanMaxTemp[i]+
      postdraws$beta[d,4,4]*flu_data_norm$MeanMaxRelHum[i]+
      postdraws$beta[d,4,5]*flu_data_norm$MeanMinRelHum[i]+
      postdraws$beta[d,4,6]*flu_data_norm$MeanMeanAbsHum[i]+
      mvn_draw[4]
  }
}


# un-normalize
for (i in 1:K) {
  predicted_parameters_100[i,4:103] = (predicted_parameters_100[i,4:103]*means_and_sds[means_and_sds$Variable == 'PeakWeek','SD']) + means_and_sds[means_and_sds$Variable == 'PeakWeek','Mean']
  predicted_parameters_100[(i+K),4:103] = (predicted_parameters_100[(i+K),4:103]*means_and_sds[means_and_sds$Variable == 'UpSlope','SD']) + means_and_sds[means_and_sds$Variable == 'UpSlope','Mean']
  predicted_parameters_100[(i+2*K),4:103] = (predicted_parameters_100[(i+2*K),4:103]*means_and_sds[means_and_sds$Variable == 'DownSlope','SD']) + means_and_sds[means_and_sds$Variable == 'DownSlope','Mean']
  predicted_parameters_100[(i+3*K),4:103] = (predicted_parameters_100[(i+3*K),4:103]*means_and_sds[means_and_sds$Variable == 'PeakVal','SD']) + means_and_sds[means_and_sds$Variable == 'PeakVal','Mean']
}


# ==============================================================================
# Plot ILI+ data, actual fit, predicted fit for mean and first 100 draws
# ==============================================================================

for (i in 1:K) {
  if (predicted_parameters$Year[i]!=flu_data_1$Year[i]) {
    print("Years are not equal!")
  }
  else if (predicted_parameters$State[i]!=flu_data_1$State[i]) {
    print("States are not equal!")
  }
  else {
    # create data frames to plot
    week <- 1:52
    LogILIplus_toplot <- data.frame(Week = week, LogILI = t(LogILIplus[LogILIplus$State==flu_data_1$State[i]&LogILIplus$Year==flu_data_1$Year[i],3:54]))
    names(LogILIplus_toplot)[2] <- "LogILIplus"
    
    wvals <- seq(from = 0, to = 52, by = 0.05)
    pred_vals <- (wvals < flu_data_1$Break1[i]) * (flu_data_1$Val1[i]+flu_data_1$Slope1[i]*wvals) + 
      (wvals >= flu_data_1$Break1[i] & wvals <= flu_data_1$Break2[i]) * (flu_data_1$Val2[i]+flu_data_1$Slope2[i]*(wvals-flu_data_1$Break1[i])) + 
      (wvals > flu_data_1$Break2[i]) * (flu_data_1$Val3[i]+flu_data_1$Slope3[i]*(wvals-flu_data_1$Break2[i]))
    pred_line_toplot <- data.frame(Wvals = wvals, PredVals = pred_vals)
    # add mean predicted fit line
    pred_line_toplot$PredValMean <- (wvals < predicted_parameters$PeakWeek[i]) * (predicted_parameters$PeakVal[i]-predicted_parameters$UpSlope[i]*(predicted_parameters$PeakWeek[i]-wvals)) + 
      (wvals >= predicted_parameters$PeakWeek[i]) * (predicted_parameters$PeakVal[i]+predicted_parameters$DownSlope[i]*(wvals-predicted_parameters$PeakWeek[i]))
    # add predicted lines for first 100 draws to data frame
    for (k in 4:103) {
      pred_line_toplot[[paste("PredVal", k, sep = "")]] <- (wvals < predicted_parameters_100[i,k]) * (predicted_parameters_100[(i+3*K),k]-predicted_parameters_100[(i+K),k]*(predicted_parameters_100[i,k]-wvals)) + 
        (wvals >= predicted_parameters_100[i,k]) * (predicted_parameters_100[(i+3*K),k]+predicted_parameters_100[(i+2*K),k]*(wvals-predicted_parameters_100[i,k]))
    }
    
    # initialize plot formatting
    ILIplus_plot <- ggplot() +
      theme_classic() +
      ggtitle(sprintf("Flu incidence in %s over the %i-%i season", flu_data_1$State[i], flu_data_1$Year[i], flu_data_1$Year[i]+1))+
      ylim(min(LogILIplus_toplot$LogILIplus)-0.5,max(LogILIplus_toplot$LogILIplus+1.5)) +
      xlab("Week") + ylab("Log(ILI+)") +
      theme(plot.title = element_text(face = 'bold', size =16, hjust = 0.5))
    
    
    # predicted lines from first 100 draws
    for (k in 4:103) {
      col_name <- paste("PredVal", k, sep = "")
      ILIplus_plot <- ILIplus_plot + geom_line(data = pred_line_toplot, aes(x=Wvals, y=!!sym(col_name)), linewidth = 1, color = 'lavenderblush3', alpha = 0.2)
    }
    
    # mean predicted fit
    ILIplus_plot <- ILIplus_plot + geom_line(data = pred_line_toplot, aes(x=Wvals, y=PredValMean), linewidth = 1, color = '#8B4789')
    
    # segmented fit line
    ILIplus_plot <- ILIplus_plot + geom_line(data = pred_line_toplot, aes(x=Wvals, y=PredVals, color = 'segmented fit'), linewidth = 1, color = 'aquamarine4')
    
    # scatter plot of ILI+ data
    ILIplus_plot <- ILIplus_plot + geom_point(data = LogILIplus_toplot, aes(x=Week,y=LogILIplus), color = 'black')
    ggsave(filename = paste(flu_data_1$State[i],flu_data_1$Year[i],".png",sep = ""),
           plot = ILIplus_plot,
           path = "~/FluForecasting/FluForecasting/plots_prediction_interval_multivariate_zero_is_min",
           width = 7, height = 5)
  }
}


# ==============================================================================
# Plot peak week vs up slope, etc with error ellipse
# ==============================================================================
mvn_pkweeks <- c(data.matrix(predicted_parameters_100[1:K,4:103]))
mvn_upslopes <- c(data.matrix(predicted_parameters_100[(K+1):(2*K),4:103]))
mvn_downslopes <- c(data.matrix(predicted_parameters_100[(2*K+1):(3*K),4:103]))
mvn_pkvals <- c(data.matrix(predicted_parameters_100[(3*K+1):(4*K),4:103]))

pkwk_vs_upslope_plot <- ggplot() +
  #geom_point(aes(x=predicted_parameters_100$V7[1:K],y=predicted_parameters_100$V7[(K+1):(2*K)]),alpha=0.3,color="blue")+
  geom_point(aes(x=mvn_pkweeks,y=mvn_upslopes),alpha=0.02,color="deepskyblue3")+
  geom_point(aes(x=flu_data_1$PeakWeek,y=flu_data_1$UpSlope),alpha=0.6)+
  stat_ellipse(aes(x=mvn_pkweeks,y=mvn_upslopes), color="deepskyblue3")+
  xlim(min(flu_data_1$PeakWeek)-7,max(flu_data_1$PeakWeek)+4) +
  ylim(min(flu_data_1$UpSlope)-0.13,max(flu_data_1$UpSlope)+0.05) +
  xlab("Peak Week") + ylab("Up Slope") +
  theme_classic()

pkwk_vs_downslope_plot <- ggplot() +
  geom_point(aes(x=mvn_pkweeks,y=mvn_downslopes),alpha=0.02,color="deepskyblue3")+
  geom_point(aes(x=flu_data_1$PeakWeek,y=flu_data_1$DownSlope),alpha=0.6)+
  stat_ellipse(aes(x=mvn_pkweeks,y=mvn_downslopes), color="deepskyblue3")+
  xlim(min(flu_data_1$PeakWeek)-7,max(flu_data_1$PeakWeek)+4) +
  ylim(min(flu_data_1$DownSlope)-0.13,max(flu_data_1$DownSlope)+0.6) +
  xlab("Peak Week") + ylab("Down Slope") +
  theme_classic()

pkwk_vs_pkval_plot <- ggplot() +
  geom_point(aes(x=mvn_pkweeks,y=mvn_pkvals),alpha=0.02,color="deepskyblue3")+
  geom_point(aes(x=flu_data_1$PeakWeek,y=flu_data_1$PeakVal),alpha=0.6)+
  stat_ellipse(aes(x=mvn_pkweeks,y=mvn_pkvals), color="deepskyblue3")+
  xlim(min(flu_data_1$PeakWeek)-7,max(flu_data_1$PeakWeek)+4) +
  ylim(min(flu_data_1$PeakVal)-0.5,max(flu_data_1$PeakVal)+0.5) +
  xlab("Peak Week") + ylab("Peak Value") +
  theme_classic()

upslope_vs_downslope_plot <- ggplot() +
  geom_point(aes(x=mvn_upslopes,y=mvn_downslopes),alpha=0.02,color="deepskyblue3")+
  geom_point(aes(x=flu_data_1$UpSlope,y=flu_data_1$DownSlope),alpha=0.6)+
  stat_ellipse(aes(x=mvn_upslopes,y=mvn_downslopes), color="deepskyblue3")+
  xlim(min(flu_data_1$UpSlope)-0.05,max(flu_data_1$UpSlope)+0.05) +
  ylim(min(flu_data_1$DownSlope)-0.1,max(flu_data_1$DownSlope)+0.6) +
  xlab("Up Slope") + ylab("Down Slope") +
  theme_classic()

upslope_vs_pkval_plot <- ggplot() +
  geom_point(aes(x=mvn_upslopes,y=mvn_pkvals),alpha=0.02,color="deepskyblue3")+
  geom_point(aes(x=flu_data_1$UpSlope,y=flu_data_1$PeakVal),alpha=0.6)+
  stat_ellipse(aes(x=mvn_upslopes,y=mvn_pkvals), color="deepskyblue3")+
  xlim(min(flu_data_1$UpSlope)-0.05,max(flu_data_1$UpSlope)+0.05) +
  ylim(min(flu_data_1$PeakVal)-0.5,max(flu_data_1$PeakVal)+0.5) +
  xlab("Up Slope") + ylab("Peak Value") +
  theme_classic()

downslope_vs_pkval_plot <- ggplot() +
  geom_point(aes(x=mvn_downslopes,y=mvn_pkvals),alpha=0.02,color="deepskyblue3")+
  geom_point(aes(x=flu_data_1$DownSlope,y=flu_data_1$PeakVal),alpha=0.6)+
  stat_ellipse(aes(x=mvn_downslopes,y=mvn_pkvals), color="deepskyblue3")+
  xlim(min(flu_data_1$DownSlope)-0.05,max(flu_data_1$DownSlope)+0.6) +
  ylim(min(flu_data_1$PeakVal)-0.5,max(flu_data_1$PeakVal)+0.5) +
  xlab("Down Slope") + ylab("Peak Value") +
  theme_classic()
