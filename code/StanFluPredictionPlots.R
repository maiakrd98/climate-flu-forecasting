# ==============================================================================
# Import libraries
# ==============================================================================

library(tidyverse) 
library(purrr)
library(rstan)
#install.packages("plotly")
library(plotly)
library(reshape2)
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

# normalize covariates
x <- array(dim=c(K,M))
x[,1] <- flu_data_norm$PopDensity
x[,2] <- flu_data_norm$Under18
#x[,3] <- flu_data_norm$Lat
x[,3] <- flu_data_norm$MeanMaxTemp
x[,4] <- flu_data_norm$MeanMaxRelHum
x[,5] <- flu_data_norm$MeanMinRelHum
x[,6] <- flu_data_norm$MeanMeanAbsHum

# normalize the things we are predicting (there's definitely a real name for that...)
y <- array(dim=c(K,4))
y[,1] <- flu_data_norm$PeakWeek
y[,2] <- flu_data_norm$UpSlope
y[,3] <- flu_data_norm$DownSlope
y[,4] <- flu_data_norm$PeakVal

# ==============================================================================
# Run linear models in stan 
# ==============================================================================
peakwk_model <- stan(file="peakwk2(2).stan",data=list(N=K, M=M, x=x, y=y[,1]), 
                     iter=1000, chains=4) 
upslope_model <- stan(file="peakwk2(2).stan",data=list(N=K, M=M, x=x, y=y[,2]), 
                     iter=1000, chains=4) 
downslope_model <- stan(file="peakwk2(2).stan",data=list(N=K, M=M, x=x, y=y[,3]), 
                      iter=1000, chains=4) 
peakval_model <- stan(file="peakwk2(2).stan",data=list(N=K, M=M, x=x, y=y[,4]), 
                      iter=1000, chains=4)

postdraws_peakwk <- extract(peakwk_model)
postdraws_upslope <- extract(upslope_model)
postdraws_downslope <- extract(downslope_model)
postdraws_peakval <- extract(peakval_model)

peakwk_plot <- plot(peakwk_model, show_density = TRUE, fill_color = "lightsteelblue1") + 
  ggtitle("Posterior distributions for peak week parameters") +
  theme(plot.title = element_text(face = 'bold', size =18, hjust = 0.6),
        axis.text=element_text(size=16))

ggsave("no_lat_peakwk_plot.png", peakwk_plot, path = "~/FluForecasting/FluForecasting/march2025poster/sections/images",
       width = 7, height = 5)

upslope_plot <- plot(upslope_model, show_density = TRUE, fill_color = "lightsteelblue1") + 
  ggtitle("Posterior distributions for up slope parameters") +
  theme(plot.title = element_text(face = 'bold', size =18, hjust = 0.6),
        axis.text=element_text(size=16))

ggsave("no_lat_upslope_plot.png", upslope_plot, path = "~/FluForecasting/FluForecasting/march2025poster/sections/images",
       width = 7, height = 5)

downslope_plot <- plot(downslope_model, show_density = TRUE, fill_color = "lightsteelblue1") + 
  ggtitle("Posterior distributions for down slope parameters") +
  theme(plot.title = element_text(face = 'bold', size =18, hjust = 0.6),
        axis.text=element_text(size=16))

ggsave("no_lat_downslope_plot.png", downslope_plot, path = "~/FluForecasting/FluForecasting/march2025poster/sections/images",
       width = 7, height = 5)

peakval_plot <- plot(peakval_model, show_density = TRUE, fill_color = "lightsteelblue1") + 
  ggtitle("Posterior distributions for peak value parameters") +
  theme(plot.title = element_text(face = 'bold', size =18, hjust = 0.6),
        axis.text=element_text(size=12))

ggsave("no_lat_peakval_plot.png", peakval_plot, path = "~/FluForecasting/FluForecasting/march2025poster/sections/images",
       width = 7, height = 5)

# plot chains to verify convergence
traceplot(peakwk_model, inc_warmup=T)
traceplot(upslope_model, inc_warmup=T)
traceplot(downslope_model, inc_warmup=T)
traceplot(peakval_model, inc_warmup=T)

# ==============================================================================
# Plot histograms of parameter estimates
# ==============================================================================

postdraws_peakwk_t <- tibble(
  alpha = postdraws_peakwk$alpha,
  b1 = postdraws_peakwk$beta[,1],      
  b2 = postdraws_peakwk$beta[,2],
  b3 = postdraws_peakwk$beta[,3],
  b4 = postdraws_peakwk$beta[,4],
  b5 = postdraws_peakwk$beta[,5],
  b6 = postdraws_peakwk$beta[,6],
  b7 = postdraws_peakwk$beta[,7],
  sigma = postdraws_peakwk$sigma
) %>% 
  pivot_longer(everything())

postdraws_peakwk_t %>% 
  ggplot(aes(x=value, y=after_stat(density))) + 
  geom_histogram(bins=20, fill="white", col="grey")+ 
  theme_classic() + 
  facet_wrap(~name, scales="free")


# ==============================================================================
# Plot scatterplots
# ==============================================================================

latvals <- seq(from=-3, to=3, by=0.1)

flu_data_norm %>% 
  ggplot(aes(x=Lat, y=PeakWeek)) + 
  geom_point(alpha=0.6) + 
  geom_line(data=tibble(x=latvals, y=mean(postdraws_peakwk$alpha) + mean(postdraws_peakwk$beta[,3])*latvals), aes(x=x, y=y)) +
  theme_classic() 
flu_data_norm %>% 
  ggplot(aes(x=Lat, y=UpSlope)) + 
  geom_point(alpha=0.6) + 
  geom_line(data=tibble(x=latvals, y=mean(postdraws_peakwk$alpha) + mean(postdraws_peakwk$beta[,3])*latvals), aes(x=x, y=y)) +
  theme_classic() 

tempvals<- seq(from=-3, to=3, by=0.1)

flu_data_norm %>% 
  ggplot(aes(x=MeanMaxTemp, y=PeakWeek)) + 
  geom_point(alpha=0.6) + 
  geom_line(data=tibble(x=tempvals, y=mean(postdraws_peakwk$alpha) + mean(postdraws_peakwk$beta[,4])*tempvals), aes(x=x, y=y)) +
  theme_classic() 

plot_ly(flu_data_norm, x = ~MeanMaxTemp, y = ~Lat, z = ~PeakWeek, type = "scatter3d", mode = "markers",
             marker = list(size = 3))

temp_lat_surface <- expand.grid(MeanMaxTemp = tempvals,Lat = latvals,KEEP.OUT.ATTRS = F)
temp_lat_surface$PeakWeek <- mean(postdraws_peakwk$alpha)+mean(postdraws_peakwk$beta[,3])*temp_lat_surface$Lat+mean(postdraws_peakwk$beta[,4])*temp_lat_surface$MeanMaxTemp
temp_lat_surface <- acast(temp_lat_surface, Lat ~ MeanMaxTemp, value.var = "PeakWeek") #y ~ x

temp_lat_plot <- plot_ly(flu_data_norm, x = ~MeanMaxTemp, y = ~Lat, z = ~PeakWeek, type = "scatter3d")

temp_lat_plot <- add_trace(p = temp_lat_plot,
                           z = temp_lat_surface,
                           x = tempvals,
                           y = latvals,
                           type = "surface")

# ==============================================================================
# Predict peak week, up slope, etc for each year and state
# ==============================================================================
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
    predicted_parameters$PeakWeek[i] = mean(postdraws_peakwk$alpha)+
      mean(postdraws_peakwk$beta[,1])*flu_data_norm$PopDensity[i]+
      mean(postdraws_peakwk$beta[,2])*flu_data_norm$Under18[i]+
      #mean(postdraws_peakwk$beta[,3])*flu_data_norm$Lat[i]+
      mean(postdraws_peakwk$beta[,3])*flu_data_norm$MeanMaxTemp[i]+
      mean(postdraws_peakwk$beta[,4])*flu_data_norm$MeanMaxRelHum[i]+
      mean(postdraws_peakwk$beta[,5])*flu_data_norm$MeanMinRelHum[i]+
      mean(postdraws_peakwk$beta[,6])*flu_data_norm$MeanMeanAbsHum[i]
    predicted_parameters$UpSlope[i] = mean(postdraws_peakwk$alpha)+
      mean(postdraws_upslope$beta[,1])*flu_data_norm$PopDensity[i]+
      mean(postdraws_upslope$beta[,2])*flu_data_norm$Under18[i]+
      #mean(postdraws_upslope$beta[,3])*flu_data_norm$Lat[i]+
      mean(postdraws_upslope$beta[,3])*flu_data_norm$MeanMaxTemp[i]+
      mean(postdraws_upslope$beta[,4])*flu_data_norm$MeanMaxRelHum[i]+
      mean(postdraws_upslope$beta[,5])*flu_data_norm$MeanMinRelHum[i]+
      mean(postdraws_upslope$beta[,6])*flu_data_norm$MeanMeanAbsHum[i]
    predicted_parameters$DownSlope[i] = mean(postdraws_peakwk$alpha)+
      mean(postdraws_downslope$beta[,1])*flu_data_norm$PopDensity[i]+
      mean(postdraws_downslope$beta[,2])*flu_data_norm$Under18[i]+
      #mean(postdraws_downslope$beta[,3])*flu_data_norm$Lat[i]+
      mean(postdraws_downslope$beta[,3])*flu_data_norm$MeanMaxTemp[i]+
      mean(postdraws_downslope$beta[,4])*flu_data_norm$MeanMaxRelHum[i]+
      mean(postdraws_downslope$beta[,5])*flu_data_norm$MeanMinRelHum[i]+
      mean(postdraws_downslope$beta[,6])*flu_data_norm$MeanMeanAbsHum[i]
    predicted_parameters$PeakVal[i] = mean(postdraws_peakwk$alpha)+
      mean(postdraws_peakval$beta[,1])*flu_data_norm$PopDensity[i]+
      mean(postdraws_peakval$beta[,2])*flu_data_norm$Under18[i]+
      #mean(postdraws_peakval$beta[,3])*flu_data_norm$Lat[i]+
      mean(postdraws_peakval$beta[,3])*flu_data_norm$MeanMaxTemp[i]+
      mean(postdraws_peakval$beta[,4])*flu_data_norm$MeanMaxRelHum[i]+
      mean(postdraws_peakval$beta[,5])*flu_data_norm$MeanMinRelHum[i]+
      mean(postdraws_peakval$beta[,6])*flu_data_norm$MeanMeanAbsHum[i]
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
# Also add error drawn from Normal(0, sigma)
# ==============================================================================


# extract first 100 draws
postdf_peakwk_100 <- tibble(
  alpha = postdraws_peakwk$alpha[1:100],
  b1 = postdraws_peakwk$beta[1:100,1],
  b2 = postdraws_peakwk$beta[1:100,2],
  b4 = postdraws_peakwk$beta[1:100,3],
  b5 = postdraws_peakwk$beta[1:100,4],
  b6 = postdraws_peakwk$beta[1:100,5],
  b7 = postdraws_peakwk$beta[1:100,6]
) 

# predict (normalized) peak week, etc. using the first 100 posterior draws
predicted_peakwk_100 <- data.frame(Year = flu_data_norm$Year,
                                   State = flu_data_norm$State)
predicted_upslope_100 <- data.frame(Year = flu_data_norm$Year,
                                   State = flu_data_norm$State)
predicted_downslope_100 <- data.frame(Year = flu_data_norm$Year,
                                   State = flu_data_norm$State)
predicted_peakval_100 <- data.frame(Year = flu_data_norm$Year,
                                   State = flu_data_norm$State)

for (i in 1:(dim(flu_data_norm)[1])) {
  if (predicted_peakwk_100$Year[i]!=flu_data_norm$Year[i]) {
    print("Years are not equal!")
  }
  else if (predicted_peakwk_100$State[i]!=flu_data_norm$State[i]) {
    print("States are not equal!")
  }
  else {
    predicted_peakwk_100[i,3:102] <- postdraws_peakwk$alpha[1:100]+
      postdraws_peakwk$beta[1:100,1]*flu_data_norm$PopDensity[i]+
      postdraws_peakwk$beta[1:100,2]*flu_data_norm$Under18[i]+
      postdraws_peakwk$beta[1:100,3]*flu_data_norm$MeanMaxTemp[i]+
      postdraws_peakwk$beta[1:100,4]*flu_data_norm$MeanMaxRelHum[i]+
      postdraws_peakwk$beta[1:100,5]*flu_data_norm$MeanMinRelHum[i]+
      postdraws_peakwk$beta[1:100,6]*flu_data_norm$MeanMeanAbsHum[i]+
      rnorm(100, mean=0, sd=postdraws_peakwk$sigma[1:100])
    
    predicted_upslope_100[i,3:102] <- postdraws_upslope$alpha[1:100]+
      postdraws_upslope$beta[1:100,1]*flu_data_norm$PopDensity[i]+
      postdraws_upslope$beta[1:100,2]*flu_data_norm$Under18[i]+
      postdraws_upslope$beta[1:100,3]*flu_data_norm$MeanMaxTemp[i]+
      postdraws_upslope$beta[1:100,4]*flu_data_norm$MeanMaxRelHum[i]+
      postdraws_upslope$beta[1:100,5]*flu_data_norm$MeanMinRelHum[i]+
      postdraws_upslope$beta[1:100,6]*flu_data_norm$MeanMeanAbsHum[i]+
      rnorm(100, mean=0, sd=postdraws_upslope$sigma[1:100])
    
    predicted_downslope_100[i,3:102] <- postdraws_downslope$alpha[1:100]+
      postdraws_downslope$beta[1:100,1]*flu_data_norm$PopDensity[i]+
      postdraws_downslope$beta[1:100,2]*flu_data_norm$Under18[i]+
      postdraws_downslope$beta[1:100,3]*flu_data_norm$MeanMaxTemp[i]+
      postdraws_downslope$beta[1:100,4]*flu_data_norm$MeanMaxRelHum[i]+
      postdraws_downslope$beta[1:100,5]*flu_data_norm$MeanMinRelHum[i]+
      postdraws_downslope$beta[1:100,6]*flu_data_norm$MeanMeanAbsHum[i]+
      rnorm(100, mean=0, sd=postdraws_downslope$sigma[1:100])
    
    predicted_peakval_100[i,3:102] <- postdraws_peakval$alpha[1:100]+
      postdraws_peakval$beta[1:100,1]*flu_data_norm$PopDensity[i]+
      postdraws_peakval$beta[1:100,2]*flu_data_norm$Under18[i]+
      postdraws_peakval$beta[1:100,3]*flu_data_norm$MeanMaxTemp[i]+
      postdraws_peakval$beta[1:100,4]*flu_data_norm$MeanMaxRelHum[i]+
      postdraws_peakval$beta[1:100,5]*flu_data_norm$MeanMinRelHum[i]+
      postdraws_peakval$beta[1:100,6]*flu_data_norm$MeanMeanAbsHum[i]+
      rnorm(100, mean=0, sd=postdraws_peakval$sigma[1:100])
  }
}

# un-normalize
for (i in 1:(dim(flu_data_norm)[1])) {
  predicted_peakwk_100[i,3:102] = (predicted_peakwk_100[i,3:102]*means_and_sds[means_and_sds$Variable == 'PeakWeek','SD']) + means_and_sds[means_and_sds$Variable == 'PeakWeek','Mean']
  predicted_upslope_100[i,3:102] = (predicted_upslope_100[i,3:102]*means_and_sds[means_and_sds$Variable == 'UpSlope','SD']) + means_and_sds[means_and_sds$Variable == 'UpSlope','Mean']
  predicted_downslope_100[i,3:102] = (predicted_downslope_100[i,3:102]*means_and_sds[means_and_sds$Variable == 'DownSlope','SD']) + means_and_sds[means_and_sds$Variable == 'DownSlope','Mean']
  predicted_peakval_100[i,3:102] = (predicted_peakval_100[i,3:102]*means_and_sds[means_and_sds$Variable == 'PeakVal','SD']) + means_and_sds[means_and_sds$Variable == 'PeakVal','Mean']
}

# print most extreme values
for (i in 1:(dim(flu_data_norm)[1])) {
  #print(paste(predicted_peakwk_100$State[i],predicted_peakwk_100$Year[i]))
  #print(paste("max peak week:", max(predicted_peakwk_100[i,3:102])))
  #print(paste("min peak week:", min(predicted_peakwk_100[i,3:102])))
  #print(paste("max up slope:", max(predicted_upslope_100[i,3:102])))
  #print(paste("min up slope:", min(predicted_upslope_100[i,3:102])))
  print(paste("max down slope:", max(predicted_downslope_100[i,3:102])))
  #print(paste("min down slope:", min(predicted_downslope_100[i,3:102])))
  #print(paste("max peak value:", max(predicted_peakval_100[i,3:102])))
  #print(paste("min peak value:", min(predicted_peakval_100[i,3:102])))
}


# ==============================================================================
# Plot ILI+ data, actual fit, predicted fit for mean and first 100 draws
# ==============================================================================
for (i in 1:(dim(flu_data_1)[1])) {
  if (predicted_parameters$Year[i]!=flu_data_1$Year[i]) {
    print("Years are not equal!")
  }
  else if (predicted_parameters$State[i]!=flu_data_1$State[i]) {
    print("States are not equal!")
  }
  else {
    #plot ILI+ data
    week <- 1:52
    plot(week, LogILIplus[LogILIplus$State==flu_data_1$State[i]&LogILIplus$Year==flu_data_1$Year[i],3:54], 
         main = sprintf("Flu incidence in %s over the %i-%i season", flu_data_1$State[i], flu_data_1$Year[i], flu_data_1$Year[i]+1),
         xlab = "Week", ylab = "Log(ILI+)")
    
    #add segmented regression model
    #plot(segmented.fitN, add=T)
    wvals <- seq(from = 0, to = 52, by = 0.05)
    lines(wvals,(wvals < flu_data_1$Break1[i]) * (flu_data_1$Val1[i]+flu_data_1$Slope1[i]*wvals) + 
            (wvals >= flu_data_1$Break1[i] & wvals <= flu_data_1$Break2[i]) * (flu_data_1$Val2[i]+flu_data_1$Slope2[i]*(wvals-flu_data_1$Break1[i])) + 
            (wvals > flu_data_1$Break2[i]) * (flu_data_1$Val3[i]+flu_data_1$Slope3[i]*(wvals-flu_data_1$Break2[i])), 
          type='l',lty = 'longdash', lwd = 3, col='aquamarine3')
    for (k in 3:102) {
      lines(wvals,(wvals < predicted_peakwk_100[i,k]) * (predicted_peakval_100[i,k]-predicted_upslope_100[i,k]*(predicted_peakwk_100[i,k]-wvals)) + 
              (wvals >= predicted_peakwk_100[i,k]) * (predicted_peakval_100[i,k]+predicted_downslope_100[i,k]*(wvals-predicted_peakwk_100[i,k])), 
            type='l', lwd = 3, col='hotpink1')
    }
    lines(wvals,(wvals < predicted_parameters$PeakWeek[i]) * (predicted_parameters$PeakVal[i]-predicted_parameters$UpSlope[i]*(predicted_parameters$PeakWeek[i]-wvals)) + 
            (wvals >= predicted_parameters$PeakWeek[i]) * (predicted_parameters$PeakVal[i]+predicted_parameters$DownSlope[i]*(wvals-predicted_parameters$PeakWeek[i])), 
          type='l', lwd = 3, col='deeppink3')
  }
}

# try ggplot #

for (i in 1:(dim(flu_data_1)[1])) {
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
    for (k in 3:102) {
      pred_line_toplot[[paste("PredVal", k, sep = "")]] <- (wvals < predicted_peakwk_100[i,k]) * (predicted_peakval_100[i,k]-predicted_upslope_100[i,k]*(predicted_peakwk_100[i,k]-wvals)) + 
              (wvals >= predicted_peakwk_100[i,k]) * (predicted_peakval_100[i,k]+predicted_downslope_100[i,k]*(wvals-predicted_peakwk_100[i,k]))
    }
    
    # for trying to make a legend
    mdf <- reshape2::melt(pred_line_toplot, id.var = "Wvals")
    
    # initialize plot formatting
    ILIplus_plot <- ggplot() +
      theme_classic() +
      ggtitle(sprintf("Flu incidence in %s over the %i-%i season", flu_data_1$State[i], flu_data_1$Year[i], flu_data_1$Year[i]+1))+
      ylim(min(LogILIplus_toplot$LogILIplus)-0.5,max(LogILIplus_toplot$LogILIplus+1.5)) +
      xlab("Week") + ylab("Log(ILI+)") +
      theme(plot.title = element_text(face = 'bold', size =16, hjust = 0.5))
    
    
    # predicted lines from first 100 draws
    for (k in 3:102) {
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
           path = "~/FluForecasting/FluForecasting/plots_prediction_interval",
           width = 7, height = 5)
  }
}

# ==============================================================================
# Cross validation homebrew
# ==============================================================================
lpds <- rep(0,K-2)
for (i in 2:(K-1)) {
  peakwk_model_loo <- stan(file="FluPredictLOO.stan",data=list(N=K, M=M, x=x, y=y[,1],i=i), 
                           iter=1000, chains=4)
  params <- extract(peakwk_model_loo)
  lpds[i-1] = mean(params$lpd)
}

# estimated loo IC, scaled because my code doesn't calculate lpd for 
# leaving out first or last data point
-2*sum(lpds)*K/(K-2)


peakwk_loo_plot <- plot(peakwk_model_loo, show_density = TRUE, fill_color = "lightsteelblue1") + 
  ggtitle("Posterior distributions for peak week parameters and LPD") +
  theme(plot.title = element_text(face = 'bold', size =18, hjust = 0.6),
        axis.text=element_text(size=16))

#plot(flu_data_norm$MeanMaxTemp,flu_data_norm$MeanMinRelHum)

# ==============================================================================
# Cross validation with loo package
# ==============================================================================
peakwk_model_lpd <- stan(file="linear_flu_model_IC.stan",data=list(N=K, M=M, x=x, y=y[,1]), 
                         iter=1000, chains=4)
upslope_model_lpd <- stan(file="linear_flu_model_IC.stan",data=list(N=K, M=M, x=x, y=y[,2]), 
                         iter=1000, chains=4)
downslope_model_lpd <- stan(file="linear_flu_model_IC.stan",data=list(N=K, M=M, x=x, y=y[,3]), 
                         iter=1000, chains=4)
peakval_model_lpd <- stan(file="linear_flu_model_IC.stan",data=list(N=K, M=M, x=x, y=y[,4]), 
                          iter=1000, chains=4)

# extract log likelihood values to use in loo package
ll_peakwk <- extract_log_lik(peakwk_model_lpd)
ll_upslope <- extract_log_lik(upslope_model_lpd)
ll_downslope <- extract_log_lik(downslope_model_lpd)
ll_peakval <- extract_log_lik(peakval_model_lpd)

# calculate information criteria
waic_peakwk <- waic(ll_peakwk)
waic_upslope <- waic(ll_upslope)
waic_downslope <- waic(ll_downslope)
waic_peakval <- waic(ll_peakval)

loo_peakwk <- loo(ll_peakwk)
loo_upslope <- loo(ll_upslope)
loo_downslope <- loo(ll_downslope)
loo_peakval <- loo(ll_peakval)

loo_peakwk$estimates[3,1]+loo_upslope$estimates[3,1]+loo_downslope$estimates[3,1]+loo_peakval$estimates[3,1]

# ==============================================================================
# Plot peak week vs up slope, etc with error ellipse
# ==============================================================================
norm_pkweeks <- c(data.matrix(predicted_peakwk_100[1:K,3:102]))
norm_upslopes <- c(data.matrix(predicted_upslope_100[1:K,3:102]))
norm_downslopes <- c(data.matrix(predicted_downslope_100[1:K,3:102]))
norm_pkvals <- c(data.matrix(predicted_peakval_100[1:K,3:102]))

pkwk_vs_upslope_plot <- ggplot() +
  geom_point(aes(x=norm_pkweeks,y=norm_upslopes),alpha=0.01,color="violetred4")+
  geom_point(aes(x=flu_data_1$PeakWeek,y=flu_data_1$UpSlope),alpha=0.6)+
  stat_ellipse(aes(x=norm_pkweeks,y=norm_upslopes), color="violetred4")+
  xlim(min(flu_data_1$PeakWeek)-7,max(flu_data_1$PeakWeek)+4) +
  ylim(min(flu_data_1$UpSlope)-0.13,max(flu_data_1$UpSlope)+0.05) +
  xlab("Peak Week") + ylab("Up Slope") +
  theme_classic()

pkwk_vs_downslope_plot <- ggplot() +
  geom_point(aes(x=norm_pkweeks,y=norm_downslopes),alpha=0.01,color="violetred4")+
  geom_point(aes(x=flu_data_1$PeakWeek,y=flu_data_1$DownSlope),alpha=0.6)+
  stat_ellipse(aes(x=norm_pkweeks,y=norm_downslopes), color="violetred4")+
  xlim(min(flu_data_1$PeakWeek)-7,max(flu_data_1$PeakWeek)+4) +
  ylim(min(flu_data_1$DownSlope)-0.13,max(flu_data_1$DownSlope)+0.6) +
  xlab("Peak Week") + ylab("Down Slope") +
  theme_classic()

pkwk_vs_pkval_plot <- ggplot() +
  geom_point(aes(x=norm_pkweeks,y=norm_pkvals),alpha=0.01,color="violetred4")+
  geom_point(aes(x=flu_data_1$PeakWeek,y=flu_data_1$PeakVal),alpha=0.6)+
  stat_ellipse(aes(x=norm_pkweeks,y=norm_pkvals), color="violetred4")+
  xlim(min(flu_data_1$PeakWeek)-7,max(flu_data_1$PeakWeek)+4) +
  ylim(min(flu_data_1$PeakVal)-0.5,max(flu_data_1$PeakVal)+0.5) +
  xlab("Peak Week") + ylab("Peak Value") +
  theme_classic()

upslope_vs_downslope_plot <- ggplot() +
  geom_point(aes(x=norm_upslopes,y=norm_downslopes),alpha=0.01,color="violetred4")+
  geom_point(aes(x=flu_data_1$UpSlope,y=flu_data_1$DownSlope),alpha=0.6)+
  stat_ellipse(aes(x=norm_upslopes,y=norm_downslopes), color="violetred4")+
  xlim(min(flu_data_1$UpSlope)-0.05,max(flu_data_1$UpSlope)+0.05) +
  ylim(min(flu_data_1$DownSlope)-0.1,max(flu_data_1$DownSlope)+0.6) +
  xlab("Up Slope") + ylab("Down Slope") +
  theme_classic()

upslope_vs_pkval_plot <- ggplot() +
  geom_point(aes(x=norm_upslopes,y=norm_pkvals),alpha=0.01,color="violetred4")+
  geom_point(aes(x=flu_data_1$UpSlope,y=flu_data_1$PeakVal),alpha=0.6)+
  stat_ellipse(aes(x=norm_upslopes,y=norm_pkvals), color="violetred4")+
  xlim(min(flu_data_1$UpSlope)-0.05,max(flu_data_1$UpSlope)+0.05) +
  ylim(min(flu_data_1$PeakVal)-0.5,max(flu_data_1$PeakVal)+0.5) +
  xlab("Up Slope") + ylab("Peak Value") +
  theme_classic()

downslope_vs_pkval_plot <- ggplot() +
  geom_point(aes(x=norm_downslopes,y=norm_pkvals),alpha=0.01,color="violetred4")+
  geom_point(aes(x=flu_data_1$DownSlope,y=flu_data_1$PeakVal),alpha=0.6)+
  stat_ellipse(aes(x=norm_downslopes,y=norm_pkvals), color="violetred4")+
  xlim(min(flu_data_1$DownSlope)-0.05,max(flu_data_1$DownSlope)+0.6) +
  ylim(min(flu_data_1$PeakVal)-0.5,max(flu_data_1$PeakVal)+0.5) +
  xlab("Down Slope") + ylab("Peak Value") +
  theme_classic()

