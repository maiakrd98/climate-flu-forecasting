# wtf is going on with temp and latitude???

# ==============================================================================
# Import libraries
# ==============================================================================

library(tidyverse) 
library(purrr)
library(rstan)
#install.packages("plotly")
library(plotly)
library(reshape2)

# ==============================================================================
# Import and prepare flu data
# ==============================================================================
flu_data <- read.csv("flu_data.csv")
LogILIplus <- read.csv("log_ILIplus.csv")

states = unique(flu_data$State)
years = unique(flu_data$Year)
S = length(states) # number of states
N = length(years)
M = 2 # number of covariates

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

# covariates
x1 <- flu_data_norm$Lat
x2 <- flu_data_norm$MeanMaxTemp

# normalize the things we are predicting (there's definitely a real name for that...)
y <- array(dim=c(K,4))
y[,1] <- flu_data_norm$PeakWeek
y[,2] <- flu_data_norm$UpSlope
y[,3] <- flu_data_norm$DownSlope
y[,4] <- flu_data_norm$PeakVal

# ==============================================================================
# Run linear models in stan 
# ==============================================================================
peakwk_fit <- stan(file="flu_predict.stan", 
                   data=list(N = K, x1 = x1, x2=x2, y = y[,1]),
                   iter=1000, chains=4)
upslope_fit <- stan(file="flu_predict.stan", 
                    data=list(N = K, x1 = x1, x2=x2, y = y[,2]),
                    iter=1000, chains=4)
downslope_fit <- stan(file="flu_predict.stan", 
                      data=list(N = K, x1 = x1, x2=x2, y = y[,3]),
                      iter=1000, chains=4)
peakval_fit <- stan(file="flu_predict.stan", 
                    data=list(N = K, x1 = x1, x2=x2, y = y[,4]),
                    iter=1000, chains=4)
plot(peakwk_fit, show_density = TRUE)
plot(upslope_fit, show_density = TRUE)
plot(downslope_fit, show_density = TRUE)
plot(peakval_fit, show_density = TRUE)

print(peakwk_fit)
print(upslope_fit)
print(downslope_fit)
print(peakval_fit)

traceplot(peakwk_fit, inc_warmup=TRUE)
traceplot(upslope_fit, inc_warmup=TRUE)
traceplot(downslope_fit, inc_warmup=TRUE)
traceplot(peakval_fit, inc_warmup=TRUE)

# Extract the posterior draws for each parameter:
postdraws_peakwk <- extract(peakwk_fit)
postdraws_upslope <- extract(upslope_fit)
postdraws_downslope <- extract(downslope_fit)
postdraws_peakval <- extract(peakval_fit)


# ==============================================================================
# Plot histograms of parameter estimates
# ==============================================================================

postdraws_peakwk_t <- tibble(
  alpha = postdraws_peakwk$alpha,
  b1 = postdraws_peakwk$b1,      
  b2 = postdraws_peakwk$b2,
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
  geom_line(data=tibble(x=latvals, y=mean(postdraws_peakwk$alpha) + mean(postdraws_peakwk$b1)*latvals), aes(x=x, y=y)) +
  theme_classic() 
flu_data_norm %>% 
  ggplot(aes(x=Lat, y=UpSlope)) + 
  geom_point(alpha=0.6) + 
  geom_line(data=tibble(x=latvals, y=mean(postdraws_peakwk$alpha) + mean(postdraws_peakwk$b1)*latvals), aes(x=x, y=y)) +
  theme_classic() 

tempvals<- seq(from=-3, to=3, by=0.1)

flu_data_norm %>% 
  ggplot(aes(x=MeanMaxTemp, y=PeakWeek)) + 
  geom_point(alpha=0.6) + 
  geom_line(data=tibble(x=tempvals, y=mean(postdraws_peakwk$alpha) + mean(postdraws_peakwk$b2)*tempvals), aes(x=x, y=y)) +
  theme_classic() 

plot_ly(flu_data_norm, x = ~MeanMaxTemp, y = ~Lat, z = ~PeakWeek, type = "scatter3d", mode = "markers",
        marker = list(size = 3))

# 3d plot for peak week
temp_lat_surface <- expand.grid(MeanMaxTemp = tempvals,Lat = latvals,KEEP.OUT.ATTRS = F)
temp_lat_surface$PeakWeek <- mean(postdraws_peakwk$alpha)+mean(postdraws_peakwk$b1)*temp_lat_surface$Lat+mean(postdraws_peakwk$b2)*temp_lat_surface$MeanMaxTemp
temp_lat_surface <- acast(temp_lat_surface, Lat ~ MeanMaxTemp, value.var = "PeakWeek") #y ~ x

temp_lat_plot <- plot_ly(flu_data_norm, x = ~MeanMaxTemp, y = ~Lat, z = ~PeakWeek, type = "scatter3d", mode = "markers",
                         marker = list(size = 3))
temp_lat_plot <- add_trace(p = temp_lat_plot,
                           z = temp_lat_surface,
                           x = tempvals,
                           y = latvals,
                           type = "surface")

# 3d plot for up slope
temp_lat_surface <- expand.grid(MeanMaxTemp = tempvals,Lat = latvals,KEEP.OUT.ATTRS = F)
temp_lat_surface$UpSlope <- mean(postdraws_upslope$alpha)+mean(postdraws_upslope$b1)*temp_lat_surface$Lat+mean(postdraws_upslope$b2)*temp_lat_surface$MeanMaxTemp
temp_lat_surface <- acast(temp_lat_surface, Lat ~ MeanMaxTemp, value.var = "UpSlope") #y ~ x

temp_lat_plot <- plot_ly(flu_data_norm, x = ~MeanMaxTemp, y = ~Lat, z = ~UpSlope, type = "scatter3d", mode = "markers",
                         marker = list(size = 3))
temp_lat_plot <- add_trace(p = temp_lat_plot,
                           z = temp_lat_surface,
                           x = tempvals,
                           y = latvals,
                           type = "surface")

# ==============================================================================
# Fit peak week to temp within a state
# ==============================================================================
for (s in 1:S) {
  xs <- flu_data_norm[flu_data_norm$State==states[s],"MeanMaxTemp"]
  ys <- flu_data_norm[flu_data_norm$State==states[s],"PeakWeek"]
  
  peakwk_fit_s <- stan(file="temp_lat_flu.stan", 
                     data=list(N = length(xs), x = xs,y = ys),
                     iter=1000, chains=4)
  print(states[s])
  print(peakwk_fit_s)
  params <- extract(peakwk_fit_s)
  #plot(xs,ys)
  #abline(mean(params$alpha),mean(params$b1))
}



# ==============================================================================
# Try regular lm 
# ==============================================================================

peakwk_lm <- lm(PeakWeek ~ Lat + MeanMaxTemp, data = flu_data_norm)
latvals <- seq(from=25, to=50, by=0.1)
tempvals <- seq(from=8, to=30, by=0.1)

temp_lat_surface <- expand.grid(Lat = latvals,MeanMaxTemp = tempvals,KEEP.OUT.ATTRS = F)
temp_lat_surface$PeakWeek <- predict.lm(peakwk_lm, newdata = temp_lat_surface)
temp_lat_surface <- acast(temp_lat_surface, MeanMaxTemp ~ Lat, value.var = "PeakWeek") #y ~ x

temp_lat_plot <- plot_ly(flu_data_1, x = ~Lat, y = ~MeanMaxTemp, z = ~PeakWeek, type = "scatter3d", mode = "markers",
                         marker = list(size = 3))
temp_lat_plot <- add_trace(p = temp_lat_plot,
                           z = temp_lat_surface,
                           x = latvals,
                           y = tempvals,
                           type = "surface")


# ==============================================================================
# Try iris data to see if there is an issue with plotting
# ==============================================================================

my_df <- iris
petal_lm <- lm(Petal.Length ~ Sepal.Length + Sepal.Width,data = my_df)

axis_x <- seq(min(my_df$Sepal.Length), max(my_df$Sepal.Length), by = 0.1)
axis_y <- seq(min(my_df$Sepal.Width), max(my_df$Sepal.Width), by = 0.1)

petal_surface <- expand.grid(Sepal.Length = axis_x,Sepal.Width = axis_y,KEEP.OUT.ATTRS = F)
petal_surface$Petal.Length <- predict.lm(petal_lm, newdata = petal_surface)
petal_surface <- acast(petal_surface, Sepal.Width ~ Sepal.Length, value.var = "Petal.Length") #y ~ x

petal_plot <- plot_ly(my_df, x = ~Sepal.Length, y = ~Sepal.Width, z = ~Petal.Length, type = "scatter3d", mode = "markers",
                         marker = list(size = 3))
petal_plot <- add_trace(p = petal_plot,
                           z = petal_surface,
                           x = axis_x,
                           y = axis_y,
                           type = "surface")


# ==============================================================================
# What if we do it in Stan, but don't normalize the data... still bad 
# ==============================================================================

# covariates
x1 <- flu_data_1$Lat
x2 <- flu_data_1$MeanMaxTemp

# the things we are predicting (there's definitely a real name for that...)
y <- array(dim=c(K,4))
y[,1] <- flu_data_1$PeakWeek
y[,2] <- flu_data_1$UpSlope
y[,3] <- flu_data_1$DownSlope
y[,4] <- flu_data_1$PeakVal

peakwk_fit <- stan(file="flu_predict.stan", 
                   data=list(N = K, x1 = x1, x2=x2, y = y[,1]),
                   iter=2000, chains=4)
upslope_fit <- stan(file="flu_predict.stan", 
                    data=list(N = K, x1 = x1, x2=x2, y = y[,2]),
                    iter=1000, chains=4)
downslope_fit <- stan(file="flu_predict.stan", 
                      data=list(N = K, x1 = x1, x2=x2, y = y[,3]),
                      iter=1000, chains=4)
peakval_fit <- stan(file="flu_predict.stan", 
                    data=list(N = K, x1 = x1, x2=x2, y = y[,4]),
                    iter=1000, chains=4)
plot(peakwk_fit)
plot(upslope_fit)
plot(downslope_fit)
plot(peakval_fit)

print(peakwk_fit)
print(upslope_fit)
print(downslope_fit)
print(peakval_fit)



# ==============================================================================
# Generate synthetic data and try to fit it 
# ==============================================================================
# true parameters that we will try to predict 
true_b1 <- 0.7
true_b2 <- 0.6
true_alpha <- 0 
true_sigma <- 1

# generate data
num = 300 # number of data points
synth_data <- tibble(lats = rnorm(num, 0, 1),
                     temps = rnorm(num, 0, 1), # these are NOT correlated for now, but they should be
                     pkwks = rep(0, num))

for (i in 1:num) {
  synth_data$pkwks[i] <- rnorm(1, true_alpha + true_b1*synth_data$lats[i] + true_b2*synth_data$temps[i], true_sigma)
}

# plot scatterplots of synthetic and actual data
pairs(synth_data)
pairs(flu_data_norm[c(17,18,11)])

# fit stan model
peakwk_fit_test <- stan(file="flu_predict.stan", 
                   data=list(N = num, x1 = synth_data$lats, x2 = synth_data$temps, y = synth_data$pkwks),
                   iter=1000, chains=4)

plot(peakwk_fit_test, show_density = TRUE)
print(peakwk_fit_test)
traceplot(peakwk_fit_test, inc_warmup = TRUE)
postdraws_peakwk <- extract(peakwk_fit_test)


# plot 3d scatter plot and predicted plane
latvals <- seq(from=-3, to=3, by=0.1)
tempvals <- seq(from=-3, to=3, by=0.1)

temp_lat_test_surface <- expand.grid(MeanMaxTemp = tempvals,Lat = latvals,KEEP.OUT.ATTRS = F)
temp_lat_test_surface$PeakWeek <- mean(postdraws_peakwk$alpha)+mean(postdraws_peakwk$b1)*temp_lat_test_surface$Lat+mean(postdraws_peakwk$b2)*temp_lat_test_surface$MeanMaxTemp
temp_lat_test_surface <- acast(temp_lat_test_surface, Lat ~ MeanMaxTemp, value.var = "PeakWeek") #y ~ x

temp_lat_test_plot <- plot_ly(synth_data, x = ~temps, y = ~lats, z = ~pkwks, type = "scatter3d", mode = "markers",
                         marker = list(size = 3))
temp_lat_test_plot <- add_trace(p = temp_lat_test_plot,
                                z = temp_lat_test_surface,
                                x = tempvals,
                                y = latvals,
                                type = "surface")

# ==============================================================================
# Generate correlated synthetic data and try to fit it 
# ==============================================================================
# true parameters that we will try to predict 
true_b1 <- 0.7
true_b2 <- 0.6
true_alpha <- 0 
true_sigma <- 1

corr <- -0.94 # correlation b/w lat and temp

# generate data
num = 300 # number of data points
synth_data <- tibble(lats = rnorm(num, 0, 1),
                     temps = rep(0, num),
                     pkwks = rep(0, num))

for (i in 1:num) {
  synth_data$temps[i] <- rnorm(1, corr*synth_data$lats[i], 0.5)
  synth_data$pkwks[i] <- rnorm(1, true_alpha + true_b1*synth_data$lats[i] + true_b2*synth_data$temps[i], true_sigma)
}

# plot scatterplots of synthetic and actual data
pairs(synth_data)
pairs(flu_data_norm[c(17,18,11)])

# fit stan model
peakwk_fit_test <- stan(file="flu_predict.stan", 
                        data=list(N = num, x1 = synth_data$lats, x2 = synth_data$temps, y = synth_data$pkwks),
                        iter=1000, chains=4)

plot(peakwk_fit_test, show_density = TRUE)
print(peakwk_fit_test)
traceplot(peakwk_fit_test, inc_warmup = TRUE)
postdraws_peakwk <- extract(peakwk_fit_test)


# plot 3d scatter plot and predicted plane
latvals <- seq(from=-3, to=3, by=0.1)
tempvals <- seq(from=-3, to=3, by=0.1)

temp_lat_test_surface <- expand.grid(MeanMaxTemp = tempvals,Lat = latvals,KEEP.OUT.ATTRS = F)
temp_lat_test_surface$PeakWeek <- mean(postdraws_peakwk$alpha)+mean(postdraws_peakwk$b1)*temp_lat_test_surface$Lat+mean(postdraws_peakwk$b2)*temp_lat_test_surface$MeanMaxTemp
temp_lat_test_surface <- acast(temp_lat_test_surface, Lat ~ MeanMaxTemp, value.var = "PeakWeek") #y ~ x

temp_lat_test_plot <- plot_ly(synth_data, x = ~temps, y = ~lats, z = ~pkwks, type = "scatter3d", mode = "markers",
                              marker = list(size = 3))
temp_lat_test_plot <- add_trace(p = temp_lat_test_plot,
                                z = temp_lat_test_surface,
                                x = tempvals,
                                y = latvals,
                                type = "surface")


# ==============================================================================
# Fit models using just temp and just latitude
# ==============================================================================

x1 <- flu_data_norm$Lat
x2 <- flu_data_norm$MeanMaxTemp
y <- array(dim=c(K,4))
y[,1] <- flu_data_norm$PeakWeek
y[,2] <- flu_data_norm$UpSlope
y[,3] <- flu_data_norm$DownSlope
y[,4] <- flu_data_norm$PeakVal

lat_peakwk_fit <- stan(file="temp_lat_flu.stan", 
                       data=list(N = K, x = x1, y = y[,1]),
                       iter=1000, chains=4)
temp_peakwk_fit <- stan(file="temp_lat_flu.stan", 
                       data=list(N = K, x = x2, y = y[,1]),
                       iter=1000, chains=4)

plot(lat_peakwk_fit, show_density = TRUE)
plot(temp_peakwk_fit, show_density = TRUE)

print(lat_peakwk_fit)
print(temp_peakwk_fit)

traceplot(lat_peakwk_fit, inc_warmup=TRUE)
traceplot(temp_peakwk_fit, inc_warmup=TRUE)

# ==============================================================================
# Fit simulated data using just temp and just latitude
# ==============================================================================

lat_peakwk_test_fit <- stan(file="temp_lat_flu.stan",
                            data=list(N = num, x = synth_data$lats, y = synth_data$pkwks),
                            iter=1000, chains=4)
temp_peakwk_test_fit <- stan(file="temp_lat_flu.stan", 
                             data=list(N = num, x = synth_data$temps, y = synth_data$pkwks),
                             iter=1000, chains=4)

plot(lat_peakwk_test_fit, show_density = TRUE)
plot(temp_peakwk_test_fit, show_density = TRUE)

print(lat_peakwk_test_fit)
print(temp_peakwk_test_fit)

traceplot(lat_peakwk_test_fit, inc_warmup=F)
traceplot(temp_peakwk_test_fit, inc_warmup=F)



