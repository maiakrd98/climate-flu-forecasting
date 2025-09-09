# First, we import and plot week vs log of flu incidence at the national level 
# for all years available

# percent visits due to ILI
ILI <- read.csv("ILINetNat.csv", header=TRUE, skip = 1) 
# percent of flu tests positive post 2015
PercentPos1 <- read.csv("PercentPosNat(2016+).csv", header=TRUE, skip = 1) 
# percent of flu tests positive pre 2015
PercentPos2 <- read.csv("PercentPosNat(pre2016).csv", header=TRUE, skip = 1) 

# Note: we will call the 2014-2015 season the 2014 season, etc. 
# That is, the n season lasts from week 40 in year n to week 39 in year n+1

# plot data by year
for (n in 2010:2023) {
  ILIn <- ILI[(ILI$YEAR==n & ILI$WEEK>=40) | (ILI$YEAR==n+1) & ILI$WEEK<40,]
  if (n<2015) {
    PercentPosn <- PercentPos2[(PercentPos2$YEAR==n & PercentPos2$WEEK>=40) | (PercentPos2$YEAR==n+1) & PercentPos2$WEEK<40,]
  }
  else{
    PercentPosn <- PercentPos1[(PercentPos1$YEAR==n & PercentPos1$WEEK>=40) | (PercentPos1$YEAR==n+1) & PercentPos1$WEEK<40,]
  }
  IncidenceN <- ILIn$X.UNWEIGHTED.ILI*PercentPosn$PERCENT.POSITIVE
  plot(1:nrow(ILIn), log(IncidenceN) ,main = sprintf("Flu incidence over the %i-%i season", n,n+1))
}