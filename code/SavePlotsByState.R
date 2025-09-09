# group plots by state

library(png)
library(grid)
library(gridExtra)
library(ggplot2)

S = length(states)

#for (s in 1:S) {
#  valid_years <- flu_data_norm[flu_data_norm$State==states[s],1] # years where that state has data
#  print(states[s])
#  print(valid_years)
#  for (n in valid_years) { 
#    assign(paste('plot',n),readPNG(paste('~/FluForecasting/FluForecasting/plots/',states[s],n,".png",sep = "")))
#  }
#  tmp <- arrangeGrob(rasterGrob(plot1),rasterGrob(plot2),ncol=1)
#}

# Alabama
plot1 <- readPNG('~/FluForecasting/FluForecasting/plots/Alabama2010.png')
plot2 <- readPNG('~/FluForecasting/FluForecasting/plots/Alabama2011.png')
plot3 <- readPNG('~/FluForecasting/FluForecasting/plots/Alabama2013.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/Alabama2014.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Alabama2015.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Alabama2016.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Alabama2017.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Alabama2018.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Alabama2019.png')

tmp <- arrangeGrob(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),ncol=3)
ggsave('Alabama.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Arizona
plot1 <- readPNG('~/FluForecasting/FluForecasting/plots/Arizona2010.png')
plot2 <- readPNG('~/FluForecasting/FluForecasting/plots/Arizona2011.png')
plot3 <- readPNG('~/FluForecasting/FluForecasting/plots/Arizona2012.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/Arizona2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Arizona2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Arizona2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Arizona2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Arizona2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Arizona2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Arizona2019.png')

tmp <- arrangeGrob(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Arizona.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# California
plot1 <- readPNG('~/FluForecasting/FluForecasting/plots/California2010.png')
plot2 <- readPNG('~/FluForecasting/FluForecasting/plots/California2011.png')
plot3 <- readPNG('~/FluForecasting/FluForecasting/plots/California2012.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/California2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/California2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/California2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/California2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/California2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/California2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/California2019.png')

tmp <- arrangeGrob(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('California.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Colorado
plot1 <- readPNG('~/FluForecasting/FluForecasting/plots/Colorado2010.png')
plot2 <- readPNG('~/FluForecasting/FluForecasting/plots/Colorado2011.png')
plot3 <- readPNG('~/FluForecasting/FluForecasting/plots/Colorado2012.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/Colorado2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Colorado2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Colorado2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Colorado2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Colorado2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Colorado2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Colorado2019.png')

tmp <- arrangeGrob(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Colorado.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Connecticut
plot2 <- readPNG('~/FluForecasting/FluForecasting/plots/Connecticut2011.png')
plot3 <- readPNG('~/FluForecasting/FluForecasting/plots/Connecticut2012.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/Connecticut2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Connecticut2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Connecticut2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Connecticut2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Connecticut2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Connecticut2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Connecticut2019.png')

tmp <- arrangeGrob(rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Connecticut.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Delaware
plot1 <- readPNG('~/FluForecasting/FluForecasting/plots/Delaware2010.png')
plot2 <- readPNG('~/FluForecasting/FluForecasting/plots/Delaware2011.png')
plot3 <- readPNG('~/FluForecasting/FluForecasting/plots/Delaware2017.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/Delaware2018.png')

tmp <- arrangeGrob(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot3),
                   rasterGrob(plot4),ncol=2)
ggsave('Delaware.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Florida
plot1 <- readPNG('~/FluForecasting/FluForecasting/plots/Florida2010.png')
plot2 <- readPNG('~/FluForecasting/FluForecasting/plots/Florida2011.png')
plot3 <- readPNG('~/FluForecasting/FluForecasting/plots/Florida2012.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/Florida2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Florida2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Florida2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Florida2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Florida2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Florida2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Florida2019.png')

tmp <- arrangeGrob(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Florida.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Georgia
plot1 <- readPNG('~/FluForecasting/FluForecasting/plots/Georgia2010.png')
plot2 <- readPNG('~/FluForecasting/FluForecasting/plots/Georgia2011.png')
plot3 <- readPNG('~/FluForecasting/FluForecasting/plots/Georgia2012.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/Georgia2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Georgia2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Georgia2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Georgia2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Georgia2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Georgia2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Georgia2019.png')

tmp <- arrangeGrob(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Georgia.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Illinois
plot1 <- readPNG('~/FluForecasting/FluForecasting/plots/Illinois2010.png')
plot2 <- readPNG('~/FluForecasting/FluForecasting/plots/Illinois2011.png')
plot3 <- readPNG('~/FluForecasting/FluForecasting/plots/Illinois2012.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/Illinois2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Illinois2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Illinois2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Illinois2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Illinois2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Illinois2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Illinois2019.png')

tmp <- arrangeGrob(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Illinois.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Indiana
plot1 <- readPNG('~/FluForecasting/FluForecasting/plots/Indiana2010.png')
plot2 <- readPNG('~/FluForecasting/FluForecasting/plots/Indiana2011.png')
plot3 <- readPNG('~/FluForecasting/FluForecasting/plots/Indiana2012.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/Indiana2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Indiana2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Indiana2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Indiana2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Indiana2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Indiana2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Indiana2019.png')

tmp <- arrangeGrob(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Indiana.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Iowa
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Iowa2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Iowa2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Iowa2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Iowa2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Iowa2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Iowa2019.png')

tmp <- arrangeGrob(rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Iowa.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Kansas
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Kansas2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Kansas2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Kansas2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Kansas2019.png')

tmp <- arrangeGrob(rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=2)
ggsave('Kansas.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Kentucky
plot1 <- readPNG('~/FluForecasting/FluForecasting/plots/Kentucky2010.png')
plot2 <- readPNG('~/FluForecasting/FluForecasting/plots/Kentucky2011.png')
plot3 <- readPNG('~/FluForecasting/FluForecasting/plots/Kentucky2012.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/Kentucky2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Kentucky2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Kentucky2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Kentucky2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Kentucky2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Kentucky2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Kentucky2019.png')

tmp <- arrangeGrob(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Kentucky.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Louisiana
plot1 <- readPNG('~/FluForecasting/FluForecasting/plots/Louisiana2010.png')
plot2 <- readPNG('~/FluForecasting/FluForecasting/plots/Louisiana2011.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Louisiana2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Louisiana2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Louisiana2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Louisiana2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Louisiana2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Louisiana2019.png')

tmp <- arrangeGrob(rasterGrob(plot1),rasterGrob(plot2),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Louisiana.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Maryland
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/Maryland2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Maryland2014.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Maryland2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Maryland2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Maryland2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Maryland2019.png')

tmp <- arrangeGrob(rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Maryland.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Massachusetts
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/Massachusetts2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Massachusetts2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Massachusetts2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Massachusetts2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Massachusetts2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Massachusetts2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Massachusetts2019.png')

tmp <- arrangeGrob(rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Massachusetts.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Michigan
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Michigan2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Michigan2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Michigan2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Michigan2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Michigan2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Michigan2019.png')

tmp <- arrangeGrob(rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Michigan.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Minnesota
plot1 <- readPNG('~/FluForecasting/FluForecasting/plots/Minnesota2010.png')
plot2 <- readPNG('~/FluForecasting/FluForecasting/plots/Minnesota2011.png')
plot3 <- readPNG('~/FluForecasting/FluForecasting/plots/Minnesota2012.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/Minnesota2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Minnesota2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Minnesota2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Minnesota2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Minnesota2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Minnesota2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Minnesota2019.png')

tmp <- arrangeGrob(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Minnesota.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Missouri
plot1 <- readPNG('~/FluForecasting/FluForecasting/plots/Missouri2010.png')
plot2 <- readPNG('~/FluForecasting/FluForecasting/plots/Missouri2011.png')
plot3 <- readPNG('~/FluForecasting/FluForecasting/plots/Missouri2012.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/Missouri2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Missouri2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Missouri2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Missouri2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Missouri2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Missouri2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Missouri2019.png')

tmp <- arrangeGrob(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Missouri.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Montana
plot3 <- readPNG('~/FluForecasting/FluForecasting/plots/Montana2012.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/Montana2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Montana2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Montana2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Montana2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Montana2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Montana2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Montana2019.png')

tmp <- arrangeGrob(rasterGrob(plot3),rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Montana.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Nebraska
plot3 <- readPNG('~/FluForecasting/FluForecasting/plots/Nebraska2012.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/Nebraska2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Nebraska2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Nebraska2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Nebraska2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Nebraska2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Nebraska2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Nebraska2019.png')

tmp <- arrangeGrob(rasterGrob(plot3),rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Nebraska.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# New Mexico
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/New Mexico2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/New Mexico2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/New Mexico2017.png')

tmp <- arrangeGrob(rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),ncol=2)
ggsave('New Mexico.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# New York
plot1 <- readPNG('~/FluForecasting/FluForecasting/plots/New York2010.png')
plot2 <- readPNG('~/FluForecasting/FluForecasting/plots/New York2011.png')
plot3 <- readPNG('~/FluForecasting/FluForecasting/plots/New York2012.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/New York2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/New York2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/New York2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/New York2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/New York2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/New York2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/New York2019.png')

tmp <- arrangeGrob(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('New York.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# North Carolina
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/North Carolina2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/North Carolina2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/North Carolina2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/North Carolina2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/North Carolina2019.png')

tmp <- arrangeGrob(rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('North Carolina.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# North Dakota
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/North Dakota2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/North Dakota2019.png')

tmp <- arrangeGrob(rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('North Dakota.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Ohio
plot1 <- readPNG('~/FluForecasting/FluForecasting/plots/Ohio2010.png')
plot2 <- readPNG('~/FluForecasting/FluForecasting/plots/Ohio2011.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/Ohio2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Ohio2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Ohio2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Ohio2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Ohio2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Ohio2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Ohio2019.png')

tmp <- arrangeGrob(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Ohio.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Oklahoma
plot1 <- readPNG('~/FluForecasting/FluForecasting/plots/Oklahoma2010.png')
plot3 <- readPNG('~/FluForecasting/FluForecasting/plots/Oklahoma2012.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/Oklahoma2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Oklahoma2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Oklahoma2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Oklahoma2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Oklahoma2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Oklahoma2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Oklahoma2019.png')

tmp <- arrangeGrob(rasterGrob(plot1),rasterGrob(plot3),rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Oklahoma.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Oregon
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Oregon2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Oregon2015.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Oregon2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Oregon2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Oregon2019.png')

tmp <- arrangeGrob(rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Oregon.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Pennsylvania
plot1 <- readPNG('~/FluForecasting/FluForecasting/plots/Pennsylvania2010.png')
plot2 <- readPNG('~/FluForecasting/FluForecasting/plots/Pennsylvania2011.png')
plot3 <- readPNG('~/FluForecasting/FluForecasting/plots/Pennsylvania2012.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/Pennsylvania2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Pennsylvania2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Pennsylvania2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Pennsylvania2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Pennsylvania2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Pennsylvania2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Pennsylvania2019.png')

tmp <- arrangeGrob(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Pennsylvania.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# South Carolina
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/South Carolina2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/South Carolina2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/South Carolina2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/South Carolina2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/South Carolina2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/South Carolina2018.png')

tmp <- arrangeGrob(rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),ncol=3)
ggsave('South Carolina.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# South Dakota
plot2 <- readPNG('~/FluForecasting/FluForecasting/plots/South Dakota2011.png')
plot3 <- readPNG('~/FluForecasting/FluForecasting/plots/South Dakota2012.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/South Dakota2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/South Dakota2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/South Dakota2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/South Dakota2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/South Dakota2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/South Dakota2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/South Dakota2019.png')

tmp <- arrangeGrob(rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('South Dakota.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Tennessee
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/Tennessee2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Tennessee2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Tennessee2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Tennessee2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Tennessee2017.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Tennessee2019.png')

tmp <- arrangeGrob(rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot10),ncol=3)
ggsave('Tennessee.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Texas
plot1 <- readPNG('~/FluForecasting/FluForecasting/plots/Texas2010.png')
plot2 <- readPNG('~/FluForecasting/FluForecasting/plots/Texas2011.png')
plot3 <- readPNG('~/FluForecasting/FluForecasting/plots/Texas2012.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/Texas2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Texas2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Texas2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Texas2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Texas2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Texas2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Texas2019.png')

tmp <- arrangeGrob(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Texas.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Texas
plot1 <- readPNG('~/FluForecasting/FluForecasting/plots/Texas2010.png')
plot2 <- readPNG('~/FluForecasting/FluForecasting/plots/Texas2011.png')
plot3 <- readPNG('~/FluForecasting/FluForecasting/plots/Texas2012.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/Texas2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Texas2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Texas2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Texas2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Texas2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Texas2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Texas2019.png')

tmp <- arrangeGrob(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Texas.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Utah
plot1 <- readPNG('~/FluForecasting/FluForecasting/plots/Utah2010.png')
plot2 <- readPNG('~/FluForecasting/FluForecasting/plots/Utah2011.png')
plot3 <- readPNG('~/FluForecasting/FluForecasting/plots/Utah2012.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/Utah2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Utah2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Utah2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Utah2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Utah2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Utah2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Utah2019.png')

tmp <- arrangeGrob(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Utah.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Vermont
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Vermont2014.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Vermont2016.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Vermont2018.png')

tmp <- arrangeGrob(rasterGrob(plot5),rasterGrob(plot7),
                   rasterGrob(plot9),ncol=2)
ggsave('Vermont.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Virginia
plot1 <- readPNG('~/FluForecasting/FluForecasting/plots/Virginia2010.png')
plot2 <- readPNG('~/FluForecasting/FluForecasting/plots/Virginia2011.png')
plot3 <- readPNG('~/FluForecasting/FluForecasting/plots/Virginia2012.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/Virginia2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Virginia2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Virginia2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Virginia2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Virginia2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Virginia2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Virginia2019.png')

tmp <- arrangeGrob(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Virginia.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Washington
plot1 <- readPNG('~/FluForecasting/FluForecasting/plots/Washington2010.png')
plot2 <- readPNG('~/FluForecasting/FluForecasting/plots/Washington2011.png')
plot3 <- readPNG('~/FluForecasting/FluForecasting/plots/Washington2012.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/Washington2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Washington2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Washington2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Washington2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Washington2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Washington2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Washington2019.png')

tmp <- arrangeGrob(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Washington.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# West Virginia
plot1 <- readPNG('~/FluForecasting/FluForecasting/plots/West Virginia2010.png')
plot2 <- readPNG('~/FluForecasting/FluForecasting/plots/West Virginia2011.png')
plot3 <- readPNG('~/FluForecasting/FluForecasting/plots/West Virginia2012.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/West Virginia2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/West Virginia2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/West Virginia2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/West Virginia2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/West Virginia2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/West Virginia2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/West Virginia2019.png')

tmp <- arrangeGrob(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('West Virginia.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")

# Wisconsin
plot1 <- readPNG('~/FluForecasting/FluForecasting/plots/Wisconsin2010.png')
plot3 <- readPNG('~/FluForecasting/FluForecasting/plots/Wisconsin2012.png')
plot4 <- readPNG('~/FluForecasting/FluForecasting/plots/Wisconsin2013.png')
plot5 <- readPNG('~/FluForecasting/FluForecasting/plots/Wisconsin2014.png')
plot6 <- readPNG('~/FluForecasting/FluForecasting/plots/Wisconsin2015.png')
plot7 <- readPNG('~/FluForecasting/FluForecasting/plots/Wisconsin2016.png')
plot8 <- readPNG('~/FluForecasting/FluForecasting/plots/Wisconsin2017.png')
plot9 <- readPNG('~/FluForecasting/FluForecasting/plots/Wisconsin2018.png')
plot10 <- readPNG('~/FluForecasting/FluForecasting/plots/Wisconsin2019.png')

tmp <- arrangeGrob(rasterGrob(plot1),rasterGrob(plot3),rasterGrob(plot4),
                   rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),rasterGrob(plot8),
                   rasterGrob(plot9),rasterGrob(plot10),ncol=3)
ggsave('Wisconsin.png',tmp,width=7,height=5,path = "~/FluForecasting/FluForecasting/plots")


