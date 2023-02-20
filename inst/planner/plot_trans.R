library(ggplot2)
library(mapdata)

#tr = read.csv('transects.csv') #->tr = path to csv file containing transects
#str = read.csv('HERAS_Strata_2020.csv') #-> path to csv file containing transect descriptions

plot_trans <- function(tr,str,xlims=c(-8,12),ylims=c(52,63)){
  ggplot()+
    borders(database='worldHires', fill="darkgrey",colour="black",xlim=xlims,ylim=ylims) +
    geom_polygon(data=str, aes(x=Lon, y=Lat, group=Stratum), fill='orange', col='blue',lwd=1,alpha=0.2)+
    geom_path(data=tr, aes(x=X, y=Y, group=Stratum), lwd=1.2)+
    xlab(expression(Longitude~(degree~E)))+
    ylab(expression(Latitude~(degree~N)))+
    coord_quickmap(xlim = xlims,ylim=ylims)+
    theme_classic()+
    theme(text=element_text(size=16))
}
