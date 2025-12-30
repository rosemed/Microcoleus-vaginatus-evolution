library(ggplot2)
library(ggpubr)
library(rgdal)
library(plotly)
library(ggsci)
library(BIOMEplot) #install_github("kunstler/BIOMEplot")
library(devtools)
library(RColorBrewer)
library(dplyr)

site<-read.table("location.txt",row.names=1,header = T,sep="\t")
site <- site %>% relocate(Lon, .after = Lat)
unique_site<-unique.data.frame(site)
names(unique_site)<-c("lat","lon")
unique_site<-na.omit(unique_site)

library(raster)
data<-raster("database/worldclim/wc2.1_10m_bio_1.tif")
data_df<- as.data.frame(rasterToPoints(data))
colnames(data_df)

wmap<-readOGR(dsn="database/natural_earth_vector/10m_physical/ne_10m_land.shp",layer="ne_10m_land")
countries<-readOGR("database/natural_earth_vector/10m_cultural/ne_10m_admin_0_countries.shp")
wmap_proj<-spTransform(wmap,CRS("+proj=longlat +datum=WGS84"))
countries_proj<-spTransform(countries,CRS("+proj=longlat +datum=WGS84"))

coordinates(unique_site)<-c("lon","lat")
proj4string(unique_site)<-CRS("+proj=longlat +datum=WGS84")
sites_proj<-spTransform(unique_site,CRS("+proj=longlat +datum=WGS84"))

wmap_proj_df<-fortify(wmap_proj)
countries_proj_df<-fortify(countries_proj)
sites_proj_df<-data.frame(sites_proj)

a=data_df$wc2.1_30s_bio_12
a=cut(a, breaks = c(0,100,300,500,1000,2000,5000,10000,15000),include.lowest = TRUE, labels = c("0-100", "100-300", "300-500", "500-1000","1000-2000","2000-5000","5000-10000",">10000"))
data_df=data.frame(data_df,group=a)
color_gradient <- rev(gray.colors(n = 8, start = 0.3, end = 0.8))
sum(is.na(data_df))

ggplot()+
  geom_polygon(data = wmap_proj_df,aes(long,lat,group=group),col="black",fill=NA, linewidth=0.2)+
  geom_tile(data =data_df, aes(x = x, y = y, fill = group)) +
  scale_fill_manual(values =color_gradient) +
  new_scale_fill() + 
  geom_scatterpie(data = data_pie,
                  aes(x=lon, y=lat, r=log2(Size+1)*2),cols = colnames(data_pie[, 4:11]),
                  alpha=0.6,color = "black", linewidth = 0.3) +
  scale_fill_manual(values = color,name='phylogroup') + 
  geom_scatterpie_legend(r=log2(data_pie$Size+1)*2, x= -150, y= -30,n=3,
                         breaks = c(3,6,10.37829),labeller = function(x) round(2^(x/3)-1))+
  coord_fixed()+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  theme_bw()+
  labs(x=NULL,y=NULL,col=NULL)+
  theme_void()
