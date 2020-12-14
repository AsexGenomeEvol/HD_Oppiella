#map of mite sampling

#install.packages("OpenStreetMap")
library(ggmap)
library(ggplot2)
?get_stamenmap


#lat <- c(47,55)
#long <- c(5,15)

lat <- c(45,57)
long <- c(3,17)

bbox <- make_bbox(long,lat,f=0.05)
#bbox <- c(left = 5, bottom = 47, right = 15, top = 55)
#a <- get_stamenmap(bbox,maptype="toner-hybrid",source="stamen", zoom=8)
#b <- get_stamenmap(bbox,maptype="toner-hybrid",source="stamen")
b <- get_stamenmap(bbox,maptype="terrain-background",source="stamen", zoom=9)
ggmap(b)


df <- read.csv("~/Dropbox/AA_projects/mites_ASD/sampling.csv")

ggmap(b)+
  geom_point(aes(x = lon, y = lat, colour=location),
  data = df, size=9, alpha = .8) +
  labs(x = "Longitude", y = "Latitude", colour="sampling location") +
  scale_colour_manual(values=c("Göttingen forest"="#AAA9AD", "Hainich forest"="#EE4000", "Kranichstein forest"="#4F94CD", "Schwäbische Alb"="#73880A")) +
  theme(axis.title = element_text(family = "Arial", color="#666666", face="plain", size=16), axis.text.x = element_text(family = "Arial", color="#666666", face="plain", size=14), axis.text.y = element_text(family = "Arial", color="#666666", face="plain", size=14)) +
  theme(panel.border = element_rect(linetype = "solid", colour = "grey", fill = NA), panel.grid.major = element_line(color = "grey", linetype = "dotted"), panel.grid.minor = element_line(colour = "grey", linetype = "dotted"), panel.background = element_blank(), axis.line = element_line(colour = "grey40")) +
  theme(legend.text= element_text(family = "Arial", color="#666666", face="plain", size=12), legend.title= element_text(family = "Arial", color="#666666", face="plain", size=14))



ggmap(b) + geom_point(data = df, aes(lon,lat, color=location),size=10,alpha=0.7) +
          labs(x = "Longitude", y = "Latitude", +
          #color = "Sampling locations") +
          scale_colour_manual(values=c("Göttingen forest"="AAA9AD", "Hainich forest"="EE4000", "Kranichstein forest"="4F94CD", "Schwäbische Alb"="73880A"))) +
          theme(axis.title = element_text(family = "Arial", color="#666666", face="plain", size=18), axis.text.x = element_text(family = "Arial", color="#666666", face="plain", size=14), axis.text.y = element_text(family = "Arial", color="#666666", face="plain", size=14)) +
          theme(panel.border = element_rect(linetype = "solid", colour = "grey", fill = NA), panel.grid.major = element_line(color = "grey", linetype = "dotted"), panel.grid.minor = element_line(colour = "grey", linetype = "dotted"), panel.background = element_blank(), axis.line = element_line(colour = "grey40"))

