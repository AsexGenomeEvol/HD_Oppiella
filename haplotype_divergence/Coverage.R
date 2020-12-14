#script for coverage histograms

library(ggplot2)
library(plyr)

#load data

#new
#bowtie2 local
genomehist <- read.table (file="~/Desktop/ME_coveragetest/new2/On_conmask_phreg.cov.rd", sep=" ", header=F)
list281 <- read.table (file="~/Desktop/ME_coveragetest/new2/On_phasedcontigs_281.header", sep=" ", header=F)
list69 <- read.table (file="~/Desktop/ME_coveragetest/new2/On_phasedcontigs_69.header", sep=" ", header=F)

buscos <- read.table (file="~/Desktop/ME_coveragetest/On_buscos_conv.cov.mn2.rd2", sep=" ", header=F)
busdupl  <- read.table (file="~/Desktop/ME_coveragetest/On_buscos_conv.cov.mn2.rd2.dupl", sep=" ", header=F)
buscomp  <- read.table (file="~/Desktop/ME_coveragetest/On_buscos_conv.cov.mn2.rd2.comp", sep=" ", header=F)

#genomebusco <- rbind(genomehist, buscos)
genomebusco <- rbind(genomehist, buscomp, busdupl)


#rename columns
colnames(genomehist) <- c("contig","mean_cov","std","region")
colnames(genomebusco) <- c("contig","mean_cov","std","region")
colnames(list281) <- c("contig")
colnames(list69) <- c("contig")
#colnames(busco_dupl) <- c("contig")
#colnames(busco_compl) <- c("contig")


#genomehist$mean_cov <- as.numeric(genomehist$mean_cov)
#genomehist$std <- as.numeric(genomehist$std)
#genomehist$region <- as.factor(genomehist$region)
#genomehist$contig <- as.character(genomehist$contig)

genomebusco$mean_cov <- as.numeric(genomebusco$mean_cov)
genomebusco$std <- as.numeric(genomebusco$std)
genomebusco$region <- as.factor(genomebusco$region)
genomebusco$contig <- as.character(genomebusco$contig)


#str(genomehist)

#subsample
genomehist281 <- subset(genomehist, contig %in% list281$contig)
genomehist69 <- subset(genomehist, contig %in% list69$contig)

#number of regions
count(genomebusco, "region")


#mean and median
meangb <- ddply(genomebusco, "region", summarise, grp.mean=mean(mean_cov))
head(meangb)
mediangb <- ddply(genomebusco, "region", summarise, grp.median=median(mean_cov))
head(mediangb)

mu69 <- ddply(genomehist69, "region", summarise, grp.mean=mean(mean_cov))
head(mu69)
mu69 <- ddply(genomehist69, "region", summarise, grp.median=median(mean_cov))
head(mu69)
mu281 <- ddply(genomehist281, "region", summarise, grp.mean=mean(mean_cov))
head(mu281)
mu281 <- ddply(genomehist281, "region", summarise, grp.median=median(mean_cov))
head(mu281)



ggplot(genomebusco, aes(x=mean_cov, fill=region)) +
  geom_histogram(position="identity", alpha=0.5, binwidth = 10)+
  geom_vline(data=mediangb, aes(xintercept=grp.median, color=region), linetype="dashed")+
  xlim(0,1000) +
  scale_fill_brewer(palette="RdBu") +
  theme(axis.title = element_text(family = "Arial", color="#666666", face="plain", size=16), axis.text.x = element_text(family = "Arial", color="#666666", face="plain", size=14), axis.text.y = element_text(family = "Arial", color="#666666", face="plain", size=14)) +
  theme(panel.border = element_rect(linetype = "solid", colour = "grey", fill = NA), panel.grid.major = element_line(color = "grey", linetype = "dotted"), panel.grid.minor = element_line(colour = "grey", linetype = "dotted"), panel.background = element_blank(), axis.line = element_line(colour = "grey40"))


ggplot(genomebusco, aes(x=mean_cov, fill=region)) +
  geom_density(position="identity", alpha=0.5, binwidth = 10)+
  #geom_vline(data=mediangb, aes(xintercept=grp.median, color=region), linetype="dashed")+
  xlim(0,1000) +
  scale_fill_brewer(palette="RdBu") +
  theme(axis.title = element_text(family = "Arial", color="#666666", face="plain", size=16), axis.text.x = element_text(family = "Arial", color="#666666", face="plain", size=14), axis.text.y = element_text(family = "Arial", color="#666666", face="plain", size=14)) +
  theme(panel.border = element_rect(linetype = "solid", colour = "grey", fill = NA), panel.grid.major = element_line(color = "grey", linetype = "dotted"), panel.grid.minor = element_line(colour = "grey", linetype = "dotted"), panel.background = element_blank(), axis.line = element_line(colour = "grey40"))


ggplot(genomebusco, aes(x=mean_cov, fill=region)) +
  geom_histogram(position="identity", alpha=0.5, binwidth = 10)+
  geom_vline(data=meangb, aes(xintercept=grp.mean, color=region), linetype="dashed") +
  xlim(0,750) #+ scale_y_continuous(trans='log10') 


