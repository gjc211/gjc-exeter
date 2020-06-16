library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(chron)
library(ggpubr)
library(ggthemes)
library(cluster)
library(rnrfa)
library(verification)
library(cclust)
library(knitr)

# 1.
characteristics <- read.csv("Characteristics.csv", stringsAsFactors=FALSE)
table1 <- summary(characteristics[,3:6])
options(knitr.table.format="html")
kable(table1) %>% cat(., file = "table1.html")
hist(characteristics$Percentage_IC, main="Percentage industrial and commercial customers")
hist(characteristics$Transformer_RATING, main="Transformer rating")
ground_monit <- filter(characteristics, TRANSFORMER_TYPE == "Grd Mtd Dist. Substation")
pole_monit <- filter(characteristics, TRANSFORMER_TYPE == "Pole Mtd Dist. Substation")
summary(ground_monit$TRANSFORMER_TYPE)
summary(pole_monit$TRANSFORMER_TYPE)
# Ground mounted more frequent with 706 stations compared to 242

# 2.
# Much fewer pole mounted stations
mean(ground_monit$Percentage_IC)
# 0.4096
mean(pole_monit$Percentage_IC)
# 0.2929
# Ground mounted with higher proportion of industrial and commerical customers
cor(ground_monit$Percentage_IC,ground_monit$Transformer_RATING)
cor(pole_monit$Percentage_IC,pole_monit$Transformer_RATING)
# Urban ground mounted indicate stronger relationship between % IC customers and size of total power generated
cor(characteristics$LV_FEEDER_COUNT,characteristics$Transformer_RATING)
# Suprisingly low relationship of feeder count and total energy supplied

# 3.
load("Autumn_2012.RData")
autumn_scaled <- Autumn_2012[-c(147:291)]
autumn_scaled <- autumn_scaled %>%
  group_by(Station) %>% 
  summarise_each(funs(mean)) 
autumn_scaled$Date <- NULL
Station <- autumn_scaled$Station 
autumn_scaled$Station <- NULL
# Data reduction needed - PCA
cov <- cov(autumn_scaled)
autm.eig <- eigen(cov)
eigenvectors <- autm.eig$vectors
eigenvalues <- autm.eig$values
sum(eigenvalues[1:3])/sum(eigenvalues)
# First 3 PCs contain 90% of info in original dataset
autumn.matrix <- as.matrix(autumn_scaled)
reduced <- autumn.matrix %*% eigenvectors[,1:3]
autumn_scaled <- data.frame(Station,reduced)
# Distance matrix and dendogram
autumn_scaled$Station <- as.character(autumn_scaled$Station)
autumn.dist <- dist(autumn_scaled[,2:4])
autumn.hclust <- hclust(autumn.dist)
plot(autumn.hclust, cex = 0.6, hang = -1)

# 4.
groups <- cutree(autumn.hclust,5)
table(groups)
autumn_scaled$Station[groups == 3]
# Two main clusters formed 
agnes <- agnes(autumn_scaled[,2:4], method = 'complete')
pltree(agnes, cex = 0.6, hang = -1, main = "Dendrogram (agnes)")
agnes.groups <- cutree(agnes,5)
table(agnes.groups)
autumn_scaled$Station[agnes.groups == 1]
# Same results as previous
# Labelling  
clusters <- cbind(Group=groups,autumn_scaled)

# 5.
autumn_scaled1 <- Autumn_2012[-c(147:291)] 
# All days
all_days <- autumn_scaled1 %>%
  group_by(Station) %>% 
  summarise_each(funs(mean)) 
all_days$Date <- NULL
all_days <- cbind(Group=groups,all_days)
cluster1_alldays <- filter(all_days,Group == 1)
cluster1_alldays$Group <- NULL
cluster1_alldays <- cluster1_alldays %>% gather(time,value,-Station)
cluster1_alldays$time <- as.numeric(cluster1_alldays$time)
newtime <- format(seq.POSIXt(as.POSIXct(Sys.Date()), as.POSIXct(Sys.Date()+1), by = "10 min"),"%H:%M", tz="GMT")
newtime = newtime[1:144]
# Cluster 1
ad1 <- ggplot(cluster1_alldays,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="Cluster 1")
# Cluster 2
cluster2_alldays <- filter(all_days,Group == 2)
cluster2_alldays$Group <- NULL
cluster2_alldays <- cluster2_alldays %>% gather(time,value,-Station)
cluster2_alldays$time <- as.numeric(cluster2_alldays$time)
ad2 <- ggplot(cluster2_alldays,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="Cluster 2")
# Cluster 3
cluster3_alldays <- filter(all_days,Group == 3)
cluster3_alldays$Group <- NULL
cluster3_alldays <- cluster3_alldays %>% gather(time,value,-Station)
cluster3_alldays$time <- as.numeric(cluster3_alldays$time)
ad3 <- ggplot(cluster3_alldays,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="Cluster 3")
# Cluster 4
cluster4_alldays <- filter(all_days,Group == 4)
cluster4_alldays$Group <- NULL
cluster4_alldays <- cluster4_alldays %>% gather(time,value,-Station)
cluster4_alldays$time <- as.numeric(cluster4_alldays$time)
ad4 <- ggplot(cluster4_alldays,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="Cluster 4")
# Cluster 5
cluster5_alldays <- filter(all_days,Group == 5)
cluster5_alldays$Group <- NULL
cluster5_alldays <- cluster5_alldays %>% gather(time,value,-Station)
cluster5_alldays$time <- as.numeric(cluster5_alldays$time)
ad5 <- ggplot(cluster5_alldays,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="Cluster 5")
figure1 <- ggarrange(ad1,ad2,ad3,ad4,ad5,ncol = 3, nrow = 2)
annotate_figure(figure1,top = text_grob("Daily Average Demand for Substations - All Days",face = "bold", size = 14))

# Weekdays
converted_dates <- dates(autumn_scaled1[,2], origin = c(month = 1,day = 1,year = 1970))
autumn_scaled1$Date <- converted_dates
autumn_scaled1$day <- weekdays(as.Date(autumn_scaled1$Date))
weekdays <- filter(autumn_scaled1, day %in% c("Monday","Tuesday","Wednesday","Thursday","Friday"))
weekdays$day <- NULL
weekdays <- weekdays %>%
  group_by(Station) %>% 
  summarise_each(funs(mean)) 
weekdays$Date <- NULL
weekdays <- cbind(Group=groups,weekdays)
cluster1_weekdays <- filter(weekdays,Group == 1)
cluster1_weekdays$Group <- NULL
cluster1_weekdays <- cluster1_weekdays %>% gather(time,value,-Station)
cluster1_weekdays$time <- as.numeric(cluster1_weekdays$time)
# Cluster 1
wd1 <- ggplot(cluster1_weekdays,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="Cluster 1")
# Cluster 2
cluster2_weekdays <- filter(weekdays,Group == 2)
cluster2_weekdays$Group <- NULL
cluster2_weekdays <- cluster2_weekdays %>% gather(time,value,-Station)
cluster2_weekdays$time <- as.numeric(cluster2_weekdays$time)
wd2 <- ggplot(cluster2_weekdays,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="Cluster 2")
# Cluster 3
cluster3_weekdays <- filter(weekdays,Group == 3)
cluster3_weekdays$Group <- NULL
cluster3_weekdays <- cluster3_weekdays %>% gather(time,value,-Station)
cluster3_weekdays$time <- as.numeric(cluster3_weekdays$time)
wd3 <- ggplot(cluster3_weekdays,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="Cluster 3")
# Cluster 4
cluster4_weekdays <- filter(weekdays,Group == 4)
cluster4_weekdays$Group <- NULL
cluster4_weekdays <- cluster4_weekdays %>% gather(time,value,-Station)
cluster4_weekdays$time <- as.numeric(cluster4_weekdays$time)
wd4 <- ggplot(cluster4_weekdays,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="Cluster 4")
# Cluster 5
cluster5_weekdays <- filter(weekdays,Group == 5)
cluster5_weekdays$Group <- NULL
cluster5_weekdays <- cluster5_weekdays %>% gather(time,value,-Station)
cluster5_weekdays$time <- as.numeric(cluster5_weekdays$time)
wd5 <- ggplot(cluster5_weekdays,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="Cluster 5")
figure2 <- ggarrange(wd1,wd2,wd3,wd4,wd5,ncol = 3, nrow = 2)
annotate_figure(figure2,top = text_grob("Daily Average Demand for Substations - Weekdays",face = "bold", size = 14))

# Saturdays
saturdays <- filter(autumn_scaled1, day == "Saturday")
saturdays$day <- NULL
saturdays <- saturdays %>%
  group_by(Station) %>% 
  summarise_each(funs(mean)) 
saturdays$Date <- NULL
clusters$Station[!(clusters$Station %in% saturdays$Station)]
# 512190 and 512202 with no values for Saturday (57th and 58th)
groups1 <- groups[-c(57,58)]
saturdays <- cbind(Group=groups1,saturdays)
cluster1_saturdays <- filter(saturdays,Group == 1)
cluster1_saturdays$Group <- NULL
cluster1_saturdays <- cluster1_saturdays %>% gather(time,value,-Station)
cluster1_saturdays$time <- as.numeric(cluster1_saturdays$time)
# Cluster 1
sd1 <- ggplot(cluster1_saturdays,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="Cluster 1")
# Cluster 2
cluster2_saturdays <- filter(saturdays,Group == 2)
cluster2_saturdays$Group <- NULL
cluster2_saturdays <- cluster2_saturdays %>% gather(time,value,-Station)
cluster2_saturdays$time <- as.numeric(cluster2_saturdays$time)
sd2 <- ggplot(cluster2_saturdays,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="Cluster 2")
# Cluster 3
cluster3_saturdays <- filter(saturdays,Group == 3)
cluster3_saturdays$Group <- NULL
cluster3_saturdays <- cluster3_saturdays %>% gather(time,value,-Station)
cluster3_saturdays$time <- as.numeric(cluster3_saturdays$time)
sd3 <- ggplot(cluster3_saturdays,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="Cluster 3")
# Cluster 4
cluster4_saturdays <- filter(saturdays,Group == 4)
cluster4_saturdays$Group <- NULL
cluster4_saturdays <- cluster4_saturdays %>% gather(time,value,-Station)
cluster4_saturdays$time <- as.numeric(cluster4_saturdays$time)
sd4 <- ggplot(cluster4_saturdays,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="Cluster 4")
# Cluster 5
cluster5_saturdays <- filter(saturdays,Group == 5)
cluster5_saturdays$Group <- NULL
cluster5_saturdays <- cluster5_saturdays %>% gather(time,value,-Station)
cluster5_saturdays$time <- as.numeric(cluster5_saturdays$time)
sd5 <- ggplot(cluster5_saturdays,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="Cluster 5")
figure3 <- ggarrange(sd1,sd2,sd3,sd4,sd5,ncol = 3, nrow = 2)
annotate_figure(figure3,top = text_grob("Daily Average Demand for Substations - Saturdays",face = "bold", size = 14))

# Sundays
sundays <- filter(autumn_scaled1, day == "Sunday")
sundays$day <- NULL
sundays <- sundays %>%
  group_by(Station) %>% 
  summarise_each(funs(mean)) 
sundays$Date <- NULL
clusters$Station[!(clusters$Station %in% sundays$Station)]
# 512190, 512202 and 536787 with no values for Sunday (57th, 58th and 249th)
groups2 <- groups[-c(57,58,249)]
sundays <- cbind(Group=groups2,sundays)
cluster1_sundays <- filter(sundays,Group == 1)
cluster1_sundays$Group <- NULL
cluster1_sundays <- cluster1_sundays %>% gather(time,value,-Station)
cluster1_sundays$time <- as.numeric(cluster1_sundays$time)
# Cluster 1
sun1 <- ggplot(cluster1_sundays,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="Cluster 1")
# Cluster 2
cluster2_sundays <- filter(sundays,Group == 2)
cluster2_sundays$Group <- NULL
cluster2_sundays <- cluster2_sundays %>% gather(time,value,-Station)
cluster2_sundays$time <- as.numeric(cluster2_sundays$time)
sun2 <- ggplot(cluster2_sundays,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="Cluster 2")
# Cluster 3
cluster3_sundays <- filter(sundays,Group == 3)
cluster3_sundays$Group <- NULL
cluster3_sundays <- cluster3_sundays %>% gather(time,value,-Station)
cluster3_sundays$time <- as.numeric(cluster3_sundays$time)
sun3 <- ggplot(cluster3_sundays,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="Cluster 3")
# Cluster 4
cluster4_sundays <- filter(sundays,Group == 4)
cluster4_sundays$Group <- NULL
cluster4_sundays <- cluster4_sundays %>% gather(time,value,-Station)
cluster4_sundays$time <- as.numeric(cluster4_sundays$time)
sun4 <- ggplot(cluster4_sundays,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="Cluster 4")
# Cluster 5
cluster5_sundays <- filter(sundays,Group == 5)
cluster5_sundays$Group <- NULL
cluster5_sundays <- cluster5_sundays %>% gather(time,value,-Station)
cluster5_sundays$time <- as.numeric(cluster5_sundays$time)
sun5 <- ggplot(cluster5_sundays,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="Cluster 5")
figure4 <- ggarrange(sun1,sun2,sun3,sun4,sun5,ncol = 3, nrow = 2)
annotate_figure(figure4,top = text_grob("Daily Average Demand for Substations - Sundays",face = "bold", size = 14))

# 6.
cluster1_stat <- unique(cluster1_alldays$Station) 
char_c1 <- filter(characteristics,SUBSTATION_NUMBER %in% c(cluster1_stat))
cluster2_stat <- unique(cluster2_alldays$Station)
char_c2 <- filter(characteristics,SUBSTATION_NUMBER %in% c(cluster2_stat))
cluster3_stat <- unique(cluster3_alldays$Station)
char_c3 <- filter(characteristics,SUBSTATION_NUMBER %in% c(cluster3_stat))
cluster4_stat <- unique(cluster4_alldays$Station)
char_c4 <- filter(characteristics,SUBSTATION_NUMBER %in% c(cluster4_stat))
cluster5_stat <- unique(cluster5_alldays$Station)
char_c5 <- filter(characteristics,SUBSTATION_NUMBER %in% c(cluster5_stat))
summary(char_c1[,3:6])
summary(char_c2[,3:6])
summary(char_c3[,3:6])
summary(char_c4[,3:6])
summary(char_c5[,3:6])

# 7.
# Cluster 1 with much fewer industrial/commercial customers 
# Cluster 1 with less power being delivered on average (rating)
# Cluster 3 much higher % IC customers and more power delivered
# Cluster 3 with fewer customers on average per substation 
c1_ground <- filter(char_c1,TRANSFORMER_TYPE == "Grd Mtd Dist. Substation")
c2_ground <- filter(char_c2,TRANSFORMER_TYPE == "Grd Mtd Dist. Substation")
c3_ground <- filter(char_c3,TRANSFORMER_TYPE == "Grd Mtd Dist. Substation")
nrow(c1_ground)/nrow(char_c1)
nrow(c2_ground)/nrow(char_c2)
nrow(c3_ground)/nrow(char_c3)
# Cluster 1 with fewer % of ground mounted (more rural)
# Clear difference in pattern between rural and urban
# What separates cluster 2 is the consistently high demand substations throughout the day
# Cluster 2 shows signs of both power demands of cluster 1 and 3
# Cluster 2 shows similar pattern to cluster 1 during weekends but shows similar pattern to cluster 3 on weekdays
# Describe the main two power demand patterns each day 
# Cluster 4 with few customers and feeders per substation, more domestic
# Cluster 5 clearly different pattern - no demand through daytime (streetlights?)

# Names: 
# C1 - Highly rural/domestic customers (low demand and many customers)
# C2 - Moderate/high industrial/commercial (Variable by substation but more uniform through day)
# C3 - Highly urban/industrial and commercial customers (high demand and few customers)
# C4 - Few customers/feeders, moderately rural/domestic 
# C5 - No demand during daytime - streetlights

# 8. 
newsub <- read_csv("NewSubstations.csv")
newsub$X1 <- NULL
#newsub$Date <- as.integer(as.Date(newsub$Date, format="%d/%m/%Y"))
#newsub$Date <- as.numeric(newsub$Date)
# Sub 511079
`511079` <- filter(newsub, Substation == "511079")
`511079`$Substation <- NULL
`511079` <- `511079` %>% gather(time,value,-Date)
`511079`$Date <- weekdays(as.Date(`511079`$Date))
`511079`$Date <- factor(`511079`$Date,levels=c("Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday"),labels=c("Weekday","Weekday","Weekday","Weekday","Weekday","Saturday","Sunday"))
`511079`$time <- as.numeric(`511079`$time)
wday079 <- filter(`511079`, Date == "Weekday")
wday079 <- wday079 %>% group_by(time) %>% summarise_each(funs(mean))
sat079 <- filter(`511079`, Date == "Saturday")
sat079 <- sat079 %>% group_by(time) %>% summarise_each(funs(mean))
sun079 <- filter(`511079`, Date == "Sunday")
sun079 <- sun079 %>% group_by(time) %>% summarise_each(funs(mean))
alldays079 <- `511079` %>% group_by(time) %>% summarise_each(funs(mean))
sub1 <- ggplot() + geom_line(alldays079,mapping=aes(time,value,colour="All"),lwd=1.1) + geom_line(wday079,mapping=aes(time,value,colour="Weekdays"),lwd=1.1) + geom_line(sat079,mapping=aes(time,value,colour="Saturdays"),lwd=1.1) + geom_line(sun079,mapping=aes(time,value,colour="Sundays"),lwd=1.1) + scale_colour_manual(name="Weekly Period",values=c(All="black",Weekdays="red",Saturdays="blue",Sundays="purple")) + theme_bw() + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="511079")
# Sub 512457
`512457` <- filter(newsub, Substation == "512457")
`512457`$Substation <- NULL
`512457` <- `512457` %>% gather(time,value,-Date)
`512457`$Date <- weekdays(as.Date(`512457`$Date))
`512457`$Date <- factor(`512457`$Date,levels=c("Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday"),labels=c("Weekday","Weekday","Weekday","Weekday","Weekday","Saturday","Sunday"))
`512457`$time <- as.numeric(`512457`$time)
wday457 <- filter(`512457`, Date == "Weekday")
wday457 <- wday457 %>% group_by(time) %>% summarise_each(funs(mean))
sat457 <- filter(`512457`, Date == "Saturday")
sat457 <- sat457 %>% group_by(time) %>% summarise_each(funs(mean))
sun457 <- filter(`512457`, Date == "Sunday")
sun457 <- sun457 %>% group_by(time) %>% summarise_each(funs(mean))
alldays457 <- `512457` %>% group_by(time) %>% summarise_each(funs(mean))
sub2 <- ggplot() + geom_line(alldays457,mapping=aes(time,value,colour="All"),lwd=1.1) + geom_line(wday457,mapping=aes(time,value,colour="Weekdays"),lwd=1.1) + geom_line(sat457,mapping=aes(time,value,colour="Saturdays"),lwd=1.1) + geom_line(sun457,mapping=aes(time,value,colour="Sundays"),lwd=1.1) + scale_colour_manual(name="Weekly Period",values=c(All="black",Weekdays="red",Saturdays="blue",Sundays="purple")) + theme_bw() + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="512457")
# Sub 532697
`532697` <- filter(newsub, Substation == "532697")
`532697`$Substation <- NULL
`532697` <- `532697` %>% gather(time,value,-Date)
`532697`$Date <- weekdays(as.Date(`532697`$Date))
`532697`$Date <- factor(`532697`$Date,levels=c("Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday"),labels=c("Weekday","Weekday","Weekday","Weekday","Weekday","Saturday","Sunday"))
`532697`$time <- as.numeric(`532697`$time)
wday697 <- filter(`532697`, Date == "Weekday")
wday697 <- wday697 %>% group_by(time) %>% summarise_each(funs(mean))
sat697 <- filter(`532697`, Date == "Saturday")
sat697 <- sat697 %>% group_by(time) %>% summarise_each(funs(mean))
sun697 <- filter(`532697`, Date == "Sunday")
sun697 <- sun697 %>% group_by(time) %>% summarise_each(funs(mean))
alldays697 <- `532697` %>% group_by(time) %>% summarise_each(funs(mean))
sub3 <- ggplot() + geom_line(alldays697,mapping=aes(time,value,colour="All"),lwd=1.1) + geom_line(wday697,mapping=aes(time,value,colour="Weekdays"),lwd=1.1) + geom_line(sat697,mapping=aes(time,value,colour="Saturdays"),lwd=1.1) + geom_line(sun697,mapping=aes(time,value,colour="Sundays"),lwd=1.1) + scale_colour_manual(name="Weekly Period",values=c(All="black",Weekdays="red",Saturdays="blue",Sundays="purple")) + theme_bw() + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="532697")
# Sub 552863
`552863` <- filter(newsub, Substation == "552863")
`552863`$Substation <- NULL
`552863` <- `552863` %>% gather(time,value,-Date)
`552863`$Date <- weekdays(as.Date(`552863`$Date))
`552863`$Date <- factor(`552863`$Date,levels=c("Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday"),labels=c("Weekday","Weekday","Weekday","Weekday","Weekday","Saturday","Sunday"))
`552863`$time <- as.numeric(`552863`$time)
wday863 <- filter(`552863`, Date == "Weekday")
wday863 <- wday863 %>% group_by(time) %>% summarise_each(funs(mean))
sat863 <- filter(`552863`, Date == "Saturday")
sat863 <- sat863 %>% group_by(time) %>% summarise_each(funs(mean))
sun863 <- filter(`552863`, Date == "Sunday")
sun863 <- sun863 %>% group_by(time) %>% summarise_each(funs(mean))
alldays863 <- `552863` %>% group_by(time) %>% summarise_each(funs(mean))
sub4 <- ggplot() + geom_line(alldays863,mapping=aes(time,value,colour="All"),lwd=1.1) + geom_line(wday863,mapping=aes(time,value,colour="Weekdays"),lwd=1.1) + geom_line(sat863,mapping=aes(time,value,colour="Saturdays"),lwd=1.1) + geom_line(sun863,mapping=aes(time,value,colour="Sundays"),lwd=1.1) + scale_colour_manual(name="Weekly Period",values=c(All="black",Weekdays="red",Saturdays="blue",Sundays="purple")) + theme_bw() + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="552863")
# Sub 563729
`563729` <- filter(newsub, Substation == "563729")
`563729`$Substation <- NULL
`563729` <- `563729` %>% gather(time,value,-Date)
`563729`$Date <- weekdays(as.Date(`563729`$Date))
`563729`$Date <- factor(`563729`$Date,levels=c("Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday"),labels=c("Weekday","Weekday","Weekday","Weekday","Weekday","Saturday","Sunday"))
`563729`$time <- as.numeric(`563729`$time)
wday729 <- filter(`563729`, Date == "Weekday")
wday729 <- wday729 %>% group_by(time) %>% summarise_each(funs(mean))
sat729 <- filter(`563729`, Date == "Saturday")
sat729 <- sat729 %>% group_by(time) %>% summarise_each(funs(mean))
sun729 <- filter(`563729`, Date == "Sunday")
sun729 <- sun729 %>% group_by(time) %>% summarise_each(funs(mean))
alldays729 <- `563729` %>% group_by(time) %>% summarise_each(funs(mean))
sub5 <- ggplot() + geom_line(alldays729,mapping=aes(time,value,colour="All"),lwd=1.1) + geom_line(wday729,mapping=aes(time,value,colour="Weekdays"),lwd=1.1) + geom_line(sat729,mapping=aes(time,value,colour="Saturdays"),lwd=1.1) + geom_line(sun729,mapping=aes(time,value,colour="Sundays"),lwd=1.1) + scale_colour_manual(name="Weekly Period",values=c(All="black",Weekdays="red",Saturdays="blue",Sundays="purple")) + theme_bw() + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="563729")

figure5 <- ggarrange(sub1,sub2,sub3,sub4,sub5,ncol = 3, nrow = 2,common.legend = TRUE,legend = "right")
annotate_figure(figure5,top = text_grob("Daily Average Demand for New Substations",face = "bold", size = 14))
# All but the second new substation (512457) show similar demand trends, that is domestic dominated
# All but the first station have higher demand for Sunday around midday compared to other times in the week 

# 9. 
# Scaling
newsub.matrix <- newsub[-c(1:2)]
i1 <- !rowSums(newsub.matrix==1)>0
newsub.matrix[i1,] <- newsub.matrix[i1,]/do.call(pmax, newsub.matrix[i1,])
newsub.scaled <- cbind(Station=newsub$Substation,Date=newsub$Date,newsub.matrix)
newsub.scaled$Station <- factor(newsub.scaled$Station,levels=c("511079","512457","532697","552863","563729"),labels=c("511079new","512457new","532697new","552863new","563729new"))
# PCA
autumn_scalednew <- Autumn_2012[-c(147:291)] 
autumn_new <- rbind(autumn_scalednew,newsub.scaled)
autumn_new <- autumn_new %>% group_by(Station) %>% summarise_each(funs(mean)) 
autumn_new$Date <- NULL
Station <- autumn_new$Station 
autumn_new$Station <- NULL
cov.new <- cov(autumn_new)
autm_new.eig <- eigen(cov.new)
eigenvectors <- autm_new.eig$vectors
eigenvalues <- autm_new.eig$values
autumn_new <- as.matrix(autumn_new)
reduced.new <- autumn_new %*% eigenvectors[,1:3]
autumn_new <- data.frame(Station,reduced.new)
autumn_new.dist <- dist(autumn_new[,2:4])
autumn_new.hclust <- hclust(autumn_new.dist)
table(cutree(autumn_new.hclust,5))
autumn_new <- cbind(groups=(cutree(autumn_new.hclust,5)),autumn_new)
newgroups <- filter(autumn_new, str_detect(Station, "new"))

# Using k-means - Euclidean nearest neighbour algorithm
library(class)
df <- data.frame(Station=autumn_new$Station,reduced.new)
dfnew <- filter(df, str_detect(Station, "new"))
success <- data.frame(groups)
success$s <- rep(0, nrow(success))
success$s[groups==1] <- 1
s <- c(success$s)
# Using 80 nearest neighbours to allocate new substations
knn(autumn_scaled[,2:4],dfnew[,2:4],k=80,groups,use.all = TRUE)
# Only change is first station now becomes cluster 2 - another step maybe to use another distance measure

# 10. 
# Would have expected all but the second new substation (512457) to be in cluster 1, as is observed.
# It however has been assigned to cluster 2 rather than cluster 3 (which would have been expected based
# on industrial/commerical trends) - although we observed similarities in cluster 2 and 3 so this 
# discrepancy can be explained by possibly limited number of clusters reducing the reliability of repeated
# clustering.

# An area for further study would be repeating the formation of clusters, but using k-means algorithm and 
# usings more clusters in order to separate out cluster 2 which show characteristics of both 1 and 3.

# 11. Can substation type help to predict daily demand trends and how does cluster structure and their trends vary over seasons?

# (1)
characteristics <- read.csv("Characteristics.csv", stringsAsFactors=FALSE)
characteristics$GRID_REFERENCE = gsub(" ", "", characteristics$GRID_REFERENCE, fixed=T)
x = osg_parse(characteristics$GRID_REFERENCE)
characteristics$east = x[[1]]
characteristics$north = x[[2]]
write.csv(characteristics, "char_clean.csv", row.names=F)
characteristics <- read.csv("char_clean.csv", stringsAsFactors=FALSE)
ggplot() + geom_point(characteristics,mapping=aes(x=east,y=north,group=TRANSFORMER_TYPE,colour=TRANSFORMER_TYPE)) + annotate(geom="text",label="Cardiff",x=3200000,y=1680000,label.size=0.35,color="black",fontface="italic") + annotate(geom="text",label="Swansea",x=2650000,y=1850000,label.size=0.35,color="black",fontface="italic") + theme_map() + theme(plot.title = element_text(face = "bold",size = 15)) + scale_colour_manual(values=c("Grd Mtd Dist. Substation" = "red","Pole Mtd Dist. Substation" = "blue")) + labs(x="Easting",y="Northing",title="Map of Substations in South Wales",colour="Transformer Type")
# Clear grouping of ground mounted transformers in the major cities of Cardiff and Swansea, with scattered pole mounted transformers in rural areas and in the valleys to the north east

# (2) Important correlations (rating and type of customer/type of customer and transformer type - between seasons)
cor(ground_monit$Transformer_RATING,ground_monit$LV_FEEDER_COUNT)
cor(pole_monit$Transformer_RATING,pole_monit$LV_FEEDER_COUNT)
cor(ground_monit$Percentage_IC,ground_monit$Transformer_RATING)
cor(pole_monit$Percentage_IC,pole_monit$Transformer_RATING)

# (3) Logistic regression with transformer type
characteristics$successes <- rep(0, nrow(characteristics))
characteristics$successes[characteristics$TRANSFORMER_TYPE == "Grd Mtd Dist. Substation" ] <- 1
glm <- glm(successes ~ Transformer_RATING + LV_FEEDER_COUNT + Percentage_IC, data = characteristics, family = binomial)
summary(glm)
glm1 <- glm(successes ~ Transformer_RATING, data = characteristics, family = binomial)
characteristics$predictions <- predict(glm1, type = "response")
glm2 <- glm(successes ~ Percentage_IC, data = characteristics, family = binomial)
characteristics$predictions1 <- predict(glm2, type = "response")
par(mfrow = c(1,2))
roc.plot(characteristics$successes,characteristics$predictions, cex.main = 0.85, main = "Transformer rating")
roc.plot(characteristics$successes,characteristics$predictions1, cex.main = 0.85, main = "Percentage of industrial/commercial customers")
mtext("ROC plots for predicting transformer type", side = 3, line = -1, outer = TRUE)
roc.area(characteristics$successes,characteristics$predictions)
roc.area(characteristics$successes,characteristics$predictions1)

# (4) Average demand pattern for ground mounted compared to pole mounted (two graphs in one figure)
autumn_scaled <- Autumn_2012[-c(147:291)]
characteristics_comp <- filter(characteristics, SUBSTATION_NUMBER %in% c(unique(autumn_scaled$Station)))
characteristics_comp$TRANSFORMER_TYPE <- as.factor(characteristics_comp$TRANSFORMER_TYPE)
ground.comp <- filter(autumn_scaled, Station %in% characteristics_comp$SUBSTATION_NUMBER[characteristics_comp$TRANSFORMER_TYPE == "Grd Mtd Dist. Substation"])
ground.means <- data.frame(time=c(1:144),means=(colMeans(ground.comp[-c(1:2)])))
pole.comp <- filter(autumn_scaled, Station %in% characteristics_comp$SUBSTATION_NUMBER[characteristics_comp$TRANSFORMER_TYPE == "Pole Mtd Dist. Substation"])
pole.means <- data.frame(time=c(1:144),means=(colMeans(pole.comp[-c(1:2)])))
ground.comp <- ground.comp %>%
  group_by(Station) %>% 
  summarise_each(funs(mean)) 
ground.comp$Date <- NULL
ground.comp <- ground.comp %>% gather(time,value,-Station)
ground.comp$time <- as.numeric(ground.comp$time)
ground.plot <- ggplot(ground.comp,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + geom_line(data=ground.means,aes(x=time,y=means),lwd=1.5) + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="Ground mounted")

pole.comp <- pole.comp %>%
  group_by(Station) %>% 
  summarise_each(funs(mean)) 
pole.comp$Date <- NULL
pole.comp <- pole.comp %>% gather(time,value,-Station)
pole.comp$time <- as.numeric(pole.comp$time)
pole.plot <- ggplot(pole.comp,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + geom_line(data=pole.means,aes(x=time,y=means),lwd=1.5) + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="Pole mounted")
figure6 <- ggarrange(ground.plot,pole.plot,ncol = 2, nrow = 1)
annotate_figure(figure6,top = text_grob("Daily Average Demand by Substation Type",face = "bold", size = 14))
# Shows no clear separation of domestic and industrial trends by substation type, but does show clear disparity within pole mounted substations
# There appear to be much fewer industrial/commercial customers where there are pole mounted substations - street lamps 
# Next step to compare by season - would expect substations to be allocated differently by clustering

# (5) Power demand changes between seasons - rural substations more variable?
load("C:/Users/greg-/AppData/Local/Temp/Temp1_Data.zip/Data/HighSummer_2012.RData")
load("C:/Users/greg-/AppData/Local/Temp/Temp1_Data.zip/Data/Summer_2012.RData")
load("C:/Users/greg-/AppData/Local/Temp/Temp1_Data.zip/Data/Winter_2012.RData")
load("C:/Users/greg-/AppData/Local/Temp/Temp1_Data.zip/Data/Spring_2013.RData")
# Summer 
Summer_2012 <- Summer_2012[-c(147:291)]
summer.means <- data.frame(time=c(1:144),means=(colMeans(Summer_2012[-c(1:2)])))
summer.scaled <- Summer_2012 %>%
  group_by(Station) %>% 
  summarise_each(funs(mean)) 
summer.scaled$Date <- NULL
summer.scaled <- summer.scaled %>% gather(time,value,-Station)
summer.scaled$time <- as.numeric(summer.scaled$time)
sum <- ggplot(summer.scaled,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + geom_line(data=summer.means,aes(x=time,y=means),lwd=1.5) + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="Summer 2012")
# High summer
HighSummer_2012 <- HighSummer_2012[-c(147:291)]
highsummer.means <- data.frame(time=c(1:144),means=(colMeans(HighSummer_2012[-c(1:2)])))
highsummer.scaled <- HighSummer_2012 %>%
  group_by(Station) %>% 
  summarise_each(funs(mean)) 
highsummer.scaled$Date <- NULL
highsummer.scaled <- highsummer.scaled %>% gather(time,value,-Station)
highsummer.scaled$time <- as.numeric(highsummer.scaled$time)
highsum <- ggplot(highsummer.scaled,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + geom_line(data=highsummer.means,aes(x=time,y=means),lwd=1.5) + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="High Summer 2012")
# Autumn
Autumn_2012 <- Autumn_2012[-c(147:291)]
autumn.means <- data.frame(time=c(1:144),means=(colMeans(Autumn_2012[-c(1:2)])))
autumn.scaled <- Autumn_2012 %>%
  group_by(Station) %>% 
  summarise_each(funs(mean)) 
autumn.scaled$Date <- NULL
autumn.scaled <- autumn.scaled %>% gather(time,value,-Station)
autumn.scaled$time <- as.numeric(autumn.scaled$time)
aut <- ggplot(autumn.scaled,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + geom_line(data=autumn.means,aes(x=time,y=means),lwd=1.5) + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="Autumn 2012")
# Winter
Winter_2012 <- Winter_2012[-c(147:291)]
winter.means <- data.frame(time=c(1:144),means=(colMeans(Winter_2012[-c(1:2)])))
winter.scaled <- Winter_2012 %>%
  group_by(Station) %>% 
  summarise_each(funs(mean)) 
winter.scaled$Date <- NULL
winter.scaled <- winter.scaled %>% gather(time,value,-Station)
winter.scaled$time <- as.numeric(winter.scaled$time)
wint <- ggplot(winter.scaled,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + geom_line(data=winter.means,aes(x=time,y=means),lwd=1.5) + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Power delivered (KW) - Scaled",title="Winter 2012")
# Spring
Spring_2013 <- Spring_2013[-c(147:291)]
spring.means <- data.frame(time=c(1:144),means=(colMeans(Spring_2013[-c(1:2)])))
spring.scaled <- Spring_2013 %>%
  group_by(Station) %>% 
  summarise_each(funs(mean)) 
spring.scaled$Date <- NULL
spring.scaled <- spring.scaled %>% gather(time,value,-Station)
spring.scaled$time <- as.numeric(spring.scaled$time)
spr <- ggplot(spring.scaled,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + geom_line(data=spring.means,aes(x=time,y=means),lwd=1.5) + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Real power delivered (KW) - Scaled",title="Spring 2013")
figure7 <- ggarrange(sum,highsum,aut,wint,spr,ncol = 3, nrow = 2)
annotate_figure(figure7,top = text_grob("Daily Average Demand By Season",face = "bold", size = 14))
# Appears to be little change in general demand trends - main difference is the average demand staying more constant through the evening in summer and spring
# but shows a marked increase at around 8pm for autumn and winter
# By substation type
# Summer
ground.summ <- filter(Summer_2012, Station %in% characteristics_comp$SUBSTATION_NUMBER[characteristics_comp$TRANSFORMER_TYPE == "Grd Mtd Dist. Substation"])
groundsumm_means <- data.frame(time=c(1:144),means=(colMeans(ground.summ[-c(1:2)])))
pole.summ <- filter(Summer_2012, Station %in% characteristics_comp$SUBSTATION_NUMBER[characteristics_comp$TRANSFORMER_TYPE == "Pole Mtd Dist. Substation"])
polesumm_means <- data.frame(time=c(1:144),means=(colMeans(pole.summ[-c(1:2)])))
ground.summ <- ground.summ %>%
  group_by(Station) %>% 
  summarise_each(funs(mean)) 
ground.summ$Date <- NULL
ground.summ <- ground.summ %>% gather(time,value,-Station)
ground.summ$time <- as.numeric(ground.summ$time)
gsumm <- ggplot(ground.summ,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + geom_line(data=groundsumm_means,aes(x=time,y=means),lwd=1.5) + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Real power delivered (KW) - Scaled",title="Ground mounted")
pole.summ <- pole.summ %>%
  group_by(Station) %>% 
  summarise_each(funs(mean)) 
pole.summ$Date <- NULL
pole.summ <- pole.summ %>% gather(time,value,-Station)
pole.summ$time <- as.numeric(pole.summ$time)
psumm <- ggplot(pole.summ,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + geom_line(data=polesumm_means,aes(x=time,y=means),lwd=1.5) + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Real power delivered (KW) - Scaled",title="Pole mounted")
figure8 <- ggarrange(gsumm,psumm,ncol = 2, nrow = 1)
annotate_figure(figure8,top = text_grob("Daily average demand by substation type - Summer",face = "bold", size = 14))
# Winter
ground.wint <- filter(Winter_2012, Station %in% characteristics_comp$SUBSTATION_NUMBER[characteristics_comp$TRANSFORMER_TYPE == "Grd Mtd Dist. Substation"])
groundwint_means <- data.frame(time=c(1:144),means=(colMeans(ground.wint[-c(1:2)])))
pole.wint <- filter(Winter_2012, Station %in% characteristics_comp$SUBSTATION_NUMBER[characteristics_comp$TRANSFORMER_TYPE == "Pole Mtd Dist. Substation"])
polewint_means <- data.frame(time=c(1:144),means=(colMeans(pole.wint[-c(1:2)])))
ground.wint <- ground.wint %>%
  group_by(Station) %>% 
  summarise_each(funs(mean)) 
ground.wint$Date <- NULL
ground.wint <- ground.wint %>% gather(time,value,-Station)
ground.wint$time <- as.numeric(ground.wint$time)
gwint <- ggplot(ground.wint,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + geom_line(data=groundwint_means,aes(x=time,y=means),lwd=1.5) + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Real power delivered (KW) - Scaled",title="Ground mounted")
pole.wint <- pole.wint %>%
  group_by(Station) %>% 
  summarise_each(funs(mean)) 
pole.wint$Date <- NULL
pole.wint <- pole.wint %>% gather(time,value,-Station)
pole.wint$time <- as.numeric(pole.wint$time)
pwint <- ggplot(pole.wint,aes(x=time,y=value)) + geom_line(alpha=.4,aes(group=Station,colour=Station)) + theme_bw() + scale_color_continuous(type = "viridis") + geom_line(data=polewint_means,aes(x=time,y=means),lwd=1.5) + theme(legend.position = "none") + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Real power delivered (KW) - Scaled",title="Pole mounted")
figure9 <- ggarrange(gwint,pwint,ncol = 2, nrow = 1)
annotate_figure(figure9,top = text_grob("Daily average demand by substation type - Winter",face = "bold", size = 14))
# Shows that both type of substation follow a similar pattern across season - e.g. street lamps active for more time due to less light in winter 
# Now can analyse whether individual substations change across seasons in terms of cluster membership

# (6) Optimal number of clusters between seasons
# Summer
Summer_2012 <- Summer_2012 %>% group_by(Station) %>% summarise_each(funs(mean))
v <- vector(length=19)
for(i in 2:20){
  kc <- kmeans(Summer_2012[-c(1:2)],i,iter.max=30)
  kv <- c(kc$tot.withinss)
  v[i] <- kv
}
yv <- v[2:20]
xclusters <- c(2:20)
df <- data.frame(x=xclusters, y=yv)
# Calculating the point of maximum curvature
out.spl <- with(df, smooth.spline(x, y, df = 3))
derivative.out <- with(df, predict(out.spl, x, deriv = 2))
derivative.out <- as.data.frame(derivative.out)
max(derivative.out1$y)
# Indicates 8 as optimal number of clusters

# Winter
Winter_2012 <- Winter_2012 %>% group_by(Station) %>% summarise_each(funs(mean))
v1 <- vector(length=19)
for(i in 2:20){
  kc1 <- kmeans(Winter_2012[-c(1:2)],i,iter.max=30)
  kv1 <- c(kc1$tot.withinss)
  v1[i] <- kv1
}
yv1 <- v1[2:20]
xclusters1 <- c(2:20)
df1 <- data.frame(x=xclusters1, y=yv1)
# Calculating the point of maximum curvature
out.spl1 <- with(df1, smooth.spline(x, y, df = 3))
derivative.out1 <- with(df1, predict(out.spl1, x, deriv = 2))
derivative.out1 <- as.data.frame(derivative.out1)
max(derivative.out1$y)
# Indicates 7 as the optimal number of clusters
ggplot() + geom_line(data=df, aes(x=xclusters, y=yv, color="Summer")) + geom_point(data=df, aes(x=xclusters, y=yv, color="Summer")) + geom_line(data=df1, aes(x=xclusters1, y=yv1, color="Winter")) + geom_point(data=df1, aes(x=xclusters1, y=yv1, color="Winter")) + geom_vline(xintercept = 8, linetype = 2, colour = "red") + geom_vline(xintercept = 7, linetype = 2, colour = "purple") + theme_bw() + labs(x="Number of clusters", y="Total within sum of squares",title="Optimal number of clusters across seasons") + scale_colour_manual(name="Season",values=c(Summer="red",Winter="purple"))  

# Using manhattan distance measure to compare results
# Summer
V <- vector(length=19)
XC <- as.matrix(Summer_2012[-c(1:2)])
for(i in 2:20){
  cclust <- cclust(XC,i,dist = "manhattan",method="kmeans")
  within <- c(sum(cclust$withinss))
  V[i] <- within
}
YV <- V[2:20]
DF <- data.frame(x=xclusters, y=YV)
# Calculating point of maximum curvature, showing optimal number of clusters to be 7
OUT.SPL <- with(DF, smooth.spline(x, y, df = 3))
DERIV.OUT <- with(DF, predict(OUT.SPL, x, deriv = 2))
DERIV.OUT <- as.data.frame(DERIV.OUT)
# Winter
V1 <- vector(length=19)
XC1 <- as.matrix(Winter_2012[-c(1:2)])
for(i in 2:20){
  cclust1 <- cclust(XC1,i,dist = "manhattan",method="kmeans")
  within1 <- c(sum(cclust1$withinss))
  V1[i] <- within1
}
YV1 <- V1[2:20]
DF1 <- data.frame(x=xclusters, y=YV1)
OUT.SPL1 <- with(DF1, smooth.spline(x, y, df = 3))
DERIV.OUT1 <- with(DF1, predict(OUT.SPL1, x, deriv = 2))
DERIV.OUT1 <- as.data.frame(DERIV.OUT1)
# Similar results are obtained with other distance measure - more variable total WSS across higher numbers of clusters with this method in particular 

# (7) Cluster membership between seasons - rural substations more variable? (summer and winter)

# Chronologically using 8 clusters for each season - summaries
# Demand profiles for domestic and industrial clusters  
# Finding common stations for consistent comparison
stations <- Reduce(intersect, list(unique(Summer_2012$Station),unique(Spring_2013$Station),unique(HighSummer_2012$Station),unique(Winter_2012$Station),unique(Autumn_2012$Station)))
Spring_2013 <- Spring_2013 %>% group_by(Station) %>% summarise_each(funs(mean)) %>% filter(Station %in% c(stations))
Summer_2012 <- Summer_2012 %>% group_by(Station) %>% summarise_each(funs(mean)) %>% filter(Station %in% c(stations))
HighSummer_2012 <- HighSummer_2012 %>% group_by(Station) %>% summarise_each(funs(mean)) %>% filter(Station %in% c(stations))
Autumn_2012 <- Autumn_2012 %>% group_by(Station) %>% summarise_each(funs(mean)) %>% filter(Station %in% c(stations))
Winter_2012 <- Winter_2012 %>% group_by(Station) %>% summarise_each(funs(mean)) %>% filter(Station %in% c(stations))
set.seed(123)
km_summer <- kmeans(Summer_2012[-c(1:2)],8,iter.max=30)
km_winter <- kmeans(Winter_2012[-c(1:2)],8,iter.max=30)
season <- data.frame(station=Summer_2012$Station,summer=km_summer$cluster,winter=km_winter$cluster)
# Format tables
table(season$summer)
table(season$winter)

# Change in cluster demand patterns
Summer_2012 <- Summer_2012 %>% group_by(Station) %>% summarise_each(funs(mean)) 
Summer_2012 <- cbind(cluster=season$summer,Summer_2012)
Summer_2012$Date <- NULL
clusters_summ <- Summer_2012 %>% group_by(cluster) %>% summarise_each(funs(mean))
clusters_summ <- clusters_summ %>% gather(time,value,-cluster)
clusters_summ$time <- as.numeric(clusters_summ$time)
clusters_summ$cluster <- as.character(clusters_summ$cluster)
a <- ggplot(clusters_summ,aes(x=time,y=value)) + geom_line(aes(group=cluster,colour=cluster),lwd=1) + theme_bw() + scale_color_discrete() + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Real power delivered (KW) - Scaled",title="Summer",colour="Cluster")

Winter_2012 <- cbind(cluster=season$winter,Winter_2012)
Winter_2012$Date <- NULL
clusters_winter <- Winter_2012 %>% group_by(cluster) %>% summarise_each(funs(mean))
clusters_winter$Station <- NULL
clusters_winter <- clusters_winter %>% gather(time,value,-cluster)
clusters_winter$time <- as.numeric(clusters_winter$time)
clusters_winter$cluster <- as.character(clusters_winter$cluster)
b <- ggplot(clusters_winter,aes(x=time,y=value)) + geom_line(aes(group=cluster,colour=cluster),lwd=1) + theme_bw() + scale_color_discrete() + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Real power delivered (KW) - Scaled",title="Winter",colour="Cluster")
figure11 <- ggarrange(a,b,ncol = 2, nrow = 1,common.legend = TRUE,legend="right")
annotate_figure(figure11,top = text_grob("Average cluster demand patterns",face = "bold", size = 14))
# Illustrates why we shouldn't look just at cluster membership as an indicator of changing patterns - due to the algorithm being random - would have to specify initial centroids
# Shows 1-2 I/C clusters in summer and winter 
# The domestic/urban patterns are much more similar in summer for each cluster  
km_winter1 <- kmeans(Winter_2012[-c(1:2)],7,iter.max=30)
Winter_2012 <- cbind(cluster1=km_winter1$cluster,Winter_2012)
Winter_2012$cluster <- NULL
clusters_winter <- Winter_2012 %>% group_by(cluster1) %>% summarise_each(funs(mean))
clusters_winter$Station <- NULL
clusters_winter <- clusters_winter %>% gather(time,value,-cluster1)
clusters_winter$time <- as.numeric(clusters_winter$time)
clusters_winter$cluster1 <- as.character(clusters_winter$cluster1)
ggplot(clusters_winter,aes(x=time,y=value)) + geom_line(aes(group=cluster1,colour=cluster1),lwd=1) + theme_bw() + scale_color_discrete() + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Real power delivered (KW) - Scaled",title="Winter with optimal number of clusters",colour="Cluster")
# When using the optimal number of clusters for winter, similar results obtained, with one cluster showing significant disparity (cluster 7 in previous graph)
# Summer show more predictable patterns for each cluster - more variable and exaggerated patterns in winter 

# Characteristics of cluster 7 for winter to assess profile with all stations for comparison
wint.cluster7_char <- filter(characteristics, SUBSTATION_NUMBER %in% c(Winter_2012$Station[Winter_2012$cluster == 7]))
summary(wint.cluster7_char[,3:6])
summary(characteristics[,3:6])
# Relatively low number of customers on average, low % IC indicating domestic, but with mostly ground mounted substations and high transformer rating

# Change in these substations in cluster 7
cluster7stations <- Winter_2012$Station[Winter_2012$cluster == 7]
cluster7sum <- filter(Summer_2012, Station %in% c(cluster7stations))
cluster7wint <- filter(Winter_2012, Station %in% c(cluster7stations))
cluster7wint$cluster <- NULL
cluster7sum$cluster <- NULL
cluster7wint <- cluster7wint %>% gather(time,value,-Station)
cluster7sum <- cluster7sum %>% gather(time,value,-Station)
cluster7sum$time <- as.numeric(cluster7sum$time)
cluster7wint$time <- as.numeric(cluster7wint$time)
cluster7sum$Station <- as.character(cluster7sum$Station)
cluster7wint$Station <- as.character(cluster7wint$Station)
c <- ggplot(cluster7sum,aes(x=time,y=value)) + geom_line(aes(group=Station,colour=Station),lwd=1) + theme_bw() + scale_color_discrete() + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Real power delivered (KW) - Scaled",title="Summer",colour="Station")
d <- ggplot(cluster7wint,aes(x=time,y=value)) + geom_line(aes(group=Station,colour=Station),lwd=1) + theme_bw() + scale_color_discrete() + scale_x_continuous(breaks = c(0,24,48,72,96,120,143), labels = c(newtime[c(1,25,49,73,97,121,144)])) + labs(x="Time",y="Real power delivered (KW) - Scaled",title="Winter",colour="Station")
figure12 <- ggarrange(c,d,ncol = 2, nrow = 1,common.legend = TRUE,legend="right")
annotate_figure(figure12,top = text_grob("Variable demand trends of stations from Cluster 7 - Winter",face = "bold", size = 14))
# Shows a significant change in demand trends of these 14 substations from Cluster 7 - Winter

