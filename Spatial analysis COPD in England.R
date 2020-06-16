install.packages("spdep")
install.packages("sf")
install.packages("CARBayes")
install.packages("rgdal")
install.packages("rgeos")
install.packages("RColorBrewer")
install.packages("ggplot2")
library(spdep)
library(sf)
library(CARBayes)
library(rgdal)
library(rgeos)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(ggthemes)

copd_ex <- read_csv("copdexpected.csv")
copd_ob <- read_csv("copdobserved.csv")

# Plotting gain in hospital admissions over time
copd_time <- data.frame(total=colSums(copd_ob[2:11]),exp=colSums(copd_ex[2:11]))
Year <- as.character(c(2001:2010))
ggplot(data=copd_time, aes(x=Year,y=total,group=1)) + geom_line(aes(),size=1.2,colour="violetred4") + geom_point(aes(),size=2) + geom_line(data=copd_time,aes(x=Year,y=exp),size=1,linetype=2) + labs(x="Year",y="Total no. of hospital admissions for COPD") + ggtitle("Observed and expected hospital admissions for COPD in England between 2001 and 2010") + theme_clean()

# Raw SMRs for years '01, '05 and '09
districts <- readOGR(dsn = '.', layer = 'englandlocalauthority')
copd_smr01 <- copd_ob$Y2001/copd_ex$E2001
copd_smr05 <- copd_ob$Y2005/copd_ex$E2005
copd_smr09 <- copd_ob$Y2009/copd_ex$E2009
df01 <- data.frame(name=copd_ex$Name, ob=copd_ob$Y2001, ex=copd_ex$E2001, smr=copd_smr01) #var=var01)
df05 <- data.frame(name=copd_ex$Name, ob=copd_ob$Y2005, ex=copd_ex$E2005, smr=copd_smr05) #var=var05)
df09 <- data.frame(name=copd_ex$Name, ob=copd_ob$Y2009, ex=copd_ex$E2009, smr=copd_smr09) #var=var09)

SMR01 <- merge(districts,df01,by = 'name')
SMR01 <- st_as_sf(SMR01)
ggplot(SMR01,aes(fill=smr)) + geom_sf(colour=NA) + 
  theme() + labs(x="Longitude",y="Latitude",fill='SMR01') +
  scale_fill_gradientn(colours=brewer.pal(9, "RdPu"),
                       breaks=c(0,0.5,1,1.5,2,2.5,3),limits=c(0,3)) + ggtitle("COPD Mapping - 2001") 

SMR05 <- merge(districts,df05,by = 'name')
SMR05 <- st_as_sf(SMR05)
ggplot(SMR05,aes(fill=smr)) + geom_sf(colour=NA) + 
  theme() + labs(x="Longitude",y="Latitude",fill='SMR05') +
  scale_fill_gradientn(colours=brewer.pal(9, "RdPu"),
                       breaks=c(0,0.5,1,1.5,2,2.5,3),limits=c(0,3)) + ggtitle("COPD Mapping - 2005")

SMR09 <- merge(districts,df09,by = 'name')
SMR09 <- st_as_sf(SMR09)
ggplot(SMR09,aes(fill=smr)) + geom_sf(colour=NA) + 
  theme() + labs(x="Longitude",y="Latitude",fill='SMR09') +
  scale_fill_gradientn(colours=brewer.pal(9, "RdPu"),
                       breaks=c(0,0.5,1,1.5,2,2.5,3),limits=c(0,3)) + ggtitle("COPD Mapping - 2009")

# Variance:
var09 <- df09$smr/df09$ex
var05 <- df05$smr/df05$ex
var01 <- df01$smr/df01$ex

ggplot(data=df01, aes(x=var01,y=smr)) + geom_point(aes(),colour="violetred4") + xlim(0,0.3) + ylim(0,2.5) + labs(x="Variance",y="Raw SMR",title="2001") + theme_clean()
ggplot(data=df05, aes(x=var05,y=smr)) + geom_point() + xlim(0,0.3) + ylim(0,2.5) + labs(x="Variance",y="Raw SMR",title="2005")
ggplot(data=df09, aes(x=var09,y=smr)) + geom_point() + xlim(0,0.3) + ylim(0,2.5) + labs(x="Variance",y="Raw SMR",title="2009")

# Smoothed SMRs
nb01 <- poly2nb(districts,row.names=rownames(districts))
mat01 <- nb2mat(nb01,style="B") 
mod01 <- S.CARleroux(formula=ob~offset(log(ex)),
                     data=df01,family="poisson",
                     W=mat01,burnin=20000,n.sample=100000,
                     thin=10,rho=1)

df01$fitted <- mod01$fitted.values
df01$smrsmooth <- df01$fitted/df01$ex
summary(df01)

# Mapping smoothed SMRs
# 2001
SSMR01 <- merge(districts,df01,by='name')
SSMR01 <- st_as_sf(SSMR01)
ggplot(SSMR01,aes(fill=smrsmooth)) + geom_sf(colour=NA) + 
  theme() + labs(x="Longitude",y="Latitude",fill='SSMR01') +
  scale_fill_gradientn(colours=brewer.pal(9, "RdPu"),
                       breaks=c(0,0.5,1,1.5,2,2.5,3),limits=c(0,3)) + ggtitle("Smoothed SMR COPD Mapping - 2001")

# 2005
mod05 <- S.CARleroux(formula=ob~offset(log(ex)),
                     data=df05,family="poisson",
                     W=mat01,burnin=20000,n.sample=100000,
                     thin=10,rho=1)

df05$fitted <- mod05$fitted.values
df05$smrsmooth <- df05$fitted/df05$ex
summary(df05)
SSMR05 <- merge(districts,df05,by='name')
SSMR05 <- st_as_sf(SSMR05)
ggplot(SSMR05,aes(fill=smrsmooth)) + geom_sf(colour=NA) + 
  theme() + labs(x="Longitude",y="Latitude",fill='SSMR05') +
  scale_fill_gradientn(colours=brewer.pal(9, "RdPu"),
                       breaks=c(0,0.5,1,1.5,2,2.5,3),limits=c(0,3)) + ggtitle("Smoothed SMR COPD Mapping - 2005")

# 2009
mod09 <- S.CARleroux(formula=ob~offset(log(ex)),
                     data=df09,family="poisson",
                     W=mat01,burnin=20000,n.sample=100000,
                     thin=10,rho=1)

df09$fitted <- mod09$fitted.values
df09$smrsmooth <- df09$fitted/df09$ex
summary(df09)
SSMR09 <- merge(districts,df09,by='name')
SSMR09 <- st_as_sf(SSMR09)
ggplot(SSMR09,aes(fill=smrsmooth)) + geom_sf(colour=NA) + 
  theme() + labs(x="Longitude",y="Latitude",fill='SSMR09') +
  scale_fill_gradientn(colours=brewer.pal(9, "RdPu"),
                       breaks=c(0,0.5,1,1.5,2,2.5,3),limits=c(0,3)) + ggtitle("Smoothed SMR COPD Mapping - 2009")

# Using scatter plot to visualise the differences
ggplot(df01,aes(smr,smrsmooth)) + geom_point(colour="blue") + geom_abline(intercept=0,slope=1,colour="red") + labs(x="Raw SMRs",y="Smooth SMRs") + ggtitle("2001")
ggplot(df05,aes(smr,smrsmooth)) + geom_point(colour="blue") + geom_abline(intercept=0,slope=1,colour="red") + labs(x="Raw SMRs",y="Smooth SMRs") + ggtitle("2005")
ggplot(df09,aes(smr,smrsmooth)) + geom_point(colour="violetred4") + geom_abline(intercept=0,slope=1) + labs(x="Raw SMRs",y="Smooth SMRs") + theme_clean() + ggtitle("2009")

# Covariance of SMRs 
cov(df01$smr,df01$smrsmooth)
cov(df05$smr,df05$smrsmooth)
cov(df09$smr,df09$smrsmooth)










