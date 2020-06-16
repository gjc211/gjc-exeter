library(tidyverse)
library(Amelia)
library(ggplot2)
library(ggpubr)
library(spdep)
library(sf)
library(CARBayes)
library(rgdal)
library(rgeos)
library(lattice)
library(R2jags)
library(runjags)
library(MCMCpack)
library(naniar)

london <- read_csv("London_Pollution.csv")
summary(london$Bloomsbury)
summary(london$Barking)
# Bloomsbury
sum(is.na(filter(london[,-3], str_detect(Date, "2000")))) 
sum(is.na(filter(london[,-3], str_detect(Date, "2001")))) 
sum(is.na(filter(london[,-3], str_detect(Date, "2002")))) 
sum(is.na(filter(london[,-3], str_detect(Date, "2003")))) 
sum(is.na(filter(london[,-3], str_detect(Date, "2004")))) 
# Barking
sum(is.na(filter(london[,-2], str_detect(Date, "2000")))) 
sum(is.na(filter(london[,-2], str_detect(Date, "2001")))) 
sum(is.na(filter(london[,-2], str_detect(Date, "2002")))) 
sum(is.na(filter(london[,-2], str_detect(Date, "2003")))) 
sum(is.na(filter(london[,-2], str_detect(Date, "2004")))) 
missmap(london[,-1])

# Both show most missing values in 2002, but Barking shows high number in 2001 and Bloomsbury in 2003 

# Plotting PM10
london1 <- london
#london1[is.na(london1)] <- 0
london1$no <- c(1:1827)
#london1$na <- rep(0,nrow(london1)) 
#london1$na[london1$Bloomsbury==0] <- 1
time <- london1$Date


g1 <- ggplot(data=london1, mapping=aes(no,Bloomsbury,group=cumsum(is.na(london1$Bloomsbury)))) + 
  geom_line() + annotate("rect",xmin=907,xmax=1237,ymin=0,ymax=Inf,alpha=0.2,fill="red") +
  annotate("rect",xmin=766,xmax=795,ymin=0,ymax=Inf,alpha=0.2,fill="red") +
  scale_x_continuous(breaks = c(0,366,731,1096,1461,1827), labels = c(time[c(1,366,731,1096,1461,1827)])) +
  theme_bw() + labs(x="Date",y="PM10",title="Bloomsbury")

g2 <- ggplot(data=london1, mapping=aes(no,Barking,group=cumsum(is.na(london1$Barking)))) + 
  geom_line() + annotate("rect",xmin=1007,xmax=1103,ymin=0,ymax=Inf,alpha=0.2,fill="red") +
  annotate("rect",xmin=496,xmax=588,ymin=0,ymax=Inf,alpha=0.2,fill="red") +
  scale_x_continuous(breaks = c(0,366,731,1096,1461,1827), labels = c(time[c(1,366,731,1096,1461,1827)])) +
  theme_bw() + labs(x="Date",y="PM10",title="Barking")

figure1 <- ggarrange(g1,g2,ncol = 2, nrow = 1)
annotate_figure(figure1,top = text_grob("PM10 Air Pollution - Significant Missing Data Highlighted",face = "bold", size = 14))

# Mapping

London <- readOGR(dsn = '.', layer = 'London')
                  
plot(London, main="Map of London")
points(530123,182014,col="red",pch=19,cex=2)  
points(548030, 183363,col="blue",pch=19,cex=2)
legend("right", legend=c("Bloomsbury", "Barking"),
       fill=c("red", "blue"))

# Bloomsbury in the west of London show a slight increasing trend in PM10 before the 
# extended period of missing data, but show a decrease thereafter
# Barking in the east of London show a similiarly slight increasing trend up to early
# 2003. Bloomsbury does show more extreme high values than Barking.
london2 <- london[1:1461,]
N <- nrow(london2)
y <- london2$Bloomsbury
dat <- list("N"=nrow(london2),"y"=london2$Bloomsbury)
#dat1 <- list("N"=1460,"y"=diff(london2$Bloomsbury))
params <- c("tau.pr","sd")
params1 <- c("y","pred")

init <- function(){
  list("mu" = rnorm(1),"tau.pr" = dnorm(1))
}

rw.mod <- function() {
  for(i in 2:N) {
    w[i] ~ dnorm(0, pow(1/sd,2))
    pred[i] <- y[i-1] + w[i]
    y[i] ~ dnorm(pred[i], tau.pr) 
    w[i] ~ dnorm(0,pow(1/sd,2))
  }
  pred[1] ~ mu
  tau.pr ~ dgamma(0.001,0.001) 
  mu ~ dnorm(0,0.01) 
  sd <- 1/sqrt(tau.pr)
}

##############################################################
rw.mod <- function() { 
for(i in 2:N) {
  w[i] ~ dnorm(0, pow(1/sd.pr,2))
  pred[i] <- y[i-1] + w[i]
  y[i] ~ dunif(pred[i], tau.pr)
}
pred[1] ~ mu # initial value
mu ~ dnorm(0, 0.01)
tau.pr ~ dgamma(0.001,0.001)
sd.pr <- 1/sqrt(tau.pr)
}
##############################################################

rw.mod <- function() {
  
sigma.obs ~ dunif(1,10) # Prior for error of observation process
sigma2.obs <- pow(sigma.obs, 2)
tau.obs <- pow(sigma.obs, -2)

sigma.proc ~ dunif(1, 10) # Prior for error of state process
sigma2.proc <- pow(sigma.proc, 2)
tau.proc <- pow(sigma.proc, -2)

y.est[1] ~ dnorm(0, 0.001) # Flat prior for initial y

#mean.r ~ dnorm(0, 0.001) # Flat Prior for mean r
for (t in 2:N){
  r[t-1] ~ dnorm(0, tau.proc)
  y.est[t] <- y.est[t-1] + r[t-1]
}

for (t in 1:N) {
  y[t] ~ dnorm(y.est[t], tau.obs)
}
} 

###############################################################

rw.jags <- jags(data=dat, parameters.to.save=params,
                        model.file=rw.mod, n.chains = 2, n.burnin=5000,
                        n.thin=1, n.iter=10000,DIC = TRUE)

densityplot(as.mcmc(rw.jags))

rw.jags1_ <- jags(data=dat, parameters.to.save=params1,
                 model.file=rw.mod, n.chains = 2, n.burnin=5000,
                 n.thin=1, n.iter=10000,DIC = TRUE)


#traceplot(rw.jags1,parameters = 'y[907:910]')

#rw.mcmc <- as.mcmc(rw.jags1)
#autocorr(rw.mcmc)

# 95% credible intervals
london2$est <- round(rw.jags1_$BUGSoutput$mean$pred,1)
df <- data.frame(rw.jags1_$BUGSoutput$summary)
df <- df[-1,]
df <- df[-c(1462:2922),]
upper <- df$X97.5.
lower <- df$X2.5.
ggplot(data=london2, mapping=aes(x=c(1:1461),y=est)) +
  scale_x_continuous(breaks = c(0,366,731,1096,1461), labels = c(time[c(1,366,731,1096,1461)])) +
  theme_bw() + labs(x="Date",y="PM10",title="Bloomsbury") + 
  geom_ribbon(london2,mapping=aes(ymin=lower,ymax=upper),fill="grey70") + geom_line() 


# Second model
set.seed(123)
rw.mod1 <- function() {
  for(i in 3:N) {
    #w[i] ~ dnorm(mean.w, pow(1/sd,2))
    pred[i] <- (2*y[i-1]) - y[i-2] + w[i]
    y[i] ~ dnorm(pred[i], tau.pr) 
    w[i] ~ dnorm(0,pow(1/sd,2))
    }
  pred[1] ~ dnorm(0,0.01)
  pred[2] ~ dnorm(0,0.01)
  tau.pr ~ dgamma(0.001,0.001) 
  #tau ~ dgamma(0.001,0.001) 
  #mean.w ~ dnorm(0,0.001)
  mu ~ dnorm(0,0.01) 
  sd <- 1/sqrt(tau.pr)
}

########
rw.mod1 <- function() {
  
  sigma.obs ~ dunif(1,10) # Prior for error of observation process
  sigma2.obs <- pow(sigma.obs, 2)
  tau.obs <- pow(sigma.obs, -2)
  
  sigma.proc ~ dunif(1, 10) # Prior for error of state process
  sigma2.proc <- pow(sigma.proc, 2)
  tau.proc <- pow(sigma.proc, -2)
  
  y.est[1] ~ dnorm(0, 0.001) # Flat prior for initial y
  y.est[2] ~ dnorm(0, 0.001)
  mean.r ~ dnorm(0, 0.001) # Flat Prior for mean r
  
  for (t in 3:N){
    r[t] ~ dnorm(mean.r, tau.proc)
    y.est[t] <- (2*y.est[t-1]) - y.est[t-2] + r[t]
  }
  
  for (t in 1:N) {
    y[t] ~ dnorm(y.est[t], tau.obs)
  }
} 



rw.jags2 <- jags(data=dat, parameters.to.save=params1,
                 model.file=rw.mod1, n.chains = 2, n.burnin=5000,
                 n.thin=1, n.iter=10000,DIC = TRUE)

rw.jags2_ <- jags(data=dat, parameters.to.save=params,
                  model.file=rw.mod1, n.chains = 2, n.burnin=5000,
                  n.thin=1, n.iter=10000,DIC = TRUE)

traceplot(rw.jags2)
traceplot(rw.jags2_)
densityplot(as.mcmc(rw.jags2_))

#london2$est1 <- round(rw.jags2$BUGSoutput$mean$y,1)
df1 <- data.frame(rw.jags2$BUGSoutput$summary)
df1 <- df1[-1,]
df1 <- df1[-c(1462:2922),]

upper1 <- df1$X97.5.
lower1 <- df1$X2.5.
ggplot(mapping=aes(x=c(1:1461),y=rw.jags2$BUGSoutput$mean$pred)) + 
  scale_x_continuous(breaks = c(0,366,731,1096,1461), labels = c(time[c(1,366,731,1096,1461)])) +
  theme_bw() + labs(x="Date",y="PM10",title="Bloomsbury") +  
  geom_ribbon(mapping=aes(ymin=lower1,ymax=upper1),fill="grey70") + geom_line() 

rw.mod1 <- function() { 
  mu1 ~ dnorm(0, 0.01)
  mu2 ~ dnorm(0, 0.01)
  tau.pro ~ dgamma(0.001,0.001)
  sd.pro <- 1/sqrt(tau.pro)
  predY[1] <- mu1 # initial value
  predY[2] <- mu2 
  for(i in 3:N) {
    w[i] ~ dnorm(0, pow(1/sd.pro,2))
    predY[i] <- 2 * y[i-1] - y[i-2] + w[i]
    y[i] ~ dnorm(predY[i], tau.pro)
  }
}

###########################################################################################################
# Specifying initial values for each missing value 
###########################################################################################################

# Models
rw.modnew <- function() {
  for(i in 2:N) {
    pred[i] ~ dnorm(pred[i-1],tau) 
    y[i] ~ dnorm(pred[i], tau.pr) 
  }
  pred[1] <- mu
  tau.pr ~ dgamma(0.001,0.001) 
  tau ~ dgamma(0.001,0.001) 
  mu ~ dnorm(0,0.01) 
  sd.pr <- 1/sqrt(tau.pr)
  sd <- 1/sqrt(tau)
}

rw.modnew2 <- function() {
  for(i in 3:N) {
    pred[i] ~ dnorm(2*pred[i-1] - pred[i-2], tau)
    y[i] ~ dnorm(pred[i], tau.pr) 
  }
  pred[1] ~ dnorm(0,0.01)
  pred[2] ~ dnorm(0,0.01)
  tau.pr ~ dgamma(0.001,0.001) 
  tau ~ dgamma(0.001,0.001) 
  sd.pr <- 1/sqrt(tau.pr)
  sd <- 1/sqrt(tau)
}

# Bloomsbury - Initial values
miss <- is.na(london2$Bloomsbury)
y.init <- rep(1,1461)
y.init[miss==FALSE] <- NA
y.init[miss==TRUE] <- 20
pred.init = rep(20, N)
y.init1 <- rep(1,1461)
y.init1[miss==FALSE] <- NA
y.init1[miss==TRUE] <- 22
pred.init1 = rep(20, N)

dat.inits <- list( "tau" = 1, "tau.pr" = 1,"y" = y.init, "pred" = pred.init) 
dat.inits1 <- list( "tau" = 1, "tau.pr" = 1,"y" = y.init1, "pred" = pred.init1) 
jags.inits <- list(dat.inits,dat.inits1)
# Barking - Initial values
miss1 <- is.na(london2$Barking)
y.init2 <- rep(1,1461)
y.init2[miss1==FALSE] <- NA
y.init2[miss1==TRUE] <- 20
pred.init2 = rep(20, N)
y.init3 <- rep(1,1461)
y.init3[miss1==FALSE] <- NA
y.init3[miss1==TRUE] <- 22
pred.init3 = rep(20, N)

dat.inits2 <- list( "tau" = 1, "tau.pr" = 1,"y" = y.init2, "pred" = pred.init2) 
dat.inits3 <- list( "tau" = 1, "tau.pr" = 1,"y" = y.init3, "pred" = pred.init3) 
jags.inits1 <- list(dat.inits2,dat.inits3)

                     


# Bloomsbury RW1
rw.jagsnew_ <- jags(data=dat, parameters.to.save=params1, inits = jags.inits,
                    model.file=rw.modnew, n.chains = 2, n.burnin=5000,
                    n.thin=1, n.iter=10000,DIC = TRUE)

# Bloomsbury RW2
rw.jagsnew2 <- jags(data=dat, parameters.to.save=params1, inits = jags.inits,
                   model.file=rw.modnew2, n.chains = 2, n.burnin=5000,
                   n.thin=1, n.iter=10000,DIC = TRUE)

rw.jagsnew2a <- jags(data=dat, parameters.to.save=c("tau","tau.pr","sd","sd.pr"), inits = jags.inits,
                    model.file=rw.modnew2, n.chains = 2, n.burnin=5000,
                    n.thin=1, n.iter=10000,DIC = TRUE)

traceplot(rw.jagsnew2a)
densityplot(as.mcmc(rw.jagsnew2a))

ggplot(mapping=aes(x=c(1:1461),y=rw.jagsnew$BUGSoutput$mean$pred)) + geom_line()                                                                                  
ggplot(mapping=aes(x=c(1:1461),y=rw.jagsnew_$BUGSoutput$mean$pred)) + geom_line()                                                                                  

# Bloomsbury
df_ <- data.frame(rw.jagsnew2$BUGSoutput$summary)
df_ <- df_[-1,]
df_ <- df_[-c(1462:2922),]
upper_ <- df_$X97.5.
lower_ <- df_$X2.5.

df1_ <- data.frame(rw.jagsnew_$BUGSoutput$summary)
df1_ <- df1_[-1,]
df1_ <- df1_[-c(1462:2922),]
upper1_ <- df1_$X97.5.
lower1_ <- df1_$X2.5.

g5 <- ggplot(mapping=aes(x=c(1:1461),y=rw.jagsnew_$BUGSoutput$mean$pred)) + geom_ribbon(mapping=aes(ymin=lower1_,ymax=upper1_),fill="grey70") +
  geom_point(mapping=aes(x=c(1:1461),y=london2$Bloomsbury),alpha=0.2) + geom_line() + theme_bw() + 
  scale_x_continuous(breaks = c(0,366,731,1096,1461), labels = c(time[c(1,366,731,1096,1461)])) + labs(x="Date",y="PM10",title="RW(1)")
                                                                                                      
g6 <- ggplot(mapping=aes(x=c(1:1461),y=rw.jagsnew2$BUGSoutput$mean$pred)) + geom_ribbon(mapping=aes(ymin=lower_,ymax=upper_),fill="grey70") +
     geom_point(mapping=aes(x=c(1:1461),y=london2$Bloomsbury),alpha=0.2) + geom_line(lwd=1) + theme_bw() + 
     scale_x_continuous(breaks = c(0,366,731,1096,1461), labels = c(time[c(1,366,731,1096,1461)])) + labs(x="Date",y="PM10",title="RW(2)")

figure3 <- ggarrange(g5,g6,ncol = 2, nrow = 1)
annotate_figure(figure3,top = text_grob("Bloomsbury",face = "bold", size = 14))

# Difference in convergence

traceplot(rw.jagsnew2) # New RW2
traceplot(rw.jags2) # Old RW2

# Smoothing effects
s1 <- ggplot(mapping=aes(x=c(1:50),y=rw.jagsnew_$BUGSoutput$mean$pred[1:50])) + geom_ribbon(mapping=aes(ymin=lower1_[1:50],ymax=upper1_[1:50]),fill="grey70") +
  geom_point(mapping=aes(x=c(1:50),y=london2$Bloomsbury[1:50]),alpha=0.2) + geom_line() + theme_bw() + labs(x="Date",y="PM10",title="RW(1)") +
  scale_x_continuous(breaks = c(0,10,20,30,40,50), labels = c(time[c(1,10,20,30,40,50)])) + labs(x="Date",y="PM10",title="RW(1)")

s2 <- ggplot(mapping=aes(x=c(1:50),y=rw.jagsnew2$BUGSoutput$mean$pred[1:50])) + geom_ribbon(mapping=aes(ymin=lower_[1:50],ymax=upper_[1:50]),fill="grey70") + 
  geom_point(mapping=aes(x=c(1:50),y=london2$Bloomsbury[1:50]),alpha=0.2) + geom_line() + theme_bw() + labs(x="Date",y="PM10",title="RW(2)") +
  scale_x_continuous(breaks = c(0,10,20,30,40,50), labels = c(time[c(1,10,20,30,40,50)])) + labs(x="Date",y="PM10",title="RW(2)")

figure3a <- ggarrange(s1,s2,ncol = 2, nrow = 1)
annotate_figure(figure3a,top = text_grob("Visualising the Smoothing Effects of Each Model",face = "bold", size = 14))

#############################################################################################
# Predicting 
#############################################################################################
# Initial values
london3 <- london[1:1468,]
#london3 <- replace_with_na(london3,list(Bloomsbury=london3$Bloomsbury[1462:1468]))
dat1 <- list("N"=nrow(london3),"y"=c(london2$Bloomsbury,NA,NA,NA,NA,NA,NA,NA))

miss_ <- is.na(london3$Bloomsbury)
y.init_ <- rep(1,1468)
y.init_[miss_==FALSE] <- NA
y.init_[miss_==TRUE] <- 20
pred.init_ = rep(20, 1468)
y.init1_ <- rep(1,1461)
y.init1_[miss_==FALSE] <- NA
y.init1_[miss_==TRUE] <- 22
pred.init1_ = rep(20, 1468)

dat.inits_ <- list( "tau" = 1, "tau.pr" = 1,"y" = y.init_, "pred" = pred.init_) 
dat.inits1_ <- list( "tau" = 1, "tau.pr" = 1,"y" = y.init1_, "pred" = pred.init1_) 
jags.inits_ <- list(dat.inits_,dat.inits1_)


rw.jags3 <- jags(data=dat1, parameters.to.save=params1, inits = jags.inits_,
                  model.file=rw.modnew, n.chains = 2, n.burnin=5000,
                  n.thin=1, n.iter=10000,DIC = TRUE)

rw.jags4 <- jags(data=dat1, parameters.to.save=params1, inits = jags.inits_,
                 model.file=rw.modnew2, n.chains = 2, n.burnin=5000,
                 n.thin=1, n.iter=10000,DIC = TRUE)

#london3 <- london[1:1468,]

time1 <- c("01/01","02/01","03/01","04/01","05/01","06/01","07/01")
  
g3 <- ggplot() + geom_line(mapping=aes(x=c(1:7),y=rw.jags3$BUGSoutput$mean$pred[1462:1468],colour="Predicted"),lty=2) + 
  geom_line(mapping=aes(x=c(1:7),y=london3$Bloomsbury[1462:1468],colour="Actual")) +
  scale_colour_manual(name=NULL,values=c(Predicted="red",Actual="blue")) +
  theme_bw() + scale_x_continuous(breaks = c(1,2,3,4,5,6,7), labels = c(time1[c(1,2,3,4,5,6,7)])) +
  labs(x="Date",y="PM10",title="RW(1)")

g4 <- ggplot() + geom_line(mapping=aes(x=c(1:7),y=rw.jags4$BUGSoutput$mean$pred[1462:1468],colour="Predicted"),lty=2) + 
  geom_line(mapping=aes(x=c(1:7),y=london3$Bloomsbury[1462:1468],colour="Actual")) +
  scale_colour_manual(name=NULL,values=c(Predicted="red",Actual="blue")) +
  theme_bw() + scale_x_continuous(breaks = c(1,2,3,4,5,6,7), labels = c(time1[c(1,2,3,4,5,6,7)])) +
  labs(x="Date",y="PM10",title="RW(2)")


figure2 <- ggarrange(g3,g4,ncol = 2, nrow = 1,common.legend=TRUE,legend="right")
annotate_figure(figure2,top = text_grob("Predictions for first week in 2004",face = "bold", size = 14))


# With  root mean squared prediction error

rw.mod_error <- function() {
  for(i in 2:N) {
    pred[i] ~ dnorm(pred[i-1],tau) 
    y[i] ~ dnorm(pred[i], tau.pr) 
    err[i] <- ((pred[i] - y[i])^2)
  }
  pred[1] ~ dnorm(0,0.01)
  tau.pr ~ dgamma(0.001,0.001) 
  tau ~ dgamma(0.001,0.001) 
  mu ~ dnorm(0,0.01) 
  sd.pr <- 1/sqrt(tau.pr)
  sd <- 1/sqrt(tau)
}


rw.mod_error2 <- function() {
  for(i in 3:N) {
    pred[i] ~ dnorm(2*pred[i-1] - pred[i-2], tau)
    y[i] ~ dnorm(pred[i], tau.pr) 
    err[i] <- ((pred[i] - y[i])^2)
  }
  pred[1] ~ dnorm(0,0.01)
  pred[2] ~ dnorm(0,0.01)
  tau.pr ~ dgamma(0.001,0.001) 
  tau ~ dgamma(0.001,0.001) 
  sd.pr <- 1/sqrt(tau.pr)
  sd <- 1/sqrt(tau)
}

rw.jags_err <- jags(data=dat1, parameters.to.save=c("err"), inits = jags.inits_,
                 model.file=rw.mod_error, n.chains = 2, n.burnin=5000,
                 n.thin=1, n.iter=10000,DIC = TRUE)

rw.jags_err2 <- jags(data=dat1, parameters.to.save=c("err"), inits = jags.inits_,
                    model.file=rw.mod_error2, n.chains = 2, n.burnin=5000,
                    n.thin=1, n.iter=10000,DIC = TRUE)


sqrt(sum(rw.jags_err$BUGSoutput$mean$err)/1468) # RW1
sqrt(sum(rw.jags_err2$BUGSoutput$mean$err)/1468) # RW2

# On a whole, the first model is a better predictor but for extended periods of missing data,
# the smoothed RW(2) model provides more reliable imputed values. 




#####################################################################################
# Barking
#####################################################################################

dat_b <- list("N"=nrow(london2),"y"=london2$Barking)

# Barking RW1
rw.jagsnew <- jags(data=dat_b, parameters.to.save=params1, inits = jags.inits1,
                 model.file=rw.modnew, n.chains = 2, n.burnin=5000,
                 n.thin=1, n.iter=10000,DIC = TRUE)

traceplot(rw.jagsnew)


# Barking RW2
rw.jagsnew1 <- jags(data=dat_b, parameters.to.save=params1,
                    model.file=rw.modnew2, n.chains = 2, n.burnin=5000, inits = jags.inits1,
                    n.thin=1, n.iter=10000,DIC = TRUE)

rw.jagsnew1a <- jags(data=dat_b, parameters.to.save=c("tau","tau.pr","sd","sd.pr"), inits = jags.inits1,
                       model.file=rw.modnew2, n.chains = 2, n.burnin=5000,
                       n.thin=1, n.iter=10000,DIC = TRUE)

densityplot(as.mcmc(rw.jagsnew1a))
traceplot(rw.jagsnew1)

df3 <- data.frame(rw.bark2$BUGSoutput$summary)
df3 <- df3[-c(1:1462),]
upper3 <- df3$X97.5.
lower3 <- df3$X2.5.
ggplot(mapping=aes(x=c(1:1461),y=rw.bark2$BUGSoutput$mean$pred)) + 
  scale_x_continuous(breaks = c(0,366,731,1096,1461), labels = c(time[c(1,366,731,1096,1461)])) +
  theme_bw() + labs(x="Date",y="PM10",title="Barking") +  
  geom_ribbon(mapping=aes(ymin=lower3,ymax=upper3),fill="grey70") + geom_line()


# Barking
df3_ <- data.frame(rw.jagsnew$BUGSoutput$summary)
df3_ <- df3_[-1,]
df3_ <- df3_[-c(1462:2922),]
upper3_ <- df3_$X97.5.
lower3_ <- df3_$X2.5.

df4_ <- data.frame(rw.jagsnew1$BUGSoutput$summary)
df4_ <- df4_[-1,]
df4_ <- df4_[-c(1462:2922),]
upper4_ <- df4_$X97.5.
lower4_ <- df4_$X2.5.

g7 <- ggplot(mapping=aes(x=c(1:1461),y=rw.jagsnew$BUGSoutput$mean$pred)) + geom_ribbon(mapping=aes(ymin=lower3_,ymax=upper3_),fill="grey70") +
  geom_point(mapping=aes(x=c(1:1461),y=london2$Barking),alpha=0.2) + geom_line() + theme_bw() + 
  scale_x_continuous(breaks = c(0,366,731,1096,1461), labels = c(time[c(1,366,731,1096,1461)])) + labs(x="Date",y="PM10",title="RW(1)")

g8 <- ggplot(mapping=aes(x=c(1:1461),y=rw.jagsnew1$BUGSoutput$mean$pred)) + geom_ribbon(mapping=aes(ymin=lower4_,ymax=upper4_),fill="grey70") +
  geom_point(mapping=aes(x=c(1:1461),y=london2$Barking),alpha=0.2) + geom_line(lwd=1) + theme_bw() + 
  scale_x_continuous(breaks = c(0,366,731,1096,1461), labels = c(time[c(1,366,731,1096,1461)])) + labs(x="Date",y="PM10",title="RW(2)")

figure4 <- ggarrange(g7,g8,ncol = 2, nrow = 1)
annotate_figure(figure4,top = text_grob("Barking",face = "bold", size = 14))

# Convergence

traceplot(rw.jagsnew) # RW1
traceplot(rw.jagsnew1) # RW2





####################################################################################
# Different priors
####################################################################################

rw.mod_prior <- function() {
  for(i in 3:N) {
    #w[i] ~ dnorm(mean.w, pow(1/sd,2))
    pred[i] ~ dnorm(2*pred[i-1] - pred[i-2], tau)
    y[i] ~ dnorm(pred[i], tau.pr) 
    #w[i] ~ dnorm(0,0.001)
  }
  pred[1] ~ dnorm(1,0.01)
  pred[2] ~ dnorm(1,0.01)
  tau.pr ~ dnorm(0.0202,(1/(0.001^2))) 
  tau ~ dnorm(3.8,2.5) 
  #mean.w ~ dnorm(0,0.001)
  #mu ~ dnorm(0,0.01) 
  sd.pr <- 1/sqrt(tau.pr)
  sd <- 1/sqrt(tau)
}

rw.mod_prior <- function() {
  for(i in 3:N) {
    #w[i] ~ dnorm(mean.w, pow(1/sd,2))
    pred[i] ~ dnorm(2*pred[i-1] - pred[i-2], tau)
    y[i] ~ dnorm(pred[i], tau.pr) 
    #w[i] ~ dnorm(0,0.001)
  }
  pred[1] ~ dnorm(1,0.01)
  pred[2] ~ dnorm(1,0.01)
  tau.pr ~ dnorm(0.0162,(1/(0.001^2))) 
  tau ~ dnorm(2.8,2.5) 
  #mean.w ~ dnorm(0,0.001)
  #mu ~ dnorm(0,0.01) 
  sd.pr <- 1/sqrt(tau.pr)
  sd <- 1/sqrt(tau)
}


rw.jags_prior <- jags(data=dat_b, parameters.to.save=params1, inits = jags.inits1,
                    model.file=rw.mod_prior, n.chains = 2, n.burnin=5000,
                    n.thin=1, n.iter=10000,DIC = TRUE)

rw.jags_prior1 <- jags(data=dat, parameters.to.save=c("tau","tau.pr","sd","sd.pr"), inits = jags.inits,
                      model.file=rw.mod_prior, n.chains = 2, n.burnin=5000,
                      n.thin=1, n.iter=10000,DIC = TRUE)

traceplot(rw.jags_prior1)
densityplot(as.mcmc(rw.jags_prior1))

rw.jags_prior1a <- update(rw.jags_prior1,n.iter = 10000)

ggplot(mapping=aes(x=c(1:1461),y=rw.jags_prior$BUGSoutput$mean$pred)) + geom_line()

traceplot(rw.jags_prior)
traceplot(rw.jagsnew2)

rw.jags_prior1b <- jags(data=dat_b, parameters.to.save=c("tau","tau.pr","sd","sd.pr"), inits = jags.inits1,
                       model.file=rw.mod_prior, n.chains = 2, n.burnin=5000,
                       n.thin=1, n.iter=10000,DIC = TRUE)

densityplot(as.mcmc(rw.jags_prior1b))

traceplot(rw.jags_prior)





