
# 1. 
ohio <- read.csv("Ohio_Data.csv")
ohio$smr <- ohio$Obs/ohio$Exp
hist(ohio$smr)

# 2.
OhioMap <-function(data, ncol=5, figmain="", digits=5, type="e",
                   lower=NULL, upper=NULL) {
  if (is.null(lower)) lower <- min(data)
  if (is.null(upper)) upper <- max(data) 
  
  if (type=="q"){p <- seq(0,1,length=ncol+1)
  br <- round(quantile(data,probs=p),2)}
  if (type=="e"){br <- round(seq(lower,upper,length=ncol+1),2)}
  shading <- gray((ncol-1):0/(ncol-1))
  data.grp <- findInterval(data,vec=br,rightmost.closed=T,all.inside=T)
  data.shad <- shading[data.grp]
  map("county", "ohio", fill=TRUE, col=data.shad)
  leg.txt<-paste("[",br[ncol],",",br[ncol+1],"]",sep="")
  for(i in (ncol-1):1){
    leg.txt<-append(leg.txt,paste("[",br[i],",",br[i+1],")",sep=""),)
  }
  leg.txt<-rev(leg.txt)
  legend(-81.4,39.4,legend=leg.txt,fill=shading,bty="n",ncol=1,cex=.8)
  title(main=figmain,cex=1.5)
  invisible()
}
library(maps)
map.text("county", "ohio")
testdat <- ohio$smr
OhioMap(testdat,ncol=8,type="e",figmain="Ohio Lung Cancer SMRs",lower=0.36,upper=1.54)

# 3.
# Relative risk is the probability of disease given exposure/probability of disease given no exposure.
# The thetas represent prior disctributions of the count data and therefore will influence formulation 
# of posterior for the count data and therefore the relative risk. This is because the prior can indicate
# the risk relative to each county based on past results. The beta naught represents the constant
# which is a part of the linear predictor.

# 4.
install.packages("R2jags", dependencies = TRUE, repos = "http://cran.us.r-project.org")
install.packages("runjags", dependencies = TRUE, repos = "http://cran.us.r-project.org")
install.packages("MCMCpack", dependencies = TRUE, repos = "http://cran.us.r-project.org")
library(lattice)
library(R2jags)
library(runjags)
library(MCMCpack)

o <- ohio$Obs
e <- ohio$Exp
N <- nrow(ohio)

o <- ohio$Obs
e <- ohio$Exp
N <- nrow(ohio)
dat <- list("o","e","N")
params <- c("beta","rr")
params1 <- c("alpha","beta")
init <- function(){
  list("beta" = rnorm(1))
}

library(R2jags)
bayes.mod <- function() {
  for(i in 1:N){
    o[i] ~ dpois(mu[i])
    theta[i] ~ dgamma(alpha,alpha)
    mu[i] <- e[i] * exp(beta) * theta[i]
    rr[i] <- exp(beta) * theta[i]
  }
  alpha ~ dgamma(1,1)
  beta ~ dunif(-100,100)
}

bayes.mod.fit <- jags(data = dat, inits = init,
                      parameters.to.save = params1, n.chains = 2, n.iter = 10000,
                      n.burnin = 5000, model.file = bayes.mod)

traceplot(bayes.mod.fit)

# 5.
# Chains appear to have converged. Discarding the first 5000 samples has appeared to be sufficient 
# to ensure appropriate estimates for the parameters.

# Converting to MCMC
mcmc.fit <- as.mcmc(bayes.mod.fit)
summary(mcmc.fit)

xyplot(mcmc.fit)
densityplot(mcmc.fit)

# Formally checking convergence 
gelman.diag(mcmc.fit)
# Here we can observe perfect convergence

# Checking for autocorrelation
autocorr.plot(mcmc.fit)
# Possibly some autocorrelation in beta due to there still being some correlation at lag=20 

bayes.mod.fitnew <- jags(data = dat, inits = init,
                      parameters.to.save = params1, n.chains = 2, n.iter = 10000,
                      n.burnin = 5000, n.thin = 20, model.file = bayes.mod)

mcmc.fitnew <- as.mcmc(bayes.mod.fitnew)
autocorr.plot(mcmc.fitnew)

# Including thinning to 10 has gone a way to reduce such autocorrelation in beta but there is still some correlation at lag=10 

# 6.

bayes.mod.fit1 <- jags(data = dat, inits = init,
                      parameters.to.save = params, n.chains = 2, n.iter = 10000,
                      n.burnin = 5000, model.file = bayes.mod)


rr <- bayes.mod.fit1$BUGSoutput$mean$rr
  
# Mapping RR
ohio$rr <- rr
OhioMap(rr,ncol=8,type="e",figmain="Ohio Lung Cancer RRs",lower=0.64,upper=1.36)

# 7.
prob.par <- c("rr.crit","beta")

prob.model <-  function() {
  # Priors
  alpha ~ dgamma(1,1)
  beta ~ dunif(-100,100)
for(i in 1:N){
    mu[i] <- e[i] * exp(beta) * theta[i]
    o[i] ~ dpois(mu[i])
    rr[i] <- exp(beta) * theta[i]
  }
  # Predictive
  #theta ~ dgamma(alpha,alpha)
for(k in 1:N){
      theta[k] ~ dgamma(alpha,alpha)
      rr.crit[k] <- ifelse(theta[k]>1.2,1,0)
  }
}



bayes.mod.prob <- jags(data = dat, inits = init,
                       parameters.to.save = prob.par, n.chains = 2, n.iter = 10000,
                       n.burnin = 5000, model.file = prob.model)

prob <- bayes.mod.prob$BUGSoutput$mean$rr.crit
ohio$prob <- prob
OhioMap(prob,ncol=8,type="e",figmain="Ohio Lung Cancer - Probability RR > 1.2",lower=0,upper=0.9)

# 8.

# Different priors for p(beta) and alpha 

                    

new.priors.mod <- function() {
  for(i in 1:N){
    o[i] ~ dpois(mu[i])
    theta[i] ~ dgamma(alpha,alpha)
    mu[i] <- e[i] * exp(beta) * theta[i]
    rr[i] <- exp(beta) * theta[i]
  }
  alpha ~ dgamma(4,0.25)
  beta ~ dnorm(-0.04,sqrt(0.001))
}

new.priors.mod.fit <- jags(data = dat, inits = init,
                      parameters.to.save = params1, n.chains = 2, n.thin = 20, n.iter = 10000,
                      n.burnin = 5000, model.file = new.priors.mod)
upd <- update(new.priors.mod.fit, n.iter=10000)

new.priors.mod.fit1 <- jags(data = dat, inits = init,
                           parameters.to.save = params, n.chains = 2, n.iter = 10000,
                           n.burnin = 5000, model.file = new.priors.mod)


rrnew <- new.priors.mod.fit1$BUGSoutput$mean$rr
ohio$rrnew <- rrnew

OhioMap(rrnew,ncol=8,type="e",figmain="Ohio Lung Cancer RR - New Priors" ,lower=0,upper=2)
OhioMap(rr,ncol=8,type="e",figmain="Ohio Lung Cancer RRs",lower=0,upper=2)


densityplot(as.mcmc(new.priors.mod.fit))
densityplot(as.mcmc(upd))

# The pattern of the distribution of alpha has changed, as it now shows an uneven distribution about 
# the mean, with two peaks either side 

# What it has shown is a greater convergence of the chains, especially considering alpha

gelman.diag(as.mcmc(bayes.mod.fitnew))

par(mfrow=c(1,2))
plot(bayes.mod.fit$BUGSoutput$sims.list$beta,bayes.mod.fit$BUGSoutput$sims.list$alpha, main="Original Priors",xlab="beta",ylab="alpha")
plot(new.priors.mod.fit$BUGSoutput$sims.list$beta,new.priors.mod.fit$BUGSoutput$sims.list$alpha, main="New Priors",xlab="beta",ylab="alpha")

