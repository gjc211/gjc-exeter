library(astsa)
library(tseries)
library(numDeriv)
library(dlm)
library(forecast)

rapid <- read.csv('moc_transports.csv')
rapid$qyyyy <- paste(rapid$year,rapid$Quarter, sep='-')
rapidmean <- tapply(rapid$Overturning_Strength,rapid$qyyyy,mean)
rapidmean.ts <- ts(as.vector(rapidmean),start=c(2004,2),frequency = 4)

# Explain theory of ARMA - apply backward shift operator
par(mfrow=c(1,2))
acf(rapidmean.ts,lag.max = 39)
pacf(rapidmean.ts,lag.max = 39)


# Clearly ARMA from the oscillating damped sine curve

rm.arma <- arma(rapidmean.ts)
rm.arma2 <- arima(rapidmean.ts,order=c(2,0,2),method="ML")
rm.arma2$aic
rm.arma3 <- arima(rapidmean.ts,order=c(3,0,2),method="ML")
rm.arma3$aic
rm.arma4 <- arima(rapidmean.ts,order=c(3,0,1),method="ML")
rm.arma4$aic


plot(rapidmean.ts)
lines((rapidmean.ts-rm.arma2$residuals),col="blue",lty=2)
lines((rapidmean.ts-rm.arma3$residuals),col="red",lty=2)
lines((rapidmean.ts-rm.arma4$residuals),col="purple",lty=2)

rapidmean.ts-rm.arma4$residuals # Optimal predictions for each quarter

# Residuals - explain white noise and iid requirements (uncorrelated, independent random variables)
tsdiag(rm.arima2)

# Explain theory of why to use ARIMA rather than ARMA (seasonality)

adf.test(rapidmean.ts)
adf.test(rapidmean.ts,k=1) # Quarterly time dependency

rm.arima <- arima(rapidmean.ts,order=c(2,1,2),method="ML")
rm.arima$aic
plot(rapidmean.ts)
lines((rapidmean.ts-rm.arima$residuals),col="blue",lty=2)

rm.arima1 <- arima(rapidmean.ts,order=c(2,1,1),method="ML")
rm.arima1$aic
lines((rapidmean.ts-rm.arima1$residuals),col="red",lty=2)

# Seasonal
rm.arima2 <- arima(rapidmean.ts,order=c(2,1,1),seasonal = list(order=c(0,1,2)), method="ML")
rm.arima2$aic
lines((rapidmean.ts-rm.arima2$residuals),col="purple",lty=2)

# Forecast

rm.fore <- forecast(rm.arima2,h=6,level=c(90))
plot(rm.fore)


#########################################
# DLM
#########################################
rm.dlm_ <- dlmModPoly(order=2) + dlmModSeas(4) # Why did I pick 3rd order? 
rm.dlm_ <- dlmModPoly(order=3,dW=c(rep(1.2,2),1)) + dlmModSeas(4) 
rm.dlm_ <- dlmModPoly(order=3,dV=3.1,dW=c(rep(1.2,2),1)) + dlmModSeas(4) 

buildFun_ <- function(x) {
  diag(W(rm.dlm_))[2:4] <- exp(x[1:3])
  V(rm.dlm_) <- exp(x[4])
  return(rm.dlm_)
}


# By MLE
(fit_ <- dlmMLE(rapidmean.ts, parm = c(-5,-20,-17,2), build = buildFun_))$conv
rm.dlm_ <- buildFun_(fit_$par)
drop(V(rm.dlm_))
diag(W(rm.dlm_))[2:3]

fit_ <- dlmMLE(rapidmean.ts, parm = c(-5,-20,-17,2), build = buildFun_, hessian = T)
rm.ll_ <- hessian(function(x) dlmLL(rapidmean.ts, buildFun_(x)), fit_$par)
all(eigen(rm.ll_, only.values = TRUE)$values > 0) # Positive definite?

#########################################
# Another for comparison 
#########################################

rm.dlm <- dlmModPoly(order=3,dV=3.1,dW=c(rep(1.2,2),1)) + dlmModSeas(4) 
buildFun <- function(x) {
  diag(W(rm.dlm))[2:4] <- exp(x[1:3])
  V(rm.dlm) <- exp(x[4])
  return(rm.dlm)
}

# By MLE
(fit <- dlmMLE(rapidmean.ts, parm = rep(0, 4), build = buildFun))$conv
rm.dlm <- buildFun(fit$par)
drop(V(rm.dlm))
diag(W(rm.dlm))[2:3]

fit <- dlmMLE(rapidmean.ts, parm = rep(0, 4), build = buildFun, hessian = T)
rm.ll <- hessian(function(x) dlmLL(rapidmean.ts, buildFun(x)), fit$par)
all(eigen(rm.ll, only.values = TRUE)$values > 0) # Not positive definite 


# smoothing estimates of states
rm.Smooth <- dlmSmooth(rapidmean.ts, mod = rm.dlm)
xs <- cbind(rapidmean.ts, dropFirst(rm.Smooth$s[,c(1,4)]))
colnames(xs) <- c("Overturning", "Trend", "Seasonal")
plot(xs, type = 'o', main = "North Atlantic Overturning Strength") 
# Separating out trend and seasonal (quarterly) cycle

# Comments on stationarity 

rm.filt <- dlmFilter(rapidmean.ts,mod=rm.dlm)
par(mfrow=c(1,1))
qqnorm(residuals(rm.filt, sd = FALSE))
qqline(residuals(rm.filt, sd = FALSE))
tsdiag(rm.filt)

rm.fore <- dlmForecast(rm.filt,nAhead=5)

sqrtR <- sapply(rm.fore$R, function(x) sqrt(x[1,1]))
pl <- rm.fore$a[,1] + qnorm(0.05, sd = sqrtR)
pu <- rm.fore$a[,1] + qnorm(0.95, sd = sqrtR)
x <- ts.union(rapidmean.ts, rm.Smooth$s[,1] ,rm.fore$a[,1], pl, pu) 
plot(x, plot.type = "single", type = 'o', pch = c(1, 0, 20, 3, 3),col = c("darkgrey", "darkgrey", "brown", "yellow", "yellow"), ylab = "North Atlantic Overturning Strength")
legend("bottomright", legend = c("Observed", "Smoothed (deseasonalized)", "Forecasted level", "90% probability limit"), bty = 'n', pch = c(1, 0, 20, 3, 3), lty = 1,col = c("darkgrey", "darkgrey", "brown", "yellow", "yellow"))

# Previous smoothing estimates of states 
rm.Smooth_ <- dlmSmooth(rapidmean.ts, mod = rm.dlm_)
xs_ <- cbind(rapidmean.ts, dropFirst(rm.Smooth_$s[,c(1,4)]))
colnames(xs_) <- c("Overturning", "Trend", "Seasonal")
plot(xs_, type = 'o', main = "North Atlantic Overturning Strength")

rm.filt_ <- dlmFilter(rapidmean.ts,mod=rm.dlm_)
par(mfrow=c(1,1))
qqnorm(residuals(rm.filt_, sd = FALSE))
qqline(residuals(rm.filt_, sd = FALSE))
tsdiag(rm.filt_)

rm.fore_ <- dlmForecast(rm.filt_,nAhead=5)

sqrtR_ <- sapply(rm.fore_$R, function(x) sqrt(x[1,1]))
pl_ <- rm.fore_$a[,1] + qnorm(0.05, sd = sqrtR_)
pu_ <- rm.fore_$a[,1] + qnorm(0.95, sd = sqrtR_)
x_ <- ts.union(rapidmean.ts, rm.Smooth_$s[,1] ,rm.fore_$a[,1], pl_, pu_) 
plot(x_, plot.type = "single", type = 'o', pch = c(1, 0, 20, 3, 3),col = c("darkgrey", "darkgrey", "brown", "yellow", "yellow"), ylab = "North Atlantic Overturning Strength")
legend("bottomright", legend = c("Observed", "Smoothed (deseasonalized)", "Forecasted level", "90% probability limit"), bty = 'n', pch = c(1, 0, 20, 3, 3), lty = 1,col = c("darkgrey", "darkgrey", "brown", "yellow", "yellow"))

# Comparing the two models by AIC
-2*(-fit$value-length(fit$par))
-2*(-fit_$value-length(fit_$par)) # New model better fit by this criteria

# MSE
mean((residuals(rm.filt, type = "raw", sd = FALSE))^2)
mean((residuals(rm.filt_, type = "raw", sd = FALSE))^2) # Old model gives slightly lower MSE


# Bayesian estimation

# Gibbs sampler
outGibbs <- dlmGibbsDIG(rapidmean.ts, dlmModPoly(2) + dlmModSeas(4), 
                        a.y = 5.9, b.y = 1000, a.theta = c(0.02,0.06), 
                        b.theta = 1000, shape.y = 12, rate.y = 1, shape.theta = 3, rate.theta = 2,
                        n.sample = 1100, ind = c(2, 3),
                        save.states = FALSE)

outGibbs1 <- dlmGibbsDIG(rapidmean.ts, dlmModPoly(3,dV=3.1,dW=c(rep(1.2,2),1)) + dlmModSeas(4), 
                        a.y = 1, b.y = 1000, a.theta = 1, 
                        b.theta = 1000,
                        n.sample = 1100, ind = c(2, 3),
                        save.states = FALSE)


# diagnostics

burn <- 100
attach(outGibbs)
attach(outGibbs1)
par(mfrow=c(2,3), mar=c(3.1,2.1,2.1,1.1))
plot(dV[-(1:burn)], type = 'l', xlab="", ylab="", main=expression(sigma^2))
plot(dW[-(1:burn),1], type = 'l', xlab="", ylab="", main=expression(sigma[beta]^2))
plot(dW[-(1:burn),2], type = 'l', xlab="", ylab="", main=expression(sigma[s]^2))
hist(dV[-(1:burn)], type = 'l', xlab="", ylab="", main=expression(sigma^2))
hist(dW[-(1:burn),1], type = 'l', xlab="", ylab="", main=expression(sigma[beta]^2))
hist(dW[-(1:burn),2], type = 'l', xlab="", ylab="", main=expression(sigma[s]^2))
# Comment on shape of posterior distribution

par(mfrow=c(1,3))
acf(sqrt(outGibbs$dV[-(1:burn)]), main="")
acf(sqrt(outGibbs$dW[-(1:burn),1]), main="Autocovariance")
acf(sqrt(outGibbs$dW[-(1:burn),2]), main="")


use <- length(dV) - burn
from <- 0.05 * use
at <- pretty(c(0,use),n=3); at <- at[at>=from]
plot(ergMean(dV[-(1:burn)], from), type="l", xaxt="n",xlab="", ylab="")
axis(1, at=at-from, labels=format(at))
plot(ergMean(dW[-(1:burn),1], from), type="l", xaxt="n",xlab="", ylab="")
axis(1, at=at-from, labels=format(at))
plot(ergMean(dW[-(1:burn),2], from), type="l", xaxt="n",xlab="", ylab="")
axis(1, at=at-from, labels=format(at))
detach()
# How stable are the ergodic means towards end?
# How non-stationary are the trace plots?


mcmcMean(cbind(outGibbs$dV[-(1:burn)], outGibbs$dW[-(1:burn), ]))

par(mfrow=c(1,3))
plot(dV[-(1:burn)],dW[-(1:burn),1],xlab=expression(sigma^2),ylab=expression(sigma[beta]^2))
plot(dV[-(1:burn)],dW[-(1:burn),2],xlab=expression(sigma^2),ylab=expression(sigma[s]^2),main="Bivariate Plots")
plot(dW[-(1:burn),1],dW[-(1:burn),2],xlab=expression(sigma[beta]^2),ylab=expression(sigma[s]^2))

quantile(dW[-burn], c(0.025, 0.975))



# RMD

```{r}
outGibbs <- dlmGibbsDIG(rapidmean.ts, dlmModPoly(2) + dlmModSeas(4), 
                        a.y = 1, b.y = 1000, a.theta = 1, 
                        b.theta = 1000,
                        n.sample = 1100, ind = c(2, 3),
                        save.states = FALSE)

```

The use of a Gibbs sampler to generate posterior estimates of the model parameters is a further method that may be used to obtain accurate results when using a more complex model. Below is an illustration of how altering the prior distributions of this Bayesian appoach can lead to more accurate results, taking the previously found optimal DLM. After observing the posterior distributions, the gamma priors for the new model were adjusted accordingly. These plots below give the ergodic means of the precision parameters for the process. 

```{r, warning=FALSE}
set.seed(123)
attach(outGibbs)
burn <- 100
par(mfrow=c(1,3), mar=c(3.1,4.1,2.1,1.1))
use <- length(dV) - burn
from <- 0.05 * use
at <- pretty(c(0,use),n=3); at <- at[at>=from]
plot(ergMean(dV[-(1:burn)], from), type="l", xaxt="n",xlab="", 
     ylab="Vague Priors", main=expression(sigma^2))
axis(1, at=at-from, labels=format(at))
plot(ergMean(dW[-(1:burn),1], from), type="l", xaxt="n",xlab="", 
     ylab="", main=expression(sigma[beta]^2))
axis(1, at=at-from, labels=format(at))
plot(ergMean(dW[-(1:burn),2], from), type="l", xaxt="n",xlab="", 
     ylab="", main=expression(sigma[s]^2))
axis(1, at=at-from, labels=format(at))
detach()
outGibbs1 <- dlmGibbsDIG(rapidmean.ts, dlmModPoly(2) + dlmModSeas(4), 
                         a.y = 5.9, b.y = 1000, a.theta = c(0.02,0.06), 
                         b.theta = 1000, shape.y = 12, rate.y = 1, shape.theta = 3, 
                         rate.theta = 2,
                         n.sample = 1100, ind = c(2, 3),
                         save.states = FALSE)
attach(outGibbs1)
use <- length(dV) - burn
from <- 0.05 * use
at <- pretty(c(0,use),n=3); at <- at[at>=from]
plot(ergMean(dV[-(1:burn)], from), 
     type="l", xaxt="n",xlab="", 
     ylab="Informative Priors", main=expression(sigma^2))
axis(1, at=at-from, labels=format(at))
plot(ergMean(dW[-(1:burn),1], from), 
     type="l", xaxt="n",xlab="", 
     ylab="", main=expression(sigma[beta]^2))
axis(1, at=at-from, labels=format(at))
plot(ergMean(dW[-(1:burn),2], from), type="l", 
     xaxt="n",xlab="", ylab="", main=expression(sigma[s]^2))
axis(1, at=at-from, labels=format(at))
mcmcMean(cbind(dV[-(1:burn)], dW[-(1:burn), ])) # Bayesian 
drop(V(rm.dlm)) # MLE
diag(W(rm.dlm))[2:3]

```

It can be observed that the process appears to have become slightly more stable at the later iterations, compared to the vague prior model. This implies that there is not much improvement in convergence of the updated model, suggesting that the most efficient method would be maximum likelihood for this model, which is not overly complex. The Gibbs sampler has resulted in a similar estimate to the maximum likelihood method for the observation precision, but the system precision estimates are inconsistent. 


