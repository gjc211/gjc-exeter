# 1.
require(geoR)
set.seed(362)
data <- read.csv('GrandBanks_Dec_1997.csv')
gdata <- as.geodata(data,coords.col=2:3,data.col=6)
dup <- dup.coords(gdata)
gdata2 <- jitterDupCoords(gdata,max=0.1,min=0.05)

par(mfrow=c(1,1))
plot(gdata2,trend="1st")
# Some outlying values of x and y coordinates 
# One of the quartiles show a more obvious non-linear pattern 

# 2. 
variog <- variog4(gdata2,trend="1st")
plot(variog, omni = TRUE)
# With the trend line added, which represents 0 direction variogram
# Indicates possible anisotropy due to 90 degree semivariogram appearing different 

# 3.
mle <- likfit(coords=gdata2$coords,data=gdata2$data, nugget = 0, ini=c(40,0.5), trend="1st", cov.model='exponential',kappa=3/2)
mle1 <- likfit(coords=gdata2$coords,data=gdata2$data, nugget = 0, ini=c(40,0.5), trend="1st", cov.model='sph',kappa=3/2)
mle2 <- likfit(coords=gdata2$coords,data=gdata2$data, nugget = 0, ini=c(40,0.5), trend="1st",cov.model='powered.exponential',kappa=3/2)
mle3 <- likfit(coords=gdata2$coords,data=gdata2$data, nugget = 0, ini=c(30,0.1), trend="1st",cov.model='circular',kappa=3/2)
mle4 <- likfit(coords=gdata2$coords,data=gdata2$data, nugget = 0, ini=c(30,0.1), trend="1st", cov.model='matern',kappa=5/2)
mle5 <- likfit(coords=gdata2$coords,data=gdata2$data, nugget = 3, ini=c(30,0.1), trend="1st", cov.model='gaussian')

### Take out bad data

baddata <- data[-c(366,664), ]
gdatanew <- as.geodata(baddata,coords.col=2:3,data.col=6)
dupnew <- dup.coords(gdatanew)
gdata2new <- jitterDupCoords(gdatanew,max=0.1,min=0.05)
mlenew <- likfit(coords=gdata2new$coords,data=gdata2new$data, nugget = 0, ini=c(40,0.5), trend="1st", cov.model='exponential',kappa=3/2)

plot(variog(gdata2,trend="1st"))
lines(mle,col="red")
lines(mle1,col="pink")
lines(mle2,col="blue")
lines(mle3,col="purple")
lines(mlenew,lty=5)
legend(0,30, legend=c("Exponential","Spherical","Powered Exponential","Circular","Exponential Without Outliers"),
       col=c("red","pink","blue","purple","black"), lty=1, cex=0.8)
legend(0,32.5,c(" "," "),title="Covariance Function",cex=1, bty='n')

# Plotting predicted values and variance
grid <- pred_grid(gdata2$coords,by = 0.25)
grid1 <- pred_grid(gdata2$coords,by = 0.5)
grid2 <- pred_grid(gdata2$coords,by = 0.35)

kc <- krige.conv(gdata2, loc = grid, krige = krige.control(obj.m = mle))
par(mfrow=c(1,2))
image(kc, loc = grid, xlab = "Coord X", ylab = "Coord Y", main="Predicted Values")
contour(kc, add = TRUE, nlevels = 10)
image(kc, val = kc$krige.var, loc = grid, xlab = "Coord X", ylab = "Coord Y", main="Prediction Variance")
contour(kc, val = kc$krige.var, add=TRUE, nlevels = 10)

# Cross-validation
xv <- xvalid(gdata2, model = mle)
par(mfcol = c(5, 2), mar = c(3, 3, 1, 0.5), mgp = c(1.5, 0.7, 0)) 
plot(xv)


# 4.

# Kriging explanation: using the variogram to interpolate data - finding the best linear unbiased predictor 
# for data point at each location based on the other data at that point in space (show forumula from L5)
# Phi is used as the prior estimator - component of all possible variograms

#head(parana)
#parana.gr <- pred_grid(parana$borders, by = 15)


ex.bayes <- krige.bayes(gdata2, coords=gdata2$coords, data=gdata2$data, loc=grid,
                        model = model.control(cov.m="exponential", kappa=2, trend.l = "1st", trend.d = "1st"),
                        prior = prior.control(phi = 1,
                                              phi.prior="exponential"))


# Plotting prediction results
op <- par(no.readonly = TRUE)
par(mfrow=c(2,2), mar=c(4,4,2.5,0.5), mgp = c(2,1,0))
image(ex.bayes, main="Predicted values")
image(ex.bayes, val="variance", main="Prediction variance")
image(ex.bayes, val= "simulation", number.col=1,
      main="Simulation from the \npredictive distribution (1)")
image(ex.bayes, val= "simulation", number.col=2,
      main="Simulation from \nthe predictive distribution (2)")

# Alternative

ex.bayes1 <- krige.bayes(gdata2, coords=gdata2$coords, data=gdata2$data, loc=grid,
                        model = model.control(cov.m="exponential", kappa=2, trend.l = "1st", trend.d = "1st"),
                        prior = prior.control(phi.prior = "rec"))
                                              
par(mfrow=c(2,2), mar=c(4,4,2.5,0.5), mgp = c(2,1,0))
image(ex.bayes1, main="predicted values")
image(ex.bayes1, val="variance", main="prediction variance")
image(ex.bayes1, val= "simulation", number.col=1,
      main="Simulation from the \npredictive distribution (1)")
image(ex.bayes1, val= "simulation", number.col=2,
      main="Simulation from \nthe predictive distribution (2)")



# Trying other models with a phi posterior distribution and nugget 

ex.bayes2 <- krige.bayes(gdata2, coords=gdata2$coords, data=gdata2$data, loc=grid1,
                         model = model.control(cov.m="gaussian", kappa=5/2, trend.l = "1st", trend.d = "1st", ),
                         prior = prior.control(phi.prior = "rec",phi.discrete = seq(0,1,l=51),
                                               tausq.rel.prior = "uni", tausq.rel.discrete = seq(3,4,l=6)))

par(mfrow=c(2,2), mar=c(4,4,2.5,0.5), mgp = c(2,1,0))
image(ex.bayes2, main="Predicted values")
image(ex.bayes2, val="variance", main="Prediction variance")
image(ex.bayes2, val= "simulation", number.col=1,
      main="Simulation from the \npredictive distribution (1)")
image(ex.bayes2, val= "simulation", number.col=2,
      main="Simulation from \nthe predictive distribution (2)")

par(mfrow=c(2,1))
hist(ex.bayes2)

###

ex.bayes3 <- krige.bayes(gdata2,coords=gdata2$coords, data=gdata2$data, loc=grid1,
                         model = model.control(cov.m="matern", kappa=5/2, trend.l = "1st", trend.d = "1st", ),
                         prior = prior.control(phi.prior = "rec",phi.discrete = seq(0,1,l=21),
                                               tausq.rel.prior = "uni", tausq.rel.discrete = seq(3,4,l=6)))

quants <- function(x){quantile(x, prob = c(0.05, 0.5, 0.95))}
lines(ex.bayes3, summ = quants, ty="l", lty=c(2,1,2), col=1)

### 

ex.bayes4 <- krige.bayes(gdata2,coords=gdata2$coords, data=gdata2$data, loc=grid1,
                         model = model.control(cov.m="powered.exponential", kappa=3/2, trend.l = "1st", trend.d = "1st", ),
                         prior = prior.control(phi.prior = "uniform",phi.discrete = seq(0,1,l=21),
                                               tausq.rel.prior = "fixed", tausq.rel = 3))

### Trying more data

ex.bayes5 <- krige.bayes(gdata2,coords=gdata2$coords, data=gdata2$data, loc=grid2,
                         model = model.control(cov.m="matern", kappa=3/2, trend.l = "1st", trend.d = "1st", ),
                         prior = prior.control(phi.prior = "uniform",phi.discrete = seq(0,1,l=21),
                                               tausq.rel.prior = "fixed", tausq.rel = 3))

###
ex.bayes7 <- krige.bayes(gdata2,coords=gdata2$coords, data=gdata2$data, loc=grid2,
                         model = model.control(cov.m="powered.exponential", kappa=3/2, trend.l = "1st", trend.d = "1st", ),
                         prior = prior.control(phi.prior = "squared.reciprocal",phi.discrete = seq(0,1,l=21),
                                               tausq.rel.prior = "uni", tausq.rel.discrete = seq(3,4,l=6)))

# Prior and posterior plots
par(mfrow=c(2,2))
plot(ex.bayes7, type = "h", tausq.rel = FALSE, col=c("red", "blue"))
title(main="Powered Exponential")
plot(ex.bayes2, type = "h", tausq.rel = FALSE, col=c("red", "blue"))
title(main="Gaussian")
plot(ex.bayes3, type = "h", tausq.rel = FALSE, col=c("red", "blue")) 
title(main="Matern 5/2")
plot(ex.bayes5, type = "h", tausq.rel = FALSE, col=c("red", "blue"))
title(main="Matern 3/2")
legend("topright", legend=c("Prior", "Posterior"),
       col=c("red", "blue"),lty=1,xpd = T, inset = c(0,-0.5))

# Explain about positive definite 
# Show how priors influence posteriors - different from normal kriging 

