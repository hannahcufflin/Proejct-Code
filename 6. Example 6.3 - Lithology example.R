### LITHOLOGY EXAMPLE  - EXAMPLE 6.3###
library(MBA)
library(fields)
library(geoR)
library(plgp)
library(spBayes)
library(mcmc)

## DATA ##
lithology = read.csv('C:/Users/hanna/OneDrive/Desktop/PROJECT/DATA/lithology data set edited.csv', header = TRUE, nrows = 85)
lith = data.frame(lithology)
lith = lith[-c(31, 45, 46,47, 71, 72, 73 ), ] # remove outlying data, far awat from other point
lith = data.frame(lith)
coords = as.matrix(lith[,c('Easting', 'ï..Northing')])
coords = coords/10000
X = as.matrix(lith[,c('Surf.Elevation', 'A.B.Elevation', 'B.C.Elevation', 'C.D.Elevation')])
X = as.data.frame(X)
Y = as.matrix(lith[,c('Thickness', 'Thickness.1', 'Thickness.2', 'Thickness.3')])

## SURFACE PLOTS ## - Figure 6.5
par(mar = c(6.1, 6.1, 5.1, 3.1))
surf1 <- mba.surf(cbind(coords,X[,1]), no.X=x.res, no.Y=y.res, h=5, m=2, extend=FALSE)$xyz.est
image.plot(surf1, xlab = 'Eastings', ylab = 'Northings' )
contour(surf1, add=T)

surf2 <- mba.surf(cbind(coords,X[,2]), no.X=x.res, no.Y=y.res, h=5, m=2, extend=FALSE)$xyz.est
image.plot(surf2, xlab = 'Eastings', ylab = 'Northings' )
contour(surf2, add=T)

surf3 <- mba.surf(cbind(coords,X[,3]), no.X=x.res, no.Y=y.res, h=5, m=2, extend=FALSE)$xyz.est
image.plot(surf3, xlab = 'Eastings', ylab = 'Northings' )
contour(surf3, add=T)

surf4 <- mba.surf(cbind(coords,X[,4]), no.X=x.res, no.Y=y.res, h=5, m=2, extend=FALSE)$xyz.est
image.plot(surf4, xlab = 'Eastings', ylab = 'Northings' )
contour(surf4, add=T)




### Mulitvariate model including surface 2 ###
X = as.matrix(lith[,c('Surf.Elevation', 'A.B.Elevation', 'B.C.Elevation', 'C.D.Elevation')])
X = as.data.frame(X)
q=4
n.samples = 1000
n.ltr = q*(q+1)/2

A.starting <- diag(c(15,3),q)[lower.tri(diag(1,q), TRUE)]
Psi.starting <- rep(2.5, q)

A.tuning <- rep(0.5,length(A.starting))
Psi.tuning <- rep(0.05, q)

n.samples <- 10000

starting <- list("phi"=rep(3/2,q),"A"=A.starting,"Psi"=Psi.starting, 'beta' = rep(0,6))
tuning <- list("phi"=rep(1,q),"A"=A.tuning, "Psi"=Psi.tuning)
priors <- list("phi.Unif"=list(rep(3/6,q), rep(3,q)), "K.IW"=list(q+1, diag(0.01,q)), "Psi.IG"=list(rep(2,q), c(0.5,0.8, 0.5, 0.5)))
### model with no covariates
m.1 <- spMvLM(list(Surf.Elevation ~ 1,A.B.Elevation ~ 1, B.C.Elevation~1, C.D.Elevation~1), coords=coords, data=X,
              starting=starting, tuning=tuning, priors=priors,
              cov.model="exponential",  n.samples=n.samples, n.report=2000)


burn.in <- 0.75*n.samples

m.1 <- spRecover(m.1, start=burn.in) #spRecover recovers regression coefficients and spatial random effects

## Table 6.3(a) values 
round(summary(m.1$p.theta.recover.samples)$quantiles[,c(1,3,5)],3)
round(summary(m.1$p.beta.recover.samples)$quantiles[,c(1,3,5)],3)

spDiag(m.1) #- calculating DIC and D
### Model metrics ###
#$DIC
#value
#bar.D       762.2116
#D.bar.Omega 549.1412
#pD          213.0704
#DIC         975.2821

#$GP
#value
#G  771.0486
#P 3034.9081
#D 3805.9568

m<- nrow(coords)

pred.X <- mkMvX(list(matrix(1,m,1),matrix(1,m,1),matrix(1,m,1),matrix(1,m,1)))


nut.pred.lith <- spPredict(m.1, start=burn.in, pred.coords=coords, pred.covars=pred.X)
indicies.1 = seq(1, 309,4)
indicies.2 = seq(2,310, 4)
indicies.3 = seq(3,311, 4)
indicies.4 = seq(4, 312, 4)
y.pred.mu.lith.surf1 <- apply(nut.pred.lith$p.y.predictive.samples[indicies.1,], 1, mean)
y.pred.mu.lith.surf2 <- apply(nut.pred.lith$p.y.predictive.samples[indicies.2,], 1, mean)
y.pred.mu.lith.surf3 = apply(nut.pred.lith$p.y.predictive.samples[indicies.3,], 1, mean)
y.pred.mu.lith.surf4 = apply(nut.pred.lith$p.y.predictive.samples[indicies.4,], 1, mean)

## Model metrics for individual surface # do not need to do this as can show the predictve sample - actual sample is zero so they are all zero - accurate predicting ## 
##Surface 1##
#RMSPE
y.difference.surf1= y.pred.mu.lith.surf1 - X[,1]
y = sum(y.difference.surf1^2)
RMSPE.surf1 = sqrt(1/78 * y)
RMSPE.surf1
# [1]0
y.difference.surf2= y.pred.mu.lith.surf2 - X[,2]
y = sum(y.difference.surf2^2)
RMSPE.surf2 = sqrt(1/78 * y)
RMSPE.surf2
#[1] 0
y.difference.surf3= y.pred.mu.lith.surf3 - X[,3]
y = sum(y.difference.surf3^2)
RMSPE.surf3 = sqrt(1/78 * y)
RMSPE.surf3
# [1] 0
y.difference.surf4= y.pred.mu.lith.surf4 - X[,4]
y = sum(y.difference.surf4^2)
RMSPE.surf4 = sqrt(1/78 * y)
RMSPE.surf4


#NSME = 0 for all as RMSPE is 0

### 95% ALCI for multivariate model 1 ###

y.pred.surf1.interval = apply(nut.pred.lith$p.y.predictive.samples[indicies.1,],1,quantile, probs = c(0.025, 0.975))
ALCI.surf1 = mean(y.pred.surf1.interval[2,]  - y.pred.surf1.interval[1,])
ALCI.surf1
# [1] 0

y.pred.surf2.interval = apply(nut.pred.lith$p.y.predictive.samples[indicies.2,],1,quantile, probs = c(0.025, 0.975))
ALCI.surf2 = mean(y.pred.surf2.interval[2,]  - y.pred.surf2.interval[1,])
ALCI.surf2
#[1] 0

y.pred.surf3.interval = apply(nut.pred.lith$p.y.predictive.samples[indicies.3,],1,quantile, probs = c(0.025, 0.975))
ALCI.surf3 = mean(y.pred.surf3.interval[2,]  - y.pred.surf3.interval[1,])
ALCI.surf3
#[1] 0

y.pred.surf4.interval = apply(nut.pred.lith$p.y.predictive.samples[indicies.4,],1,quantile, probs = c(0.025, 0.975))
ALCI.surf4 = mean(y.pred.surf4.interval[2,]  - y.pred.surf4.interval[1,])
ALCI.surf4
#[1] 0





## Multivariate model not including surface 2  ##
X = as.data.frame(X)
X2 = X[,-2]
q=3
n.samples = 1000
n.ltr = q*(q+1)/2

A.starting <- diag(c(15,15,3),q)[lower.tri(diag(1,q), TRUE)]
Psi.starting <- rep(2.5, q)

A.tuning <- rep(0.5,length(A.starting))
Psi.tuning <- rep(0.05, q)

n.samples <- 10000

starting <- list("phi"=rep(3/2,q),"A"=A.starting,"Psi"=Psi.starting, 'beta' = rep(0,6))
tuning <- list("phi"=rep(1,q),"A"=A.tuning, "Psi"=Psi.tuning)
priors <- list("phi.Unif"=list(rep(3/6,q), rep(3,q)), "K.IW"=list(q+1, diag(0.01,q)), "Psi.IG"=list(rep(2,q), c(0.5,0.5, 0.5)))
### model with no covariates
m.2 <- spMvLM(list(Surf.Elevation ~ 1, B.C.Elevation~1, C.D.Elevation~1), coords=coords, data=X2,
              starting=starting, tuning=tuning, priors=priors,
              cov.model="exponential",  n.samples=n.samples, n.report=2000)
burn.in <- 0.75*n.samples


m.2 <- spRecover(m.2, start=burn.in) #spRecover recovers regression coefficients and spatial random effects

## Table 6.3(b) values
round(summary(m.2$p.theta.recover.samples)$quantiles[,c(1,3,5)],3)
round(summary(m.2$p.beta.recover.samples)$quantiles[,c(1,3,5)],3)

## Model metrics for the entire model ##
spDiag(m.2)

#$DIC
#value
#bar.D       306.4815
#D.bar.Omega 135.7760
#pD          170.7055
#DIC         477.1871

#$GP
#value
#G  633.0506
#P 1587.4955
#D 2220.5461


m<- nrow(coords)

pred.X <- mkMvX(list(matrix(1,m,1),matrix(1,m,1),matrix(1,m,1)))


nut.pred.lith <- spPredict(m.2, start=burn.in, pred.coords=coords, pred.covars=pred.X)
indicies.1 = seq(1, 232,3)
indicies.2 = seq(2,233, 3)
indicies.3 = seq(3,234, 3)

y.pred.mu.lith.surf1 <- apply(nut.pred.lith$p.y.predictive.samples[indicies.1,], 1, mean)
y.pred.mu.lith.surf2 <- apply(nut.pred.lith$p.y.predictive.samples[indicies.2,], 1, mean)
y.pred.mu.lith.surf3 = apply(nut.pred.lith$p.y.predictive.samples[indicies.3,], 1, mean)

# Values in table 6.4
y.difference.surf1= y.pred.mu.lith.surf1 - X2[,1]
y = sum(y.difference.surf1^2)
RMSPE.surf1 = sqrt(1/78 * y)
RMSPE.surf1
#[1] 0

y.difference.surf2= y.pred.mu.lith.surf2 - X2[,2]
y = sum(y.difference.surf2^2)
RMSPE.surf2 = sqrt(1/78 * y)
RMSPE.surf2
# [1] 0

y.difference.surf3= y.pred.mu.lith.surf3 - X2[,3]
y = sum(y.difference.surf3^2)
RMSPE.surf3 = sqrt(1/78 * y)
RMSPE.surf3
# [1] 0

###NSME values are 1 based on values of RMSPE

y.pred.surf1.interval = apply(nut.pred.lith$p.y.predictive.samples[indicies.1,],1,quantile, probs = c(0.025, 0.975))
ALCI.surf1 = mean(y.pred.surf1.interval[2,]  - y.pred.surf1.interval[1,])
ALCI.surf1
# [1] 0

y.pred.surf2.interval = apply(nut.pred.lith$p.y.predictive.samples[indicies.2,],1,quantile, probs = c(0.025, 0.975))
ALCI.surf2 = mean(y.pred.surf2.interval[2,]  - y.pred.surf2.interval[1,])
ALCI.surf2
#[1] 0

y.pred.surf3.interval = apply(nut.pred.lith$p.y.predictive.samples[indicies.3,],1,quantile, probs = c(0.025, 0.975))
ALCI.surf3 = mean(y.pred.surf3.interval[2,]  - y.pred.surf3.interval[1,])
ALCI.surf3
#[1] 0



