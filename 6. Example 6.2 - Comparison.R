# Comparing data
library(spBayes)
library(MBA)
library(geoR)
library(fields)
library(plgp)
library(mcmc)
library(SpatialExtremes)
# Examlple 1 - Colorado data - actually the second example
colorado = read.table('C:/Users/hanna/OneDrive/Desktop/PROJECT/DATA/coloradolmc.dat', header = TRUE)
colorado$alt = c(2131,2773, 2564, 1806, 2660, 1551, 2338, 904, 3394, 1105, 1203, 2692, 2917, 2524, 1397, 2023, 2495, 2310, 1509, 1941, 3198, 1852, 2469, 3132, 2691, 2149, 2144, 2380, 1333, 2603, 1521, 1024, 2341, 2911, 1846, 1821, 2109, 3150, 3097, 2444, 1596, 1908, 2212, 1687, 1592, 1619, 2302, 1426, 1911, 2592 )
colorado$alt = colorado$alt/1000
coords = as.matrix(colorado[, c('Lat' ,'Lon')])
Ppt = colorado$Ppt
Temp = colorado$Temp
alt = colorado$alt

p = 2 # 
#no of columns in the design matrix = coefficient and alt

#Attempt at posterior samples

n.samples <- 10000
m.1 <- spLM(Temp~alt, data=colorado, coords=coords, starting=list("phi"=3/200,"sigma.sq"=0.08,
                                                                  "tau.sq"=0.02), tuning=list("phi"=0.1, "sigma.sq"=0.05, "tau.sq"=0.05),
            priors=list("phi.Unif"=c(3/1500, 3/50), "sigma.sq.IG"=c(2, 0.08),"tau.sq.IG"=c(2, 0.02)), cov.model="exponential",n.samples=n.samples)


m.2 <- spLM(Ppt~alt, data=colorado, coords=coords, starting=list("phi"=3/200,"sigma.sq"=0.08,
                                                                 "tau.sq"=0.02, 'beta' = 3), tuning=list("phi"=0.1, "sigma.sq"=0.05, "tau.sq"=0.05),
            priors=list("phi.Unif"=c(3/1500, 3/50), "sigma.sq.IG"=c(2, 0.08),"tau.sq.IG"=c(2, 0.02)), cov.model="exponential",n.samples=n.samples)

burn.in <- 0.75*n.samples

m.1 <- spRecover(m.1, start=burn.in) #spRecover recovers regression coefficients and spatial random effects

m.2 <- spRecover(m.2, start=burn.in)



m <- nrow(coords)

pred.X <- mkMvX(list(matrix(c(rep(1,50), colorado$alt), m, 2)))

nut.pred.temp <- spPredict(m.1, start=burn.in, thin=10, pred.coords=coords, pred.covars=pred.X)


nut.pred.ppt <- spPredict(m.2, start=burn.in, thin=10, pred.coords=coords, pred.covars=pred.X)

### DIC and posterior predictive loss for univariate model###

spDiag(m.1) # for temperature

#$DIC
#value
#bar.D       -154.65360
#D.bar.Omega -174.83579
#pD            20.18219
#DIC         -134.47141

#$GP for posterior predictive loss
#value
#G 0.4809984
#P 1.2705786
#D 1.7515769

spDiag(m.2) # For precipitation

#$DIC
#value
#bar.D       223.512668
#D.bar.Omega 220.432596
#pD            3.080072
#DIC         226.592740

#$GP
#value
#G 1504.065
#P 1597.452
#D 3101.517



### RMSPE for univariate ###
y.pred.mean.temp = apply(nut.pred.temp$p.y.predictive.samples, 1, mean)
y.observed.temp = Temp
y.difference.temp = y.pred.mean.temp - y.observed.temp
y.diff.sum.temp = sum(y.difference.temp^2)
RMSPE.temp = sqrt(1/50 * y.diff.sum.temp)
RMSPE.temp
# [1] 3.968234e-17

y.pred.mean.ppt = apply(nut.pred.ppt$p.y.predictive.samples, 1, mean)
y.observed.ppt = Ppt
y.difference.ppt = y.pred.mean.ppt - y.observed.ppt
y.diff.sum.ppt = sum(y.difference.ppt^2)
RMSPE.ppt = sqrt(1/50 * y.diff.sum.ppt)
RMSPE.ppt
# [1] 1.53837e-16

### NSME ###
y.observed.temp.diff = (Temp - mean(Temp))^2
y.observed.temp.diff.sum = sum(y.observed.temp.diff)
NSME.temp = 1 - (y.diff.sum.temp /y.observed.temp.diff.sum)
NSME.temp

#[1] 1


y.observed.ppt.diff = (Ppt - mean(Ppt))^2
y.observed.ppt.diff.sum = sum(y.observed.ppt.diff)
NSME.ppt = 1 - (y.diff.sum.ppt /y.observed.ppt.diff.sum)
NSME.ppt
#[1] 1


### 95% ALCI for univariate model  ####

nut.pred.temp.interval = apply(nut.pred.temp$p.y.predictive.samples, 1, quantile, probs = c (0.025, 0.975))
Y.temp  = nut.pred.temp.interval[2,]-nut.pred.temp.interval[1,]
ALCI.temp = mean(Y.temp)
ALCI.temp
# [1] 2.479683e-15

nut.pred.ppt.interval = apply(nut.pred.ppt$p.y.predictive.samples, 1, quantile, probs = c (0.025, 0.975))
Y.ppt = nut.pred.ppt.interval[2,]-nut.pred.ppt.interval[1,]
ALCI.ppt = mean(Y.ppt)
ALCI.ppt
# 4.396483e-15




colorado = read.table('C:/Users/hanna/OneDrive/Desktop/PROJECT/DATA/coloradolmc.dat', header = TRUE)
colorado$alt = c(2131,2773, 2564, 1806, 2660, 1551, 2338, 904, 3394, 1105, 1203, 2692, 2917, 2524, 1397, 2023, 2495, 2310, 1509, 1941, 3198, 1852, 2469, 3132, 2691, 2149, 2144, 2380, 1333, 2603, 1521, 1024, 2341, 2911, 1846, 1821, 2109, 3150, 3097, 2444, 1596, 1908, 2212, 1687, 1592, 1619, 2302, 1426, 1911, 2592 )
colorado$alt = colorado$alt/1000
coords = as.matrix(colorado[, c('Lat' ,'Lon')])
weather = colorado[, c('Temp', 'Ppt', 'alt')]



q=2
n.samples = 1000
n.ltr = q*(q+1)/2

A.starting <- diag(0.1,q)[lower.tri(diag(1,q), TRUE)]
Psi.starting <- rep(0.5, q)

A.tuning <- rep(0.0005,length(A.starting))
Psi.tuning <- rep(0.0005, q)

n.samples <- 10000

starting <- list("phi"=rep(3/20,q),"A"=A.starting,"Psi"=Psi.starting)
tuning <- list("phi"=rep(0.1,q),"A"=A.tuning, "Psi"=Psi.tuning)
priors <- list("phi.Unif"=list(rep(3/60,q), rep(3/10,q)), "K.IW"=list(q+1, diag(0.001,q)), "Psi.IG"=list(rep(2,q), c(0.05,0.08)))

m.1 <- spMvLM(list(Temp~alt,Ppt~alt), coords=coords, data=weather,
              starting=starting, tuning=tuning, priors=priors,
              cov.model="exponential",  n.samples=n.samples, n.report=2000)

m.1 <- spRecover(m.1, start=burn.in, thin=10)

m<- nrow(coords)

pred.X <- mkMvX(list(matrix(c(rep(1,50), colorado$alt), m, 2),matrix(c(rep(1,50), colorado$alt), m, 2) ))


nut.pred.colorado <- spPredict(m.1, start=burn.in, thin=10, pred.coords=coords, pred.covars=pred.X)
even_indicies = seq(2,100,2)
odd_indicies = seq(1,99,2)
y.pred.temp <- nut.pred.colorado$p.y.predictive.samples[odd_indicies,]
y.pred.ppt <- nut.pred.colorado$p.y.predictive.samples[even_indicies,]


### values in Table 6.4 
### DIC and predictive posterior loss ###
spDiag(m.1)
#value
#bar.D       72.83685
#D.bar.Omega 50.62251
#pD          22.21434
#DIC         95.05119

#$GP
#value
#G 1533.392
#P 1485.699
#D 3019.091


### RMSPE for multivariate ###
y.pred.mean.temp = apply(y.pred.temp, 1, mean)
y.observed.temp = Temp
y.difference.temp = y.pred.mean.temp - y.observed.temp
y = sum(y.difference.temp^2)
RMSPE.temp = sqrt(1/50 * y)
RMSPE.temp
# [1] 4.567036e-17

y.pred.mean.ppt = apply(y.pred.ppt, 1, mean)
y.observed.ppt = Ppt
y.difference.ppt = y.pred.mean.ppt - y.observed.ppt
y = sum(y.difference.ppt^2)
RMSPE.ppt = sqrt(1/50 * y)
RMSPE.ppt
#[1] 1.691041e-16


### NSME ###
y.difference.temp = y.pred.mean.temp - y.observed.temp
y = sum(y.difference.temp^2)
y.difference.mean.temp = sum((Temp - mean(Temp))^2)
NSME.temp = 1 - (y/y.difference.mean.temp)
NSME.temp
#[1] 1

y.difference.ppt = y.pred.mean.ppt - y.observed.ppt
y = sum(y.difference.ppt^2)
y.difference.mean.ppt = sum((Ppt- mean(Ppt))^2)
NSME.ppt = 1 - (y/y.difference.mean.ppt)
NSME.ppt
#[1] 1


### 95% ALCI for multivariate model ###

y.pred.ppt.interval = apply(y.pred.ppt,1,quantile, probs = c(0.025, 0.975))
ALCI.ppt = mean(y.pred.ppt.interval[2,]  - y.pred.ppt.interval[1,])
ALCI.ppt
# [1] 6.208367e-15

y.pred.temp.interval = apply(y.pred.temp,1,quantile, probs = c(0.025, 0.975))
ALCI.temp = mean(y.pred.temp.interval[2,]  - y.pred.temp.interval[1,])
ALCI.temp
#[1] 2.885747e-15








