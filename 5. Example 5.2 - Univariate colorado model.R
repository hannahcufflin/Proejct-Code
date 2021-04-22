library(spBayes)
library(MBA)
library(geoR)
library(fields)
# Examlple 1 - Colorado data - actually the second example
colorado = read.table('C:/Users/hanna/OneDrive/Desktop/PROJECT/DATA/coloradolmc.dat', header = TRUE)
colorado$alt = c(2131,2773, 2564, 1806, 2660, 1551, 2338, 904, 3394, 1105, 1203, 2692, 2917, 2524, 1397, 2023, 2495, 2310, 1509, 1941, 3198, 1852, 2469, 3132, 2691, 2149, 2144, 2380, 1333, 2603, 1521, 1024, 2341, 2911, 1846, 1821, 2109, 3150, 3097, 2444, 1596, 1908, 2212, 1687, 1592, 1619, 2302, 1426, 1911, 2592 )
colorado$alt = colorado$alt/1000
coords = as.matrix(colorado[, c('Lat' ,'Lon')])
Ppt = colorado$Ppt
Temp = colorado$Temp
alt = colorado$alt

# Surface plot of temperaute precipitation and altitude

x.res <- 100; y.res <- 100
par(mar = c(6.1, 6.1, 5.1, 3.1))
# Figure 5.2(a)
surf1 <- mba.surf(cbind(coords,Ppt ), no.X=x.res, no.Y=y.res, h=5, m=2, extend=FALSE)$xyz.est
image.plot(surf1, xaxs = "r", yaxs = "r", xlab="Latitude", ylab="Longitude")
points(coords)
# Figure 5.2(b)
surf2 <- mba.surf(cbind(coords, Temp), no.X=x.res, no.Y=y.res, h=5, m=2, extend=FALSE)$xyz.est
image.plot(surf2, xaxs = "r", yaxs = "r", xlab="Latitude", ylab="Longitude")
points(coords)
# Figure 5.2(c)
surf3 <- mba.surf(cbind(coords,alt), no.X=x.res, no.Y=y.res, h=5, m=2, extend=FALSE)$xyz.est
image.plot(surf3, xaxs = "r", yaxs = "r", xlab="Latitude", ylab="Longitude")
points(coords)

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


# Table 5.2(a)
round(summary(m.1$p.theta.recover.samples)$quantiles[,c(1,3,5)],2)
round(summary(m.1$p.beta.recover.samples)$quantiles[,c(1,3,5)],2)

# Table 5.2(b)
round(summary(m.2$p.theta.recover.samples)$quantiles[,c(1,3,5)],2)
round(summary(m.2$p.beta.recover.samples)$quantiles[,c(1,3,5)],2)




## The posterior samples of the regression coefficients and the spatial effects can then be obtained as
beta.samples1 = m.1$p.beta.recover.samples
w.samples1 = m.1$p.w.recover.samples
beta.samples2 = m.2$p.beta.recover.samples
w.samples2 = m.2$p.w.recover.samples



## Obtain trace plots for regression coefficients
par(mfrow=c(1,1))
par(mar = c(5.1, 4.1, 4.1, 2.1))
# Figure 5.3(a) and Figure 5.3(b)
plot(beta.samples1, auto.layout=FALSE, density=FALSE, main = c('', '')) #coefficient trace plot
#Figure 5.3(c) and Figure 5.3(d)
plot(beta.samples2, auto.layout=FALSE, density=FALSE, main = c('', '')) #coefficient trace plot 
 

## Obtain posterior means and sd's of spatial residuals for each location
w.hat.mu.temp <- apply(w.samples1,1,mean)
w.hat.sd.temp <- apply(w.samples1,1,sd)
w.hat.mu.ppt <- apply(w.samples2,1,mean)
w.hat.sd.ppt <- apply(w.samples2,1,sd)


par(mfrow=c(1,1))
par(mar = c(6.1, 6.1, 5.1, 3.1))

#Figure 5.4(a)
surf <- mba.surf(cbind(coords, w.hat.mu.temp), no.X=x.res, no.Y=y.res, extend=FALSE)$xyz.est
z.lim <- range(surf[[3]], na.rm=TRUE)
image.plot(surf, xaxs = "r", yaxs = "r", zlim=z.lim, xlab = 'Latitude', ylab = 'Longitude')
points(coords) # mean spatial effects 

#Figure 5.4(b)
surf <- mba.surf(cbind(coords, w.hat.mu.ppt), no.X=x.res, no.Y=y.res, extend=FALSE)$xyz.est
z.lim <- range(surf[[3]], na.rm=TRUE)
image.plot(surf, xaxs = "r", yaxs = "r", zlim=z.lim, xlab = 'Latitude', ylab = 'Longitude')
points(coords) #samples of the spatial surface error

#plot predicted value surface plots will be useful to do for multivariate - prediction at each observed site

m <- nrow(coords)

pred.X <- mkMvX(list(matrix(c(rep(1,50), colorado$alt), m, 2)))

nut.pred.temp <- spPredict(m.1, start=burn.in, thin=10, pred.coords=coords, pred.covars=pred.X)
y.pred.mu.temp <- apply(nut.pred.temp$p.y.predictive.samples, 1, mean)

#Figure 5.5(a)
surf.pred.temp <- mba.surf(cbind(coords,y.pred.mu.temp), no.X=x.res, no.Y=y.res, h=5, m=2, extend=FALSE)$xyz.est
image.plot(surf.pred.temp, xaxs = "r", yaxs = "r", xlab="Latitude", ylab="Longitude")
points(coords, cex = 0.9) # predicted values of temperature


nut.pred.ppt <- spPredict(m.2, start=burn.in, thin=10, pred.coords=coords, pred.covars=pred.X)
y.pred.mu.ppt <- apply(nut.pred.ppt$p.y.predictive.samples, 1, mean)
# Figure 5.5(b)
surf.pred.ppt <- mba.surf(cbind(coords,y.pred.mu.ppt), no.X=x.res, no.Y=y.res, h=5, m=2, extend=FALSE)$xyz.est
image.plot(surf.pred.ppt, xaxs = "r", yaxs = "r", xlab="Latitude", ylab="Longitude")
points(coords, cex = 0.9) #  predicted values of precipitation









