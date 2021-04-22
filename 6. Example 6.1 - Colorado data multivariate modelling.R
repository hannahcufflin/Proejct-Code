### Example 6.1 - Multivariate modelling of Colorado Data ###
library(MBA)
library(fields)
library(geoR)
library(plgp)
library(spBayes)
library(mcmc)
colorado = read.table('C:/Users/hanna/OneDrive/Desktop/PROJECT/DATA/coloradolmc.dat', header = TRUE)
colorado$alt = c(2131,2773, 2564, 1806, 2660, 1551, 2338, 904, 3394, 1105, 1203, 2692, 2917, 2524, 1397, 2023, 2495, 2310, 1509, 1941, 3198, 1852, 2469, 3132, 2691, 2149, 2144, 2380, 1333, 2603, 1521, 1024, 2341, 2911, 1846, 1821, 2109, 3150, 3097, 2444, 1596, 1908, 2212, 1687, 1592, 1619, 2302, 1426, 1911, 2592 )
colorado$alt = colorado$alt/1000
coords = as.matrix(colorado[, c('Lat' ,'Lon')])
weather = colorado[, c('Temp', 'Ppt', 'alt')]

## Surface plots of Temperature and Precipitation - already shown in Example 5.2 ###
par(mfrow = c(1,1))
par(mar = c(5.5,5.5,2,2))
surf1 = mba.surf(cbind(coords, data = weather[1]), no.X = 100, no.Y = 100)$xyz.est
image.plot(surf1, xlab = 'latitude', ylab = 'longitude' )
points(coords)
surf2 = mba.surf(cbind(coords, data = weather[2]), no.X = 100, no.Y = 100)$xyz.est
image.plot(surf2, xlab = 'latitude', ylab = 'longitude' )
points(coords)


### Model set-up ###
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

#Multivariate model#
m.1 <- spMvLM(list(Temp~alt,Ppt~alt), coords=coords, data=weather,
              starting=starting, tuning=tuning, priors=priors,
              cov.model="exponential",  n.samples=n.samples, n.report=2000)



burn.in <- 0.75*n.samples

m.1 <- spRecover(m.1, start=burn.in) #spRecover recovers regression coefficients and spatial random effects

#Posterior summaries of the spatial parameters - Table 6.1
round(summary(m.1$p.theta.recover.samples)$quantiles[,c(1,3,5)],2)
#Posterior summaries of the beta coefficients - Table 6.1
round(summary(m.1$p.beta.recover.samples)$quantiles[,c(1,3,5)],2)



beta.samples = m.1$p.beta.recover.samples
theta.samples = m.1$p.theta.recover.samples
# Creating variables for the components of D
theta.samples.D1 = theta.samples[,4]
theta.samples.D2 = theta.samples[,5]




# trace plot of D[1,1] and D[2,2] - Figure 6.1
par(mfrow=c(1,1))
par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(theta.samples.D1, auto.layout=FALSE, density=FALSE, main = c('', ''), ylim = c(0.005,0.030)) #coefficient trace plot 
plot(theta.samples.D2, auto.layout=FALSE, density=FALSE, main = c('', ''), ylim = c(15,45)) #coefficient trace plot 

# Residuals of the linear model 
Temp.resids <- resid(lm(Temp~alt, data=weather))
Ppt.resids <- resid(lm(Ppt~alt, data=weather))

#The mean spatial effects of both variable
w <- rowMeans(m.1$p.w.recover.samples)
w.Temp <- w[seq(1,length(w),q)]
w.Ppt <- w[seq(2,length(w),q)]

res <- 100
par(mfrow=c(2,2))
par(mfrow = c(1,1))
par(mar = c(6.1, 6.1, 5.1, 3.1))
surf.r <- mba.surf(cbind(coords,Temp.resids), no.X=res, no.Y=res, extend=FALSE)$xyz.est
surf.w <- mba.surf(cbind(coords,w.Temp), no.X=res, no.Y=res, extend=FALSE)$xyz.est


# Surface plots of the mean spatial effects and lm residual errors for Temperature - Figure 6.2 (a) and (b)
image.plot(surf.r, zlim=c(-0.5,0.6), xlab = 'Latitude', ylab = 'Longitude') # Temperature lm residuals
points(coords)
image.plot(surf.w,zlim = c(-0.6,0.55), xlab = 'Latitude', ylab = 'Longitude') # temp error residuals 
points(coords) 

surf.r <- mba.surf(cbind(coords,Ppt.resids), no.X=res, no.Y=res, extend=FALSE)$xyz.est
surf.w <- mba.surf(cbind(coords,w.Ppt), no.X=res, no.Y=res, extend=FALSE)$xyz.est
z.lim <- range(c(surf.r[["z"]], surf.w[["z"]]), na.rm=TRUE)

# Surface plots of the mean spatial effects and lm residual errors for Precipitation - Figure 6.2 (c) and (d)
image.plot(surf.r, zlim=z.lim, xlab = 'Latitude', ylab = 'Longitude') #Precipitation spatial residuals
points(coords) # error plot of non spatial residual 
image.plot(surf.w, zlim=c(-0.07,0.07), xlab = 'Latitude', ylab = 'Longitude') #Precipitation error residuals
points(coords) # error plot of spatial residuals


# Predicting for data variables
m<- nrow(coords)

pred.X <- mkMvX(list(matrix(c(rep(1,50), colorado$alt), m, 2),matrix(c(rep(1,50), colorado$alt), m, 2) ))


nut.pred.colorado <- spPredict(m.1, start=burn.in, thin=10, pred.coords=coords, pred.covars=pred.X)
even_indicies = seq(2,100,2)
odd_indicies = seq(1,99,2)
y.pred.mu.colorado.temp <- apply(nut.pred.colorado$p.y.predictive.samples[odd_indicies,], 1, mean)
y.pred.mu.colorado.ppt <- apply(nut.pred.colorado$p.y.predictive.samples[even_indicies,], 1, mean)

# Surface plot of mean sampled predicted values for temperature - Figure 6.3(a) 
surf.pred.colorado.temp <- mba.surf(cbind(coords,y.pred.mu.colorado.temp), no.X=100, no.Y=100, h=5, m=2, extend=FALSE)$xyz.est
z.lim <- range(surf[[3]], na.rm=TRUE)
image.plot(surf.pred.colorado.temp ,xaxs = "r", yaxs = "r", xlab="Latitude", ylab="Longitude") # predicted values plot 
points(coords, cex = 0.9)

# Surface plot of mean sampled predicted values for precipitation - Figure 6.3(b) 
surf.pred.colorado.ppt <- mba.surf(cbind(coords,y.pred.mu.colorado.ppt), no.X=100, no.Y=100, h=5, m=2, extend=FALSE)$xyz.est
z.lim <- range(surf[[3]], na.rm=TRUE)
image.plot(surf.pred.colorado.ppt ,xaxs = "r", yaxs = "r", xlab="Latitude", ylab='Longitude') # predicted values plot 
points(coords, cex =0.9)


### PREDICTION AT THE NEW COORDINATES ###

### additional prediction coordinates ###
# (41, -109,2506) (41, -102, 1063) (41, -103, 1319) (40, -104, 1455) (39, -108, 3058) (39,-106,2806 ) (38, -108, 2766) (38, -106, 2313) (38, -102, 1018) (37, -102, 1101) # 10 new prediction coordinates  

lat.new = c(41, 41, 41, 40, 39, 39, 38, 38, 38, 37)
long.new = c(-109, -102, -103, -104, -108,-106, -108, -106, -102, -102)
alt.new = c(2506, 1063, 1319, 1455,  3058, 2806, 2766, 2313, 1018, 1101)/1000
colorado.new = as.data.frame(cbind((lat.new), (long.new), (alt.new)), col.names = names('Lat', 'Long', 'Alt'))
colnames(colorado.new) = c('Lat', 'Long', 'Alt')
pred.coords = colorado.new[,c('Lat', 'Long')]

m <- nrow(pred.coords)


pred.X <- mkMvX(list(matrix(c(rep(1,10), colorado.new$Alt), m, 2),matrix(c(rep(1,10), colorado.new$Alt), m, 2) ))


nut.pred <- spPredict(m.1, start=burn.in, thin=10, pred.coords=pred.coords, pred.covars=pred.X)

y.pred.mu <- apply(nut.pred$p.y.predictive.samples, 1, mean)

odd_indicies = seq(1,19,by=2)
even_indicies = seq(2,20,by=2)

y.pred.mu.temp = y.pred.mu[odd_indicies]
y.pred.mu.ppt = y.pred.mu[even_indicies]

total.lat = c(colorado$Lat, colorado.new$Lat)
total.long = c(colorado$Lon, colorado.new$Long)
total.coords = as.data.frame(cbind(total.lat, total.long))
total.temp = c(colorado$Temp, y.pred.mu.temp)
total.ppt = c(colorado$Ppt, y.pred.mu.ppt)

# Surface plot of new values + observed values for Temperature - Figure 6.4(a)
surf <- mba.surf(cbind(total.coords,total.temp), no.X= 100, no.Y=100, extend=FALSE)$xyz.est
image.plot(surf, xlab = 'Latitude', ylab = 'Longitude')
points(coords, cex=0.9)
points(pred.coords, pch= 19, cex = 0.9)

# Surface plot of new values + observed values for Precipitation - Figure 6.4(b)
surf <- mba.surf(cbind(total.coords,total.ppt), no.X= 100, no.Y=100, extend=FALSE)$xyz.est
image.plot(surf, xlab = 'Latitude', ylab = 'Longitude')
points(coords, cex=0.9)
points(pred.coords, pch= 19, cex = 0.9)


