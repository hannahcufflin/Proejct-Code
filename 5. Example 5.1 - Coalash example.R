### Example 5.1 + surface plot Figure 2.3 
coal.ash = read.table('C:/Users/hanna/OneDrive/Desktop/PROJECT/DATA/coal.ash.txt', header = TRUE)
x.res=100
y.res=100
ash = coal.ash$coal
coords1 = cbind(coal.ash$x, coal.ash$y)
surf = mba.surf(cbind(coords1, ash), no.X=x.res, no.Y=y.res,n=1,m=1,h=5, extend=FALSE)
surf = surf$xyz.est
xr = range(coal.ash$x) + c(-0.01, 0.01)*diff(range(coal.ash$x))
yr = range(coal.ash$y) + c(-0.01, 0.01)*diff(range(coal.ash$y))

#Figure 2.3 
image.plot(surf, xlim = xr, ylim = yr, xlab = 'x', ylab = 'y')
contour(surf, add=T)

p = 1
n.samples <- 10000
m.1.exp <- spLM(coal~1, data=coal.ash, coords=coords1, starting=list("phi"=0.01,"sigma.sq"=0.01,
                                                                     "tau.sq"=0.01), tuning=list("phi"=0.1, "sigma.sq"=0.05, "tau.sq"=0.05),
                priors=list("phi.Unif"=c(3/1500, 3/50), "sigma.sq.IG"=c(3, 0.08),"tau.sq.IG"=c(2, 0.02)), cov.model="exponential",n.samples=n.samples)

m.1.sph <- spLM(coal~1, data=coal.ash, coords=coords1, starting=list("phi"=0.01,"sigma.sq"=0.01,
                                                                     "tau.sq"=0.01), tuning=list("phi"=0.1, "sigma.sq"=0.05, "tau.sq"=0.05),
                priors=list("phi.Unif"=c(3/1500, 3/50), "sigma.sq.IG"=c(3, 0.08),"tau.sq.IG"=c(2, 0.02)), cov.model="spherical",n.samples=n.samples)

m.1.mat <- spLM(coal~1, data=coal.ash, coords=coords1, starting=list("phi"=0.01,"sigma.sq"=0.01,
                                                                     "tau.sq"=0.01, 'nu' = 0.01), tuning=list("phi"=0.1, "sigma.sq"=0.05, "tau.sq"=0.05, 'nu' = 0.1),
                priors=list("phi.Unif"=c(3/1500, 3/50), "sigma.sq.IG"=c(3, 0.08),"tau.sq.IG"=c(2, 0.02), 'nu.Unif' = c(3/1500, 3/50)), cov.model="matern",n.samples=n.samples)

burn.in <- 0.75*n.samples


m.1.exp <- spRecover(m.1.exp, start=burn.in)
m.1.sph <- spRecover(m.1.sph, start=burn.in)
m.1.mat <- spRecover(m.1.mat, start=burn.in)
# Table 5.1(a)
round(summary(m.1.exp$p.theta.recover.samples)$quantiles[,c(1,3,5)],2)
round(summary(m.1.exp$p.beta.recover.samples)$quantiles[c(1,3,5)],2)
beta.samples1.exp = m.1.exp$p.beta.recover.samples
theta.samples1.exp = m.1.exp$p.theta.recover.samples
w.samples1.exp = m.1.exp$p.w.recover.samples
w.hat.mu.coal.exp <- apply(w.samples1.exp,1,mean)

#Table 5.1(b)
round(summary(m.1.sph$p.theta.recover.samples)$quantiles[,c(1,3,5)],2)
round(summary(m.1.sph$p.beta.recover.samples)$quantiles[c(1,3,5)],2)
beta.samples1.sph = m.1$p.beta.recover.samples
theta.samples1.sph = m.1$p.theta.recover.samples
w.samples1.sph = m.1$p.w.recover.samples
w.hat.mu.coal.sph <- apply(w.samples1.sph,1,mean)

#Table 5.1(c)
round(summary(m.1.mat$p.theta.recover.samples)$quantiles[,c(1,3,5)],2)
round(summary(m.1.mat$p.beta.recover.samples)$quantiles[c(1,3,5)],2)
beta.samples1.mat = m.1$p.beta.recover.samples
theta.samples1.mat = m.1$p.theta.recover.samples
w.samples1.mat = m.1$p.w.recover.samples
w.hat.mu.coal.mat <- apply(w.samples1.mat,1,mean)


# Spherical model

m.1 = m.1.sph


m <- nrow(coords1)

pred.X <- mkMvX(list(matrix(1,m,1)))


par(mfrow=c(1,1))
par(mar = c(6.1, 6.1, 3.1, 5.1))

nut.pred.coal <- spPredict(m.1, start=burn.in, thin=10, pred.coords=coords1, pred.covars=pred.X)
y.pred.mu.coal <- apply(nut.pred.coal$p.y.predictive.samples, 1, mean)
surf.pred.coal <- mba.surf(cbind(coords1,y.pred.mu.coal), no.X=x.res, no.Y=y.res, h=5, m=2, extend=FALSE)$xyz.est
z.lim <- range(surf[[3]], na.rm=TRUE)
image.plot(surf.pred.coal,zlim = c(7.1, 13.9) ,xaxs = "r", yaxs = "r", xlab="x", ylab="y")
contour(surf.pred.coal, add=T)


surf <- mba.surf(cbind(coords1, w.hat.mu.coal.sph), no.X=x.res, no.Y=y.res, extend=FALSE)$xyz.est
z.lim <- range(surf[[3]], na.rm=TRUE)
image.plot(surf, xaxs = "r", yaxs = "r", zlim=z.lim +c(-1,0), xlab = 'x', ylab = 'y')
contour(surf, add=T)


