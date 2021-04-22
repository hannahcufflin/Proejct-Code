# Gaussian process realizations with different covariance functions
library(mvtnorm)
library(plgp)
par(mar=c(5.1,4.1,4.1,2.1))
# generating values

n=100
X = matrix(seq(0,5, length = n), ncol = 1)
D = sqrt(distance(X))
eps = sqrt(.Machine$double.eps)

# plotting for different values of l:

#squared exponential - Figure 3.3(a)


se.1 = exp(-(D^2/(2*(2^2)))) # l = 2
se.2 = exp(-(D^2/(2*(0.5^2)))) # l = 0.5
se.3 = exp(-(D^2/(2*(0.0005^2)))) # l tends to zero (0.0005)
Y1 = rmvnorm(1, sigma = se.1)
Y2 = rmvnorm(1, sigma = se.2)
Y3 = rmvnorm(1, sigma = se.3)
Y = cbind(t(Y1), t(Y2), t(Y3))
matplot(X, Y, type = 'l', ylab = 'Y', lty = c(1,1,1)) # black = l = 2, red = l = 0.5   
matplot(X, t(Y2), type = 'l', xlim = c(0,5), ylim = c(-2,2))



#Matern Covariance functions - Figure 3.3(b)

l.1 = 2
MT.1 = exp(-D/l.1) 
MT.2 = (1+ sqrt(3)*D/(l.1))*exp(-sqrt(3)*D/(l.1))
MT.3 = (1+ sqrt(5)*D/(l.1) + (5*(D^2)/(3*l.1^2)))*exp(-sqrt(5)*D/(l.1))

Y1 = rmvnorm(1, sigma = MT.1)
Y2 = rmvnorm(1, sigma = MT.2)
Y3 = rmvnorm(1, sigma = MT.3)
Y = cbind(t(Y1), t(Y2), t(Y3))
matplot(X, Y, type = 'l', ylab = 'Y', lty = c(1,1,1)) # black = 1/2, red = 3/2

matplot(X, t(Y2), type = 'l', xlim = c(0,5), ylim = c(-2,2))



