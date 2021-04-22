### Modellling the semi-variograms ###

### Exponential ###
tau = 0.2
sigma = 0.5 
phi = 2

x = seq(from = 0,to = 4, by = 0.00001)
y.exp = tau + sigma*(1 - exp(-phi*x))
plot(x, y.exp, type = 'l', xlab = 'h', ylab = 'Semivariance', xlim = c(0,3))
abline(h = 0.7, col = 'green')
abline(h = 0.2, col = 'red')


### Spherical ###
x.phi = seq(from = 0, to = 1/phi, by = 0.00001)
x.phi.2 = seq(from = 1/phi, to = 4, by = 0.00001)
y.sph.1 = tau + sigma
y.sph.2 = tau + sigma*((3*phi*x.phi)/2 - 0.5*(phi*x.phi)**3)
y.sph = rep(y.sph.1, 350001)
x = c(x.phi, x.phi.2)
y = c(y.sph.2, y.sph)
plot(x,y, type = 'l', xlab = 'h', ylab = 'Semivariance')
abline(h = 0.2, col = 'red')
#abline(h = 0.71, col = 'green')
points(x = c(0,1/phi), y = c(0.7, 0.7), type = 'l', col = 'green')
abline(v = 1/phi, col = 'blue')

