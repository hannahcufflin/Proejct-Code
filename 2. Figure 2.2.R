coal.ash = read.table('C:/Users/hanna/OneDrive/Desktop/PROJECT/DATA/coal.ash.txt', header = TRUE)
plot(coal.ash$x, coal.ash$y)

#Figure 5.2(a)
library(geoR)
library(plgp)
par(mar=c(4.5,4.5,1.5,1.5))
coords1 = as.matrix(coal.ash[,c('x','y')])
max.dist =  0.5 * max(distance(coords1)
)
bins = 1000
vario.coal = variog(coords = coords1, data = coal.ash$coal, uvec = (seq(0, max.dist, length = bins)) )
fit.coal = variofit(vario.coal, ini.cov.pars = c(2,1), cov.model = 'exponential', minimisation.function = 'optim', weights = 'equal')
plot(vario.coal, xlim = c(0,25))
lines(fit.coal, col = 'red')
abline(h= fit.coal$nugget, col = 'blue')
abline(h= fit.coal$cov.pars[1] +fit.coal$nugget, col = 'green')
abline(v = fit.coal$cov.pars[2], col = 'orange')

#Values of tau = 0.9325994

#Figure 5.2(b)
vario.coal = variog(coords = coords1, data = coal.ash$coal, uvec = (seq(0, max.dist, length = bins)) )
fit.coal = variofit(vario.coal, ini.cov.pars = c(2,1), cov.model = 'gaussian', minimisation.function = 'optim', weights = 'equal')
plot(vario.coal, xlim = c(0,25))
lines(fit.coal, col = 'red')
abline(h= fit.coal$nugget, col = 'blue')
abline(h= fit.coal$cov.pars[1] +fit.coal$nugget, col = 'green')
abline(v = fit.coal$cov.pars[2], col = 'orange')


#value of tau = 0.5563717