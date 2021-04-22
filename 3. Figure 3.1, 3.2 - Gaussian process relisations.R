## Chapter 3, Figure 3.1 and 3.2 - sampling from the Gaussian processes
par(mar=c(4.5,4.5,2,2))
n=100
X = matrix(seq(-2, 2, length=100), ncol=1)
Sigma = matrix(0,100,100)
for (i in 1:100){
  for (j in 1:100){
    Sigma[i,j] = (abs(X[i]-X[j]))^2
  }
}
Sigma1 = exp(-Sigma/2) # squared exponential covairance function l = 1
## Figure 3.1 - 5 different relaisation from the prior 
Z = rmvnorm(5, mean = rep(0,100), sigma=Sigma1) # squared exponential function
matplot (X,t(Z), type = 'l',lty = c(1,1,1),lwd = 1, ylab = 'f(x)', xlab = 'x')
# Repeat for two differnt plots 

par(mar=c(4.5,4.5,2.5,2.5))


# set the observed data
n = 5 # Repeat for n=8 to get Figure 3.2(b)
X = matrix(seq(0,2*pi,length=n), ncol=1)
Y = cos(X) # Observed data is from the cos function 
D = distance(X) #n by n matrix
Sigma = matrix(0,n,n)
for (i in 1:n){
  for (j in 1:n){
    Sigma[i,j] = (abs(X[i]-X[j]))^2
  }
}
Sigma = exp(-0.5*Sigma/2)

eps <- sqrt(.Machine$double.eps) # finding covariance matrix for the test data
XX <- matrix(seq(-0.5, 2*pi+0.5, length=100), ncol=1)
DXX <- distance(XX)
SXX = matrix(0,100,100)
for (i in 1:100){
  for (j in 1:100){
    SXX[i,j] = (abs(XX[i]-XX[j]))^2
  }
}
SXX = exp(-0.5*SXX/2) #for noisy observations need to add the variance

#covariance matrix for the text and training data
DX <- distance(XX, X) #100 by 5 marix
SX = matrix(0,100,n)
for (i in 1:100){
  for (j in 1:n){
    SX[i,j] = (abs(XX[i]-X[j]))^2
  }
}
SX = exp(-0.5*SX/2)

Si <- solve(Sigma)
mup <- SX %*% Si %*% Y #100 by 1 this is the mean of the posterior 
Sigmap <- SXX - SX %*% Si %*% t(SX) # 100X100
YY <- rmvnorm(100, mup, Sigmap)

q1 <- mup + qnorm(0.05, 0, sqrt(diag(Sigmap)))
q2 <- mup + qnorm(0.95, 0, sqrt(diag(Sigmap)))
## Figure 3.2(a),(b)
matplot(XX, t(YY), type="l", col="gray", lty=1, xlab="x", ylab="f(x)")
points(X, Y, pch=20, cex=2) # the observed data points 
lines(XX, mup, lwd=2) # the black line - the mean line
lines(XX, cos(XX), col="blue") # the blue line (what the observed data should be)
lines(XX, q1, lwd=2, lty=2, col=2) # lower dotted red line
lines(XX, q2, lwd=2, lty=2, col=2)# upper dotted red line (95% quantile)















