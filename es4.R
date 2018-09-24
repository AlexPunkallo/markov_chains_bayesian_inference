##### c)

# Observations
y=c(4,5,4,1,0,4,3,4,0,6,3,3,4,0,2,6,3,3,5,4,5,3,1,4,
    4,1,5,5,3,4,2,5,2,2,3,4,2,1,3,2,1,1,1,1,1,3,0,0,
    1,0,1,1,0,0,3,1,0,3,2,2,0,1,1,1,0,1,0,1,0,0,0,2,
    1,0,0,0,1,1,0,2,2,3,1,1,2,1,1,1,1,2,4,2,0,0,0,1,
    4,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0)

n.obs = length(y)
n.iter=10000

stat.y.firstperiod=cumsum(y)[-length(y)]
stat.y.secondperiod=sum(y)-stat.y.firstperiod

n.iter=10000
lambda=rep(NA,n.iter+1)
phi=rep(NA,n.iter+1)
m=rep(NA,n.iter+1)
m.support=seq(1,n.obs-1)

# Starting Values
lambda.start=lambda[1]=5
phistart=phi[1]=9
m.start=m[1]=50

# Hyperparameters
ip.alpha=0.001
ip.beta=0.001
ip.a=0.001
ip.b=0.001

# Gibbs Algorithm
for(gibbs in 1:n.iter){
  lambda[gibbs+1]=rgamma(1,shape=ip.alpha+sum(y[1:m[gibbs]]),
                         rate=ip.beta+m[gibbs])
  phi[gibbs+1]=rgamma(1,shape=ip.a+sum(y[(m[gibbs]+1):n.obs]),
                      rate=ip.b+n.obs-m[gibbs])
  logci = - lambda[gibbs+1]*m.support + stat.y.firstperiod*log(lambda[gibbs+1]) +
    log(phi[gibbs+1])*stat.y.secondperiod + phi[gibbs+1]*m.support
  m.full.conditional.nn=exp(logci-max(logci))
  m[gibbs+1]=sample(x=m.support,size=1,prob=m.full.conditional.nn)
}

# Burn-in procedure
lambda=lambda[-(1:1001)]
phi=phi[-(1:1001)]
m=m[-(1:1001)]

#####  d)

mean(lambda)
mean(phi)
mean(m)


par(mfrow=c(3,2))
plot(lambda,type="l",main="Trace plot Lambda", xlab="gibbs", 
     ylab="lambda")
abline(h=mean(lambda),col='blue')
acf(lambda, main='ACF lambda')
plot(phi,type="l",main="Trace plot Phi", xlab="gibbs",
     ylab="phi ")
abline(h=mean(phi),col='blue')
acf(phi, main='ACF phi')
plot(m,type="l",main="Trace plot m", xlab="gibbs",
     ylab="m")
abline(h=mean(m),col='blue')
acf(m, main='ACF m')

par(mfrow=c(1,3))
hist(lambda,freq = F,col='yellow',main='Marginal posterior d. of lambda')
x1 = seq(min(lambda), max(lambda), length=length(lambda))
y1 = dgamma(x1, ip.alpha+sum(y[1:mean(m)]), ip.beta+mean(m))
lines(x1, y1, col="blue", lwd=2)
hist(phi,freq = F,col='orange',main='Marginal posterior d. of phi')
x2 = seq(min(phi), max(phi), length=length(phi))
y2 = dgamma(x2, sum(y[(mean(m)+1):n.obs])+ip.b, n.obs-mean(m)+ip.b)
lines(x2, y2, col="blue", lwd=2)
plot(table(m),type="h",main='Histogram of m',
     col='red',ylab="Density")
#library(coda)

######  e) 

plot(y,main = "Number of fatal accidents in UK coal mining sites", 
     xlab = "year", ylab = "number of accidents", type = "h")
abline(v = (round(mean(m))), col = "blue",lwd=2)

cat('There has been a change in the year: ', 1851 + (round(mean(m))))
# quindi ? il 1890 (1851+40)

######  f)

accid_reduct = (lambda-phi)/lambda*100
mean(accid_reduct)

hist(accid_reduct, freq = F,col='green',main='Histogram of the reduction of the rate
     of accidents')

cat('The expected reduction is: ', mean(accid_reduct),'%')

