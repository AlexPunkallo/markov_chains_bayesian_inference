#####  c)

# Dugongs data

# Length
x = c( 1.0, 1.5, 1.5, 1.5, 2.5, 4.0, 5.0, 5.0, 7.0, 8.0, 8.5,
       9.0, 9.5, 9.5, 10.0, 12.0, 12.0, 13.0, 13.0, 14.5, 15.5, 
       15.5, 16.5, 17.0, 22.5, 29.0, 31.5)

# Age
Y = c(1.80, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47, 2.19, 
      2.26, 2.40, 2.39, 2.41, 2.50, 2.32, 2.32, 2.43, 2.47, 2.56, 
      2.65, 2.47, 2.64, 2.56, 2.70, 2.72, 2.57)
N = 27

n = length(x)

#Hyperparameters
sigma.a = 100
sigma.b = 100
a = 0.001
b = 0.001

set.seed(1234)
library(pscl)
library(MCMCpack)

# Full alpha conditional
full.cond.alpha = function(beta, gamma, tau){
  # Alpha mu
  mu = (sigma.a^2 * sum(beta*gamma^x + Y)) / (tau + n*sigma.a^2 )
  # Alpha var
  var = (tau*sigma.a^2) / (tau + n*sigma.a^2)
  a.constr = 0
  while(a.constr < 1){
    a.constr = rnorm(1, mu, sqrt(var))
  }
  return(a.constr)
}

# Full beta conditional
full.cond.beta = function(alpha, gamma, tau){
  # Beta mu
  mu = (sigma.b^2*sum(alpha*gamma^x - Y*gamma^x))/(tau + sigma.b^2*sum(gamma^(2*x)))
  # Beta var
  var = (tau*sigma.b^2) / (tau + sigma.b^2*sum(gamma^(2*x)))
  b.constr = 0
  while(b.constr < 1){
    b.constr = rnorm(1, mu, sqrt(var))
  }
  return(b.constr)
}

# Full tau conditional
full.cond.tau = function(alpha, beta, gamma){
  # Tau Shape
  shape = n/2 + a
  # Tau rate
  rate = b + (sum((beta*gamma^x + Y - alpha)^2))/2
  tau = rigamma(1, shape, rate)
  return(tau)
}

# Full gamma conditional
full.cond.gamma = function(alpha, beta, tau, gamma){
  arg = -1/(2*tau) * sum((beta*gamma^x + Y - alpha)^2)
  return(exp(arg))
}


# Metropolis-within-Gibbs algorithm
metro_fun = function(alpha.old, beta.old, gamma.old, tau.old, n.sim){
  
  for (gibbs in 1:n.sim){
    alpha.old[gibbs + 1] = full.cond.alpha(beta.old[gibbs], gamma.old[gibbs], tau.old[gibbs])
    beta.old[gibbs + 1] = full.cond.beta(alpha.old[gibbs + 1], gamma.old[gibbs], tau.old[gibbs])
    gamma.cand = runif(1, 0, 1)
    gamma.prob.old = full.cond.gamma(alpha.old[gibbs + 1], beta.old[gibbs + 1], tau.old[gibbs],
                                     gamma.old[gibbs])
    gamma.prob.cand = full.cond.gamma(alpha.old[gibbs + 1], beta.old[gibbs + 1], tau.old[gibbs],
                                      gamma.cand)
    if (gamma.prob.cand/gamma.prob.old<1) {
      p = gamma.prob.cand/gamma.prob.old
    } else {
      p = 1
    }
    gamma.old[gibbs + 1] = sample(c(gamma.cand, gamma.old[gibbs]), size = 1, prob = c(p, 1-p))
    tau.old[gibbs + 1] = full.cond.tau(alpha.old[gibbs + 1], beta.old[gibbs + 1], gamma.old[gibbs + 1])
  }
  return(list(alpha.old,beta.old,gamma.old,tau.old))
}

n.sim = 13000
a = 0.1

alpha.old = rep(NA,n.sim+1)
beta.old = rep(NA,n.sim+1)
gamma.old = rep(NA,n.sim+1)
tau.old = rep(NA, n.sim+1)

alpha.old[1] = 2
beta.old[1] = 1
gamma.old[1] = 0.5
tau.old[1] = 1

mf =  metro_fun(alpha.old,beta.old,gamma.old,tau.old,n.sim)

# Burn-in procedure
alpha = mf[[1]][-(1:1000)]
beta = mf[[2]][-(1:1000)]
gamma = mf[[3]][-(1:1000)]
tau = mf[[4]][-(1:1000)]

alpha.hat = mean(alpha)
beta.hat = mean(beta)
gamma.hat = mean(gamma)
tau.hat = mean(tau)

######  d)

par(mfcol=c(2,2))
plot(seq(1, length(alpha)), alpha, type = "l", xlab = "t", ylab = expression(alpha)
     ,main='Alpha')
plot(seq(1, length(beta)), beta, type = "l", xlab = "t", ylab = expression(beta)
     ,main='Beta')
plot(seq(1, length(gamma)), gamma, type = "l", xlab = "t", ylab = expression(gamma)
     ,main='Gamma')
plot(seq(1, length(tau)), tau, type = "l", xlab = "t", ylab = expression(tau^2)
     ,main='Tau^2')

#####  e)

par(mfcol=c(2,2))
plot(cumsum(alpha)/1:length(alpha), type="l", ylab = "alpha")
plot(cumsum(beta)/1:length(beta), type="l", ylab = "beta")
plot(cumsum(gamma)/1:length(gamma), type="l", ylab = "gamma")
plot(cumsum(tau)/1:length(tau), type="l", ylab = "tau2")

#####  f)

library(batchmeans)
app.err.alpha = bm(alpha)$se
app.err.beta = bm(beta)$se
app.err.gamma = bm(gamma)$se
app.err.tau = bm(tau)$se


#####  g)

pu.alpha = sqrt(var(alpha))/mean(alpha)
pu.alpha = sqrt(var(beta))/mean(beta)
pu.alpha = sqrt(var(gamma))/mean(gamma)
pu.alpha = sqrt(var(tau))/mean(tau)

#####  h)

correlation.par = cor(data.frame(alpha, beta, gamma, tau))

par(mfcol=c(1,1))
df.corr = data.frame(alpha,beta,gamma,tau)
library(corrplot)
correl <- cor(df.corr)
corrplot.mixed(correl)

#####  i)

pos.pred.distr = function(x){
  mu = alpha - beta*gamma^x
  arg = rnorm(length(alpha), mu, tau)
  return(arg)
}

new_dug1 = pos.pred.distr(20)
new_dug_mean1 = mean(new_dug1)


######  j)

new_dug2 = pos.pred.distr(30)
new_dug_mean2 = mean(new_dug2)


######  k)

app.err1.dug1 = bm(new_dug1)$se
app.err2.dug2 = bm(new_dug2)$se
pos.unc1.dug1 = sqrt(var(new_dug1))/mean(new_dug1)
pos.unc2.dug2 = sqrt(var(new_dug2))/mean(new_dug2)
