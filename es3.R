#####  a)

set.seed(123)

# Transition Matrix
mpt<-matrix(c(0,0.5,0.5,
              5/8,1/8,1/4,
              2/3,1/3,0),
            nrow=3,byrow=T)
S=c(1,2,3)
# Initial State at initial time t=0
X0<-1

# Number of simulations
nsim <- 1000
# Markov Chain simulation function
mc_fun = function(mat, nsim, states, startingstate){
  chain<-rep(NA,nsim+1)
  # First state of the chain
  chain[1]<-startingstate
  # Simulation
  for(t in 1:nsim){
    chain[t+1]<-sample(states,size=1,prob=mat[chain[t],])
  }
  return(chain)
}
mc_sim1 = mc_fun(mpt,nsim,S,X0)
plot(mc_sim1, ylim=c(0,4), ylab = "Chain State", xlab = "Time",
     main = "Simulated Markov Chain")

#####  b)

# Relative frequencies
rel_freq_chain1 = table(mc_sim1)/(nsim+1)
cat('Relative Frequencies of the chain with X0 = 1: ', rel_freq_chain1)

plot(rel_freq_chain1,xlab="States",ylab="Relative Frequency", 
     main = "Simulated Markov Chain")

#####  c)

rep_mc_fun<- function(nrepetition, mat, nsim, states, 
                      startingstate){
  final.state = rep(NA, nrepetition)
  for (r in 1:nrepetition){
    chain_temp = rep(NA, nsim+1)
    chain_temp[1] = startingstate
    for(t in 1:nsim){
      chain_temp[t+1] = sample(S, size=1, prob=mat[chain_temp[t],])
    }
    final.state[r] = chain_temp[nsim+1]
  }
  return(final.state)
}

nrepetit = 500
rep_mc_sim1 = rep_mc_fun(nrepetit,mpt,nsim,S,X0)
rel_freq_sim1 = table(rep_mc_sim1)/(nrepetit+1)
cat('Relative frequencies of the simulation with X0 = 1: ', rel_freq_sim1)

plot(rel_freq_sim1,xlab = 'Chain State',ylab = 'Relative Frequencies',
     main='Empirical Relative Frequencies')

#####  d)

st_fun <- function(matrix1, matrix2){
  pi = solve(matrix1,matrix2)
  return(pi)
}

m_1 = matrix(c(-1,5/8,2/3,1/2,-0.875,1/3,1,1,1),nrow=3,byrow = T)
m_2 = matrix(data=c(0,0,1), nrow=3, ncol=1, byrow=FALSE)
pi = st_fun(m_1,m_2)
cat('Stationary distribution: ', pi)

#####  e)

Distr_plots <- read.table(
  header=TRUE, text='Distributions        State Probability
  1   Stationary_d       1      0.3917526
  2   Stationary_d       2      0.3298969
  3   Stationary_d       3      0.2783505
  4   Simulation_1       1      0.386
  5   Simulation_1       2      0.348
  6   Simulation_1       3      0.267
  7   Simulation_2       1      0.356
  8   Simulation_2       2      0.332
  9   Simulation_2       3      0.312')

library(ggplot2)
ggplot(Distr_plots, aes(State, Probability, fill = Distributions)) + 
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1")
# scale_fill_manual(values=c("#39c922", "#f4f4f4", "#ff0000"))

#####  f)

set.seed(123)

X0 <- 2

mc_sim2 = mc_fun(mpt,nsim,S,X0)
rel_freq_chain2 = table(mc_sim2)/(nsim+1)
cat('Relative frequencies of the chain with X0 = 2: ', rel_freq_chain2)

rep_mc_sim2 = rep_mc_fun(nrepetit,mpt,nsim,S,X0)
rel_freq_sim2 = table(rep_mc_sim2)/(nrepetit+1)
cat('Relative frequencies of the simulation with X0 = 2: ', rel_freq_sim2)



cat('Relative frequencies of the chain with X0 = 1: ', rel_freq_chain1)
cat('Relative frequencies of the chain with X0 = 2: ', rel_freq_chain2)
cat('Relative frequencies of the simulation with X0 = 1: ', rel_freq_sim1)
cat('Relative frequencies of the simulation with X0 = 2: ', rel_freq_sim2)
cat('Stationary distribution: ', pi)