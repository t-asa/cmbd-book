#----------------------------------------------------------#
# Model functions (for parameter fit)
#----------------------------------------------------------#

# Q-learning
func_qlearning <- function(param, data, prior = NULL)
{
  
  alpha <- param[1]
  beta <- param[2]
  c <- data$choice
  r <- data$reward
  T <- length(c)
  
  p1 <- numeric(T) 

  # set Q values (#option x T)
  Q <- matrix(numeric(2*T), nrow=2, ncol=T)
  
  # initialize log-likelihood
  ll <- 0
  
  for (t in 1:T) {
    
    # choosing prob 1
    p1[t] <- 1/(1+exp(-beta * (Q[1,t]-Q[2,t])))
    
    ll <- ll + (c[t]==1) * log(p1[t]) +  (c[t]==2) * log(1-p1[t])
    
    # update values 
    if (t < T) {
      
      Q[c[t],t+1] <- Q[c[t],t] + alpha * (r[t] - Q[c[t],t] ) 
      
      # for unchosen option
      Q[3-c[t],t+1] <- Q[3-c[t],t]
    }
  }
  
  # log prior density 
  if (is.null(prior)) {
    lprior <- 0
  } else {
    lprior <- dbeta(alpha,prior$alpha_a, prior$alpha_b,log = T) + 
      dgamma(beta,shape=prior$beta_shape, scale=prior$beta_scale,log = T) 
  }
  
  return(list(negll = -ll - lprior,Q = Q, p1 = p1))
}

# win-stay lose-shift
func_wsls <- function(param, data, prior = NULL)
{
  
  epsilon <- param[1] # error rate
  c <- data$choice
  r <- data$reward
  T <- length(c)
  
  p1 <- numeric(T) 
  
  # initialize log-likelihood
  ll <- 0
  
  for (t in 1:T) {
    
    # choosing prob 1
    if (t == 1) {
      p1[t] <- 0.5
    } else {
      if (r[t-1]==1)
        p1[t] <- (c[t-1]==1) * (1-epsilon) + (c[t-1]==2) * epsilon
      else 
        p1[t] <- (c[t-1]==1) * (epsilon) + (c[t-1]==2) * (1-epsilon)
    }
    
    ll <- ll + (c[t]==1) * log(p1[t]) +  (c[t]==2) * log(1-p1[t])
  }
  
  # log prior density 
  if (is.null(prior)) {
    lprior <- 0
  } else {
    lprior <- dbeta(epsilon, prior$epsilon_a, prior$epsilon_b, log = T)
  }
  
  return(list(negll = -ll - lprior, p1 = p1))
}

priorList <- list(
  list(alpha_a = 2, alpha_b = 2, beta_shape = 2, beta_scale = 3), # q-learning
  list(epsilon_a = 0.1, epsilon_b = 4)  # winstay-lose shift
)

# biased random choice
func_random <- function(param, data, prior = NULL)
{
  
  p1 <- param[1] # choice prob
  c <- data$choice
  T <- length(c)
  
  # initialize log-likelihood
  ll <- 0
  
  ll <- sum(c == 1) * log(p1) + sum(c == 2) * log(1-p1) 
  
  
  # log prior density 
  if (is.null(prior)) {
    lprior <- 0
  } else {
    lprior <- dbeta(p1, prior$p1_a, prior$p1_b, log = T)
  }
  
  return(list(negll = -ll - lprior, p1 = p1))
}

priorList <- list(
  list(alpha_a = 2, alpha_b = 2, beta_shape = 2, beta_scale = 3), # q-learning
  list(epsilon_a = 0.1, epsilon_b = 4),  # winstay-lose shift
  list(p1_a = 2, p1_b = 2)  # biased choice
)

modelfunctions <- c(func_qlearning, func_wsls, func_random)

nParamList <- c(2,1,1)

# # lower bound
# lblist <- list() 
# # upper bound
# ublist <- list()

nModel <- length(nParamList)

# # test
# source("parameter_fit_functions.R")
# paramfitML(modelfunctions, nParamList)

