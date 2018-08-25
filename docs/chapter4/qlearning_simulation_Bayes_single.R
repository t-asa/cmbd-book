# clear 
rm(list=ls())
graphics.off()

library(tidyverse)
library("gridExtra")
require(rstan)  

set.seed(141)

source("model_functions.R")
source("parameter_fit_functions.R")

#----------------------------------------------------------#
# Simulation for generating samples
#----------------------------------------------------------#

# number of simulation trials
T <- 80

# initialize Q values (#option x T)
Q <- matrix(numeric(2*T), nrow=2, ncol=T)

c <- numeric(T)
r <- numeric(T) 
p1 <- numeric(T) 

alpha <- 0.3     # learning rate
beta <- 2.0      # Inverse temperature

# reward probability for each option
pr <- c(0.7,0.3)

for (t in 1:T) {
  
  # choosing prob 1
  p1[t] <- 1/(1+exp(-beta*(Q[1,t]-Q[2,t])))
  
  if (runif(1,0,1) < p1[t]) {
    # choose option 1
    c[t] <- 1
    r[t] <- as.numeric(runif(1,0,1) < pr[1])
  } else {
    # choose option 2
    c[t] <- 2
    r[t] <- as.numeric(runif(1,0,1) < pr[2])
  }
  
  # update values 
  if (t < T) {
    
    Q[c[t],t+1] <- Q[c[t],t] + alpha * (r[t] - Q[c[t],t] ) 
    
    # for unchosen option
    Q[3-c[t],t+1] <- Q[3-c[t],t]
  }
}

data <- list(reward=r,choice=c)

#----------------------------------------------------------#
# Parameter fit ML 
#----------------------------------------------------------#
cat("----------- ML -----------\n")
resultML <- paramfitML(modelfunctions, data, nParamList)

#----------------------------------------------------------#
# Parameter fit MAP 
#----------------------------------------------------------#
cat("----------- MAP -----------\n")
resultMAP <- paramfitMAP(modelfunctions, data, nParamList, priorList)

#----------------------------------------------------------#
# Parameter fit Bayes (MCMC)
#----------------------------------------------------------#

# Put the data into a list.
dataList = list(    
  c = c,
  r = r,
  T = T
)

# run MCMC
stanFit <- stan(file='qlearning_single_subject.stan', 
                data=dataList, iter=20000, thin=10, chains=3)

# plot prior------------------ #
x11(width = 14*2, height = 12*2)

galpha_prior <- ggplot(data.frame(alpha = c(0, 1)), 
                       aes(x = alpha)) + 
  stat_function(fun = dbeta, args=list(2, 2),linetype = 1, size = 1.5) +
  ylab('density') +
  xlab('alpha') 

gbeta_prior <- ggplot(data.frame(beta = c(0, 20)), 
                      aes(x = beta)) + 
  stat_function(fun = dgamma, args=list(shape=2, scale=3),linetype = 1, size = 1.5) +
  ylab('density') +
  xlab('beta') 

# MCMCサンプルの抽出
df_post <- data.frame(rstan::extract(stanFit,"alpha"),
                       rstan::extract(stanFit,"beta"))

galpha_posterior <- ggplot(df_post,
            aes(x = alpha)) + 
  geom_histogram(aes(y = ..density..)) + 
  geom_density(size=1,linetype=1) +
  geom_vline(xintercept = resultMAP$paramlist[[1]][1], 
             size=1,linetype=2) + # MAP estimates
  geom_vline(xintercept = resultML$paramlist[[1]][1], 
             size=0.5,linetype=3) + # ML estimates
  ylab('density') +
  xlab('alpha') 

gbeta_posterior <- ggplot(df_post,
            aes(x = beta)) + 
  geom_histogram(aes(y = ..density..)) + 
  geom_density(size=1,linetype=1) + 
  geom_vline(xintercept = resultMAP$paramlist[[1]][2], 
             size=1,linetype=2) + # MAP estimates
  geom_vline(xintercept = resultML$paramlist[[1]][2], 
             size=0.5,linetype=3) + # ML estimates
  ylab('density') +
  xlab('beta') 

grid.arrange(galpha_prior, gbeta_prior, 
             galpha_posterior, gbeta_posterior, ncol = 2)

# g <- arrangeGrob(galpha_prior, gbeta_prior, 
#                 galpha_posterior, gbeta_posterior, ncol = 2)
# ggsave(file="./figs/posterior_single_qlearning.eps", g) 

