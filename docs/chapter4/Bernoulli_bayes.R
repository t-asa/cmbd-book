# clear 
rm(list=ls())
graphics.off()

library(tidyverse)
library(gridExtra)

ggplot() + theme_set(theme_bw(base_size = 18,base_family="Arial")) 

# parameter for prior 
a <- 5
b <- 5

data <- c(1,1,1)

likelihood <- function(mu, data) {
  C1 <- sum(data) # count for data=1
  C0 <- sum(1-data) # count for data=1
  return(exp(C1 * log(mu) + C0*log(1-mu)))
}

prior <- function(mu, a, b) {
  return( dbeta(mu, a, b, log = F) )
}

# log posterior (unnormalized)
post <- function(mu, data, a, b) {
  
  C1 <- sum(data) # count for data=1
  C0 <- sum(1-data) # count for data=1
  
  return( dbeta(mu, a+C1, b+C0, log = F) )
}

x11()
# log likelihood
gll <- ggplot(data.frame(mu = c(0, 1)), aes(x=mu)) +
  stat_function(fun = likelihood, 
                args = list(data=data), linetype = 1, size = 1.5) +
  ylab('pdf') +
  xlab('mu') +
  ggtitle("likelihood") +
  scale_y_continuous(breaks=c(0, 1, 2), labels = c(0, 1, 2)) + 
  coord_cartesian(ylim = c(0,1))

# log prior
glprior <- ggplot(data.frame(mu = c(0, 1)), aes(x=mu)) +
  stat_function(fun = prior, 
                args = list(a=a, b=b), linetype = 1, size = 1.5) +
  ylab('pdf') +
  xlab('mu') +
  ggtitle("prior")  + 
  scale_y_continuous(breaks=c(0, 1, 2, 3), labels = c(0, 1, 2, 3)) + 
  coord_cartesian(ylim = c(0,3))

# log posterior
glp <- ggplot(data.frame(mu = c(0, 1)), aes(x=mu)) +
  stat_function(fun = post, 
                args = list(data=data, a=a, b=b), linetype = 1, size = 1.5) +
  ylab('pdf') +
  xlab('mu') +
  ggtitle("poterior") + 
  scale_y_continuous(breaks=c(0, 1, 2, 3), labels = c(0, 1, 2, 3)) +
  coord_cartesian(ylim = c(0,3))

grid.arrange(gll, glprior, glp, nrow=3) 

# 図を保存する場合は以下を実行
# g <- arrangeGrob(gll, glp,glp, nrow=2) #generates g
# ggsave(file="./figs/bernoulli_bayes.eps", g) #saves g

