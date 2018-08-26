# メモリ，グラフのクリア
rm(list=ls())
graphics.off()

# 最尤推定，MAP推定，最適化のためRsolnpを読み込み
library(Rsolnp)

# RStanの読み込み
require(rstan)

# 以下を実行するとStanのサンプリングを並列化で高速にできる
# rstan_options(auto_write=TRUE)
# options(mc.cores=parallel::detectCores())

nseed <- 182
set.seed(nseed)

flgCompileModel <- 1
flgHB <- 1

# Simulation parameters:
N <- 40             # Number of subjct
T <-  80            # Number of trials for each subject

nParam <- 2
paramName <- c("alpha", "beta")
collist <- rep("black",5) 
pchlist <- c(1, 4, 0)

nModel <- length(smodels)

# Model function (for parameter fit)
#----------------------------------------------------------#
# Model function (for individual parameter fit)
############## ---------------------------------------------#
func_qlearning <- function(param, choice, reward, prior = NULL)
{
  T <- length(choice)
  alpha <- param[1]
  beta <- param[2]
  c <- choice
  r <- reward
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

#----------------------------------------------------------#
# Model function for fixed effect, pooled group parameter fit
# ---------------------------------------------#
func_qlearning_FE <- function(param, choice, reward)
{
  
  N <- dim(choice)[1]
  T <- dim(choice)[2]
  alpha <- param[1]
  beta <- param[2]
  c <- choice
  r <- reward
  
  # initialize log-likelihood
  ll <- 0
  
  for (idxsub in 1:N) {
    # set Q values (#option x T)
    Q <- matrix(numeric(2*T), nrow=2, ncol=T)
    p1 <- numeric(T) 
    
    for (t in 1:T) {
      
      # choosing prob 1
      p1[t] <- 1/(1+exp(-beta * (Q[1,t]-Q[2,t])))
      
      ll <- ll + (c[idxsub, t]==1) * log(p1[t]) +  (c[idxsub, t]==2) * log(1-p1[t])
      
      # update values 
      if (t < T) {
        
        Q[c[idxsub, t],t+1] <- Q[c[idxsub, t],t] + alpha * (r[idxsub, t] - Q[c[idxsub, t],t] ) 
        
        # for unchosen option
        Q[3-c[idxsub, t],t+1] <- Q[3-c[idxsub, t],t]
      }
    }
  }
  return(-ll)
}

#----------------------------------------------------------#
# objective function to be minimized (for parameter fit) 
# just pass the negative-likelihood
#----------------------------------------------------------#
func_minimize <- function(param, choice, reward, prior = NULL)
{
  ret = func_qlearning(param, choice, reward, prior)
  
  # return negative log-likelihood
  return(ret$negll)
}

#----------------------------------------------------------#
# Simulation for generating samples
#----------------------------------------------------------#

c <- matrix(numeric(N * T), nrow = N, ncol = T)
r <- matrix(numeric(N * T), nrow = N, ncol = T)
p1 <- numeric(T)

alpha <- rnorm(N, 0.4, 0.2)     # learning rate
beta <- rnorm(N, 3.0, 1.0)     # inverse temperature
alpha[alpha < 0.01] <- 0.01
alpha[alpha > 0.99] <- 0.99
beta[beta < 0.0] <- 0.01

true_param <- list(alpha = alpha, beta = beta)

# reward probability for each option
p <- c(0.7, 0.3)

for (i in 1:N) {
  # initialize Q values (#option x T)
  Q <- matrix(numeric(2 * T), nrow = 2, ncol = T)
  for (t in 1:T) {
    # choosing prob 1
    p1[t] <- 1 / (1 + exp(-beta[i] * (Q[1, t] - Q[2, t])))
    
    if (runif(1, 0, 1) < p1[t]) {
      # choose option 1
      c[i, t] <- 1
      r[i, t] <- as.numeric(runif(1, 0, 1) < p[1])
    } else {
      # choose option 2
      c[i, t] <- 2
      r[i, t] <- as.numeric(runif(1, 0, 1) < p[2])
    }
    
    # update values
    if (t < T) {
      Q[c[i, t], t + 1] <- Q[c[i, t], t] + alpha[i] * (r[i, t] - Q[c[i, t], t])
      
      # for unchosen option
      Q[3 - c[i, t], t + 1] <- Q[3 - c[i, t], t]
    }
  }
}

# Load the data:
dataList = list(# Put the information into a list.
  c = c,
  r = r,
  N = N,
  T = T)

#------------------------------------------------------------------------#
# single subject MLE
#------------------------------------------------------------------------#
paramest <- matrix(nrow = N, ncol = nParam)
paraminvhessian <- matrix(nrow = N, ncol = nParam)
for (idxsub in 1:N) {
  fval <- Inf
  for (idxopt in 1:10) {
    initparam <- runif(nParam, 0, 1.0)
    
    # solnpで最適化
    res <- solnp(
      initparam,
      fun = func_minimize,
      LB = c(0, 0),
      UB = c(1, 20),
      control = list(trace = 0),
      choice = c[idxsub, ],
      reward = r[idxsub, ]
    )
    
    if (fval > res$values[length(res$values)]) {
      fval <- res$values[length(res$values)]
      paramest[idxsub, ] <- res$pars
      
      if (det(res$hessian) > 0.01) {
        paraminvhessian[idxsub, ] <- diag(solve(res$hessian))
      } else {
        paraminvhessian[idxsub, ] <-
          rep(1, nParam) # diag(solve(res$hessian))
      }
      
      cat("\n solnp: ", res$pars, " fval: ", fval)
    }
  }
}
SS_param <- list(paramest[, 1], paramest[, 2])

#------------------------------------------------------------------------#
# single subject MAP
#------------------------------------------------------------------------#
prior <-
  list(
    alpha_a = 2,
    alpha_b = 2,
    beta_shape = 2,
    beta_scale = 3
  )
paramest <- matrix(nrow = N, ncol = nParam)
paraminvhessian <- matrix(nrow = N, ncol = nParam)
for (idxsub in 1:N) {
  fval <- Inf
  for (idxopt in 1:10) {
    initparam <- runif(nParam, 0, 1.0)
    
    res <- solnp(
      initparam,
      fun = func_minimize,
      LB = c(0, 0),
      UB = c(1, 20),
      control = list(trace = 0),
      choice = c[idxsub, ],
      reward = r[idxsub, ],
      prior = prior
    )
    
    if (fval > res$values[length(res$values)]) {
      fval <- res$values[length(res$values)]
      paramest[idxsub, ] <- res$pars
      
      if (0) {
      if (det(res$hessian) > 0.01) {
        paraminvhessian[idxsub, ] <- diag(solve(res$hessian))
      } else {
        paraminvhessian[idxsub, ] <-
          rep(1, nParam) # diag(solve(res$hessian))
      }
      }
      
      cat("\n solnp: ", res$pars, " fval: ", fval)
    }
  }
}
SSMAP_param <- list(paramest[, 1], paramest[, 2])

#------------------------------------------------------------------------#
# pooled, fixed effect MLE
#------------------------------------------------------------------------#
cat("\n fixed effect MLE----------")
fval <- Inf
for (idxopt in 1:10) {
  # fixed effects
  initparam <- runif(nParam, 0, 1.0)
  
  res <- solnp(
    initparam,
    fun = func_qlearning_FE,
    LB = c(0, 0),
    UB = c(1, 20),
    control = list(trace = 0),
    choice = c,
    reward = r
  )
  
  if (fval > res$values[length(res$values)]) {
    fval <- res$values[length(res$values)]
    paramest <- res$pars
    cat("\n solnp: ", res$pars , " fval: ", res$values[length(res$values)])
    
  }
  
  FE_param <- list(paramest[1], paramest[2])
}

#------------------------------------------------------------------------#
# 階層ベイズモデルの推定
#------------------------------------------------------------------------#

# initize chains from MLE estimates -------------------- #
initf <- function(chain_id = 1, SS_param) {
  logit_alpha <-
    log(pmin(SS_param[[1]], 0.9) / (1 - pmin(SS_param[[1]], 0.9)))
  logit_beta <-
    log(pmin(SS_param[[2]] / 20, 0.9) / (1 - pmin(SS_param[[2]] / 20, 0.9)))
  # on doing
  mean(logit_alpha)
  mu_p_alpha = mean(logit_alpha)
  sigma_p_alpha = min(sd(logit_alpha), 1.4)
  mu_p_beta = mean(logit_beta)
  sigma_p_beta = min(sd(logit_beta), 1.4)
  list(
    alpha = SS_param[[1]],
    beta = SS_param[[2]],
    mu_p_alpha = mu_p_alpha,
    sigma_p_alpha = sigma_p_alpha,
    mu_p_beta = mu_p_beta,
    sigma_p_beta = sigma_p_beta,
    eta_alpha = (logit_alpha - mu_p_alpha) / sigma_p_alpha,
    eta_beta = (logit_beta - mu_p_beta) / sigma_p_beta
  )
}
# generate a list of lists to specify initial values
initsList <- lapply(1:nChains,
                    function(id)
                      initf(chain_id = id, SS_param = SS_param))

# Get MC sample of posterior:
cat("\n HB ----------")

# run MCMC
stanFit <- stan(file='model_qlearning_group.stan', 
                data=dataList, iter=10000, thin=1, 
                chains=3,
                init = initsList)

#------------------------------------------------------------------------#
# Evaluating results
#------------------------------------------------------------------------#

# 個人レベルパラメータのMCMCサンプルを抽出
param <- c("alpha","beta") # 分析対象とする個人レベルパラメータ
ms <- rstan::extract(stanFit, param)
HB_MAP_param <- list()

# 参加者ごとのMAP推定値を求める
for (p in 1:length(param)){
  HB_MAP_param[[p]] <- array(0,N)
  
  for (i in 1:N) {
    dres <- density( ms[[p]][,i] )
    HB_MAP_param[[p]][i] <- dres$x[which.max(dres$y)]
  }
}

# ------------------------------------ #
# correlation plot SS ML vs HB
# ------------------------------------ #
listxylim <- list(c(0,1), c(0,7))
x11(width = 7*nParam,height = 8)
par(mfcol=c(1,nParam))
par(oma=c(1.5,2,2,1))
par(mar = c(4, 4, 4.0, 2))

for (idxparam in 1:nParam) {
  
  COR_SS[idxparam] <- cor(x=true_param[[idxparam]], 
                          y = SS_param[[idxparam]])
  COR_HB[idxparam] <- cor(x=true_param[[idxparam]], 
                          y = HB_MAP_param[[idxparam]])
  
  # x11()
  maxse <- max(sqrt(paraminvhessian[,idxparam]))
  minse <- min(sqrt(paraminvhessian[,idxparam]))
  meanse <- mean(sqrt(paraminvhessian[,idxparam]))
  
  plot(true_param[[idxparam]], SS_param[[idxparam]], 
       pch = 16,
       col = "gray35",
       cex = ( (sqrt(paraminvhessian[,idxparam]) - minse) / (maxse - minse) * 
                 3 + 0.5),
       main = paste(paramName[idxparam]), 
       las = 1, 
       xlab = "True parameter",
       ylab = "Estimates",
       xlim = listxylim[[idxparam]], 
       ylim = listxylim[[idxparam]], 
       cex.lab  = 1.5,
       cex.axis = 1.5, 
       cex.main = 1.5)
  
  # points(x = true_theta, y=theta_EB, pch = pchlist[4], col = collist[4])
  
  points(x = true_param[[idxparam]], y= SSMAP_param[[idxparam]], 
         pch = pchlist[2])
  
  # plot fixed effect parameter
  abline(h = FE_param[[idxparam]], lty=5)
  
  points(x = true_param[[idxparam]], y= HB_MAP_param[[idxparam]], pch = pchlist[3])
  
  for (idx in 1:N) {
    arrows(true_param[[idxparam]][idx], SS_param[[idxparam]][idx], 
           true_param[[idxparam]][idx], HB_MAP_param[[idxparam]][idx], 
           col="black", 
           length = 0.1)
  }
  
  abline(a = 0, b = 1, lty = "dashed")
  # legend("topleft", 
  #        legend=c(
  #          sprintf('Single-subject MLE, r = %.3f', 
  #                  COR_SS[idxparam]), 
  #          sprintf('Hierarchical Bayes, r = %.3f', 
  #                  COR_HB[idxparam])
  #        ),
  #        bty = "n",
  #        inset = c(0.02, 0.02),pch = c(pchlist[1],  pchlist[3]))
  # 
}

dev.copy2eps(file="./figs/SS_ML_vs_HB_MAP.eps")
