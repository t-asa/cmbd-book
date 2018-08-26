# メモリ，グラフのクリア
rm(list=ls())
graphics.off()

# 最尤推定，MAP推定，最適化のためRsolnpを読み込み
library(Rsolnp)

# RStanの読み込み
require(rstan)

# hessianをnumDerivで求める場合
require(numDeriv)

# 以下を実行するとStanのサンプリングを並列化して実行できる
# rstan_options(auto_write=TRUE)
# options(mc.cores=parallel::detectCores())

nseed <- 182
set.seed(nseed)

# 参加者数
N <- 40 

# 各参加者の試行数
T <-  80   

nParam <- 2
paramName <- c("alpha", "beta")
collist <- rep("black",5) 
pchlist <- c(1, 4, 0)

#----------------------------------------------------------#
# 最尤推定のためのモデル関数の設定
#----------------------------------------------------------#

# 個人レベルのQ学習モデル
func_qlearning <- function(param, choice, reward, prior = NULL)
{
  T <- length(choice)
  alpha <- param[1]
  beta <- param[2]
  c <- choice
  r <- reward
  pA <- numeric(T) 
  
  # set Q values (#option x T)
  Q <- matrix(numeric(2*T), nrow=2, ncol=T)
  
  # initialize log-likelihood
  ll <- 0
  
  for (t in 1:T) {
    
    # choosing prob 1
    pA[t] <- 1/(1+exp(-beta * (Q[1,t]-Q[2,t])))
    
    ll <- ll + (c[t]==1) * log(pA[t]) +  (c[t]==2) * log(1-pA[t])
    
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
  
  return(list(negll = -ll - lprior,Q = Q, pA = pA))
}

# 固定効果のQ学習モデル
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
    pA <- numeric(T) 
    
    for (t in 1:T) {
      
      # choosing prob 1
      pA[t] <- 1/(1+exp(-beta * (Q[1,t]-Q[2,t])))
      
      ll <- ll + (c[idxsub, t]==1) * log(pA[t]) +  (c[idxsub, t]==2) * log(1-pA[t])
      
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


func_minimize <- function(param, choice, reward, prior = NULL)
{
  ret <- func_qlearning(param, choice, reward, prior)
  
  # return negative log-likelihood
  return(ret$negll)
}

#----------------------------------------------------------#
# Q学習モデルのシミュレーションによるデータ生成
#----------------------------------------------------------#
c <- matrix(numeric(N * T), nrow = N, ncol = T)
r <- matrix(numeric(N * T), nrow = N, ncol = T)
pA <- numeric(T)

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
    pA[t] <- 1 / (1 + exp(-beta[i] * (Q[1, t] - Q[2, t])))
    
    if (runif(1, 0, 1) < pA[t]) {
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
# 個人レベルの最尤推定
#------------------------------------------------------------------------#
paramest <- matrix(nrow = N, ncol = nParam)
param_se  <- matrix(nrow = N, ncol = nParam)
for (idxsub in 1:N) {
  # 参加者ごとに...
  
  fval <- Inf
  for (idxopt in 1:10) {
    
    # 初期値を一様乱数から生成
    initparam <- runif(nParam, 0, 1.0)
    
    # solnpで負の対数尤度を最小にするパラメータを求める
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
      
      # numDerivパッケージのhessianを使った方がヘッセ行列は精度よく求まる
      # require(numDeriv)
      hess <- hessian(func=func_minimize,
                      x = res$pars,
                      choice = c[idxsub, ],
                      reward = r[idxsub, ])
      
      # hess <- res$hessian
      param_se[idxsub, ] <- sqrt(diag(solve(hess)))
      
      cat("\n solnp: ", res$pars, " fval: ", fval)
    }
  }
}
SS_param <- list(paramest[, 1], paramest[, 2])

#------------------------------------------------------------------------#
# 個人レベルのMAP推定
#------------------------------------------------------------------------#

# 事前分布のパラメータの設定
prior <-
  list(
    alpha_a = 2,
    alpha_b = 2,
    beta_shape = 2,
    beta_scale = 3
  )

paramest <- matrix(nrow = N, ncol = nParam)

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

nChains <- 3

# MCMCの初期値を個人レベルの最尤推定値に近くなるよう設定する
initf <- function(chain_id = 1, SS_param) {
  logit_alpha <-
    log(pmin(SS_param[[1]], 0.9) / (1 - pmin(SS_param[[1]], 0.9)))
  logit_beta <-
    log(pmin(SS_param[[2]] / 20, 0.9) / (1 - pmin(SS_param[[2]] / 20, 0.9)))
  mu_p_alpha <- mean(logit_alpha)
  sigma_p_alpha = min(sd(logit_alpha), 1.4)
  mu_p_beta <- mean(logit_beta)
  sigma_p_beta <- min(sd(logit_beta), 1.4)
  
  # リストにまとめて返す
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
# チェインごとの初期値を生成 (この例では全て同じ初期値)
initsList <- lapply(1:nChains,
                    function(id)
                      initf(chain_id = id, SS_param = SS_param))

# run MCMC
stanFit <- stan(file = 'model_qlearning_group.stan', 
                data = dataList, iter = 10000, 
                thin = 1, 
                chains = nChains, 
                warmup = 2000,
                init = initsList)

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

#------------------------------------------------------------------------#
# パラメータの真値と各推定値の相関図
#------------------------------------------------------------------------#

# 各パラメータのx軸，y軸共通の範囲を設定
listxylim <- list(c(0,1.05), c(0,7))

# 描画ウィンドウの生成
x11(width = 7*nParam, height = 8)

# 描画パラメータ
par(mfcol = c(1,nParam))
par(oma = c(1.7,2,2,1))
par(mar = c(4, 4, 4.0, 2))

for (idxparam in 1:nParam) {
  
  maxse <- max(param_se[,idxparam])
  minse <- min(param_se[,idxparam])
  meanse <- mean(param_se[,idxparam])
  sdse <- sd(param_se[,idxparam])
  
  plot(true_param[[idxparam]], SS_param[[idxparam]], 
       pch = 16,
       col = "gray55",
       cex = ( (param_se[,idxparam]) / sdse * 0.7 + 0.8),
       main = paste(paramName[idxparam]), 
       las = 1, 
       xlab = "パラメータの真値",
       ylab = "推定値",
       xlim = listxylim[[idxparam]], 
       ylim = listxylim[[idxparam]], 
       cex.lab  = 1.5,
       cex.axis = 1.5, 
       cex.main = 1.5)
  
  # points(x = true_theta, y=theta_EB, pch = pchlist[4], col = collist[4])
  
  points(x = true_param[[idxparam]], y= SSMAP_param[[idxparam]], 
         pch = pchlist[2])
  
  # 固定効果分析の推定値を水平線で表示
  abline(h = FE_param[[idxparam]], lty=5)
  
  points(x = true_param[[idxparam]], y= HB_MAP_param[[idxparam]], pch = pchlist[3])
  
  for (idx in 1:N) {
    arrows(true_param[[idxparam]][idx], SS_param[[idxparam]][idx], 
           true_param[[idxparam]][idx], HB_MAP_param[[idxparam]][idx], 
           col="black", 
           length = 0.1)
  }
  
  abline(a = 0, b = 1, lty = "dashed")
}

# 図を保存するときは以下を実行
# dev.copy2eps(file="./figs/SS_ML_vs_HB_MAP.eps")
