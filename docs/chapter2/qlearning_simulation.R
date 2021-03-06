# メモリ，グラフのクリア
rm(list=ls())
graphics.off()

# 描画のためのライブラリ読み込み
library(tidyverse)
library(gridExtra)

# 乱数のシードの設定 
set.seed(141)

#----------------------------------------------------------#
# Q学習モデルのシミュレーションによるデータ生成
#----------------------------------------------------------#

# 試行数
T <- 80

# Q 値の初期化( 選択肢の数x T)
Q <- matrix(numeric(2*T), nrow=2, ncol=T)

c <- numeric(T)
r <- numeric(T) 
p1 <- numeric(T) 

alpha <- 0.3     # 学習率
beta <- 2.0      # 逆温度

# それぞれの選択肢の報酬確率
pr <- c(0.7,0.3)

for (t in 1:T) {
  
  # ソフトマックスで選択肢A の選択確率を決定する
  p1[t] <- 1/(1+exp(-beta*(Q[1,t]-Q[2,t])))
  
  if (runif(1,0,1) < p1[t]) {
    # Aを選択
    c[t] <- 1
    r[t] <- as.numeric(runif(1,0,1) < pr[1])
  } else {
    # Bを選択
    c[t] <- 2
    r[t] <- as.numeric(runif(1,0,1) < pr[2])
  }
  
  # 行動価値の更新
  if (t < T) {
    
    Q[c[t],t+1] <- Q[c[t],t] + alpha * (r[t] - Q[c[t],t] ) 
    
    # 選択肢していない行動の価値はそのままの値を次の時刻に引き継ぐ。
    # 3-c でc=1 なら2, c=2 なら1, というように
    # 逆側の選択肢のインデックスが求まる。
    Q[3-c[t],t+1] <- Q[3-c[t],t]
  }
}

#----------------------------------------------------------#
# フィットするモデルの設定
#----------------------------------------------------------#

# Q-learning
func_qlearning <- function(param, choice, reward)
{
  T <- length(choice)
  alpha <- param[1]
  beta <- param[2]
  
  # 短い名前の変数に持ち替える
  c <- choice
  r <- reward
  
  p1 <- numeric(T) 

  # Q 値の初期化( 選択肢の数x T)
  Q <- matrix(numeric(2*T), nrow=2, ncol=T)
  
  # 対数尤度を格納する変数
  ll <- 0
  
  for (t in 1:T) {
    
    # ソフトマックスで選択肢A の選択確率を決定する
    p1[t] <- 1/(1+exp(-beta * (Q[1,t]-Q[2,t])))
    
    # 試行tの対数尤度は実際の選択がA (c=1) であれば log(p1[t]), 
    # B (c=2) であればlog(1 - p1[t]) となる
    ll <- ll + (c[t]==1) * log(p1[t]) +  (c[t]==2) * log(1-p1[t])
    
    # 行動価値の更新
    if (t < T) {
      
      Q[c[t],t+1] <- Q[c[t],t] + alpha * (r[t] - Q[c[t],t] ) 
      
      # 選択肢していない行動の価値はそのままの値を次の時刻に引き継ぐ
      Q[3-c[t],t+1] <- Q[3-c[t],t]
    }
  }
  return(list(negll = -ll,Q = Q, p1 = p1))
}

# 最適化により最小化する負の対数尤度を返す関数
func_minimize <- function(modelfunc, param, choice, reward)
{
  ret <- modelfunc(param, choice, reward)
  
  # 負の対数尤度のみ返す
  return(ret$negll)
}

#----------------------------------------------------------#
# 最適化による最尤推定の実行
#----------------------------------------------------------#

# 負の対数尤度の最小値を格納する変数 (最初は無限大にしておく)  
fvalmin <- Inf

for (idx in 1:10) {
  
  # 初期値を一様乱数から決める
  initparam <- runif(2, 0, 1.0)
  
  # 最適化の実行
  res <- optim(initparam, func_minimize,
              hessian = TRUE, 
              modelfunc = func_qlearning, 
              choice=c, reward=r)
  
  # 今までの解より負の対数尤度が小さかったらその結果を採用する
  if (res$value < fvalmin) {
    paramest <- res$par
    fvalmin <- res$value
  }
}

print(sprintf("alpha - True value: %.2f, Estimated value: %.2f", alpha, paramest[1]))
print(sprintf("beta  - True value: %.2f, Estimated value: %.2f", beta, paramest[2]))
print(sprintf("Model 1: log-likelihood: %.2f, AIC: %.2f", -fvalmin, 2*fvalmin + 2*2))

# 求めた最尤推定値をもとに，行動価値や選択確率のP(a=A)を改めて計算する
ret <- func_qlearning(paramest, choice=c, reward=r)
Qest <- ret$Q
p1est <- ret$p1


#----------------------------------------------------------#
# 結果の描画
#----------------------------------------------------------#

# プロットする最大試行数 (50にすると本文の図のように50試行までをプロット)
maxtrial <- Inf
ggplot() + theme_set(theme_bw(base_size = 18,base_family="Arial")) 

x11()
gQ <- list()

df <- data.frame(trials = 1:T, 
                 Q1est = Qest[1,],
                 Q2est = Qest[2,],
                 Q1 = Q[1,],
                 Q2 = Q[2,],
                 c = c,
                 r = as.factor(r),
                 p1 = p1,
                 p1est = p1est)

dfplot <- df %>% filter(trials <= maxtrial)

# Q(A) の図の作成
idxc <- 1
g_qvalues <- ggplot(dfplot, aes(x = trials, y = Q1)) +
  geom_line(aes(y = Q1), linetype = 1, size=1.2) +
  geom_line(aes(y = Q1est), linetype = 5, size=1.0) +
  geom_point(data = dfplot %>% filter(c==idxc & r == 1), 
             aes(x = trials, y = 1.12), shape = 25, size = 1.5) + 
  scale_y_continuous(breaks=c(0, 0.5, 1.0), labels = c(0,0.5,1)) +
  geom_linerange(data = dfplot %>% filter(c==idxc), 
                 aes(
                   x = trials, 
                   ymin = 1.0, 
                   ymax = 1.06), 
                 size=1) + 
  theme(legend.position = "none") + 
  xlab("試行") +
  ylab("Q(A)") 

gQ[[idxc]] <- g_qvalues

# Q(B) の図の作成
idxc <- 2
g_qvalues <- ggplot(dfplot, aes(x = trials, y = Q2)) +
  geom_line(aes(y = Q2), linetype = 1, size=1.2) +
  geom_line(aes(y = Q2est), linetype = 5, size=1.0) +
  geom_point(data = dfplot %>% filter(c==idxc & r == 1), 
             aes(x = trials, y = 1.12), shape = 25, size = 1.5) + 
  geom_linerange(data = dfplot %>% filter(c==idxc), 
                 aes(
                   x = trials, 
                   ymin = 1.0, 
                   ymax = 1.06), 
                 size=1) + 
  theme(legend.position = "none") + 
  scale_y_continuous(breaks=c(0, 0.5, 1.0), labels = c(0,0.5,1)) +
  xlab("試行") +
  ylab("Q(B)")

gQ[[idxc]] <- g_qvalues

# P(a=A) の図の作成
g_p1 <- ggplot(dfplot, aes(x = trials, y = p1)) + 
  xlab("試行") + 
  ylab("P(a = A)") +
  geom_line(aes(y = p1), linetype = 1, size=1.2) +
  geom_line(aes(y = p1est), linetype = 2, size=1.0) +
  geom_point(data = dfplot %>% filter(c==1 & r == 1), 
             aes(x = trials, y = 1.12), shape = 25, size = 1.5) + 
  geom_point(data = dfplot %>% filter(c==2 & r == 1), 
             aes(x = trials, y = -0.12), shape = 2, size = 1.5) + 
  theme(legend.position = "none") + 
  geom_linerange(data = dfplot %>% filter(c==1), 
                 aes(
                   x = trials, 
                   ymin = 1.0, 
                   ymax = 1.05), 
                 size=1) + 
  scale_y_continuous(breaks=c(0, 0.5, 1.0), labels = c(0,0.5,1)) +
  geom_linerange(data = dfplot %>% filter(c==2), 
                 aes(
                   x = trials, 
                   ymin = -0.05, 
                   ymax = 0.0), 
                 size=1) +
  geom_path(size = 1.2) 

grid.arrange(gQ[[1]], gQ[[2]],g_p1, nrow=3) 

# 図を保存する場合は以下を実行
# g <- arrangeGrob(gQ[[1]], gQ[[2]], g_p1, nrow=3) 
# ggsave(file="./figs/qlarning_fit.eps", g) 



