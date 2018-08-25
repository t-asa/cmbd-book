
#----------------------------------------------------------#
# Plot results
#----------------------------------------------------------#

library(ggplot2)
library(gridExtra)

ggplot() + theme_set(theme_bw(base_size = 18,base_family="Arial")) 

x11()

gQ <- list()

maxtrial <- 50
trials_plot <- 1:maxtrial

df <- data.frame(trials = 1:T, 
                 Q1est = Qest[1,],
                 Q2est = Qest[2,],
                 Q1 = Q[1,],
                 Q2 = Q[2,],
                 c = c,
                 r = as.factor(r),
                 p1 = p1,
                 p1est_wsls = p1est_wsls,
                 p1est = p1est, 
                 p1est_random = p1est_random)

dfplot <- df %>% filter(trials <= 50)


idxc <- 1
g_qvalues <- ggplot(dfplot, aes(x = trials, y = Q1)) +
  # geom_path(size = 1.2) + 
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
  #theme_bw() + 
  theme(legend.position = "none") + 
  xlab("Trial") +
  ylab("QA") 

gQ[[idxc]] <- g_qvalues

idxc <- 2
g_qvalues <- ggplot(dfplot, aes(x = trials, y = Q2)) +
  # geom_path(size = 1.2) + 
  # geom_line(aes(y = Q2, linetype = "true"), size=1.2) +
  # geom_line(aes(y = Q2est, linetype = "fit"), size=1.2) +
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
  #theme_bw() + 
  theme(legend.position = "none") + 
  scale_y_continuous(breaks=c(0, 0.5, 1.0), labels = c(0,0.5,1)) +
  xlab("Trial") +
  ylab("QB")

gQ[[idxc]] <- g_qvalues

g_p1 <- ggplot(dfplot, aes(x = trials, y = p1)) + 
  xlab("Trial") + 
  ylab("Prob. choosing 1") +
  geom_line(aes(y = p1), linetype = 1, size=1.2) +
  geom_line(aes(y = p1est), linetype = 2, size=1.0) +
  #geom_line(aes(y = p1est_wsls), linetype = 2, size=1.0, color="red") +
  #geom_line(aes(y = p1est_random), linetype = 2, size=1.0, color="green") +
  geom_point(data = dfplot %>% filter(c==1 & r == 1), 
             aes(x = trials, y = 1.12), shape = 25, size = 1.5) + 
  geom_point(data = dfplot %>% filter(c==2 & r == 1), 
             aes(x = trials, y = -0.12), shape = 2, size = 1.5) + 
  #theme_bw() + 
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

g <- arrangeGrob(gQ[[1]], gQ[[2]], g_p1, nrow=3) #generates g
ggsave(file="./figs/qlarning_fit.eps", g) #saves g


#----------------------#
# likelihood comparison
#----------------------#
x11(width = 7, height = 3)

g_p1 <- ggplot(df, aes(x = trials, y = p1est)) + 
  xlab("Trial") + 
  ylab("Prob. choosing 1") +
  # geom_line(aes(y = p1), linetype = 1, size=1.2) +
  geom_line(aes(y = p1est), linetype = 2, size=1.0) +
  geom_line(aes(y = p1est_wsls), linetype = 1, size=1, color="gray44") +
  geom_line(aes(y = p1est_random), linetype = 3, size=1, color="gray55") +
  geom_point(data = df %>% filter(c==1 & r == 1), 
             aes(x = trials, y = 1.12), shape = 25, size = 1.5) + 
  geom_point(data = df %>% filter(c==2 & r == 1), 
             aes(x = trials, y = -0.12), shape = 2, size = 1.5) + 
  #theme_bw() + 
  theme(legend.position = "none") + 
  geom_linerange(data = df %>% filter(c==1), 
                 aes(
                   x = trials, 
                   ymin = 1.0, 
                   ymax = 1.05), 
                 size=1) + 
  scale_y_continuous(breaks=c(0, 0.5, 1.0), labels = c(0,0.5,1)) +
  geom_linerange(data = df %>% filter(c==2), 
                 aes(
                   x = trials, 
                   ymin = -0.05, 
                   ymax = 0.0), 
                 size=1) +
  geom_path(size = 1.2) 
print(g_p1)
ggsave(file="./figs/qlarning_ll_comparison.eps", g_p1) 
