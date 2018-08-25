
# Contour map ================================= # 

alpha <- seq(0.01,1.0, by = 0.01)
beta <- seq(0.1,4.0, by = 0.1)

llval <- matrix(NA, nrow = length(alpha), ncol = length(beta))

for (idxa in 1:length(alpha)) {
  for (idxb in 1:length(beta)) {
    llval[idxa,idxb] <- - func_minimize(c(alpha[idxa], beta[idxb]), modelfunc = func_qlearning, 
                                        choice=c, reward=r)
  }
}

llval <- pmax(llval,-100)

x11()

image(alpha, beta, llval,
      xlim = c(min(alpha),max(alpha)),
      ylim = c(min(beta),max(beta)),
      xlab = expression(
        paste("alpha")
      ),
      ylab=expression(
        paste("beta")
      ),
      #col = heat.colors(50),
      col = gray.colors(50),
      main = sprintf("Log likelihood"),
      las = 1)


contour(alpha, beta, llval, nlevels = 20,
        # add = T,
        xlim = c(min(alpha),max(alpha)),
        ylim = c(min(beta),max(beta)),
        xlab = expression(
          paste("alpha")
        ),
        ylab=expression(
          paste("beta")
        ),
        main = sprintf("Log likelihood"),
        las = 1,
        # asp=1,
        labcex = 1.1,
        add = T)

# text indicating the estimate
points(paramest[1], paramest[2], pch = "*")
text(paramest[1] + 0.07, paramest[2] - 0.1, sprintf("\n LL = %.2f", - ret$negll ))

x11()
persp(alpha, beta, llval)


require(plot3D)  

x11()
#plot3d(alpha, beta, llval)
plot3D::persp3D(alpha, beta, llval, 
                # lighting = TRUE, 
                phi = 25,
                theta = 40, 
                # alpha = 0.2, 
                # image = list(col = grey (seq(0, 1, length.out = 100))), 
                colvar = llval,
                col = "gray", # ltheta = 120, shade = 0.75,
                ticktype = "detailed", 
                xlab = "alpha",
                ylab = "beta",
                zlab = "Log-likelihood",
                zlim = c(-90,-40), 
                axes=TRUE, 
                contour = list(col = "black", nlevels = 25, side = c("zmin", "z"))
                )
