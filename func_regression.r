pol.reg <- function(V, kp = 3){
  V <- 1/V
  k <- dim(V)[1]
  n <- dim(V)[2]
  #Y <- (-(k/2):((k/2)-1)+0.5)
  #X <- (-(n/2):((n/2)-1)+0.5)
  Y <- eig.val2.mean[1:k]
  X <- eig.val1.mean[1:n]
  Y <- as.vector(sapply(Y, function(x) rep(x,n)))
  X <- as.vector(t(sapply(X, function(x) rep(x,k))))
  Z <- as.vector(t(V))
  pol <- polym(X,Y,degree = kp)
  model <- lm(Z ~ pol, x = T) 
  #fit.val <- 
  #plot3d(X,Y,Z)
}