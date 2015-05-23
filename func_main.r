require(pracma)
require(abind)
require(rgl)

pnorm <- function(x, p = 2){
  return( sum(x^p)^(1/p) )
}

normalize <- function(A, bycol = T, p = 2){
  if (bycol){
     return(apply(X = A, MARGIN = 2, FUN = function(x){
       return(x/pnorm(x, p = p))
     }))
  }else{
    return( t( normalize( t(A), bycol=T, p = p ) ) )
  }
}

h <- function(x){
  return(Conj(t(x)))
}

diag2 <- function(M){
  return(diag(apply(M, MARGIN = 1,FUN = rev)))
}

diag.max <- function(M){
  return(apply(M, MARGIN = 2,FUN = max))
}


sample.semicircle <- function(nsamples = 1, eig = F){
  if (eig == F){
    # http://stats.stackexchange.com/questions/12843/generating-random-samples-from-a-custom-distribution
    # sample random numbers from p.d.f: p(x) = (2/pi)*sqrt(1-x^2)
    x <- runif(nsamples)
    f <- function(x,u) (x*sqrt(1 - x^2) + asin(x) + pi/2 )/pi - u
    my.uniroot <- function(x) uniroot(f, c(-1, 1), tol = 1e-4, u = x)$root;
    r <- vapply(x, my.uniroot, numeric(1))
    return(r)
  }else{
    # sample eigenvalues of GUE matrix
    gue <- generate.GUE(n = nsamples, compl = F, norm = T)
    return(eigen(x = gue, symmetric = T, only.values = T)$values)
  }
}

generate.bernoulli <- function(n = 1000, norm = T){
  M <- matrix(numeric(n^2), ncol = n, nrow = n)
  M[upper.tri(M)] <- sample(c(-1,1), (n*(n-1))/2, replace = T)
  M[lower.tri(M)] <- M[upper.tri(M)]
  diag(M) <- sample(c(-1,1), n, replace = T)
  if (norm) M = M/sqrt(2*n)
  return(M)
}

generate.GUE <- function(n = 1000, compl = F, norm = T){
  if (compl){
    M <- matrix( data = complex(real = rnorm(n^2, sd = 0.5),
            imaginary = rnorm(n^2, sd = 0.5)),
            ncol = n,
            nrow = n);
    diag(M) <- rnorm(n)
  }else{
    M <- matrix(data = rnorm(n^2), ncol = n, nrow = n)
    #diag(M) <- rnorm(n,sd = 10)
  }
  M = .5*(M + h(M))
  if (norm) M = M/sqrt(2*n)
  return(M)
}

generate.GUE.diag <- function(n = 1000, compl = F, norm = T){
  # naozaj potrebujem ?
  M <- diag(rnorm(n))
  if (norm) M = M/sqrt(2*n)
  return(M)
}

generate.iid <- function(n = 100, norm = T){
  # Generate random iid (n,n) matrix with Gaussian dist.
  M <- matrix(data = rnorm(n^2), ncol = n, nrow = n)
  if (norm) M = M/sqrt(2*n)
  return(M)
}

generate.wishart <- function(n = 100, k = 150){
  X <- matrix(data = rnorm(n*k), ncol = k, nrow = n)
  return(X %*% t(X))
}

generate.Q.haar <- function(n, k, random = T){
  # k - orthogonal vectors using Haar measure (R^n)
  M <- generate.GUE(n = n,compl = F,norm = T)
  return(eigen(M, symmetric = T, only.values = F)$vectors[,1:k])
}

generate.Q.stan <- function(n, k, random = T){
  # generate k - standard orthogonal vectors
  if (random){
    P <- diag(rep(x = 1,n))[,sample(x = 1:n,size = k)]
  }else{
    #P <- diag(rep(x = 1,n))[,-((n-k+1):n)]
    P <- diag(rep(x = 1,n))[,(1:k)]
  }
  return(P)
}

generate.Q.GUE <- function(n, k, G){
  # G is GUE
  #G <- generate.GUE(n = n, compl = F, norm = T)
  Q <- generate.Q.haar(n = n, k = k)
  return(G %*% Q)
}

generate.P.haar <- function(n, k){
  Q <- generate.Q.haar(n = n, k = k)
  P <- Q %*% h(Q)
  return(P)
}

generate.P.stan <- function(n, k, random = T){
  # generate projection on random standard subspace
  Q <- generate.Q.stan(n = n, k = k, random = T)
  P <- Q %*% h(Q)
  return(P)
}

generate.even.spectrum <- function(n, a, b){
  # generates matrix with spectrum evenly distributed in [a,b] and orthonormal eigvectors
  eigs <- seq(from = a, to = b,length.out = n)
  Q <- generate.Q.haar(n = n,k = 0)
  #Q <- diag(rep(1,n))
  return(h(Q) %*% diag(eigs) %*% Q)
}

create.sum.matrix <- function(n,k,q){
  # generate pair (left/right) of helping matrices
  # For matrix (q*k x q*n) V :
  #       V2 := L %*% V %*% R
  pR <- c(rep(c(rep(1/q,q), rep(0,q*n)), n-1), rep(1/q,q) )
  pR <- matrix(data = pR,nrow = q*n, ncol = n, byrow = F)
  pL <- c(rep(c(rep(1/q,q), rep(0,q*k)), k-1), rep(1/q,q) )
  pL <- t(matrix(data = pL,nrow = q*k, ncol = k, byrow = F))
  return(list(L = pL, R = pR))
}

draw <- function(no){
  q <- matrix(numeric(no*3*3),nrow = 3)
  P <- generate.P.stan(3,1)
  G <- generate.GUE(n = 3, compl = F, norm = T)
  for (i in 1:no){
    Q <- generate.Q.haar(n = 3, k = 3)
    q[,c(i*3-2,i*3-1,i*3)] <- P %*% G %*% Q
  }
  plot3d(t(q))
}

