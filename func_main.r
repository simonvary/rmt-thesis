require(pracma)
require(abind)

pnorm <- function(x, p = 2){
	# p-norm of a vectors x
	return( sum(x^p)^(1/p) )
}

normalize <- function(A, bycol = T, p = 2){
	# Normalize columns/rows of A, using p-norm
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

diag.max <- function(M){
	# Extract max. values in every column
	return(apply(M, MARGIN = 2,FUN = max))
}


sample.semicircle <- function(nsamples = 1, eig = F){
	# Sample from semicircle distribution
	# 	eig == F -> n independent samples (without eig. repulsion)
	#	eig == T -> n dependent random eigenvalues of random GUE matrix (with eig. repulsion)
	if (eig == F){
		# http://stats.stackexchange.com/questions/12843/generating-random-samples-from-a-custom-distribution
		# sample independent random numbers from pdf: p(x) = (2/pi)*sqrt(1-x^2)
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

fun.semicircle <- function(x, m=1){
	# semicircle function with support [-m^2, m^2]
	return((2/(pi*m^2))*sqrt(m^2-x^2))
}

fun.marcenko.pastur <- function(x,q){
	# Marchenko-Pastur function with parameter q
	a = (1-sqrt(q))^2
	b = (1+sqrt(q))^2
	return(1/(2*pi*x*q)*sqrt((b-x)*(x-a)) )
}

generate.GUE <- function(n = 1000, compl = F, norm = T){
	# Generate GUE matrix:
	# 	 compl == T/F -> GUE/GOE
	#	 norm == T 	  -> spectrum with support [-1,1]
	if (compl){
		M <- matrix( data = complex(real = rnorm(n^2, sd = 0.5), imaginary = rnorm(n^2, sd = 0.5)),
			     ncol = n,
			     nrow = n);
		diag(M) <- rnorm(n)
	}else{
		M <- matrix(data = rnorm(n^2), ncol = n, nrow = n)
	}
	M = .5*(M + h(M))
	if (norm) M = M/(sqrt(2*n))
	return(M)
}

generate.GUE.diag <- function(n = 1000, compl = F, norm = T){
	# naozaj potrebujem ? asi vymazat
	M <- diag(rnorm(n))
	if (norm) M = M/sqrt(2*n)
	return(M)
}

generate.wigner.unif <- function(n = 100, norm = T){
	# Generate Wigner matrix with unifrom off-diagonal and Bernoulli diagonal
	M <- matrix(data = runif(n^2, -sqrt(12)/2, sqrt(12)/2), ncol = n, nrow = n)
	diag(M) <- sample(c(-30,30), size = n, replace = T)
	M = .5*(M + h(M))
	if (norm) M = M/sqrt(2*n)
}

generate.wigner.bernoulli <- function(n = 100, norm = T){
	# Generate Wigner Bernoulli ensemble
	# 	 norm == T -> spectrum with support [-1,1]
	M <- matrix(data = sample(c(-1,1),size = n^2, replace = T), ncol = n, nrow = n)
	M = .5*(M + h(M))
	if (norm) M = M/sqrt(2*n)
}

generate.iid <- function(n = 100, norm = T){
	# Generate random iid (n,n) matrix with Gaussian dist.
	M <- matrix(data = rnorm(n^2), ncol = n, nrow = n)
	if (norm) M <- M/sqrt(2*n)
	return(M)
}

generate.wishart <- function(n = 100, m = 150, norm = T){
	# Generate (n,n) Wishart matrix with q = m/n
	X <- matrix(data = rnorm(n*m), ncol = m, nrow = n)
	M <- X %*% t(X)
	if (norm) M <- M/m
	return(M)
}

generate.Q.haar <- function(n, k, random = T){
	# Generate (n,k) matrix whose columns are orthogonal and Haar distributed
	M <- generate.GUE(n = n,compl = F,norm = T)
	return(eigen(M, symmetric = T, only.values = F)$vectors[,1:k])
}

generate.Q.stan <- function(n, k, random = T){
	# Generate (n,k) matrix whose columns are orthogonal, standard basis vectors
	# 	  random == T/F -> random basis vectors/ first k basis vectors
	if (random){
		P <- diag(rep(x = 1,n))[,sample(x = 1:n,size = k)]
	}else{
		P <- diag(rep(x = 1,n))[,(1:k)]
	}
	return(P)
}

generate.P.haar <- function(n, k){
	# Generate projection matrix on k orthogonal, Haar distributed vectors
	Q <- generate.Q.haar(n = n, k = k)
	P <- Q %*% h(Q)
	return(P)
}

generate.P.stan <- function(n, k, random = T){
	# Generate projection matrix on random k orthogonal vectors from standard basis
	Q <- generate.Q.stan(n = n, k = k, random = T)
	P <- Q %*% h(Q)
	return(P)
}

generate.even.spectrum <- function(n, a, b){
	# generates matrix with evenly distributed spectrum in [a,b] and orthonormal eigvectors
	eigs <- seq(from = a, to = b,length.out = n)
	Q <- generate.Q.haar(n = n,k = 0)
	return(h(Q) %*% diag(eigs) %*% Q)
}

round.max <- function(q, q.max){
	# round q with upper bound q.max
	return(min(round(q), q.max))
}

