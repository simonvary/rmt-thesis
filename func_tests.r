V.symm.test <- function(V, test){
	# V is assumed to be n x k x N3d matrix 
	n = dim(V)[1]
	k = dim(V)[2]
	N = dim(V)[3]
	q = n/k
	p.vals <- matrix(numeric(n*k), nrow = n)
	for (i in (1:n)){
		for (j in (1:k)){
			print(paste('Zacinam: ',i,',',j,' test.'))
			p.vals[i,j] <- test(abs(V[i,j,])^2, abs(V[round.max(j*q,n),round.max(i/q,k),])^2 )$p.value
			if (p.vals[i,j] < 0.05){
				print(p.vals[i,j])	
			}
		}
	}
	return(p.vals)
}

V.symm.mean <- function(V){
	n = dim(V)[1]
	k = dim(V)[2]
	q = n/k
	V.diff <- matrix(numeric(n*k), nrow = n)
	for (i in (1:n)){
		for (j in (1:k)){
			V.diff[i,j] = V[i,j]-V.mean[round.max(j*q,n),round.max(i/q,k)]
		}
	}
	return(V.diff)
}

N = 1000
n = 100
k = 100
x = matrix(numeric(n*k*N), nrow = n)
B <- generate.Q.haar(n,n)
for (i in 1:N){
	Q <- generate.Q.haar(n,k)
	#e <- eigen(Q,symmetric = F,only.values = T)$values
	#x[(i*2-1):(i*2)] <- Re(e)
	#y[(i*2-1):(i*2)] <- Im(e)
	for (j in (1:n)){
		x[j,(i*k-k+1):(i*k)] <- Q[j,] %*% B[j,]	
	}
}
