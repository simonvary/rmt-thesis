test.q.symmetric <- function(V, test, q.ave =F){
	# V is assumed to be n x k x N matrix 
	# S <- V[,c(70:1),] , for other diag.
	n = dim(V)[1]
	k = dim(V)[2]
	N = dim(V)[3]
	q = n/k
	no.fail <- 0
	p.vals <- matrix(numeric(n*k), nrow = n)
	for (i in (1:n)){
		for (j in (1:k)){
			X <- V[i,j,]
			i2 <- j*q
			j2 <- i/q
			if (q.ave){
				#Nefunguje dobre ! nepouzivat
				#norm.q <- (i2-i2.low)*(j2-j2.low)^2 + (i2-i2.low)*(j2.upp-j2)^2 + (i2.upp-i2)*(j2-j2.low)^2 + (i2.upp-i2)*(j2.upp-j2)^2
				i2.upp <- ceiling(i2)
				i2.low <- floor(i2)
				j2.upp <- ceiling(j2)
				j2.low <- floor(j2)
				Y <- c( V[i2.low, j2.low, (1:round( ((i2-i2.low)*(j2-j2.low))*N ))], 
					    V[i2.low, j2.upp, (1:round( ((i2-i2.low)*(j2.upp-j2))*N ))],
					    V[i2.upp, j2.low, (1:round( ((i2.upp-i2)*(j2-j2.low))*N ))],
					    V[i2.upp, j2.upp, (1:round( ((i2.upp-i2)*(j2.upp-j2))*N ))]
				)
			}else{
				Y <- V[round.max(i2,n), round.max(j2,k),]	
			}
			p.vals[i,j] <- test(X,Y)$p.value
			if (p.vals[i,j] < 0.05){
				no.fail = no.fail + 1
				print(paste('Failol: ',i,',',j,' test. Failnutych: ',no.fail,'/',n*k))
				print(no.fail/(i*k-k + j))
				print(p.vals[i,j])	
			}
		}
	}
	return(p.vals)
}

test.normality <- function(V, test, msg =F){
	# V is assumed to be n x k x N matrix 
	# S <- V[,c(70:1),] , for other diag.
	n = dim(V)[1]
	k = dim(V)[2]
	N = dim(V)[3]
	q = n/k
	no.fail <- 0
	p.vals <- matrix(numeric(n*k), nrow = n)
	for (i in (1:n)){
		for (j in (1:k)){
			X <- V[i,j,1:3000]
			i2 <- j*q
			j2 <- i/q
			print(paste('Robim: ',i,',',j))
			p.vals[i,j] <- test(X)$p.value
			if ((p.vals[i,j] < 0.05) & (msg)){
				no.fail = no.fail + 1
				print(paste('Failol: ',i,',',j,' test. Failnutych: ',no.fail,'/',n*k))
				print(no.fail/(i*k-k + j))
				print(p.vals[i,j])	
			}
		}
	}
	return(p.vals)
}

test.id <- function(V1, V2, no, test){
	# V1, V2 are assumed to be n x k x N matrix 
	n = dim(V1)[1]
	k = dim(V1)[2]
	N = dim(V1)[3]
	no.fail <- 0
	p.vals <- matrix(numeric(n*k), nrow = n)
	for (i in (1:n)){
		for (j in (1:k)){
			X <- V1[i,j,1:no]
			Y <- V2[i,j,1:no]	
			p.vals[i,j] <- test(X,Y)$p.value
			if (p.vals[i,j] < 0.05){
				no.fail = no.fail + 1
				print(paste('Failol: ',i,',',j,' test. Failnutych: ',no.fail,'/',n*k))
				print(no.fail/(i*k-k + j))
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

MLE.value <- function(X, nbins=50){
	hi <- hist(X, breaks = nbins, plot = F)
	return(hi$mids[which.max(hi$counts)])
}

# N = 1000
# n = 100
# k = 100
# x = matrix(numeric(n*k*N), nrow = n)
# B <- generate.Q.haar(n,n)
# for (i in 1:N){
# 	Q <- generate.Q.haar(n,k)
# 	#e <- eigen(Q,symmetric = F,only.values = T)$values
# 	#x[(i*2-1):(i*2)] <- Re(e)
# 	#y[(i*2-1):(i*2)] <- Im(e)
# 	for (j in (1:n)){
# 		x[j,(i*k-k+1):(i*k)] <- Q[j,] %*% B[j,]	
# 	}
# }

# k = 150
# s<-run.detail(N = 6000,m = 1, n = 200,k = k, goe.stan)
# MLE.value(s$V[1,1,])
# MLE.value(s$V[200,k,])
