goe.stan <- function(n, k, q){
	# Generate matrices
	g1 <- generate.GUE(n = n*q, compl = F, norm = T)
	p  <- generate.Q.stan(n = n*q, k = k*q, random = T)
	g2 <- h(p) %*% g1 %*% p
	
	# Spectral decomposition
	eig1 <- eigen(x = g1, symmetric = T)
	eig2 <- eigen(x = g2, symmetric = T)
	vec1 <- h(p) %*% eig1$vector
	vec1 <- normalize(vec1)
	vec2 <- eig2$vectors
	
	# Extract info
	pV <- h(vec2) %*% vec1
	pPPT <- p %*% h(p)
	p.diag.max <- diag.max(abs(pV))
	return(list(V = pV, 
		    PPT = pPPT,
		    eig.val1  = eig1$values, 
		    eig.val2 = eig2$values , 
		    diag = p.diag.max))
}

run.q.span <- function(N, n, k, q0, dq, q1, exper){
	qs <- seq(q0, q1, dq)
	len <- length(qs)
	
	mean.Vs <- vector("list", len)
	mean.V2s <- vector("list", len)   # priemer kvadratur
	mean.PPTs <- vector("list", len)  # Vnutorny nasobok matic
	mean.PPT2s <- vector("list", len)  # Vnutorny nasobok matic
	mean.eigs1 <- vector("list", len)
	mean.eigs2 <- vector("list", len)
	mean.diags <- vector("list", len)
	
	for (q in qs){
		print(paste(q0,'/',q,'/',q1))
		# Inicializovat polia:
		V <- vector("list", N)
		PPT <- vector("list", N)
		eig.val1 <- vector("list", N)
		eig.val2 <- vector("list", N)
		diag <- vector("list", N)
		
		for (i in 1:N){
			print(paste('zacinam: ',i,'/',N,' -> q:',q))
			tmp <- exper(n, k, q)
			V[[i]] <- tmp$V
			PPT[[i]] <- tmp$PPT
			eig.val1[[i]] <- tmp$eig.val1
			eig.val2[[i]] <- tmp$eig.val2
			diag[[i]] <- tmp$diag
		}
		V <- abind(V,along = 3)
		PPT <- abind(PPT,along = 3)
		eig.val1 <- abind(eig.val1, along = 2)
		eig.val2 <- abind(eig.val2, along = 2)
		diag <- abind(diag, along = 2)
		print(paste('Skoncene iteracie. Zacinam robit priemery pre q: ', q))
		
		iter <- (q-q0)/dq + 1   # kolkate q?
		
		mean.Vs[[iter]] <- apply(abs(V), c(1,2), mean)
		mean.V2s[[iter]] <- apply(V^2, c(1,2), mean)
		mean.PPTs[[iter]] <- apply(abs(PPT), c(1,2), mean)
		mean.PPT2s[[iter]] <- apply(PPT^2, c(1,2), mean)
		mean.eigs1[[iter]] <- apply(eig.val1, 1, mean)
		mean.eigs2[[iter]] <- apply(eig.val2, 1, mean)
		mean.diags[[iter]] <- apply(diag, 1, mean)
	}
	print('KONIEC')
	return(list(    mean.Vs =  mean.Vs, 
			mean.V2s = mean.V2s, 
			mean.PPTs = mean.PPTs, 
			mean.PPT2s = mean.PPT2s, 
			mean.eigs1 = mean.eigs1, 
			mean.eigs2 = mean.eigs2, 
			mean.diags = mean.diags))
}

run.k.span <- function(N, n, k0, dk, k1, exper){
	ks <- seq(k0, k1, dk)
	len <- length(ks)
	
	mean.Vs <- vector("list", len)
	mean.V2s <- vector("list", len)   # priemer kvadratur
	mean.PPTs <- vector("list", len)  # Vnutorny nasobok matic
	mean.PPT2s <- vector("list", len)  # Vnutorny nasobok matic
	mean.eigs1 <- vector("list", len)
	mean.eigs2 <- vector("list", len)
	mean.diags <- vector("list", len)
	
	for (k in ks){
		print(paste(k0,'/',k,'/',k1))
		# Inicializovat polia:
		V <- vector("list", N)
		PPT <- vector("list", N)
		eig.val1 <- vector("list", N)
		eig.val2 <- vector("list", N)
		diag <- vector("list", N)
		
		for (i in 1:N){
			print(paste('zacinam: ',i,'/',N,' -> k:',k))
			tmp <- exper(n, k, 1)
			V[[i]] <- tmp$V
			PPT[[i]] <- tmp$PPT
			eig.val1[[i]] <- tmp$eig.val1
			eig.val2[[i]] <- tmp$eig.val2
			diag[[i]] <- tmp$diag
		}
		V <- abind(V,along = 3)
		PPT <- abind(PPT,along = 3)
		eig.val1 <- abind(eig.val1, along = 2)
		eig.val2 <- abind(eig.val2, along = 2)
		diag <- abind(diag, along = 2)
		print(paste('Skoncene iteracie. Zacinam robit priemery pre k: ', k))
		
		iter <- (k-k0)/dk + 1   # kolkate k?
		
		mean.Vs[[iter]] <- apply(abs(V), c(1,2), mean)
		mean.V2s[[iter]] <- apply(V^2, c(1,2), mean)
		mean.PPTs[[iter]] <- apply(abs(PPT), c(1,2), mean)
		mean.PPT2s[[iter]] <- apply(PPT^2, c(1,2), mean)
		mean.eigs1[[iter]] <- apply(eig.val1, 1, mean)
		mean.eigs2[[iter]] <- apply(eig.val2, 1, mean)
		mean.diags[[iter]] <- apply(diag, 1, mean)
	}
	print('KONIEC')
	return(list(    mean.Vs =  mean.Vs, 
			mean.V2s = mean.V2s, 
			mean.PPTs = mean.PPTs, 
			mean.PPT2s = mean.PPT2s, 
			mean.eigs1 = mean.eigs1, 
			mean.eigs2 = mean.eigs2, 
			mean.diags = mean.diags))
}