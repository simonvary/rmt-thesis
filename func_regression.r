rescale.coordinates <- function(V,X,Y){
	# rescale persp. matrix V based on the X,Y coordinates
	n <- dim(V)[1]
	# n = length(Y)
	k <- dim(V)[2]
	# k = length(X)
	
	pX <- as.vector(t(sapply(X, function(x) rep(x,k))))
	pY <- as.vector(sapply(Y, function(x) rep(x,n)))
	pZ <- as.vector(V)
	return(list(X = pX, Y = pY, Z = pZ))
}

nlm.reg <- function(X,Y,Z,kp){
	pol <- polym(X,Y,degree = kp)
	lm.model <- lm(Z ~ pol, x = T)
	coeff <- lm.model$coefficients
	coeff[summary(lm.model)$coefficient[,4] > 0.05] <- 0
	
	np <- dim(pol)[2]	# number of moments in pol
	f0 <- function(a){
		return( sum( (Z - pol%*%as.vector(a[1:np]) + 1/(a[np+1]*X^2-a[np+2]*Y^2)
			   )^2 ) )
	}
	
	f1 <- function(pX,pY,a){
		pol1 <- polym(pX,pY,degree = kp)
		np1 <- dim(pol1)[2]	# number of moments in pol1
		return( pol1%*%as.vector(a[1:np1]) #+ 
			#	1/(a[np1+1]*pX^2-a[np1+2]*pY^2))
		)
	}
	v <- nlm(f = f0, p = c(coeff,1,1), steptol = 1e-8)
	pfitted.Z <- f1(X, Y, v$estimate)
	return(list(solution = v, fitted.Z = pfitted.Z))
}

pol.reg <- function(X,Y,Z, kp = 3){
	pol <- polym(X,Y,degree = kp)
	model <- lm(Z ~ pol, x = T) 
	#fit.val <- 
	#plot3d(X,Y,Z)
}