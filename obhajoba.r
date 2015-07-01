# OTAZKA 1:
v1 <- run.q.span(N = 200, n = 100, k = 70, 1, 1, 5, goe.stan)

n = 100; k = 70;
for (q in (1:5)){
	matica1 <- v1$mean.V2s[[q]]
	matica2 <- ( sqrt(70/100)*v1$mean.eigs1[[q]]%*%t(rep(1,k*q)) - t(v1$mean.eigs2[[q]]%*%t(rep(1,n*q))))^2
	persp(z=matica1*matica2,
	      ticktype='detailed',
	      col = surf.colors(matica1*matica2,topo.colors(10)),theta = -45,phi = 10,
	      main = paste("porovnanie odhadu pre GOE (n =",n*q,", k = ",k*q,")"),
	      xlab = "",
	      ylab = "",
	      zlab = ""
	      #zlab = "V^2/odhad"
	      #zlab = "V^2 / (L_j - sqrt(k/n)*L_i)^2"
	)
}

v2 <- run.k.span(N = 200, n = 200, k0 = 100, dk = 10, k1 = 200, goe.stan)
n = 200; k0 = 100; dk = 10; k1 = 200
i = 1
for (k in seq(from = k0, to = k1, by = dk) ){
	matica1 <- v2$mean.V2s[[i]]
	matica2 <- ( sqrt(k/n)*v2$mean.eigs1[[i]]%*%t(rep(1,k)) - t(v2$mean.eigs2[[i]]%*%t(rep(1,n))))^2
	persp(z=matica1*matica2,
	      ticktype='detailed',
	      col = surf.colors(matica1*matica2,topo.colors(10)),theta = -45,phi = 10,
	      main = paste("porovnanie odhadu pre GOE (n =",n,", k = ",k,")"),
	      xlab = "",
	      ylab = "",
	      zlab = ""
	      #zlab = "V^2/odhad"
	      #zlab = "V^2 / (L_j - sqrt(k/n)*L_i)^2"
	)
	i = i+1
}
# OBRAZKY NA OTAZKU 2:
v <- run.detail(N = 1000, n = 100, m = 150, k = 70, wish.stan)
g <- run.detail(N = 1000, n = 100, k = 70, goe.stan)
persp3d(x = sort(g$mean.eig.val1), y = sort(g$mean.eig.val2), z=g$mean.V[100:1,70:1])

par(mfrow=c(1,2))
persp(z=v$mean.V[100:1,70:1],
      zlim = c(0,0.65),ticktype='detailed',
      col = surf.colors(v$mean.V[100:1,70:1],topo.colors(10)),
      main = expression(paste('|V|', " Wishart (n = 100, m = 150, k = 70)")),
      xlab = "<0,1>",
      ylab = "<0,1>"
)

persp(x = sort(v$mean.eig.val1), y = sort(v$mean.eig.val2), z=v$mean.V[100:1,70:1],
      zlim = c(0,0.65),ticktype='detailed',
      col = surf.colors(v$mean.V[100:1,70:1],topo.colors(10)),
      main = expression(paste('|V|', " Wishart (n = 100, m = 150, k = 70) - natiahnute")),
      xlab = expression(lambda_i),
      ylab = expression(lambda_j)
)

dev.off()

par(mfrow=c(1,2))
persp(z=g$mean.V[100:1,70:1],
      zlim = c(0,0.65),ticktype='detailed',
      col = surf.colors(g$mean.V[100:1,70:1],topo.colors(10)),
      main = expression(paste('|V|', " GOE (n = 100, k = 70)")),
      xlab = "<0,1>",
      ylab = "<0,1>"
)

persp(x = sort(g$mean.eig.val1), y = sort(g$mean.eig.val2), z=g$mean.V[100:1,70:1],
      zlim = c(0,0.65),ticktype='detailed',
      col = surf.colors(g$mean.V[100:1,70:1],topo.colors(10)),
      main = expression(paste('|V|', " GOE (n = 100, k = 70) - natiahnute")),
      xlab = expression(lambda_i),
      ylab = expression(lambda_j)
)
dev.off()

par(mfrow=c(1,2))
persp(x = sort(g$mean.eig.val1), y = sort(g$mean.eig.val2), z=g$mean.V[100:1,70:1]^(-2),
      ticktype='detailed',
      col = surf.colors(g$mean.V[100:1,70:1]^(-2),topo.colors(10)),
      main = expression(paste('|V|'^bold("-2"), " GOE (n = 100, k = 70) - natiahnute")),
      xlab = expression(lambda_i),
      ylab = expression(lambda_j)
)

persp(x = sort(v$mean.eig.val1), y = sort(v$mean.eig.val2), z=v$mean.V[100:1,70:1]^(-2),
       ticktype='detailed',
       col = surf.colors(v$mean.V[100:1,70:1]^(-2),topo.colors(10)),
       main = expression(paste('|V|'^bold("-2"), " Wishart (n = 100, m = 150, k = 70) - natiahnute")),
       xlab = expression(lambda_i),
       ylab = expression(lambda_j)
)
dev.off()

# USEKNUTA SRANDA
V = v$mean.V[100:1,70:1]^(-2)
V[V>5000] = 5000
persp(x = sort(v$mean.eig.val1), y = sort(v$mean.eig.val2), z=V,
      ticktype='detailed',
      zlim = c(0,5000),
      col = surf.colors(V,topo.colors(10)),
      theta = -45,
      main = expression(paste('|V|'^bold("-2"), " Wishart (n = 100, m = 150, k = 70) - natiahnute")),
      xlab = expression(lambda_i),
      ylab = expression(lambda_j)
)
