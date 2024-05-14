

## Density multivariate Skew-normal to contour graphic

denMSN <-function(y, mu, Sigma, shape){

  m<-nrow(y)
  p<-ncol(y)
  ver<-matrix(0,m,1)
  mu<-matrix(mu,p,1)

  for (i in 1:m){
    y1=matrix(y[i,],p,1)
    ver[i] <- dmvESN(x = as.vector(y1),mu = as.vector(mu),Sigma = Sigma,lambda = shape,tau = 0)
  }
  return(ver)
}


denMN <-function(y, mu, Sigma){

  m<-nrow(y)
  p<-ncol(y)
  ver<-matrix(0,m,1)
  mu<-matrix(mu,p,1)

  for (i in 1:m){
    y1=matrix(y[i,],p,1)
    ver[i] <- mvtnorm::dmvnorm(x = as.vector(y1),mean = as.vector(mu),sigma = Sigma)
  }
  return(ver)
}


#denMSN(y = test1$y,mu = mu[[1]],Sigma = Sigma[[1]],shape = shape[[1]])
#denMN(y = test1$y,mu = mu[[1]],Sigma = Sigma[[1]])


denMT <- function(y, mu, Sigma, nu){
  #y: deve ser uma matrix onde cada linha tem um vetor de dados multivariados de dimensao ncol(y) = p. nrow(y) = tamanho da amostra
  #mu: deve ser do tipo vetor de mesma dimensao igual a ncol(y) = p
  #Sigma: Matrix p x p
  m<-nrow(y)
  p<-ncol(y)
  ver<-matrix(0,m,1)
  mu<-matrix(mu,p,1)

  for(i in 1:m){
    y1=y[i,]
    ver[i] <- (gamma((p+nu)/2)/(gamma(nu/2)*pi^(p/2)))*nu^(-p/2)*det(Sigma)^(-1/2)*(1 + mahalanobis(y1, mu, Sigma)/nu)^(-(p+nu)/2)
  }
  return(ver)
}

#nu = 3
#denMT(y = test1$y,mu = mu[[1]],Sigma = Sigma[[1]],nu = nu)

## Density mixture multivariate Skew-normal to contour graphic

mixedMSN <- function(y, pii, mu, Sigma, shape){
  #y: Matrix de dados m x p
  #pii : deve ser do tipo vetor de dimensao g
  #mu: deve ser do tipo list com g entradas. Cada entrada do list deve ser um vetor de dimensao p
  #Sigma: deve ser do tipo list com g entradas. Cada entrada do list deve ser uma matriz p x p
  g <- length(pii)
  dens <- 0
  for (j in 1:g) dens <- dens + pii[j]*denMSN(y = y, mu = mu[[j]], Sigma = Sigma[[j]], shape = shape[[j]])
  return(dens)
}


mixedMN <- function(y, pii, mu, Sigma){
  #y: Matrix de dados m x p
  #pii : deve ser do tipo vetor de dimensao g
  #mu: deve ser do tipo list com g entradas. Cada entrada do list deve ser um vetor de dimensao p
  #Sigma: deve ser do tipo list com g entradas. Cada entrada do list deve ser uma matriz p x p
  g <- length(pii)
  dens <- 0
  for (j in 1:g) dens <- dens + pii[j]*denMN(y = y, mu = mu[[j]], Sigma = Sigma[[j]])
  return(dens)
}

## Density mixture multivariate Skew-normal to contour graphic

mixedMT <- function(y, pii, mu, Sigma, nu){
  #y: Matrix de dados m x p
  #pii : deve ser do tipo vetor de dimensao g
  #mu: deve ser do tipo list com g entradas. Cada entrada do list deve ser um vetor de dimensao p
  #Sigma: deve ser do tipo list com g entradas. Cada entrada do list deve ser uma matriz p x p
  g <- length(pii)
  dens <- 0
  for (j in 1:g) dens <- dens + pii[j]*denMT(y = y, mu = mu[[j]], Sigma = Sigma[[j]], nu = nu[j])
  return(dens)
}

#set.seed(20)
#test2 = rmMSN(n = 100,pii = pii,mu = mu,Sigma = Sigma,shape = shape,percen = 0.2,each = FALSE)

# par(mfrow = c(1,2))
# hist(test1$y,breaks=30,probability = T)
# hist(test2$y,breaks=30,probability = T)
# mixedMSN(y = test1$y, pii = pii, mu = mu, Sigma = Sigma , shape = shape)


## Graphic function (contour plot) auxiliary function ##

seq_r = function(n,a,r){
  x = c()
  x[1] = a
  for(i in 2:n){
    x[i] = x[i-1] + (r - (i-2))
  }
  return(x)
}
