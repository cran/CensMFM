
## Information matrix ##
##--- FUNCTION TO IM SIGMA ---##
Der.F <-function(M){

m1<-dim(M)[1]
m2<-dim(M)[2]
d<-list()
for(h in 1:m1){
  d[[h]]<-list()
  for(k in 1:(m2+1-h)){
    d[[h]][[k]]<-matrix(0,m1,m2)
    if(k==1){d[[h]][[k]][h,h]<-1}
    else{
      d[[h]][[k]][h,h+(k-1)]<-d[[h]][[k]][h+(k-1),h]<-1}
  }
}
return(d=d)
}

##--- FUNCTION TO IM LAMBDA ---##
Der.R<-function(R){
  r1 <- length(R)
  d <- list()
  for(k in 1:r1){
    d[[k]] <- matrix(0,r1,r1)
    di <- R
    di[k] <- 2*R[k]
    d[[k]][k,] = di
    d[[k]][,k] = di
  }
  return(d=d)
}

##--- LAMBDA ---##
VecLam = function(x,a){
  out = rep(0,length(x))
  out[a] = 1
  return(out)
}

##------------------------------##
sqrtm <- function(A)
{
  if(length(A)==1)
    Asqrt=sqrt(A)
  else{
    sva <- svd(A)
    if (min(sva$d)>=0){
      Asqrt <- sva$u%*%diag(sqrt(sva$d))%*%t(sva$v)  # svd e decomposi??o espectral
    }else{
      stop("Matrix square root is not defined")
    }
  }
  return(as.matrix(Asqrt))
}

#
invmills = function(x,mu=0,sd=1){
  z = (x-mu)/sd
  if(z < -1e4){
    return(-z/sd)
  }else{
    return(exp(dnorm(x,mu,sd,log = TRUE) - pnorm(q = x,mean = mu,sd = sd,log.p = TRUE)))
  }
}

#############################################################################
### Derivada Sigma respecto alpha_jk
################################################################################

deriv.sigma <- function(Sigma,k,p){
  k     <- as.numeric(k)
  Sigma <- as.matrix(Sigma)
  p     <- dim(Sigma)[2]
  if (dim(Sigma)[1] != dim(Sigma)[2]) stop("Sigma is not square matrix\n")
  if (k > (p+1)*p/2)  stop("k out of bounds\n")
  e    <- rep(0,(p+1)*p/2)
  e[k] <- 1
  dev  <- matrix(0,ncol = p, nrow = p)
  dev[lower.tri(dev, diag = TRUE)] <- e
  dev  <- dev + t(dev)
  if(sum(diag(dev)) > 1) dev <- dev/2
  return(dev)
}

################################################################################
### label matrix de Info
################################################################################

nombres <- function(g,p, order = "TRUE", uni.Gama = "FALSE"){

  NAME <- piiN <- c()
  if (uni.Gama == "FALSE"){
    if (order == "TRUE") {
      for (i in 1:g) {
        SigmaN <- muN <- shapeN <- c()

        piiN <- paste("pii_", i, sep = "")

        for (k in 1:p) {
          muN <- c(muN, paste("mu", i, "_", k, sep = ""))
        }

        l <- m <- 1
        for (k in 1:((p + 1) * p/2)) {
          Vis <- FALSE
          SigmaN <- c(SigmaN, paste("Sigma", i, "_", l,
                                    m, sep = ""))
          if (((l * m - p * floor((l * m)/p)) == 0) &&
              (l != m)) {
            l <- l + 1
            m <- l
            Vis <- TRUE
          }
          if (!Vis)
            m <- m + 1
        } # fim do seg for k

        NAME <- c(NAME, muN, shapeN, SigmaN, piiN)

      } # fim do seg for i
      return(NAME[-length(NAME)])
    }

    if (order == "FALSE"){
      for (i in 1:g) {
        SigmaN <- muN <- shapeN <- c()

        for (k in 1:p) {
          muN <- c(muN, paste("mu", i, "_", k, sep = ""))
        }

        l <- m <- 1
        for (k in 1:((p + 1) * p/2)) {
          Vis <- FALSE
          SigmaN <- c(SigmaN, paste("Sigma", i, "_", l,
                                    m, sep = ""))
          if (((l * m - p * floor((l * m)/p)) == 0) &&
              (l != m)) {
            l <- l + 1
            m <- l
            Vis <- TRUE
          }
          if (!Vis)
            m <- m + 1
        } # fim do seg for k



        NAME <- c(NAME, muN, shapeN, SigmaN, piiN)

      } # fim do seg for i

      for (w in 1:(g - 1)) piiN <- c(piiN, paste("pii_", w, sep = ""))
      return(c(NAME, piiN))
    }


  } #fim uni.Gama = FALSE


  if (uni.Gama == "TRUE"){
    if (order == "TRUE") {

      for (i in 1:g) {
        SigmaN <- muN <- shapeN <- c()
        piiN <- paste("pii_", i, sep = "")

        for (k in 1:p)  muN <- c(muN, paste("mu", i, "_", k, sep = ""))
        NAME <- c(NAME, muN, piiN)
      } # fim do seg for i
      NAME <- NAME[-length(NAME)]

      l <- m <- 1
      for (k in 1:((p + 1) * p/2)) {
        Vis <- FALSE
        SigmaN <- c(SigmaN, paste("Sigma_", l, m, sep = ""))
        if (((l * m - p * floor((l * m)/p)) == 0) && (l != m)) {
          l <- l + 1
          m <- l
          Vis <- TRUE
        }
        if (!Vis)
          m <- m + 1
      } # fim do seg for k
      NAME <- c(NAME, SigmaN)
      return(NAME)
    }

    if (order == "FALSE") {

      for (i in 1:g) {
        SigmaN <- muN <- shapeN <- c()
        piiN <- paste("pii_", i, sep = "")

        for (k in 1:p)  muN <- c(muN, paste("mu", i, "_", k, sep = ""))
        NAME <- c(NAME, muN, piiN)
      } # fim do seg for i
      NAME <- NAME[-length(NAME)]

      l <- m <- 1
      for (k in 1:((p + 1) * p/2)) {
        Vis <- FALSE
        SigmaN <- c(SigmaN, paste("Sigma_", l, m, sep = ""))
        if (((l * m - p * floor((l * m)/p)) == 0) && (l != m)) {
          l <- l + 1
          m <- l
          Vis <- TRUE
        }
        if (!Vis)
          m <- m + 1
      } # fim do seg for k
      NAME <- c(NAME, SigmaN)
      return(NAME)
    }


  } #fim uni.Gama = TRUE

}

################################################################################
### Raiz Cuadrada de una Matriz
################################################################################

matrix.sqrt <- function(A) {
  sva <- svd(A)
  if (min(sva$d)>=0)
    Asqrt <- t(sva$v %*% (t(sva$u) * sqrt(sva$d)))
  else
    stop("Matrix square root is not defined")
  return(Asqrt)
}


################################################################################
### Densidades Normal e Misturas
################################################################################
#library(mvtnorm)

dmvNCens<-function(cc, LI, LS, y, mu, Sigma){

  m <- nrow(y)
  p <- ncol(y)

  ver <- matrix(0,m,1)
  mu  <- matrix(mu,p,1)

  for (j in 1:m ){
    cc1=matrix(cc[j,],p,1)
    y1=matrix(y[j,],p,1)
    LI1=matrix(LI[j,],p,1)
    LS1=matrix(LS[j,],p,1)

    if(sum(cc1)==0){
      ver[j]<-mvtnorm::dmvnorm(as.vector(y1), as.vector(mu), Sigma)
    }
    if(sum(cc1)>0){
      if(sum(cc1)==p){
        #ver[j]<-  sadmvn(LI1, LS1, as.vector(mu), Sigma) #pmnorm(as.vector(y1),as.vector(mu),Sigma)
        ver[j] = pmvnorm(lower = c(LI1),upper = c(LS1),mean = as.vector(mu),sigma = Sigma)
      }
      else{
        muc<- mu[cc1==1]+Sigma[cc1==1,cc1==0]%*%solve(Sigma[cc1==0,cc1==0])%*%(y1[cc1==0]-mu[cc1==0])
        Sc<- Sigma[cc1==1,cc1==1]-Sigma[cc1==1,cc1==0]%*%solve(Sigma[cc1==0,cc1==0])%*%Sigma[cc1==0,cc1==1]
        #ver[j]<-dmnorm(as.vector(y1[cc1==0]),as.vector(mu[cc1==0]),Sigma[cc1==0,cc1==0])*(pmvnorm(lower = c(LI1[cc1==1]), upper = c(LS1[cc1==1]), mean = as.vector(muc), sigma = Sc))#(pmnorm(as.vector(y1[cc1==1]),as.vector(muc),Sc))
        ver[j]<-dmvnorm(x = as.vector(y1[cc1==0]),mean = as.vector(mu[cc1==0]),sigma = Sigma[cc1==0,cc1==0,drop=FALSE])*(pmvnorm(lower = c(LI1[cc1==1]), upper = c(LS1[cc1==1]), mean = as.vector(muc), sigma = Sc))#(pmnorm(as.vector(y1[cc1==1]),as.vector(muc),Sc))
        }
    }
  }

  return(ver)

}


### Mixturas

d.mixedNCens <- function(cc, LI, LS, y, pi1, mu, Sigma){
  #y: Matrix de dados m x p
  #pi1: deve ser do tipo vetor de dimensao g
  #mu: deve ser do tipo list com g entradas. Cada entrada do list deve ser um vetor de dimensao p
  #Sigma: deve ser do tipo list com g entradas. Cada entrada do list deve ser uma matriz p x p

  g    <- length(pi1)
  dens <- 0
  for (j in 1:g) dens <- dens + pi1[j]*dmvNCens(cc, LI, LS, y , mu[[j]], Sigma[[j]])
  return(dens)
}


################################################################################
### Densidades SN Misturas
################################################################################

VerCensSN<-function(cc, LI, LS, y, mu, Sigma, lambda){
  m<-nrow(y)
  p<-ncol(y)
  ver<-matrix(0,m,1)
  mu<-matrix(mu,p,1)

  for (j in 1:m){
    cc1=matrix(cc[j,],p,1)
    LI1=matrix(LI[j,],p,1)
    LS1=matrix(LS[j,],p,1)

    y1=matrix(y[j,],p,1)

    if(sum(cc1)==0){
      ver[j] <- MomTrunc::dmvESN(x = as.vector(y1),mu = as.vector(mu),Sigma = Sigma,lambda = lambda,tau = 0)
      }
    if(sum(cc1)>0){
      if(sum(cc1)==p){
        ver[j]<- MomTrunc::pmvESN(lower = c(LI1), upper = c(LS1), mu = as.vector(mu),Sigma = Sigma,lambda = lambda,tau = 0)
      }
      else{
        muc<- mu[cc1==1] + Sigma[cc1==1,cc1==0]%*%solve(Sigma[cc1==0,cc1==0])%*%(y1[cc1==0]-mu[cc1==0])
        Sc<- Sigma[cc1==1,cc1==1]-Sigma[cc1==1,cc1==0]%*%solve(Sigma[cc1==0,cc1==0])%*%Sigma[cc1==0,cc1==1]
        varphi = solve(sqrtm(Sigma))%*%lambda
        tilvarphi.o = varphi[cc1==0] + solve(Sigma[cc1==0,cc1==0])%*%Sigma[cc1==0,cc1==1]%*%varphi[cc1==1]
        c.oc = as.numeric(1/sqrt(1 + t(varphi[cc1==1])%*%Sc%*%varphi[cc1==1]))
        tau.co = as.numeric(t(tilvarphi.o)%*%(y1[cc1==0]-mu[cc1==0]))

        ver[j] <- MomTrunc::dmvESN(x = as.vector(y1[cc1==0]),mu = as.vector(mu[cc1==0]),Sigma = Sigma[cc1==0,cc1==0],lambda = c.oc*sqrtm(Sigma[cc1==0,cc1==0])%*%tilvarphi.o,tau = 0)*(MomTrunc::pmvESN(lower = c(LI1[cc1==1]), upper = c(LS1[cc1==1]), mu = as.vector(muc),Sigma = Sc,lambda = sqrtm(Sc)%*%varphi[cc1==1],tau = tau.co))#(pmnorm(as.vector(y1[cc1==1]),as.vector(muc),Sc))
       }
    }
  }
  return(ver)
}



d.mixedSNCens <- function(cc, LI, LS, y, pi1, mu, Sigma, shape){
  #y: Matrix de dados m x p
  #pi1: deve ser do tipo vetor de dimensao g
  #mu: deve ser do tipo list com g entradas. Cada entrada do list deve ser um vetor de dimensao p
  #Sigma: deve ser do tipo list com g entradas. Cada entrada do list deve ser uma matriz p x p

  g    <- length(pi1)
  dens <- 0
  for (j in 1:g) dens <- dens + pi1[j]*VerCensSN(cc, LI, LS, y, mu[[j]], Sigma[[j]], shape[[j]])
  return(dens)
}


################################################################################
### Densidades T Student e Misturas
################################################################################

dmvTCens<-function(cc, LI, LS, y, mu, Sigma, nu){

  # GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)

  m <- nrow(y)
  p <- ncol(y)

  ver<-matrix(0,m,1)

  mu<-matrix(mu,p,1)

  #    nu=ceiling(nu)             ###???

  for(j in 1:m ){
    cc1=matrix(cc[j,],p,1)
    y1=matrix(y[j,],p,1)
    LI1=matrix(LI[j,],p,1)
    LS1=matrix(LS[j,],p,1)

    if(sum(cc1)==0){
      #ver[j]<- dmt(as.vector(y1),as.vector(mu),Sigma,df=nu)
      ver[j]<- dmvt(x = as.vector(y1),delta = as.vector(mu),sigma=Sigma,df=nu,log = FALSE)
    }
    if(sum(cc1)>0){

      if(sum(cc1)==p){

        #ver[j]<-sadmvt(df=nu, as.vector(LI1),as.vector(LS1),as.vector(mu),Sigma)
        #ver[j]<-sadmvt(df=nu, as.vector(LI1),as.vector(LS1),as.vector(mu),Sigma)

        ver[j]<-prob_opt(lower = as.vector(LI1),upper = as.vector(LS1),mean = as.vector(mu),sigma = Sigma,nu = nu)
        #ver[j]<-pmvt(df=nu, lower = as.vector(LI1),upper = as.vector(LS1),delta = as.vector(mu),sigma = Sigma,algorithm = GenzBretz(maxpts = 2000*p, abseps = 1e-06, releps = 0))

        }

      else {
        nu1<-(nu+length(cc1[cc1==0]))
        muc<-mu[cc1==1]+Sigma[cc1==1,cc1==0]%*%solve(Sigma[cc1==0,cc1==0])%*%(y1[cc1==0]-mu[cc1==0])
        Sc <-Sigma[cc1==1,cc1==1]-Sigma[cc1==1,cc1==0]%*%solve(Sigma[cc1==0,cc1==0])%*%Sigma[cc1==0,cc1==1]
        Qy1<-t(y1[cc1==0]-mu[cc1==0])%*%solve(Sigma[cc1==0,cc1==0])%*%(y1[cc1==0]-mu[cc1==0])
        auxcte<-as.numeric((nu+Qy1)/(nu+length(cc1[cc1==0])))
        Sc22<-auxcte*Sc
        muUi<-muc
        SigmaUi<-Sc22
        SigmaUi<-(SigmaUi+t(SigmaUi))/2
        Sigma[cc1==0,cc1==0]<-(Sigma[cc1==0,cc1==0]+t(Sigma[cc1==0,cc1==0]))/2
        #ver[j,]<-dmt(as.vector(y1[cc1==0]),as.vector(mu[cc1==0]),as.matrix(Sigma[cc1==0,cc1==0]),df=nu)*
         # sadmvt(df=nu1,as.vector(LI1[cc1==1]),as.vector(LS1[cc1==1]),as.vector(muUi),SigmaUi)
        #ver[j,]<-dmvt(x = as.vector(y1[cc1==0]),delta = c(mu[cc1==0]),sigma = as.matrix(Sigma[cc1==0,cc1==0]),df=nu,log = FALSE)*
         # pmvt(df=nu1,lower = as.vector(LI1[cc1==1]),upper = as.vector(LS1[cc1==1]),delta = as.vector(muUi),sigma = SigmaUi,algorithm = GenzBretz(maxpts = 2000*p, abseps = 1e-06, releps = 0))
        ver[j,]<-dmvt(x = as.vector(y1[cc1==0]),delta = c(mu[cc1==0]),sigma = as.matrix(Sigma[cc1==0,cc1==0]),df=nu,log = FALSE)*
          prob_opt(lower = as.vector(LI1[cc1==1]),upper = as.vector(LS1[cc1==1]),mean = as.vector(muUi),sigma = SigmaUi,nu=nu1)
        }

    }
  }

  return(ver)
}


### Mixturas

d.mixedTCens <- function(cc, LI,LS, y, pi1, mu, Sigma,nu){
  #y: Matrix de dados m x p
  #pi1: deve ser do tipo vetor de dimensao g
  #mu: deve ser do tipo list com g entradas. Cada entrada do list deve ser um vetor de dimensao p
  #Sigma: deve ser do tipo list com g entradas. Cada entrada do list deve ser uma matriz p x p

  g <- length(pi1)
  dens <- 0
  for (j in 1:g) dens <- dens + pi1[j]*dmvTCens(cc, LI, LS, y , mu[[j]], Sigma[[j]], nu[j])
  return(dens)
}




################################################################################
##
################################################################################

skewness = function (x, na.rm = FALSE){
    if (is.matrix(x))
      apply(x, 2, skewness, na.rm = na.rm)
    else if (is.vector(x)) {
      if (na.rm) x <- x[!is.na(x)]
      n <- length(x)
      (sum((x-mean(x))^3)/n)/(sum((x-mean(x))^2)/n)^(3/2)
    }
    else if (is.data.frame(x))
      sapply(x, skewness, na.rm = na.rm)
    else skewness(as.vector(x), na.rm = na.rm)
  }
