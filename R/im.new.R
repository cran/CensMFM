## Information matrix function ##

imm.msnc.new = function(y,model,cc,LI,LS,family,uni.Gama){

  #model = test_fit

  ## Some variables
  n <- nrow(y)
  p <- ncol(y)
  g <- length(model$pii)

  if(family=="SN"){

    if(uni.Gama) stop("Sorry. The information matrix cannot be calculated when the uni.Gama was used!")

    ro <- g*sum(c(p,p,1,p*(1 + p)/2)) # mu, lambda, pi, sigmas

  ## Creating the information matrix
  MI <- matrix(0,ro - 1,ro - 1)
  esc <- matrix(nrow = n, ncol = ro - 1)

  ## Calculation

  for(i in 1:g) model$Sigma[[i]] <- model$Sigma[[i]]%*%model$Sigma[[i]]

  # Essential variable
  delta <- Delta <- Gamma <- list()
  mu <- model$mu
  Sigma <- model$Sigma
  shape <- model$shape
  pii <- model$pii
  Zij <- model$Zij

  for(i in 1:n){
    si <- matrix(nrow = ro - 1 ,ncol = 1)
    cc1 <- matrix(cc[i,],p,1)
    LI1<-matrix(LI[i,],p,1)
    LS1<-matrix(LS[i,],p,1)
    y1  <- matrix(y[i,],p,1)

    for(j in 1:g){
      si.mu <- si.lambda <- rep(0,dim(y)[2])
      si.sigma <- rep(0,p*(p + 1)/2)
      si.pii <-  rep(0,g-1)

      SSj   <- sqrtm(Sigma[[j]])
      iSSj  <- solve(SSj)
      delta[[j]] <- shape[[j]]/as.numeric(sqrt(1 + t(shape[[j]])%*%shape[[j]]))
      Delta[[j]] <- SSj%*%delta[[j]]
      Gamma[[j]] <- Sigma[[j]] - Delta[[j]]%*%t(Delta[[j]])
      M2     <- 1/(1+sum(shape[[j]]^2))
      varphi <- iSSj%*%shape[[j]]

      if(sum(cc1)==0){
        #aux values
        aux1 = t(varphi)%*%(y1-mu[[j]])
        aux2 = M2*t(Delta[[j]])%*%solve(Gamma[[j]])%*%(y1-mu[[j]])
        #output
        T1.y = as.numeric(aux2 + sqrt(M2)*invmills(aux1))
        T2.y = as.numeric(aux2^2 + sqrt(M2)*aux2*invmills(aux1) + M2)

        y.hat  = matrix(y1,p,1)
        yy.hat =  y1%*%t(y1)
        t1.hat = T1.y
        t2.hat = T2.y
        ty.hat = T1.y*y.hat

      }else{
        if(sum(cc1)==p){
          eta    = sqrt(2/pi)/sqrt(1+sum(shape[[j]]^2))

          aux3   = meanvarTMD(lower = c(LI1),upper = c(LS1),mu = mu[[j]],Sigma = Sigma[[j]],lambda = shape[[j]],dist = "SN")
          y.hat  = aux3$mean
          yy.hat =  aux3$EYY

          #print(y.hat)

          w0.hat   = onlymeanTMD(lower = c(LI1),upper = c(LS1),mu = c(mu[[j]]),Sigma = Gamma[[j]],dist = "normal")

          # Ltemp  = as.numeric(pmvnorm(lower = c(LI1),upper = c(LS1),mean = c(mu[[j]]),sigma = Gamma[[j]],algorithm = GenzBretz(maxpts = 25000))[1])
          # LLtemp = pmvESN(lower = c(LI1),upper = c(LS1),mu = as.vector(mu[[j]]),Sigma = Sigma[[j]],lambda = shape[[j]],tau = 0,algorithm = GenzBretz(maxpts = 25000))
          # gamma  = eta*Ltemp/LLtemp

          Ltemp  = prob_opt(lower = c(LI1),upper = c(LS1),mean = c(mu[[j]]),sigma = Gamma[[j]],uselog2 = TRUE)
          LLtemp = MomTrunc::pmvSN(lower = c(LI1),upper = c(LS1),mu = as.vector(mu[[j]]),Sigma = Sigma[[j]],lambda = shape[[j]],log2 = TRUE)
          gamma  = eta*2^(Ltemp - LLtemp)

          val       = c(t(shape[[j]])%*%iSSj%*%(y.hat - mu[[j]]))
          gamma.ap  = invmills(val)
          boolgamma = is.na(gamma) |  gamma <= 0 | is.infinite(gamma) | abs(gamma) > 100*abs(gamma.ap)
          gamma    = ifelse(boolgamma,gamma.ap,gamma)

          #output
          t1.hat = as.numeric(M2*t(Delta[[j]])%*%solve(Gamma[[j]])%*%(y.hat-mu[[j]]) + sqrt(M2)*gamma)
          t2.hat = as.numeric(M2^2*t(Delta[[j]])%*%solve(Gamma[[j]])%*%(yy.hat - 2*y.hat%*%t(mu[[j]]) + mu[[j]]%*%t(mu[[j]]))%*%solve(Gamma[[j]])%*%Delta[[j]] + M2 +
                                gamma*M2^(3/2)*t(Delta[[j]])%*%solve(Gamma[[j]])%*%(w0.hat - mu[[j]]))
          ty.hat = as.numeric(M2*(yy.hat - y.hat%*%t(mu[[j]]))%*%solve(Gamma[[j]])%*%Delta[[j]] + gamma*sqrt(M2)*w0.hat)
        }else{
          ###################################################################################################################
          #new conditional parameters
          muc         = mu[[j]][cc1==1] + Sigma[[j]][cc1==1,cc1==0]%*%solve(Sigma[[j]][cc1==0,cc1==0])%*%(y1[cc1==0]-mu[[j]][cc1==0])
          Sc          = Sigma[[j]][cc1==1,cc1==1]-Sigma[[j]][cc1==1,cc1==0]%*%solve(Sigma[[j]][cc1==0,cc1==0])%*%Sigma[[j]][cc1==0,cc1==1]
          Sc          = (Sc+t(Sc))/2
          tilvarphi.o = varphi[cc1==0] + solve(Sigma[[j]][cc1==0,cc1==0])%*%Sigma[[j]][cc1==0,cc1==1]%*%varphi[cc1==1]
          tau.co      = as.numeric(t(tilvarphi.o)%*%(y1[cc1==0]-mu[[j]][cc1==0]))
          lambda.co   = sqrtm(Sc)%*%varphi[cc1==1]
          tautil.cc   = tau.co/sqrt(1+sum(lambda.co^2))

          #auxiliar conditional parameters
          SS.cc     = sqrtm(Sc)
          iSS.cc    = solve(SS.cc)
          varphi.cc = iSS.cc%*%lambda.co
          Delta.cc  = SS.cc%*%lambda.co/sqrt(1+sum(lambda.co^2))
          Gamma.cc  = Sc - Delta.cc%*%t(Delta.cc)
          mub.cc    = tautil.cc*Delta.cc
          eta.cc    = invmills(tau.co,0,sqrt(1+sum(lambda.co^2)))

          #first and second conditional moments

          aux3   = meanvarTMD(lower = c(LI1[cc1==1]),upper = c(LS1[cc1==1]),mu = as.vector(muc),Sigma = Sc,lambda = lambda.co,tau = tau.co,dist = "ESN")
          w1.hat = aux3$mean
          w2.hat =  aux3$EYY

          #ratio

          # Ltemp  = as.numeric(pmvnorm(c(LI1[cc1==1]),c(LS1[cc1==1]),mean = c(muc - mub.cc),sigma = Gamma.cc,algorithm = GenzBretz(maxpts = 25000))[1])
          # LLtemp = pmvESN(lower = LI1[cc1==1],upper = LS1[cc1==1],mu = as.vector(muc),Sigma = Sc,lambda = lambda.co,tau = tau.co,algorithm = GenzBretz(maxpts = 25000))
          # gamma.cc = eta.cc*Ltemp/LLtemp

          Ltemp  = prob_opt(lower = c(LI1[cc1==1]),upper = c(LS1[cc1==1]),c(muc - mub.cc),sigma = Gamma.cc,uselog2 = TRUE)
          LLtemp = MomTrunc::pmvESN(lower = LI1[cc1==1],upper = LS1[cc1==1],mu = as.vector(muc),Sigma = Sc,lambda = lambda.co,tau = tau.co,log2 = TRUE)
          gamma.cc  = eta.cc*2^(Ltemp - LLtemp)

          # ratio approximation

          val         = c(tau.co + t(lambda.co)%*%iSS.cc%*%(w1.hat - muc))
          gamma.cc.ap = invmills(val)
          boolgamma   = is.na(gamma.cc) |  gamma.cc <= 0 | is.infinite(gamma.cc) | abs(gamma.cc) > 100*abs(gamma.cc.ap)
          gamma.cc    = ifelse(boolgamma,gamma.cc.ap,gamma.cc)

          w0.hat = onlymeanTMD(lower = LI1[cc1==1],upper = LS1[cc1==1],mu = c(muc - mub.cc),Sigma = Gamma.cc,dist = "normal")

          y0.hat = matrix(y1,p,1)
          y0.hat[cc1==1] = w0.hat
          aux4 = gamma.cc*M2^(3/2)*t(Delta[[j]])%*%solve(Gamma[[j]])%*%(y0.hat - mu[[j]])
          aux5 = gamma.cc*sqrt(M2)*y0.hat

          y.hat  = matrix(y1,p,1)
          y.hat[cc1==1] = w1.hat
          yy.hat =  y1%*%t(y1)
          yy.hat[cc1==0,cc1==1] = y1[cc1==0]%*%t(w1.hat)
          yy.hat[cc1==1,cc1==0] = w1.hat%*%t(y1[cc1==0])
          yy.hat[cc1==1,cc1==1] = w2.hat

          t1.hat = as.numeric(M2*t(Delta[[j]])%*%solve(Gamma[[j]])%*%(y.hat-mu[[j]]) + sqrt(M2)*gamma.cc)
          t2.hat = as.numeric(M2^2*t(Delta[[j]])%*%solve(Gamma[[j]])%*%(yy.hat - 2*y.hat%*%t(mu[[j]]) + mu[[j]]%*%t(mu[[j]]))%*%solve(Gamma[[j]])%*%Delta[[j]] + M2 + aux4)
          ty.hat = as.numeric(M2*(yy.hat - y.hat%*%t(mu[[j]]))%*%solve(Gamma[[j]])%*%Delta[[j]] + aux5)
        }
      }# end else sum(cc1)==0

      # Information matrix

      Gamma.inv = solve(Gamma[[j]])

      # beta
      si.mu <- Zij[i,j]*(Gamma.inv%*%(y.hat - mu[[j]] - Delta[[j]]*t1.hat))

      #Sigma
      D_F = Der.F(SSj)
      md2<-dim(SSj)[1]
      kont <- 0
      for(i1 in 1:md2){
        for(i2 in 1:(md2+1-i1)){
          kont <- kont+1
          der.F <- D_F[[i1]][[i2]]
          dev.gamma <- (SSj%*%der.F + der.F%*%SSj - der.F%*%delta[[j]]%*%t(delta[[j]])%*%SSj - SSj%*%delta[[j]]%*%t(delta[[j]])%*%der.F)
          sigma.aux <- Gamma.inv%*%dev.gamma%*%Gamma.inv
          si.sigma[kont] <- -0.5*Zij[i,j]*(sum(diag(Gamma.inv%*%dev.gamma)) - sum(diag(yy.hat%*%sigma.aux)) + t(y.hat)%*%(sigma.aux)%*%mu[[j]]
                                           - t(ty.hat)%*%(Gamma.inv%*%der.F%*%delta[[j]] - sigma.aux%*%Delta[[j]])  + t(mu[[j]])%*%sigma.aux%*%y.hat - t(mu[[j]])%*%sigma.aux%*%mu[[j]]
                                           + t1.hat*t(mu[[j]])%*%(Gamma.inv%*%der.F%*%delta[[j]] - sigma.aux%*%Delta[[j]]) - (t(delta[[j]])%*%der.F%*%Gamma.inv - t(Delta[[j]])%*%sigma.aux)%*%ty.hat
                                           + (t(delta[[j]])%*%der.F%*%Gamma.inv - t(Delta[[j]])%*%sigma.aux)%*%mu[[j]]*t1.hat + t2.hat*(t(delta[[j]])%*%der.F%*%Gamma.inv%*%Delta[[j]]
                                                                                                                                        - t(Delta[[j]])%*%sigma.aux%*%Delta[[j]] + t(Delta[[j]])%*%Gamma.inv%*%der.F%*%delta[[j]]))
        }
      }


      # lambda
      lamb.aux = as.vector(1 + t(shape[[j]])%*%shape[[j]])
      Lamb.M = shape[[j]]%*%t(shape[[j]])
      D_R = Der.R(shape[[j]]) # Derivada de lambda%*%t(lambda) em relação a lambda_jr
      md2 <- length(shape[[j]])
      #kont <- 0
      for(i1 in 1:md2){
        #i1 = 1
        #for(i2 in 1:(md2+1-i1)){
        #kont <- kont+1
        der.R <- D_R[[i1]]
        lamb1 <- SSj%*%((der.R*lamb.aux - 2*Lamb.M*shape[[j]][i1])/(lamb.aux)^2)%*%SSj
        lamb2 <- Gamma.inv%*%lamb1%*%Gamma.inv
        lamb3 <- SSj%*%((VecLam(shape[[j]],i1)*lamb.aux - shape[[j]][i1]*shape[[j]])/(lamb.aux)^(3/2)  )
        lamb4 <- (( t(VecLam(shape[[j]],i1))*lamb.aux - shape[[j]][i1]*t(shape[[j]]))/(lamb.aux)^(3/2))%*%SSj
        si.lambda[i1] <- 0.5*Zij[i,j]*sum(diag(Gamma.inv%*%lamb1)) - 0.5*Zij[i,j]*(sum(diag(yy.hat%*%lamb2)) - t(y.hat)%*%lamb2%*%mu[[j]]
                                                                                   - t(ty.hat)%*%(lamb2%*%Delta[[j]]  + Gamma.inv%*%lamb3) - t(mu[[j]])%*%lamb2%*%y.hat +  t(mu[[j]])%*%lamb2%*%mu[[j]]
                                                                                   + t1.hat*t(mu[[j]])%*%(lamb2%*%Delta[[j]] + Gamma.inv%*%lamb3) - (lamb4%*%Gamma.inv + t(Delta[[j]])%*%lamb2)%*%ty.hat
                                                                                   + t1.hat*(lamb4%*%Gamma.inv + t(Delta[[j]])%*%lamb2)%*%mu[[j]] + t2.hat*(lamb4%*%Gamma.inv%*%Delta[[j]] + t(Delta[[j]])%*%lamb2%*%Delta[[j]] + t(Delta[[j]])%*%Gamma.inv%*%lamb3) )
        #}
      }

      si[seq((j - 1)*(ro/g - 1) + 1,j*(ro/g) - j),] <- c(si.mu,si.sigma,si.lambda)
      #S <- c(S,si.mu,si.sigma,si.lambda)

    }# end for j

    if(g > 1){
      for(j in 1:(g-1)){ si.pii[j] <- (Zij[i,j]/as.numeric(pii[j])) - (Zij[i,g]/as.numeric(pii[g]))
      si[ro - (g - j),] <- si.pii[j]
      }
    }

    #esc[i,] <- si
    MI <- MI + si%*%t(si)
  }# end for i
  NAME <- piiN <- c()
  for(i in 1:g){
    SigmaN <- muN <- shapeN <- c()
    for (k in 1:p){
      muN <- c(muN,paste("mu",i,"_",k,sep=""))
      shapeN <- c(shapeN,paste("shape",i,"_",k,sep=""))
    }
    l <- m <- 1
    for(k in 1:((p+1)*p/2)) {
      Vis <- FALSE
      SigmaN <- c(SigmaN,paste("alpha",i,"_",l,m,sep=""))
      if(((l*m - p*floor((l*m)/p)) == 0) && (l != m)) {
        l <- l+1
        m <- l
        Vis <- TRUE
      }
      if(!Vis) m <- m+1
    }
    NAME <- c(NAME,muN,shapeN,SigmaN)
  }

  if(g >1){
    for(i in 1:(g-1)) piiN <- c(piiN, paste("pii_",i,sep=""))
    NAME <- c(NAME,piiN)
  }
  dimnames(MI)[[1]] <- NAME
  dimnames(MI)[[2]] <- NAME

  MI <- (MI + t(MI))/2
  IM = sqrt(diag(solve(MI)))
  #return(list(MI=IM,esc = esc))
  return(list(MI=IM))

  }# end family SN

#-- Family Normal
if(family=="Normal"){

  if(uni.Gama == FALSE) ro <- g*(((p+1)*p/2)+p+((g-1)/g))
  if(uni.Gama == TRUE)  ro <- g*(((p+1)*p/(2*g))+p+((g-1)/g))

  ## Creating the information matrix
  MI      <- matrix(0, ncol = ro, nrow = ro)
  alpha.g <- matrix(0, ncol = (p+1)*p/2, nrow = n)
  si.pi   <- si.mu <- alpha <- alpha.2 <- si <- si.g <- sipi <- NULL

  ## Calculation
  # Essential variable
  mu <- model$mu
  Sigma <- model$Sigma
  nu <- model$nu
  pii <- model$pii
  Zij <- model$Zij

  tuyi  <- matrix(0,n,p)
  tui   <- matrix(0,n,1)
  tuyyi <- list()

  for (j in 1:g){

    mus    <- mu[[j]]
    Sigmas <- Sigma[[j]]

    for (i in 1:n ){

      cc1 <- matrix(cc[i,],p,1)
      LI1<-matrix(LI[i,],p,1)
      LS1<-matrix(LS[i,],p,1)
      y1  <- matrix(y[i,],p,1)

      if(sum(cc1)==0){

        uy  <- matrix(y1,p,1)
        uyy <- y1%*%t(y1)

      }

      if(sum(cc1)>=1){

        if(sum(cc1)==p){
          muc<-mus
          Sc<-Sigmas
          aux <- meanvarTMD(lower = c(LI1), upper = c(LS1), mu = as.vector(muc), Sigma = Sc, dist = "normal")
          uy <- aux$mean
          uyy <- aux$EYY
          #aux<-  mtmvnorm(as.vector(muc), Sc, c(LI1), c(LS1))#MomemNT(muc,Sc,y1)
          #uy<-aux$tmean
          #uyy<- aux$tvar+(uy)%*%t(uy)#aux$Eyy
        }

        else {
          muc = mus[cc1==1]+Sigmas[cc1==1,cc1==0]%*%solve(Sigmas[cc1==0,cc1==0])%*%(y1[cc1==0]-mus[cc1==0])
          Sc <-Sigmas[cc1==1,cc1==1]-Sigmas[cc1==1,cc1==0]%*%solve(Sigmas[cc1==0,cc1==0])%*%Sigmas[cc1==0,cc1==1]
          Sc<-(Sc+t(Sc))/2
          #aux <- mtmvnorm(as.vector(muc), Sc, LI1[cc1==1], LS1[cc1==1])#MomemNT(muc,Sc,y1[cc1==1])
          aux <- meanvarTMD(lower = c(LI1[cc1==1]), upper = c(LS1[cc1==1]), mu = as.vector(muc), Sigma = Sc, dist = "normal")
          uy <- matrix(y1,p,1)
          #uy[cc1==1]<- aux$tmean
          uy[cc1==1]<- aux$mean
          uyy<-matrix(0,p,p)
          #uyy[cc1==1,cc1==1]<-aux$tvar
          #uyy<- uyy+uy%*%t(uy)
          uyy[cc1==1,cc1==1] <- aux$EYY
          }

      }

      tuyi[i,]   <- t((Zij[i,j])*uy)
      tui[i,]    <- Zij[i,j]
      tuyyi[[i]] <- (Zij[i,j])*(uyy)

      if(uni.Gama == FALSE){
        si.mu <- solve(Sigma[[j]]) %*% (tuyi[i,] - (tui[i]*mu[[j]]))
        si.pi <- (Zij[i,j]/as.numeric(pii[j])) - (Zij[i,g]/as.numeric(pii[g]))
        for(t in 1:((p+1)*p/2)) {
          alpha[t] <- -1/2*sum(diag( ( (Zij[i,j]*solve(Sigma[[j]])) %*% deriv.sigma(Sigma[[j]],t,p) ) -
                                       ((tuyyi[[i]] - mu[[j]]%*%t(tuyi[i,]) - tuyi[i,]%*%t(mu[[j]]) + ((tui[i]*mu[[j]])%*%t(mu[[j]]))) %*%
                                          solve(Sigma[[j]]) %*% deriv.sigma(Sigma[[j]],t,p) %*% solve(Sigma[[j]])) ))

        }
        si <- rbind(si,c(si.mu, alpha, si.pi))
      } else {
        si.mu <- solve(Sigma[[j]]) %*% (tuyi[i,] - (tui[i]*mu[[j]]))
        si.pi <- (Zij[i,j]/as.numeric(pii[j])) - (Zij[i,g]/as.numeric(pii[g]))
        for(t in 1:((p+1)*p/2)) {
          alpha[t] <- -1/2*sum(diag( ( (Zij[i,j]*solve(Sigma[[j]])) %*% deriv.sigma(Sigma[[j]],t,p) ) -
                                       ((tuyyi[[i]] - mu[[j]]%*%t(tuyi[i,]) - tuyi[i,]%*%t(mu[[j]]) + ((tui[i]*mu[[j]])%*%t(mu[[j]]))) %*%
                                          solve(Sigma[[j]]) %*% deriv.sigma(Sigma[[j]],t,p) %*% solve(Sigma[[j]])) ))
        }
        alpha.2 <- rbind(alpha.2,alpha)
        si <- rbind(si,c(si.mu,si.pi))
      }

    } ## end for de i

    if(uni.Gama == TRUE) alpha.g <- alpha.g + alpha.2
    si.g    <- cbind(si.g, si)
    si      <- NULL
    alpha.2 <- NULL

  }   ##FIM DO FOR  de g

  if(uni.Gama == FALSE){
    si.g <- si.g[,-dim(si.g)[2]]
    if(g==3) si.g <- cbind(si.g[,-c((p+(p+1)*p/2 + 1),(p+(p+1)*p/2 + 1)*2)],si.g[,c((p+(p+1)*p/2 + 1),(p+(p+1)*p/2 + 1)*2)])
    if(g==2) si.g <- cbind(si.g[,-(p+(p+1)*p/2 + 1)],si.g[,(p+(p+1)*p/2 + 1)])
    for(v in 1:n) MI = MI + si.g[v,]%*%t(si.g[v,])
    colnames(MI) <- rownames(MI) <- nombres(g, p, order = "FALSE")
    #round(MI); print(paste("det(MI) =",det(MI),sep=" "))
    IM <- sqrt(diag(solve(MI)))
  } else {
    si.g <- si.g[,-dim(si.g)[2]]
    si.g <- cbind(si.g,alpha.g)
    for(v in 1:n) MI = MI + si.g[v,]%*%t(si.g[v,])
    colnames(MI) <- rownames(MI) <- nombres(g, p, uni.Gama = "TRUE")
    #round(MI); print(paste("det(MI) =",det(MI),sep=" "))
    IM <- sqrt(diag(solve(MI)))
  }
  return(list(MI=IM))
  } #end family Normal

  #-- Family t
  if(family=="t"){

    if(uni.Gama == FALSE) ro <- g*(((p+1)*p/2)+p+((g-1)/g))
    if(uni.Gama == TRUE)  ro <- g*(((p+1)*p/(2*g))+p+((g-1)/g))

    ## Creating the information matrix
    MI      <- matrix(0, ncol = ro, nrow = ro)
    alpha.g <- matrix(0, ncol = (p+1)*p/2, nrow = n)
    si.pi   <- si.mu <- alpha <- alpha.2 <- si <- si.g <- sipi <- NULL

    ## Calculation
    # Essential variable
    mu <- model$mu
    Sigma <- model$Sigma
    nu <- model$nu
    pii <- model$pii
    Zij <- model$Zij

    tuyi  <- matrix(0,n,p)
    tui   <- matrix(0,n,1)
    tuyyi <- list()

    for (j in 1:g){

      mus    <- mu[[j]]
      Sigmas <- Sigma[[j]]

      for (i in 1:n ){

        cc1 <- matrix(cc[i,],p,1)
        y1  <- matrix(y[i,],p,1)
        LI1<-matrix(LI[i,],p,1)
        LS1<-matrix(LS[i,],p,1)

        dm  <- t(y1-mus)%*%solve(Sigmas)%*%(y1-mus)
        cdm <- as.numeric((nu[j]+p)/(nu[j]+dm))

        if(sum(cc1)==0){

          tuy  <- (matrix(y1,p,1))*cdm
          tuyy <- (y1%*%t(y1))*cdm
          tu   <- cdm

          #ver[j,]<- dmvt(as.vector(y1),as.vector(mu),as.matrix(Sigma),df=nu,log=FALSE)
        }

        if(sum(cc1)>0){

          if(sum(cc1) == p){

            muUi     <- mus
            SigmaUi  <- Sigmas
            SigmaUiA <- SigmaUi*nu[j]/(nu[j]+2)
            auxupper <- y1-muUi
            auxU1    <- pmvt(df=nu[j]+2, lower = as.vector(LI1), upper = as.vector(LS1), delta = as.vector(muUi),sigma = SigmaUiA,type="shitfed")
            auxU2    <- pmvt(df=nu[j], lower = as.vector(LI1), upper = as.vector(LS1), delta = as.vector(muUi), sigma= SigmaUi,type="shitfed")
            #auxU1    <- sadmvt(df=nu+2, as.vector(LI1), as.vector(LS1), as.vector(muUi),SigmaUiA)
            #auxU2    <- sadmvt(df=nu, as.vector(LI1), as.vector(LS1), as.vector(muUi), SigmaUi)
            #MoMT     <- Mtmvt(muUi,SigmaUiA,nu+2,rep(-Inf,p),y1)
            MoMT     <- meanvarTMD(as.vector(LI1),as.vector(LS1),as.vector(muUi),SigmaUiA,dist = "t",nu = nu[j]+2)
            U0       <- as.numeric(auxU1/auxU2)
            U1       <- auxU1/auxU2*MoMT$mean
            U2       <- auxU1/auxU2*MoMT$EYY

            tuy  <- U1
            tuyy <- U2
            tu   <- U0

            #   ver[j,]<-pmvt(upper = c(auxupper), sigma = SigmaUi, df = nu,algorithm = GB) [1]


          }

          else {

            PsiA <- Sigmas*nu[j]/(nu[j]+2)
            nu1  <- (nu[j]+length(cc1[cc1==0]))

            muc  <- mus[cc1==1]+Sigmas[cc1==1,cc1==0]%*%solve(Sigmas[cc1==0,cc1==0])%*%(y1[cc1==0]-mus[cc1==0])
            Sc   <- Sigmas[cc1==1,cc1==1]-Sigmas[cc1==1,cc1==0]%*%solve(Sigmas[cc1==0,cc1==0])%*%Sigmas[cc1==0,cc1==1]
            Sc   <- (Sc+t(Sc))/2
            ScA  <- nu[j]/(nu[j]+2)*Sc

            Qy1  <- t(y1[cc1==0]-mus[cc1==0])%*%solve(Sigmas[cc1==0,cc1==0])%*%(y1[cc1==0]-mus[cc1==0])
            Qy2  <- t(y1[cc1==0]-mus[cc1==0])%*%solve(PsiA[cc1==0,cc1==0])%*%(y1[cc1==0]-mus[cc1==0])

            auxcte  <- as.numeric((nu[j]+Qy1)/(nu[j]+length(cc1[cc1==0])))
            auxcte1 <- as.numeric((nu[j]+2+Qy2)/(nu[j]+2+length(cc1[cc1==0])))

            Sc22 <- auxcte*Sc

            muUi    <- muc
            SigmaUi <- Sc22

            SigmaUiA <- auxcte1*ScA
            auxupper <- y1[cc1==1]-muUi


            auxU1    <- pmvt(df=nu1+2, lower = as.vector(LI1[cc1==1]), upper = as.vector(LS1[cc1==1]), delta = as.vector(muUi),sigma = SigmaUiA,type="shitfed")
            auxU2    <- pmvt(df=nu1, lower = as.vector(LI1[cc1==1]), upper = as.vector(LS1[cc1==1]), delta = as.vector(muUi), sigma= SigmaUi,type="shitfed")
            #auxU1<-sadmvt(df=nu1+2, as.vector(LI1[cc1==1]), as.vector(LS1[cc1==1]), as.vector(muUi), SigmaUiA)
            #auxU2 <- sadmvt(df=nu1, as.vector(LI1[cc1==1]), as.vector(LS1[cc1==1]), as.vector(muUi), SigmaUi)
            #MoMT <- Mtmvt(muUi,SigmaUiA,nu1+2,rep(-Inf,length(cc1[cc1==1])),y1[cc1==1])
            MoMT <- meanvarTMD(as.vector(LI1[cc1==1]),as.vector(LS1[cc1==1]),muUi,SigmaUiA,dist = "t",nu = nu1+2)

            U0 <- as.numeric(auxU1/auxU2)/auxcte
            U1 <- (U0)*(MoMT$mean)
            U2 <- (U0)*(MoMT$EYY)

            Auxtuy <- (matrix(y1,p,1))

            tuy <- Auxtuy*U0
            tuy[cc1==1]<- U1

            tuyy <- (Auxtuy%*%t(Auxtuy))

            AAx <- tuyy[cc1==0,cc1==0]*U0
            ABx <- Auxtuy[cc1==0]%*%t(U1)
            BAx <- t(ABx)
            BBx <- U2

            tuyy[cc1==0,cc1==0] <- AAx
            tuyy[cc1==0,cc1==1] <- ABx
            tuyy[cc1==1,cc1==0] <- BAx
            tuyy[cc1==1,cc1==1] <- BBx


            tu <- U0
          }

        }

        tuyi[i,]   <- t((Zij[i,j])*tuy)
        tui[i]     <- (Zij[i,j])*tu
        tuyyi[[i]] <- (Zij[i,j])*(tuyy)

        if(uni.Gama == FALSE){
          si.mu <- solve(Sigma[[j]]) %*% (tuyi[i,] - (tui[i]*mu[[j]]))
          si.pi <- (Zij[i,j]/as.numeric(pii[j])) - (Zij[i,g]/as.numeric(pii[g]))
          for(t in 1:((p+1)*p/2)) {
            alpha[t] <- -1/2*sum(diag( ( (Zij[i,j]*solve(Sigma[[j]])) %*% deriv.sigma(Sigma[[j]],t,p) ) -
                                         ((tuyyi[[i]] - mu[[j]]%*%t(tuyi[i,]) - tuyi[i,]%*%t(mu[[j]]) + ((tui[i]*mu[[j]])%*%t(mu[[j]]))) %*%
                                            solve(Sigma[[j]]) %*% deriv.sigma(Sigma[[j]],t,p) %*% solve(Sigma[[j]])) ))

          }
          si <- rbind(si,c(si.mu, alpha, si.pi))
        } else {
          si.mu <- solve(Sigma[[j]]) %*% (tuyi[i,] - (tui[i]*mu[[j]]))
          si.pi <- (Zij[i,j]/as.numeric(pii[j])) - (Zij[i,g]/as.numeric(pii[g]))
          for(t in 1:((p+1)*p/2)) {
            alpha[t] <- -1/2*sum(diag( ( (Zij[i,j]*solve(Sigma[[j]])) %*% deriv.sigma(Sigma[[j]],t,p) ) -
                                         ((tuyyi[[i]] - mu[[j]]%*%t(tuyi[i,]) - tuyi[i,]%*%t(mu[[j]]) + ((tui[i]*mu[[j]])%*%t(mu[[j]]))) %*%
                                            solve(Sigma[[j]]) %*% deriv.sigma(Sigma[[j]],t,p) %*% solve(Sigma[[j]])) ))
          }
          alpha.2 <- rbind(alpha.2,alpha)
          si <- rbind(si,c(si.mu,si.pi))
        }

      } ## end for de i

      if(uni.Gama == TRUE) alpha.g <- alpha.g + alpha.2
      si.g    <- cbind(si.g, si)
      si      <- NULL
      alpha.2 <- NULL

    }   ##FIM DO FOR  de g

    if(uni.Gama == FALSE){
      si.g <- si.g[,-dim(si.g)[2]]
      if(g==3) si.g <- cbind(si.g[,-c((p+(p+1)*p/2 + 1),(p+(p+1)*p/2 + 1)*2)],si.g[,c((p+(p+1)*p/2 + 1),(p+(p+1)*p/2 + 1)*2)])
      if(g==2) si.g <- cbind(si.g[,-(p+(p+1)*p/2 + 1)],si.g[,(p+(p+1)*p/2 + 1)])
      for(v in 1:n) MI = MI + si.g[v,]%*%t(si.g[v,])
      colnames(MI) <- rownames(MI) <- nombres(g, p, order = "FALSE")
      #round(MI); print(paste("det(MI) =",det(MI),sep=" "))
     IM <- sqrt(diag(solve(MI)))
    } else {
      si.g <- si.g[,-dim(si.g)[2]]
      si.g <- cbind(si.g,alpha.g)
      for(v in 1:n) MI = MI + si.g[v,]%*%t(si.g[v,])
      colnames(MI) <- rownames(MI) <- nombres(g, p, uni.Gama = "TRUE")
      #round(MI); print(paste("det(MI) =",det(MI),sep=" "))
      IM <- sqrt(diag(solve(MI)))
    }
    return(list(MI=IM))
  } #end family t


}# end function
