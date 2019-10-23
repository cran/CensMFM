## EM function ##


fit.EM.MSNC <- function(cc, LI, LS, y, mu=NULL, Sigma = NULL, shape=NULL, nu = NULL, pii = NULL, g = 3, get.init = TRUE,
                        criteria = TRUE, group = FALSE, family = "SN", error = 0.00001,
                        iter.max = 350, uni.Gama = FALSE, obs.prob= FALSE, kmeans.param = NULL, cal.im = FALSE)
{ # begin function

  # Some Validation #
  y <- as.matrix(y)
  dimnames(y) <- NULL
  if(!is.matrix(y)) stop("The response is not in a matrix format\n")

  if(is.list(family)==TRUE) {
    if(family[[1]] != family[[2]]) stop("Diferent family\n")
    family = family[[1]]
  }

  if(!is.matrix(y)) stop("The response is not in a matrix format\n")
  if((family != "SN") && (family != "Normal") && (family != "t")) stop(paste("Family",family,"not recognized.",sep=" "))
  if((length(g) == 0) && ((length(mu)==0) || (length(Sigma)==0) || (length(pii)==0)))  stop("The model is not specified correctly.\n")

  if((length(g)!= 0) && (g < 1)) stop("g must be greater than 0.\n")

  # Some variables #
  n   <- nrow(y)
  p   <- ncol(y)

  # Initial values #
  if(get.init == TRUE){

    if(length(g) == 0) stop("g is not specified correctly.\n")

    k.iter.max <- 50
    n.start    <- 1
    algorithm  <- "Hartigan-Wong"

    if(is.element(NA, y)){

      y1 = na.omit(y)
      if(g > 1){

        init <- kmeans(y1,g,k.iter.max,n.start,algorithm, trace = FALSE)
        pii  <- init$size/dim(y1)[1] # sum(pii)==0.5 tem que ser 1
        mu   <- Sigma <- shape <-list()

        for (s in 1:g){
          mu[[s]] <- as.vector(init$centers[s,])
          Sigma[[s]] <- var(y1[init$cluster == s,])
          shape[[s]] <- sign(apply((y1[init$cluster == s, ] - matrix(rep(mu[[s]], nrow(y1[init$cluster == s, ])), nrow = nrow(y1[init$cluster == s, ]), ncol=p, byrow = TRUE))^3, 2, sum))
          dimnames(Sigma[[s]]) <- NULL
          names(mu[[s]])       <- NULL
          names(shape[[s]]) <- NULL
        }
      }else{
        mu  <- Sigma <- shape<-list()
        pii <- 1
        mu[[1]]    <- as.vector(colMeans(y1))
        Sigma[[1]] <- var(y1)
        shape[[1]]<-matrix(3*sign(apply(y1,2,skewness)),p,1)
        dimnames(Sigma[[1]]) <- NULL
        names(mu[[1]])       <- NULL
        names(shape[[1]])   <- NULL
      }
    }else{
      if(g > 1){

        init <- kmeans(y,g,k.iter.max,n.start,algorithm, trace = FALSE)
        pii  <- init$size/dim(y)[1] # sum(pii)==0.5 tem que ser 1
        mu   <- Sigma <- shape <-list()

        for (s in 1:g){
          mu[[s]] <- as.vector(init$centers[s,])
          Sigma[[s]] <- var(y[init$cluster == s,])
          shape[[s]] <- sign(apply((y[init$cluster == s, ] - matrix(rep(mu[[s]], nrow(y[init$cluster == s, ])), nrow = nrow(y[init$cluster == s, ]), ncol=p, byrow = TRUE))^3, 2, sum))
          dimnames(Sigma[[s]]) <- NULL
          names(mu[[s]])       <- NULL
          names(shape[[s]]) <- NULL
        }
      }else {
        mu  <- Sigma <- shape<-list()
        pii <- 1
        mu[[1]]    <- as.vector(colMeans(y))
        Sigma[[1]] <- var(y)
        shape[[1]]<-matrix(3*sign(apply(y,2,skewness)),p,1)
        dimnames(Sigma[[1]]) <- NULL
        names(mu[[1]])       <- NULL
        names(shape[[1]])   <- NULL
      }
    }
  }# end get.initial true

  if(get.init == FALSE){

    if(length(g) == 0) stop("g is not specified correctly.\n")

    if(g > 1){
      for (j in 1:g){
        dimnames(Sigma[[j]]) <- NULL
        names(mu[[j]])       <- NULL
        names(shape[[j]])   <- NULL
      }
    }else{
      dimnames(Sigma[[1]]) <- NULL
      names(mu[[1]])       <- NULL
      names(shape[[1]])   <- NULL
    }
  }# end get.initial false

  # sort the data by pii weight
  if(family == "SN"){
    if( sqrt(sum((order(pii)-g:1)*(order(pii)-g:1))) != 0 ){
      musi  <- mu
      Sigsi <- Sigma
      shapesi<-shape
      for(l in 1:g){
        mu[[l]]    <- musi[[(order(pii, decreasing = TRUE)[l])]]
        Sigma[[l]] <- Sigsi[[(order(pii, decreasing = TRUE)[l])]]
        shape[[l]] <- shapesi[[(order(pii, decreasing = TRUE)[l])]]
      }
      pii <- pii[order(pii, decreasing = TRUE)]
    }
  } else {
    if( sqrt(sum((order(pii)-g:1)*(order(pii)-g:1))) != 0 ){
      musi  <- mu
      Sigsi <- Sigma
      for(l in 1:g){
        mu[[l]]    <- musi[[(order(pii, decreasing = TRUE)[l])]]
        Sigma[[l]] <- Sigsi[[(order(pii, decreasing = TRUE)[l])]]
      }
      pii <- pii[order(pii, decreasing = TRUE)]
    }
  }


  # begin SN family
  if (family == "SN"){

    start.time  <- Sys.time()

    delta <- Delta <- Gamma <- list()

    for (k in 1:g){
      delta[[k]] <- shape[[k]] / as.numeric(sqrt(1 + t(shape[[k]])%*%shape[[k]]))
      Delta[[k]] <- as.vector(sqrtm(Sigma[[k]])%*%delta[[k]])
      Gamma[[k]] <- Sigma[[k]] - Delta[[k]]%*%t(Delta[[k]])

    }

    if(uni.Gama){
      Gamma.uni <- Gamma[[1]]
      if(g > 1) for(k in 2:g) Gamma.uni <- Gamma.uni + Gamma[[k]]
      Gamma.uni <- Gamma.uni/g
      Gamma.uni <- (Gamma.uni + t(Gamma.uni))/2
      for(k in 1:g){Gamma[[k]] <- Gamma.uni
      Sigma[[k]] <- Gamma[[k]] + Delta[[k]]%*%t(Delta[[k]])
      shape[[k]] <- (solve(sqrtm(Sigma[[k]]))%*%Delta[[k]])/as.numeric((1 - t(Delta[[k]])%*%solve(Sigma[[k]])%*%Delta[[k]]) )^(1/2)
      }
    }

    criterio <- 1
    count    <- 0
    lkante   <- 1
    yest<-matrix(nrow = n,ncol = p)

    while((criterio > error) && (count <= iter.max)){

      count <- count + 1
      #print(count)

      Zij <- matrix(0, n, g) # cluster variable
      lik <- matrix(0,n,g)

      # begin groups (j)
      for(j in 1:g){

        # sums for M-step
        soma1 <- matrix(0,p,1) # mu
        soma2 <- matrix(0,p,1) # Delta
        soma3 <- 0 # E(T^2)
        soma4 <- matrix(0,p,p) # Gamma

        # sort the data by pii weight
        if( sqrt(sum((order(pii)-g:1)*(order(pii)-g:1))) != 0 ){
          musi  <- mu
          Sigsi <- Sigma
          shapesi <- shape
          for(l in 1:g){
            mu[[l]]    <- musi[[(order(pii, decreasing = TRUE)[l])]]
            Sigma[[l]] <- Sigsi[[(order(pii, decreasing = TRUE)[l])]]
            shape[[l]] <- shapesi[[(order(pii, decreasing = TRUE)[l])]]
          }
          pii <- pii[order(pii, decreasing = TRUE)]
        }

        d1 <- VerCensSN(cc, LI, LS, y, mu[[j]], Sigma[[j]], shape[[j]] )
        if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin

        d2 <-  d.mixedSNCens(cc, LI, LS, y, pii, mu, Sigma, shape)
        if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin

        Zij[,j] <- (d1*pii[j])/(d2) # M-step Zij

        pii[j] <- (1/n)*sum(Zij[,j]) # M-step pii

        SSj   <- sqrtm(Sigma[[j]])
        iSSj  <- solve(SSj)
        M2     <- 1/(1+sum(shape[[j]]^2))
        varphi <- iSSj%*%shape[[j]]

        for(i in 1:n){# begin ind. (i)
          #print(i)
          cc1 <- matrix(cc[i,],p,1)
          LI1<-matrix(LI[i,],p,1)
          LS1<-matrix(LS[i,],p,1)
          y1  <- matrix(y[i,],p,1)

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
              Sigma[[j]] = (Sigma[[j]] + t(Sigma[[j]]))/2
              aux3   = meanvarTMD(lower = c(LI1),upper = c(LS1),mu = mu[[j]],Sigma = Sigma[[j]],lambda = shape[[j]],dist = "SN")
              y.hat  = aux3$mean
              yy.hat =  aux3$EYY

              #print(y.hat)

              w0.hat   = onlymeanTMD(lower = c(LI1),upper = c(LS1),mu = c(mu[[j]]),Sigma = Gamma[[j]],dist = "normal")

              Ltemp  = as.numeric(pmvnorm(lower = c(LI1),upper = c(LS1),mean = c(mu[[j]]),sigma = Gamma[[j]],algorithm = GenzBretz(maxpts = 25000))[1])
              LLtemp = pmvESN(lower = c(LI1),upper = c(LS1),mu = as.vector(mu[[j]]),Sigma = Sigma[[j]],lambda = shape[[j]],tau = 0,algorithm = GenzBretz(maxpts = 25000))
              gamma  = eta*Ltemp/LLtemp

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
              #print(Sc)
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

              Ltemp  = as.numeric(pmvnorm(c(LI1[cc1==1]),c(LS1[cc1==1]),mean = c(muc - mub.cc),sigma = Gamma.cc,algorithm = GenzBretz(maxpts = 25000))[1])
              LLtemp = pmvESN(lower = LI1[cc1==1],upper = LS1[cc1==1],mu = as.vector(muc),Sigma = Sc,lambda = lambda.co,tau = tau.co,algorithm = GenzBretz(maxpts = 25000))
              gamma.cc = eta.cc*Ltemp/LLtemp

              # ratio approximation

              val         = c(tau.co + t(lambda.co)%*%iSS.cc%*%(w1.hat - muc))
              #print(val)
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

          # M-step
          soma1 <- soma1  + (Zij[i,j])*(y.hat - t1.hat*Delta[[j]])
          soma2 <- soma2  + (Zij[i,j])*(ty.hat - t1.hat*mu[[j]])
          soma3 <- soma3  + (Zij[i,j])*t2.hat
          soma4 <- soma4  + (Zij[i,j])*(yy.hat + mu[[j]]%*%t(mu[[j]]) + t2.hat*Delta[[j]]%*%t(Delta[[j]]) - (mu[[j]]%*%t(y.hat) + ty.hat%*%t(Delta[[j]]) - t1.hat*Delta[[j]]%*%t(mu[[j]]) + (y.hat%*%t(mu[[j]]) + Delta[[j]]%*%t(ty.hat) - t1.hat*mu[[j]]%*%t(Delta[[j]]))))

          yest[i,] = t(y.hat)

        }# end ind. (i)

        mu[[j]] <- soma1/sum(Zij[,j])
        Delta[[j]] <- soma2/soma3
        Gamma[[j]]<-  soma4/sum(Zij[,j])
        Gamma[[j]] <- (Gamma[[j]]+t(Gamma[[j]]) )/2

        #recovering
        Sigma[[j]]  = Gamma[[j]] + Delta[[j]]%*%t(Delta[[j]])
        SSj   = sqrtm(Sigma[[j]])
        iSSj  = solve(SSj)
        shape[[j]] = iSSj%*%Delta[[j]]/as.numeric(sqrt(1 - t(Delta[[j]])%*%solve(Sigma[[j]])%*%Delta[[j]]))

        pii[g] <- 1 - (sum(pii) - pii[g])

      }# end groups (j)


      lk <- sum(log(d.mixedSNCens(cc, LI, LS, y, pii, mu, Sigma,shape) ))
      (lk)
      criterio <- abs((lk/lkante-1))
      lkante <- lk

    }# end while

    end.time <- Sys.time()
    time.taken <- end.time - start.time

    if(criteria == TRUE){
      if(uni.Gama){
        d <- g*2*p + length(Sigma[[1]][upper.tri(Sigma[[1]], diag = TRUE)]) + (g-1) #mu + shape + Sigma + pi
      }
      else{
        d <- g*(2*p + length(Sigma[[1]][upper.tri(Sigma[[1]], diag = TRUE)]) ) + (g-1) #mu + shape + Sigma + pi
      }
      aic <- -2*lk + 2*d
      bic <- -2*lk + log(n)*d
      edc <- -2*lk + 0.2*sqrt(n)*d
      out <- list(mu = mu, Sigma = Sigma, shape = shape, Gamma = Gamma, pii = pii, Zij = Zij, yest = yest, logLik = lk, aic = aic, bic = bic, edc = edc, iter = count, n = nrow(y), group = apply(Zij, 1, which.max),time = time.taken, family = family)
    }

    if(criteria == FALSE) out <- list(mu = mu, Sigma = Sigma, shape = shape, Gamma = Gamma, pii = pii, Zij = Zij, yest = yest, iter = count, n = nrow(y), group = apply(Zij, 1, which.max),time = time.taken, family = family)

    for(i in 1:length(out$Sigma)) out$Sigma[[i]] <- sqrtm(out$Sigma[[i]])

    if(cal.im){
      out$MI = imm.msnc.new(y,out,cc,LI,LS,family,uni.Gama)
    }

    return(out)
  } # end family SN

  # Begin family Normal
  if (family == "Normal"){

    start.time  <- Sys.time()

    if(uni.Gama){
      Sigma.uni <- Sigma[[1]]
      if(g > 1) for(k in 2:g) Sigma.uni <- Sigma.uni + Sigma[[k]]
      Sigma.uni <- Sigma.uni / g
      for(k in 1:g) Sigma[[k]] <- Sigma.uni
    }

    criterio <- 1
    count    <- 0
    lkante   <- 1
    yest <- matrix(nrow = n,ncol = p)

    while((criterio > error) && (count <= iter.max)){

      count <- count + 1
      #print(count)

      Zij   <- matrix(0, n, g)


      for (j in 1:g){

        #Ordenar por tamanho das componentes de pii
        if( sqrt(sum((order(pii)-g:1)*(order(pii)-g:1))) != 0 ){
          musi  <- mu
          Sigsi <- Sigma
          for(l in 1:g){
            mu[[l]]    <- musi[[(order(pii, decreasing = TRUE)[l])]]
            Sigma[[l]] <- Sigsi[[(order(pii, decreasing = TRUE)[l])]]
          }
          pii <- pii[order(pii, decreasing = TRUE)]
        }

        soma1 <- matrix(0,p,1)
        soma2 <- matrix(0,p,p)

        d1 <- dmvNCens(cc, LI, LS, y, mu[[j]], Sigma[[j]])
        if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin

        d2 <-  d.mixedNCens(cc, LI, LS, y, pii, mu, Sigma)
        if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin

        Zij[,j] <- d1*pii[j] / d2

        ### M-step: atualizar mu, Sigma and pii ###

        pii[j] <- (1/n)*sum(Zij[,j])


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
              muc<-mu[[j]]
              Sc<-Sigma[[j]]
              aux3   = meanvarTMD(lower = c(LI1),upper = c(LS1),mu = mu[[j]],Sigma = Sigma[[j]],dist = "normal")
              uy  = aux3$mean
              uyy =  aux3$EYY

            }

            else {
              muc = mu[[j]][cc1==1]+Sigma[[j]][cc1==1,cc1==0]%*%solve(Sigma[[j]][cc1==0,cc1==0])%*%(y1[cc1==0]-mu[[j]][cc1==0])
              Sc <-Sigma[[j]][cc1==1,cc1==1]-Sigma[[j]][cc1==1,cc1==0]%*%solve(Sigma[[j]][cc1==0,cc1==0])%*%Sigma[[j]][cc1==0,cc1==1]
              Sc<-(Sc+t(Sc))/2
              aux3   = meanvarTMD(lower = c(LI1[cc1==1]),upper = c(LS1[cc1==1]),mu = as.vector(muc),Sigma = Sc,dist = "normal")
              w1.hat = aux3$mean
              w2.hat =  aux3$EYY
              uy  = matrix(y1,p,1)
              uy[cc1==1] = w1.hat
              uyy =  y1%*%t(y1)
              uyy[cc1==0,cc1==1] = y1[cc1==0]%*%t(w1.hat)
              uyy[cc1==1,cc1==0] = w1.hat%*%t(y1[cc1==0])
              uyy[cc1==1,cc1==1] = w2.hat
            }

          }

          soma1<- soma1 +  (Zij[i,j])*uy
          soma2<- soma2 +  (Zij[i,j])*(uyy-mu[[j]]%*%t(uy)-(uy)%*%t(mu[[j]])+mu[[j]]%*%t(mu[[j]]))

          yest[i,] <- t(uy)


        } ## end for de i

        mu[[j]]    <- soma1 / sum(Zij[,j])
        Sigma[[j]] <- soma2 / sum(Zij[,j])
        Sigma[[j]] <- (Sigma[[j]]+t(Sigma[[j]]) )/2


        if(uni.Gama == TRUE){
          GS <- 0
          for (w in 1:g) GS <- GS+kronecker(Zij[,w],Sigma[[w]])
          Sigma.uni <- t(rowSums(array(t(GS),dim=c(p,p,n)),dims=2))/n
          for (w in 1:g) Sigma[[w]] <- Sigma.uni
        }


        pii[g] <- 1 - (sum(pii) - pii[g])

        zero.pos <- NULL
        zero.pos <- which(pii == 0)
        if(length(zero.pos) != 0){
          pii[zero.pos] <- 1e-10
          pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
        }


      }   ##FIM DO FOR  de g

      lk <- sum(log(d.mixedNCens(cc, LI, LS, y, pii, mu, Sigma) ))
      #print(lk)
      criterio <- abs((lk/lkante-1))
      lkante <- lk

    } # fim do while
    end.time <- Sys.time()
    time.taken <- end.time - start.time

    if(criteria == TRUE){
      if(uni.Gama){
        d <- g*p + length(Sigma[[1]][upper.tri(Sigma[[1]], diag = TRUE)]) + (g-1) #mu + Sigma + pi
      }
      else{
        d <- g*(p + length(Sigma[[1]][upper.tri(Sigma[[1]], diag = TRUE)]) ) + (g-1) #mu + Sigma + pi
      }
      aic <- -2*lk + 2*d
      bic <- -2*lk + log(n)*d
      edc <- -2*lk + 0.2*sqrt(n)*d
      out <- list(mu = mu, Sigma = Sigma, pii = pii, Zij = Zij, yest = yest, logLik = lk, aic = aic, bic = bic, edc = edc, iter = count, n = nrow(y), group = apply(Zij, 1, which.max),time = time.taken, family = family)
    }

    if(criteria == FALSE) out <- list(mu = mu, Sigma = Sigma, pii = pii, Zij = Zij, yest = yest, iter = count, n = nrow(y), group = apply(Zij, 1, which.max),time = time.taken, family = family)

    if(cal.im){
      out$MI = imm.msnc.new(y,out,cc,LI,LS,family,uni.Gama)
    }

    return(out)

  } # fim family Normal

  if (family == "t"){

    start.time  <- Sys.time()

    if(uni.Gama){
      Sigma.uni <- Sigma[[1]]
      if(g > 1) for(k in 2:g) Sigma.uni <- Sigma.uni + Sigma[[k]]
      Sigma.uni <- Sigma.uni / g
      for(k in 1:g) Sigma[[k]] <- Sigma.uni
    }

    criterio <- 1
    count    <- 0
    lkante   <- 1
    yest <- matrix(nrow = n,ncol = p)

    while((criterio > error) && (count <= iter.max)){

      count <- count + 1
      #print(count)

      Zij   <- matrix(0, n, g)
      tuyi  <- matrix(0,n,p)
      tui   <- matrix(0,n,1)
      tuyyi <- list()

      for (j in 1:g){
        #Ordenar por tamanho das componentes de pii
        if( sqrt(sum((order(pii)-g:1)*(order(pii)-g:1))) != 0 ){
          musi  <- mu
          Sigsi <- Sigma
          for(l in 1:g){
            mu[[l]]    <- musi[[(order(pii, decreasing = TRUE)[l])]]
            Sigma[[l]] <- Sigsi[[(order(pii, decreasing = TRUE)[l])]]
          }
          pii <- pii[order(pii, decreasing = TRUE)]
        }

        soma1 <- matrix(0,p,1)
        soma2 <- 0
        soma3 <- matrix(0,p,p)

        d1 <- dmvTCens(cc, LI, LS, y, mu[[j]], Sigma[[j]],nu = nu[j])
        if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin # de donde pertenece

        d2 <-  d.mixedTCens(cc, LI, LS, y, pii, mu, Sigma, nu)
        if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin # de donde pertenece

        Zij[,j] <- d1*pii[j] / d2

        ### M-step: atualizar mu, Sigma and pii ###

        pii[j] <- (1/n)*sum(Zij[,j])
        #######
        mus    <- mu[[j]]
        Sigmas <- Sigma[[j]]
        #######

        for (i in 1:n){
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

          }

          if(sum(cc1)>0){

            if(sum(cc1) == p){

              muUi     <- mus
              SigmaUi  <- Sigmas
              SigmaUiA <- SigmaUi*nu[j]/(nu[j]+2)
              auxupper <- y1-muUi
              auxU1    <- sadmvt(df=nu[j]+2, as.vector(LI1), as.vector(LS1), as.vector(muUi),SigmaUiA)
              auxU2    <- sadmvt(df=nu[j], as.vector(LI1), as.vector(LS1), as.vector(muUi), SigmaUi)
              #auxU1    <- pmvt(df=nu[j]+2, lower = as.vector(LI1), upper = as.vector(LS1), delta = as.vector(muUi),sigma = SigmaUiA)
              #auxU2    <- pmvt(df=nu[j], lower = as.vector(LI1), upper = as.vector(LS1), delta = as.vector(muUi), sigma= SigmaUi)
              MoMT     <- meanvarTMD(as.vector(LI1),as.vector(LS1),as.vector(muUi),SigmaUiA,dist = "t",nu = nu[j]+2)
              U0       <- as.numeric(auxU1/auxU2)
              U1       <- auxU1/auxU2*MoMT$mean
              U2       <- auxU1/auxU2*MoMT$EYY

              tuy  <- U1
              tuyy <- U2
              tu   <- U0

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

              auxU1 <- sadmvt(df=nu1+2, as.vector(LI1[cc1==1]), as.vector(LS1[cc1==1]), as.vector(muUi), SigmaUiA)
              auxU2 <- sadmvt(df=nu1, as.vector(LI1[cc1==1]), as.vector(LS1[cc1==1]), as.vector(muUi), SigmaUi)
              #auxU1    <- pmvt(df=nu1+2, lower = as.vector(LI1[cc1==1]), upper = as.vector(LS1[cc1==1]), delta = as.vector(muUi),sigma = SigmaUiA)
              #auxU2    <- pmvt(df=nu1, lower = as.vector(LI1[cc1==1]), upper = as.vector(LS1[cc1==1]), delta = as.vector(muUi), sigma= SigmaUi)

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

          soma1<- soma1 + (Zij[i,j])*tuy
          soma2<- soma2 + (Zij[i,j])*tu
          soma3<- soma3 + (Zij[i,j])*(tuyy-(tuy)%*%t(mus)-(mus)%*%t(tuy)+ tu*mus%*%t(mus))

          yest[i,] <- t(tuy)

        } ## end for de i

        mu[[j]]    <- soma1 / soma2
        Sigma[[j]] <- soma3 / sum(Zij[,j])
        Sigma[[j]] <- (Sigma[[j]]+t(Sigma[[j]]))/2

        if(uni.Gama == TRUE){
          GS <- 0
          for (w in 1:g) GS <- GS + kronecker(Zij[,w],Sigma[[w]])
          Sigma.uni <- t(rowSums(array(t(GS),dim=c(p,p,n)),dims=2))/n
          for (w in 1:g)  Sigma[[w]] <- Sigma.uni
        }


        pii[g] <- 1 - (sum(pii) - pii[g])

        zero.pos <- NULL
        zero.pos <- which(pii == 0)
        if(length(zero.pos) != 0){
          pii[zero.pos] <- 1e-10
          pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
        }

      }   ##FIM DO FOR  de g

      lk <- sum(log(d.mixedTCens(cc,LI,LS, y, pii, mu, Sigma, nu) ))
      #print(lk)
      criterio <- abs((lk/lkante-1))
      lkante <- lk

    } # fim do while
    end.time        <- Sys.time()
    time.taken      <- end.time - start.time

    if(criteria == TRUE){
      if(uni.Gama){
        d <- g*p + length(Sigma[[1]][upper.tri(Sigma[[1]], diag = TRUE)]) + (g-1) #mu + Sigma + pi
      }
      else{
        d <- g*(p + length(Sigma[[1]][upper.tri(Sigma[[1]], diag = TRUE)]) ) + (g-1) #mu + Sigma + pi
      }
      aic <- -2*lk + 2*d
      bic <- -2*lk + log(n)*d
      edc <- -2*lk + 0.2*sqrt(n)*d
      out <- list(mu = mu, Sigma = Sigma, pii = pii, nu = nu, Zij = Zij, yest = yest, logLik = lk, aic = aic, bic = bic, edc = edc, iter = count, n = nrow(y), group = apply(Zij, 1, which.max),time = time.taken,family = family)
    }

    if(criteria == FALSE) out <- list(mu = mu, Sigma = Sigma, pii = pii, nu = nu, Zij = Zij, yest = yest, iter = count, n = nrow(y), group = apply(Zij, 1, which.max),time = time.taken,family = family)

    if(cal.im){
      out$MI = imm.msnc.new(y,out,cc,LI,LS,family,uni.Gama)
    }

    return(out)

  } # fim family t



} # end EM function
