fit.FMMSNC <-
function(cc, LI, LS, y, mu=NULL, Sigma = NULL, shape=NULL, pii = NULL, nu = NULL,
g = NULL, get.init = TRUE, criteria = TRUE, family = "SN", error = 0.00001,
iter.max = 350, uni.Gama = FALSE, kmeans.param = NULL, cal.im = FALSE)
{
  #Running the algorithm
out <- fit.EM.MSNC(cc=cc, LI=LI, LS=LS, y=y, mu=mu, Sigma=Sigma, shape=shape,pii=pii,
nu=nu, g=g, get.init = get.init,criteria = criteria, family = family,error = error,
iter.max = iter.max, uni.Gama = uni.Gama, kmeans.param = kmeans.param, cal.im = cal.im)

  #show result
  message('\n')
  message('-------------------------------------------------------------')
  message('  Finite mixture modeling of censored and missing data using ')
  message('         the multivariate skew-normal distribution           ')
  message('-------------------------------------------------------------')
  message('\n')
  message('-----------')
  message('  Details  ')
  message('-----------')
  #message('\n')
  message('Iterations = ',out$iter)
  #message('\n')
  message("Processing time = ",out$time,units(out$time))
  #message('\n')
  message('Family = ',family)
  #message('\n')
  message('Components = ',g)
 # message('\n')
  message('Row = ',nrow(y),", ","Columns= ",ncol(y))
  message('\n')
  message('-------------')
  message('  Estimates  ')
  message('-------------')
  #message('\n')
  #Mu vector
  namesrowMedj <- matrix("",ncol(y),g)
  ttable       <- data.frame(out$mu)
  for(i in 1:g){for(j in 1:ncol(y)){namesrowMedj[j,i]   <- paste("mu",i,j,sep="")}}
  muM <- matrix(0,ncol(y),1)
  for(j in 1:g)
  {
    muM <- as.matrix(ttable[,j])
    rownames(muM) <- c(namesrowMedj[,j])
    colnames(muM) <- paste("mu",j,sep="")
    print(muM)
    message('\n')
  }
  #message('\n')

  if(family=="SN")
  {
    for(k in 1:g)
    {
      a <- out$Sigma[[k]]
      a <- a[lower.tri(a, diag=TRUE)]
      b <- matrix(0, ncol(y), ncol(y))
      b[lower.tri(b, diag=TRUE)] <- a
      b <- t(b)
      b <- round(b,digits = 4)

      namesrowSigmas   <- c()
      for(j in 1:ncol(y)){namesrowSigmas[j] <- paste("alpha",j,sep="")}
      rownames(b) <- colnames(b) <- namesrowSigmas
      message("F",k)
      #message('\n')
      print(b)
      message('\n')
    }

    #Shape vector
    namesrowShj <- matrix("",ncol(y),g)
    ttable       <- data.frame(out$shape)
    for(i in 1:g){for(j in 1:ncol(y)){namesrowShj[j,i]   <- paste("shape",i,j,sep="")}}
    ShM <- matrix(0,ncol(y),1)
    for(j in 1:g)
    {
      ShM <- as.matrix(ttable[,j])
      rownames(ShM) <- c(namesrowShj[,j])
      colnames(ShM) <- paste("Shape",j,sep="")
      print(ShM)
      message('\n')
    }
    #message('\n')
  }else{
    #Sigma matirx
    for(k in 1:g)
    {
      a <- out$Sigma[[k]]
      a <- a[lower.tri(a, diag=TRUE)]
      b <- matrix(0, ncol(y), ncol(y))
      b[lower.tri(b, diag=TRUE)] <- a
      b <- t(b)
      b <- round(b,digits = 4)

      namesrowSigmas   <- c()
      for(j in 1:ncol(y)){namesrowSigmas[j] <- paste("sigma",j,sep="")}
      rownames(b) <- colnames(b) <- namesrowSigmas
      message("Sigma",k)
      #message('\n')
      print(b)
      message('\n')
    }
  }

  #Pii vector
  if(g != 1){
    namesrowPij   <- matrix("",ncol = 1,nrow = g-1)
    PiM <- matrix(out$pii[-g],ncol = 1,nrow = g-1)
    for(j in 1:(g-1)) {namesrowPij[j,] <- paste("pii_",j,sep="")}
    rownames(PiM) <- namesrowPij
    colnames(PiM) <- "Pi"
    print(PiM)
   message('\n')
  }else{
    namesrowPij   <- matrix("",ncol = 1,nrow = g)
    PiM <- matrix(out$pii,ncol = 1,nrow = g)
    for(j in 1:g) {namesrowPij[j,] <- paste("pii_",j,sep="")}
    rownames(PiM) <- namesrowPij
    colnames(PiM) <- "Pi"
    print(PiM)
    message("\n")
  }

  if(family=="t"){
    namesrownu   <- matrix("",ncol = 1,nrow = g)
    Nu <- matrix(out$nu,ncol = 1,nrow = g)
    for(j in 1:g) {namesrownu[j,] <- paste("nu_",j,sep="")}
    rownames(Nu) <- namesrownu
    colnames(Nu) <- "Nu"
    print(Nu)
   message('\n')
  }


  if(cal.im==TRUE){
  message('---------------------------')
  message('  Stantard Error (SE)  ')
  message('---------------------------')
 # message('\n')
  SE = out$MI
  MM       <- data.frame(SE)
  colnames(MM) <- c("SE")
  print(MM)
  message('\n')
  }

  if(criteria==TRUE){
  message('----------------------------')
  message('  Model selection criteria  ')
  message('----------------------------')
  #message('\n')
  critFin           <- c(out$logLik, out$aic, out$bic, out$edc)
  critFin           <- round(t(as.matrix(critFin)),digits=3)
  dimnames(critFin) <- list(c("Value"),c("Loglik", "AIC", "BIC","EDC"))
  print(critFin)
  }
  #message('\n')
  # message('-----------')
  # message('  Details  ')
  # message('-----------')
  # #message('\n')
  # message('Iterations =',out$iter)
  # #message('\n')
  # message("Processing time =",out$time,units(out$time))
  # #message('\n')

  if(family=="SN"){
    if(criteria==TRUE){
  res            <- list(group = out$group, family = out$family, yest = out$yest, SE = out$MI, iter = out$iter, mu=out$mu, Sigma=out$Sigma, Gamma = out$Gamma, shape = out$shape, pii=out$pii, loglik=out$logLik, aic=out$aic, bic=out$bic, edc=out$edc, time = out$time)
  obj.out        <- list(res = res)
  class(obj.out) <-  "FM-MSNC"
  return(obj.out)
    }else{
      res            <- list(group =out$group, yest = out$yest, family = out$family, SE = out$MI, iter = out$iter, mu=out$mu, Sigma=out$Sigma, Gamma = out$Gamma, shape = out$shape, pii=out$pii, time = out$time)
      obj.out        <- list(res = res)
      class(obj.out) <-  "FM-MSNC"
      return(obj.out)
    }
  }
  if(family=="Normal"){
    if(criteria==TRUE){
    res            <- list(group =out$group, yest = out$yest, family = out$family, SE = out$MI, iter = out$iter, mu=out$mu, Sigma=out$Sigma, pii=out$pii, loglik=out$logLik, aic=out$aic, bic=out$bic, edc=out$edc, time = out$time)
    obj.out        <- list(res = res)
    class(obj.out) <-  "FM-MNC"
    return(obj.out)
    }else{
      res            <- list(group =out$group, yest = out$yest, family = out$family, SE = out$MI, iter = out$iter, mu=out$mu, Sigma=out$Sigma, pii=out$pii, time = out$time)
      obj.out        <- list(res = res)
      class(obj.out) <-  "FM-MNC"
      return(obj.out)
    }
  }
  if(family=="t"){
    if(criteria==TRUE){
      res            <- list(group =out$group, yest = out$yest, family = out$family, SE = out$MI, iter = out$iter, mu=out$mu, nu = out$nu, Sigma=out$Sigma, pii=out$pii, loglik=out$logLik, aic=out$aic, bic=out$bic, edc=out$edc, time = out$time)
      obj.out        <- list(res = res)
      class(obj.out) <-  "FM-MtC"
      return(obj.out)
    }else{
      res            <- list(group =out$group, yest = out$yest, family = out$family, SE = out$MI, iter = out$iter, mu=out$mu, nu = out$nu, Sigma=out$Sigma, pii=out$pii, time = out$time)
      obj.out        <- list(res = res)
      class(obj.out) <-  "FM-MtC"
      return(obj.out)
    }
  }
}
