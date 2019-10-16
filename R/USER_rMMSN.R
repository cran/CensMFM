rMMSN <-
function(n = NULL, mu = NULL, Sigma = NULL, shape = NULL, percent = NULL,
each = FALSE, pii = NULL ,family = "SN"){
  #n: Number of generating data;
  #mu: Must be a list with g entries, where each entries must be a vector with dimension p;
  #Sigma: Must be a list with g entries, where each entries must be a mztrix with dimension pxp.
  #shape: Must be a list with g entries, where each entries must be a vector with dimension p;
  #percent: percenttage of censored data in each group or data as a whole.
  #each: If  TRUE indicates that data should be censored in each group,
  #so percent must be a vector of dimension p, if FALSE indicates that data should be
  #censored in the whole set, and percent must be a vector of dimension 1;
  #pii: Must be a vector with dimension g;
  #family: Normal, SN

  if((family != "Normal") && (family != "SN")) stop(paste("Family",family,"not recognized.\n",sep=" "))

  if(sum(pii) != 1) stop("Sum of pii diferent of one")
  if(each==TRUE){
    if(length(percent) != length(mu)) stop("Censure in a component")
    if(sum(percent) > length(mu)) stop("Sum censure greather that one")
  }else{
    if(length(percent) != 1) stop("Censure in a component")
  }

  w <- rmultinom(n, size = 1, prob = pii)
  G <- rowSums(w)
  z <- h <- NULL
  p <- length(mu[[1]])


  if(sum(percent) != 0){

    if(each==TRUE){
      cutof <- matrix(nrow = 1,ncol = length(pii))
      for (r in 1:length(G)){
        y <- matrix(0,nrow=G[r],ncol=p)
        if(family == "Normal") y <- rMSN(n = G[r], mu = mu[[r]], Sigma = Sigma[[r]], shape = rep(0, length(mu[[r]]))) # gerando y
        if(family == "SN")     y <- rMSN(n = G[r], mu = mu[[r]], Sigma = Sigma[[r]], shape = shape[[r]]) # gerando y

        m   <- dim(y)[1]
        aa  <- sort(y,decreasing = FALSE)
        aa1 <- matrix(t(aa),m*p,1)
        cutof[,r] <- aa1[ceiling(percent[r]*m*p)]
        cc    <- (y < cutof[,r])+0
        y[cc==1] <- cutof[,r]
        z <- rbind(z,y)
        h <- rbind(h,cc)

      }
    }
    if(each == FALSE){
      for (r in 1:length(G)){
        y <- matrix(0,nrow=G[r],ncol=p)
        if(family == "Normal") y <- rMSN(n = G[r], mu = mu[[r]], Sigma = Sigma[[r]], shape = rep(0, length(mu[[r]]))) # gerando y
        if(family == "SN")     y <- rMSN(n = G[r], mu = mu[[r]], Sigma = Sigma[[r]], shape = shape[[r]]) # gerando y
        z <- rbind(z,y)
      }
      m   <- dim(z)[1]
      aa  <- sort(z,decreasing = FALSE)
      aa1 <- matrix(t(aa),m*p,1)
      cutof <- aa1[ceiling(percent*m*p)]
      cc    <- (z < cutof)+0
      z[cc==1] <- cutof
      h <- rbind(h,cc)
    }
    G = rep(1:length(G),G)
    return(list(y = z, cc = h, G = G,cutoff = cutof))
  }else{
    for (r in 1:length(G)){
      y <- matrix(0,nrow=G[r],ncol=p)
      if(family == "Normal") y <- rMSN(n = G[r], mu = mu[[r]], Sigma = Sigma[[r]], shape = rep(0, length(mu[[r]]))) # gerando y
      if(family == "SN")     y <- rMSN(n = G[r], mu = mu[[r]], Sigma = Sigma[[r]], shape = shape[[r]]) # gerando y
      z <- rbind(z,y)
    }
    G = rep(1:length(G),G)
    return(list(y = z,G = G))
  }
}
