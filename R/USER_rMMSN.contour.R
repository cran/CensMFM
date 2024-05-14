rMMSN.contour <-
function(model = NULL, y = NULL, mu=NULL, Sigma=NULL, shape=NULL,nu = NULL ,
pii=NULL, Zij = NULL, contour = FALSE, hist.Bin = 30, contour.Bin = 10, slice=100, col.names = NULL,
length.x = c(0.5,0.5),length.y = c(0.5,0.5),family = "SN"){

  dens = g = x = x1 = x2 = NULL

  # This function only serves to graph the contours of finite mixtures of bivariate SMSN !!!
  # model: is an object resultant from fit.EM.MSNC
  # slice is a parameter to the contour function
  if((family != "Normal") && (family != "SN") && (family != "t")) stop(paste("Family",family,"not recognized.\n",sep=" "))
  if(is.null(model)){
    if(!all(lengths(mu)) | !all(lengths(Sigma)) | !all(lengths(shape))) stop("Parameters have not conforming size")
    if(length(pii) != length(unique(Zij))) stop("Number of groups have not conforming size")
    if(sum(pii) != 1) stop("Sum of pii diferent of one")
  }else{
    y <- model$res$yest
    Zij <- model$res$group
    mu <- model$res$mu
    if(family == "SN"){
      for(i in 1:length(model$res$Sigma)) model$res$Sigma[[i]] <- model$res$Sigma[[i]]%*%model$res$Sigma[[i]]
      Sigma <- model$res$Sigma
      shape <- model$res$shape
    }else{
      Sigma <- model$res$Sigma
      }
    pii <- model$res$pii
    #nu <- model$res$nu
  }

  if(!is.null(col.names) && length(col.names) != ncol(y)) stop("col.names must be the same length than ncol y \n")
  if(!is.null(col.names) && length(col.names) == ncol(y)) colnames(y) <- col.names

  p = dim(y)[2]

  if(contour == TRUE && p != 2) stop("The contour function is only appropriate for the bivariate analysis.\n")

  graphc = list()

  cont1 = 1
  cont2 = 1

  aux.grap1 = seq_r(p,1,p)
  aux = seq(1:max(aux.grap1))
  aux.grap2 = setdiff(aux,aux.grap1)

  for(i1 in 1:p){

    for(i2 in i1:p){

      if(contour == TRUE & p == 2){
        lim.mim.x <- min(y[,i1]) - length.x[1]
        lim.max.x <- max(y[,i1]) + length.x[2]
        lim.mim.y <- min(y[,i2]) - length.y[1]
        lim.max.y <- max(y[,i2]) + length.y[2]
        d1.SN <- expand.grid("x1" = seq(lim.mim.x, lim.max.x,length.out = slice), "x2" = seq(lim.mim.y, lim.max.y,length.out = slice))
        cc = matrix(0,ncol =2, nrow = dim(d1.SN)[1])
        LI = cc
        LS = cc
        if(family=="SN") d1.SN$dens  <-  mixedMSN(y = as.matrix(d1.SN), pii = pii, mu = mu, Sigma = Sigma , shape = shape)
        if(family=="Normal") d1.SN$dens  <-  mixedMN(y = as.matrix(d1.SN), pii = pii, mu = mu, Sigma = Sigma)
        if(family=="t") d1.SN$dens  <-  mixedMT(y = as.matrix(d1.SN), pii = pii, mu = mu, Sigma = Sigma , nu = nu)

        # mat2 = as.matrix(d1.SN)
        #
        # matx1 = aggregate(dens ~ x1, mat2, sum)
        # matx2 = aggregate(dens ~ x2, mat2, sum)

        # plot(matx1[,1], matx1[,2],type = "l")
        # plot(matx2[,1], matx2[,2],type = "l")
        #
        # dfx1 <- data.frame(x1 = matx1[,1], dens = matx1[,2])
        # dfx2 <- data.frame(x2 = matx2[,1], dens = matx2[,2])

        df <- data.frame(x = y[,i1], y = y[,i2],g = Zij)

        if(i1 == i2){
          ## Hist
          df.hist <- data.frame(x = y[,cont1],g = as.factor(Zij))
          if(is.null(colnames(y))){
            df.hist$title = paste("y", cont1, sep = "")
            commonTheme = list(labs(x="", y=""))
          }else{
            df.hist$title = paste(colnames(y)[cont1])
            commonTheme = list(labs(x="", y=""))}

          graphc[[aux.grap1[cont1]]] = ggplot(df.hist, aes(x=x,color = g)) +
            geom_histogram(fill="gray",alpha=0.5, position="identity", bins = hist.Bin) + #,aes(y=..density..)
            commonTheme+
            # geom_line(data = dfx1,aes(x = x1,y = dens),colour = "gray50") +
            theme_bw(base_size = 12)+
            theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank()) +
            facet_grid(. ~ title)


          cont1 = cont1 + 1
        }else{
          if(is.null(colnames(y))){
            commonTheme = list(labs(x=paste("y", i1, sep = ""), y=paste("y", i2, sep = "")))
          }else{
            commonTheme = list(labs(x=paste(colnames(y)[i1]), y=paste(colnames(y)[i2])))
          }
          graphc[[aux.grap2[cont2]]] = ggplot(df, aes(x = x, y = y)) +
            geom_point(aes(shape=factor(g), color=factor(g)),size = 1) +
            geom_contour(data = d1.SN, aes(x1, x2, z = dens),colour = "gray50",bins = contour.Bin)+
            commonTheme +
            theme_bw(base_size = 12) +
            coord_cartesian(xlim = c(lim.mim.x, lim.max.x), ylim = c(lim.mim.y, lim.max.y)) +
            theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank())
          cont2 = cont2 + 1
        }
      }else{
        lim.mim.x <- min(y[,i1]) - length.x[1]
        lim.max.x <- max(y[,i1]) + length.x[2]
        lim.mim.y <- min(y[,i2]) - length.y[1]
        lim.max.y <- max(y[,i2]) + length.y[2]

        # create a ggplot2 scatterplot
        df <- data.frame(x = y[,i1], y = y[,i2],g = Zij)

        if(i1 == i2){
          ## Hist
          df.hist <- data.frame(x = y[,cont1],g = as.factor(Zij))
          if(is.null(colnames(y))){
            df.hist$title = paste("y", cont1, sep = "")
            commonTheme = list(labs(x="", y=""))
          }else{
            df.hist$title = paste(colnames(y)[cont1])
            commonTheme = list(labs(x="", y=""))}
          graphc[[aux.grap1[cont1]]] = ggplot(df.hist, aes(x=x,color = g)) +
            geom_histogram(fill="gray",alpha=0.5, position="identity",bins = hist.Bin) +
            commonTheme+
            theme_bw(base_size = 12)+
            theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank()) +
            facet_grid(. ~ title)

          cont1 = cont1 + 1
        }else{
          if(is.null(colnames(y))){
            commonTheme = list(labs(x=paste("y", i1, sep = ""), y=paste("y", i2, sep = "")))
          }else{
            commonTheme = list(labs(x=paste(colnames(y)[i1]), y=paste(colnames(y)[i2])))
          }
          graphc[[aux.grap2[cont2]]] = ggplot(df, aes(x = x, y = y)) +
            geom_point(aes(shape=factor(g), color=factor(g)),size = 1) +
            commonTheme +
            theme_bw(base_size = 12) +
            coord_cartesian(xlim = c(lim.mim.x, lim.max.x), ylim = c(lim.mim.y, lim.max.y)) +
            theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank())
          cont2 = cont2 + 1
        }
      }
    }#i2
  }#i1

  aux.grap2

  Ma <- matrix(1:p^2, p, p,byrow=FALSE)
  lower.tri(Ma)
  Ma[upper.tri(Ma)] <- NA
  return(
  graphcs = grid.arrange(
    grobs = graphc,
    widths = rep(p,p),
    heights = rep(p,p),
    layout_matrix = Ma
  )
  )

}
