ellipseYS <- function(center, shape, radius) {
  segments<-51
  if (! (is.vector(center) && 2==length(center))) stop("center must be a vector of length 2")
  if (! (is.matrix(shape) && all(2==dim(shape)))) stop("shape must be a 2 by 2 matrix")
  if (max(abs(shape - t(shape)))/max(abs(shape)) > 1e-10) stop("shape must be a symmetric matrix")
  angles <- (0:segments)*2*pi/segments 
  unit.circle <- cbind(cos(angles), sin(angles)) 
  #	ellipse <- t(center + radius*t(unit.circle %*% chol(shape,pivot=TRUE))) 
  Q <- chol(shape, pivot=TRUE)
  order <- order(attr(Q, "pivot"))
  ellipse <- t( center + radius*t( unit.circle %*% Q[,order]))
  colnames(ellipse) <- c("x", "y")
  ellipse
}

aPCoA<-function (formula,data,maincov,drawEllipse=TRUE,drawCenter=TRUE) 
{ 
  Terms <- terms(formula, data = data)
  lhs <- formula[[2]]
  lhs <- eval(lhs, data, parent.frame())
  formula[[2]] <- NULL
  rhs.frame <- model.frame(formula, data, drop.unused.levels = TRUE)

  rhs <- model.matrix(formula, rhs.frame)

  grps <- attr(rhs, "assign")
  qrhs <- qr(rhs)
  rhs <- rhs[, qrhs$pivot, drop = FALSE]
  rhs <- rhs[, 1:qrhs$rank, drop = FALSE]
  grps <- grps[qrhs$pivot][1:qrhs$rank]
  u.grps <- unique(grps)
  nterms <- length(u.grps) - 1
  if (nterms < 1) 
    stop("right-hand-side of formula has no usable terms")
    dmat <- as.matrix(lhs^2)

  X<-rhs
  y<-lhs

  X<-X[rownames(dmat),]
  X<-as.matrix(X[,-1],nrow=nrow(X))

  H<-X%*%solve(t(X)%*%X)%*%t(X)
  
  A<--1/2*((diag(nrow(H))-H)%*%as.matrix(y)^2%*%(diag(nrow(H))-H))
  J<-diag(nrow(X))-matrix(rep(1/(nrow(X)),length(A)),nrow=nrow(A))
  E<-J%*%A%*%J

  rownames(E)<-rownames(data)
  colnames(E)<-rownames(data)
  eigenE<-eigen(E)$vectors
  eigenvalue<-eigen(E)$values
  tempvector<-as.character(data[,as.character(substitute(maincov))])
  color1<-distinctColorPalette(length(unique(tempvector)))
  names(color1)<-unique(tempvector)
  rownames(eigenE)<-rownames(data)
  centernames<-unique(tempvector)
  centers<-centernames

  for(i in 1:length(centernames)){
    centers[i]<-rownames(pam(as.matrix(y)[tempvector==centernames[i],
                                     tempvector==centernames[i]], 1)$medoids)
  }
  newcenters<-centers
  for(i in 1:length(centernames)){
    newcenters[i]<-rownames(pam(E[tempvector==centernames[i],
                                  tempvector==centernames[i]], 1)$medoids)
  }

  origpcoa<-pcoa(y)


if(drawEllipse){

  xMaxPlot<-NULL
  xMinPlot<-NULL
  yMaxPlot<-NULL
  yMinPlot<-NULL
  
  for(i in 1:length(unique(tempvector))){
    whichmainCov<-unique(tempvector)[i]
    X<-cbind(origpcoa$vectors[tempvector==whichmainCov,1],
             origpcoa$vectors[tempvector==whichmainCov,2])
    dfn <- 2
    dfd <- nrow(X) - 1
    radius<-sqrt(dfn * qf(0.95, dfn, dfd ))
    temp<- ellipseYS(center = colMeans(X), shape = cov(X),radius =radius)
    xMaxPlot<-max(c(xMaxPlot,temp[,1]))
    xMinPlot<-min(c(xMinPlot,temp[,1]))
    yMaxPlot<-max(c(yMaxPlot,temp[,2]))
    yMinPlot<-min(c(yMinPlot,temp[,2]))
  }
  plot(origpcoa$vectors[,1],origpcoa$vectors[,2],
       xlim=c(xMinPlot,xMaxPlot),ylim = c(yMinPlot,yMaxPlot),
       main="Original PCoA colored \n by the main covariate",
       pch=19,col=color1[tempvector],cex=2,
       xlab=paste("1st Coordinate",round(origpcoa$values$Relative_eig[1]*100,2),"%"),
       ylab=paste("2nd Coordinate",round(origpcoa$values$Relative_eig[2]*100,2),"%"))

    dataEllipse(origpcoa$vectors[,1],origpcoa$vectors[,2],
                data[,as.character(substitute(maincov))],lwd=3,
                levels = 0.95,add=TRUE,plot.points = FALSE,
                group.labels = NULL,
                ellipse.label=NULL,center.pch = NULL,
                col=color1)
}else{
  plot(origpcoa$vectors[,1],origpcoa$vectors[,2],
       main="Original PCoA colored \n by the main covariate",
       pch=19,col=color1[tempvector],cex=2,
       xlab=paste("1st Coordinate",round(origpcoa$values$Relative_eig[1]*100,2),"%"),
       ylab=paste("2nd Coordinate",round(origpcoa$values$Relative_eig[2]*100,2),"%"))
}
  if(drawCenter){
    for(i in 1:length(centernames)){
      for(j in rownames(data)[tempvector==centernames[i]]){
        segments(origpcoa$vectors[centers[i],1],origpcoa$vectors[centers[i],2],
                 origpcoa$vectors[j,1],origpcoa$vectors[j,2],col=color1[i],lwd=3)
      } 
    }
  }
  legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,
         unique(tempvector),pch=19,col=color1)
  

  if(drawEllipse){
    xMaxPlot<-NULL
    xMinPlot<-NULL
    yMaxPlot<-NULL
    yMinPlot<-NULL
    for(i in 1:length(unique(tempvector))){
      whichmainCov<-unique(tempvector)[i]
      datamatrix<-cbind(eigenE[,1]*eigenvalue[1]^(1/2),eigenE[,2]*eigenvalue[2]^(1/2))
      X<-cbind(datamatrix[tempvector==whichmainCov,1],
               datamatrix[tempvector==whichmainCov,2])
      dfn <- 2
      dfd <- nrow(X) - 1
      radius<-sqrt(dfn * qf(0.95, dfn, dfd ))
      temp<- ellipseYS(center = colMeans(X), shape = cov(X),radius =radius)
      xMaxPlot<-max(c(xMaxPlot,temp[,1]))
      xMinPlot<-min(c(xMinPlot,temp[,1]))
      yMaxPlot<-max(c(yMaxPlot,temp[,2]))
      yMinPlot<-min(c(yMinPlot,temp[,2]))
    }
    plot(eigenE[,1]*eigenvalue[1]^(1/2),eigenE[,2]*eigenvalue[2]^(1/2),
         xlim=c(xMinPlot,xMaxPlot),ylim = c(yMinPlot,yMaxPlot),
         main="Covariate Adjusted PCoA colored \n by the main covariate",
         pch=19,col=color1[tempvector],cex=2,
         xlab=paste("1st Coordinate",round(eigenvalue[1]/sum(eigenvalue)*100,2),"%"),
         ylab=paste("2nd Coordinate",round(eigenvalue[2]/sum(eigenvalue)*100,2),"%"))
    dataEllipse(eigenE[,1]*eigenvalue[1]^(1/2),eigenE[,2]*eigenvalue[2]^(1/2),
                factor(tempvector),lwd=3,
                levels = 0.95,add=TRUE,plot.points = FALSE,group.labels = NULL,
                ellipse.label=NULL,center.pch = NULL,
                col=color1)
  }else{
    plot(eigenE[,1]*eigenvalue[1]^(1/2),eigenE[,2]*eigenvalue[2]^(1/2),
         main="Covariate Adjusted PCoA colored \n by the main covariate",
         pch=19,col=color1[tempvector],cex=2,
         xlab=paste("1st Coordinate",round(eigenvalue[1]/sum(eigenvalue)*100,2),"%"),
         ylab=paste("2nd Coordinate",round(eigenvalue[2]/sum(eigenvalue)*100,2),"%"))
    }
  if(drawCenter){
    for(i in 1:length(centernames)){
      for(j in rownames(data)[tempvector==centernames[i]]){
        segments(eigenE[newcenters[i],1]*eigenvalue[1]^(1/2),eigenE[newcenters[i],2]*eigenvalue[2]^(1/2),
                 eigenE[j,1]*eigenvalue[1]^(1/2),eigenE[j,2]*eigenvalue[2]^(1/2),
                 col=color1[i],lwd=3)
      } 
    }
  }
  plotMatrix<-eigenE*matrix(rep(eigenvalue^(1/2),each=nrow(eigenE)),nrow=nrow(eigenE))
  plotMatrix<-plotMatrix[,!is.na(apply(plotMatrix,2,sum))]
  result<-list(plotMatrix=plotMatrix)
  result
}

