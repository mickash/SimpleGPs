#' Plot Gaussian Process
#' 
#' Plot Gaussian Process
#' @param model The Gaussian process to model
#' @param ci The confidence interval level, in terms of standard deviation
#' @param steps The number of points to calculate
#' @param top Whether the top confidence interval should be shown
#' @param bot Whether the bottom confidence interval should be shown
#' @param samples The number of function samples to display
#' @param ... Other arguments to plot or rgl::plot3d
#' @return Invisible NULL
#' @export
plot.gp = function (model,ci=2,steps=50,top=T,bot=T,samples=0,residuals=T,xLim=NULL,yLim=NULL,...) {
  if (model$dim == 2) {
    plot(cbind(model$rX,model$rY),xlim=xLim,ylim=yLim,...)
    xmin=min(c(model$rX,xLim))
    xmax=max(c(model$rX,xLim))
    xlen=xmax-xmin
    x=seq(xmin-.1*xlen,xmax+.1*xlen,xlen/steps)
    p=predict(model,matrix(x,ncol=1))
    lines(x,p[,1])
    if (residuals) {
      print("Plotting residuals")
      segments(model$rX,model$rY,model$rX,predict(model,matrix(model$rX,ncol=1))[,1],col="red")
    }
    if (top) lines(x,p[,1]+ci*sqrt(p[,2]),col="green")
    if (bot) lines(x,p[,1]-ci*sqrt(p[,2]),col="green")
    if (samples>0) {
      for (i in 1:samples) {
        points(x,sample(model,matrix(x,ncol=1)),type="l",col=i)
      }
    }
  }
  else if (model$dim==3) {
    rgl::plot3d(cbind(model$rX,model$rY),...)
    x1min=min(model$rX[,1])
    x1max=max(model$rX[,1])
    x1len=x1max-x1min
    x1seq=seq(x1min-.1*x1len,x1max+.1*x1len,x1len/steps)
    x2min=min(model$rX[,2])
    x2max=max(model$rX[,2])
    x2len=x2max-x2min
    x2seq=seq(x2min-.1*x2len,x2max+.1*x2len,x2len/steps)
    x=cbind(rep(x1seq,each=length(x2seq)),x2seq)
    p=predict(model,x)
    rgl::surface3d(
      matrix(x[,1],ncol=length(x2seq)),
      matrix(x[,2],ncol=length(x1seq)),
      matrix(p[,1],nrow=length(x1seq)),
      col="blue",alpha=.6)
    if (residuals) {
      segs=gdata::interleave(cbind(model$rX,model$rY),cbind(model$rX,predict(model,model$rX)[,1]))
      rgl::segments3d(segs,col="red")
    }
    if (top) 
      rgl::surface3d(
        matrix(x[,1],ncol=length(x2seq)),
        matrix(x[,2],ncol=length(x1seq)),
        matrix(p[,1]+ci*sqrt(p[,2]),nrow=length(x1seq)),
        col="green",alpha=.6)
    if (bot) 
      rgl::surface3d(
        matrix(x[,1],ncol=length(x2seq)),
        matrix(x[,2],ncol=length(x1seq)),
        matrix(p[,1]-ci*sqrt(p[,2]),nrow=length(x1seq)),
        col="green",alpha=.6)
    if (samples>0) {
      for (i in 1:samples) {
        rgl::surface3d(
          matrix(x[,1],ncol=length(x2seq)),
          matrix(x[,2],ncol=length(x1seq)),
          sample(model,x),
          col=i,alpha=.6)
      }
    }
  }
  else stop("Cannot plot Gaussian process in space with higher than three dimensions.")
  invisible(NULL)
}

#' Plot a categorical Gaussian Process
#' @param model The gpc model
#' @param steps The number of points to calculate
#' @param ci The confidence interval level, in terms of standard deviation
#' @param top Whether the top confidence interval should be shown
#' @param bot Whether the bottom confidence interval should be shown
#' @param ... Other arguments to plot or rgl::plot3d
#' @return Invisible NULL
#' @export
plot.gpc = function (model,...) {
  if (is.null(dim(model$X)) || ncol(model$X)==1) {
    plot(model$X[,1],model$Y,col=model$Y+2,...)    
    est=colMeans(1/(1+exp(-model$mult*model$f)))
    abline(h=.5)
    points(model$X[model$unknownIndices,],est[model$unknownIndices],col="blue")
  }
  else if (ncol(model$X)==2) {
    rgl::plot3d(model$X[,1],model$X[,2],model$Y,col=model$Y+2,...)    
    x1min=min(model$X[,1])
    x1max=max(model$X[,1])
    x1len=x1max-x1min
    x1seq=seq(x1min-.1*x1len,x1max+.1*x1len,x1len/steps)
    x2min=min(model$X[,2])
    x2max=max(model$X[,2])
    x2len=x2max-x2min
    x2seq=seq(x2min-.1*x2len,x2max+.1*x2len,x2len/steps)
    x=cbind(rep(x1seq,each=length(x2seq)),x2seq)
    p=predict(model,x)
    rgl::surface3d(
      matrix(x[,1],ncol=length(x2seq)),
      matrix(x[,2],ncol=length(x1seq)),
      matrix(p[,1],nrow=length(x1seq)),
      col="blue",alpha=.6)
    if (top) 
      rgl::surface3d(
        matrix(x[,1],ncol=length(x2seq)),
        matrix(x[,2],ncol=length(x1seq)),
        matrix(1/(1+exp(-(p[,2]+ci*sqrt(p[,3])))),nrow=length(x1seq)),
        col="blue",alpha=.6)
    if (bot) 
      rgl::surface3d(
        matrix(x[,1],ncol=length(x2seq)),
        matrix(x[,2],ncol=length(x1seq)),
        matrix(1/(1+exp(-(p[,2]-ci*sqrt(p[,3])))),nrow=length(x1seq)),
        col="blue",alpha=.6)    
  }
  else if (ncol(model$X)==3) {
    plot(model$X,model$Y)    
  }
  invisible(NULL) 
}