#' Demo
#' 
#' Demo example of 2d and 3d Gaussian processes. Note that hyper-parameters
#' are not fit - this is a visual taster.
#' @param dim Dimensionality (2 or 3)
#' @param samples Samples to perform (for 2d)
#' @param n Number of training points. Default set in function depending on dimensionality.
#' @param seed Random seed to use. Default set in function depending on dimensionality.
#' @param laplace_param Laplacian kernel parameter. Default set in function depending on dimensionality.
#' @param lambda Lambda parameter. Default set in function depending on dimensionality.
#' @param steps Steps for plotting. Default set in function depending on dimensionality.
#' @param xLim X limits for plotting. Default set in function depending on dimensionality.
#' @param yLim Y limits for plotting. Default set in function depending on dimensionality.
#' @param zLim Z limits for 3d plotting. Default set in function depending on dimensionality.
#' @export
demo = function (
  dim=2,
  samples=2,
  n=NA,
  seed=NA,
  laplace_param=NA,
  lambda=NA,
  steps=NA,
  xLim=NA,
  yLim=NA,
  zLim=NA
) {
  # Deal with defaults for 2d and 3d
  if (dim==2) {
    if (is.na(n))
      n=39
    if (is.na(seed))
      seed=3
    if (is.na(laplace_param))
      laplace_param=.2
    if (is.na(lambda))
      lambda=.15
    if (is.na(steps))
      steps=1500
    if (length(xLim)!=2 || is.na(xLim))
      xLim=c(0,1)
    if (length(yLim)!=2 || is.na(yLim))
      yLim=c(-6,6)
  }
  else if (dim==3) {
    if (is.na(n))
      n=100
    if (is.na(seed))
      seed=110
    if (is.na(laplace_param))
      laplace_param=.8
    if (is.na(lambda))
      lambda=.5
    if (is.na(steps))
      steps=100
    if (length(xLim)!=2 || is.na(xLim))
      xLim=c(0,1)
    if (length(yLim)!=2 || is.na(yLim))
      yLim=c(-3,3)
    if (length(zLim)!=2 || is.na(zLim))
      zLim=c(-10,10)
  }
  data=getTestData(n,dim,seed)
  model=gp(data,makeLaplacianKernel(laplace_param),lambda=lambda)
  if (dim==2)
    plot(model,steps=steps,xLim=xLim,yLim=yLim,xlab="Time",ylab="Target",samples=samples)
  else if (dim==3)
    plot(model,steps=steps,xlim=xLim,ylim=yLim,zlim=zLim,xlab="Time",ylab="Exogenous",zlab="Target",samples=samples)
  return (model)
}
getTestData = function (n=10,dim=2,seed=NA,sd=1) {
  if (!is.na(seed))
    set.seed(seed)
  x=getTestX(n,dim-1)
  coef=runif(dim-1)
  y=x%*%coef + rnorm(n,0,sd)
  cbind(x,y)
}
getTestX=function(n,dim,times=NA,do_ts=FALSE) {
  if (length(times)==1 && is.na(times))
    times=runif(n)
  if (dim>1) {
    # Get Brownian motion for other features on times
    if (do_ts)
      x_=sapply(1:(dim-1),function(i)ts(times))
    else {
      coef=runif(dim-1)
      x_=sapply(coef,function(i)times*i + rnorm(n))
    }
    x=cbind(times,x_)
  }
  else {
    x=matrix(times,ncol=1)
  }
  return (x)
}
ts=function(times,t1=1,sd=1) {
  r=runif(2,-1,1)
  y0=r[1]
  y1=r[2]
  baseGP=gp(matrix(c(0,t1,y0,y1),ncol=2))
  sample(baseGP,matrix(times,ncol=1))
}


# Unfinished demos

nepal=function() {
  center = c(28,85);
  zoom = 8
  myMap <- RgoogleMaps::GetMap(center=center, zoom=zoom,destfile = "gmap.png");
  #  RgoogleMaps::PlotOnStaticMap(myMap)  
}
uppsala = function () {
  center = c(59.85,17.66);
  zoom = 12
  myMap <- RgoogleMaps::GetMap(center=center, zoom=zoom,destfile = "uppsala.png");
  RgoogleMaps::PlotOnStaticMap(myMap)  
}
uppRiot=function(){
  set.seed(10)
  centerRiot=mvtnorm::rmvnorm(30,c(-40,30),diag(160,2))
  valRiot=mvtnorm::rmvnorm(5,c(-120,-170),diag(160,2))
  nonEvents=cbind(runif(100,-190,180),runif(100,-220,200))
  uppsala()
  points(centerRiot,col="red")
  points(valRiot,col="red")
  points(nonEvents,col="blue")
}
uppData=function(includeVal=F){
  set.seed(10)
  centerRiot=mvtnorm::rmvnorm(30,c(0,0),diag(20,2))
  valRiot=mvtnorm::rmvnorm(5,c(-35,-45),diag(40,2))
  nonEvents=cbind(runif(100,-100,80),runif(100,-80,80))
  nonEvents=nonEvents[-which(nonEvents[,1]<(-50)&(nonEvents[,2]>25|nonEvents[,2]<(-20))),]
  nonEvents=nonEvents[-which(nonEvents[,1]>50),]
  nonEvents=nonEvents[-which(nonEvents[,1]>25&nonEvents[,2]>25),]
  if (includeVal)
    rbind(
      cbind(centerRiot,1),
      cbind(valRiot,1),
      cbind(nonEvents,0)
    )
  else 
    rbind(
      cbind(centerRiot,1),
      cbind(nonEvents,0)
    )
}
upp3d=function(includeVal=F) {
  d=uppData(includeVal)
  rgl::open3d()
  rgl::plot3d(cbind(d[,1:2],0),col=d[,3]+1,zlim=c(0,1.5),xlim=c(-100,100),ylim=c(-100,100))
  rgl::show2d(filename="uppsala.png")
}
map3d <- function(map,...){
  if(length(map$tiles)!=1){
    stop("multiple tiles not implemented")
  }
  
  nx = map$tiles[[1]]$xres
  ny = map$tiles[[1]]$yres
  
  xmin = map$tiles[[1]]$bbox$p1[1]
  xmax = map$tiles[[1]]$bbox$p2[1]
  ymin = map$tiles[[1]]$bbox$p1[2]
  ymax = map$tiles[[1]]$bbox$p2[2]
  
  xc = seq(xmin,xmax,len=ny)
  yc = seq(ymin,ymax,len=nx)
  colours = matrix(map$tiles[[1]]$colorData,ny,nx)
  m = matrix(0,ny,nx)
  surface3d(xc,yc,m,col=colours,...)
  
}












