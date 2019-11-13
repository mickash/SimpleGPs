#' Create a Gaussian Process
#' 
#' Create a Gaussian Process
#' @param X The X data as a matrix, or able to be case to a matrix
#' @param k The kernel function
#' @param lambda The lambda parameter
#' @param processMean The mean of the process. If NA the mean of the Y data is used.
#' @param scl Whether the data should be scaled and centered
#' @param Y The Y data as a vector. If NULL, the last column of X is used as Y.
#' @return A gp object
#' @export
gp<-function(
  X,
  k=makeLaplacianKernel(1),
  lambda=0,
  processMean=0,
  scl=T,
  Y=NULL  
){
  # Deal with X and Y
  if (is.data.frame(X)) X=as.matrix(X)
  if (is.null(Y)) {
    Y=X[,ncol(X)]
    X=X[,-ncol(X),drop=FALSE]
  }
  
  # Scale Data
  sX=NULL
  sY=NULL
  if (scl) {
    sX=scale(X)
    sY=scale(Y)      
  }
  else {
    sX=X
    sY=Y
  }
  # Create kernel matrix
  K=createKernelMatrix(as.matrix(sX),k)
  # Add noise
  K=K+diag(lambda,nrow(K))
  # Calculate Mean
  M=NA
  if (is.na(processMean)) {
    if (scl)
      M=(mean(as.vector(Y))-attr(model$Y,"scaled:center"))/attr(model$Y,"scaled:scale")  
    else
      M=mean(as.vector(Y))
  } else {
    if (scl)
      M=(processMean-attr(sY,"scaled:center"))/attr(sY,"scaled:scale")
    else
      M=processMean
  }
  # Return GP
  out=list(
    M=M,
    K=K,
    K_=solve(K),
    rX=X,
    rY=Y,
    X=sX,
    Y=sY,
    k=k,
    scaled=scl,
    lambda=lambda,
    dim=ncol(sX)+1
  )
  class(out)="gp"
  return(out)
}

#' Create a kernel matrix
#' 
#' Create a kernel matrix from k(x,x)
#' @param x The data points to use
#' @param k The kernel function to use
#' @return The kernel matrix
#' @keywords internal
createKernelMatrix<-function(x,k){
  if (is.null(nrow(x))) x=matrix(x,ncol=1)
  as.matrix(sapply(1:nrow(x),function(i)sapply(1:nrow(x),function(j)k(x[i,],x[j,]))))
}
#' Create a kernel matrix 2
#' 
#' Create a kernel matrix from k(x,X)
#' @param x The first set of data points to use
#' @param X The second set of data points to use
#' @param k The kernel function to use
#' @return The kernel matrix
#' @keywords internal
createKernelMatrix2<-function(x,X,k){
  if (is.null(nrow(x))) x=matrix(x,ncol=1)
  if (is.null(nrow(X))) X=matrix(X,ncol=1)
  t(sapply(1:nrow(x),function(i)sapply(1:nrow(X),function(j)k(x[i,],X[j,]))))
} 
