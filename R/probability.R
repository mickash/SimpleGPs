#' Probability of Y data using an empty GP
#' 
#' Probability of Y data using an empty GP
#' @param X The X data
#' @param k The kernal function
#' @param lambda The lambda parameter
#' @param processMean The process mean
#' @param scl Whether the data should be scaled and centered
#' @param Y The Y data as a vector. If NULL, the last column of X is used as Y.
#' @return A scalar giving the probability of the data
#' @keywords internal
emptyProbability = function (
  X,
  k=makeLaplacianKernel(1),
  lambda=0,
  processMean=0,
  scl=T,
  Y=NULL  
  ) {
  # Deal with X and Y
  if (is.data.frame(X)) X=as.matrix(X)
  if (is.null(Y)) {
    Y=X[,ncol(X)]
    X=X[,-ncol(X)]
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
  
  # Calculate probability of Y data
  mvtnorm::dmvnorm(Y,rep(0,length(Y)),K)  
}
