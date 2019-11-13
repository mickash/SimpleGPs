#' Fit a single hyper-parameter kernel's parameter and lambda parameter
#' 
#' Fit a single hyper-parameter kernel's parameter and lambda parameter
#' @param X Input data, or input and target variables, with Y the last column
#' @param lambda A vector of lambda parameters to search over
#' @param sigma A vector of sigma kernel parameters to search over
#' @param folds The number of folds to perform. If 0, all but one validation is performed.
#' @param processMean The process mean of the models
#' @param scl Whether the data should be scaled and centered
#' @param Y The target variables, if not included in X. If NULL it is assumed they are included in X.
#' @return The cross validation likelihood grid for the parameter combinations.
#' @export
fit1 = function (
  X,
  lambda,
  sigma,
  kernelFunc=makeLaplacianKernel,
  folds=0,
  processMean=0,
  scl=T,
  Y=NULL,
  display=F
) {
  # Grid search over sigma (kernel parameter) and lambda.
  outer(lambda,sigma,function(l,s){
    if (folds==0)
      apply(cbind(l,s),1,function(p) abo.gp(X,kernelFunc(p[2]),p[1],processMean,scl,Y,display))
    else
      apply(cbind(l,s),1,function(p) kfolds.gp(X,folds,makeLaplacianKernel(p[2]),p[1],processMean,scl,Y))
  })
}

#' Perform a grid search to fit a single hyper-parameter kernel's parameter and lambda parameter
#' 
#' Perform a grid search to fit a single hyper-parameter kernel's parameter and lambda parameter
#' @param X Input data, or input and target variables, with Y the last column
#' @param maxIter The maximum number of iterations to continue the search for
#' @param sensitivity Required convergence difference to terminate on
#' @param folds The number of folds to perform. If 0, all but one validation is performed.
#' @param processMean The process mean of the models
#' @param scl Whether the data should be scaled and centered
#' @param Y The target variables, if not included in X. If NULL it is assumed they are included in X.
#' @return The model built from the best found hyper-parameters
#' @export
gridSearch1 = function (
    X,
    kernelFunc=makeLaplacianKernel,
    maxIter=100,
    sensitivity=.001,
    folds=0,
    processMean=0,
    scl=T,
    Y=NULL,
    display=F,
    detailed=F,
    sigmaMax=1000,
    lambdaMax=1000,
    sigmaMin=.0001,
    lambdaMin=.0001
) {
  makeSigmaSteps=function(center,step) {
      out=exp(c(center-2*step,center-step,center,center+step,center+2*step))
      out[which(out>sigmaMax)]=sigmaMax
      out[which(out<sigmaMin)]=sigmaMin
      return (out)
  }
  makeLambdaSteps=function(center,step) {
    out=exp(c(center-2*step,center-step,center,center+step,center+2*step))
      out[which(out>lambdaMax)]=lambdaMax      
      out[which(out<lambdaMin)]=lambdaMin      
      return (out)
  }
  
  stepLambda=1
  centerLambda=0
  stepSigma=1
  centerSigma=0
  found=NULL
  for (i in 1:maxIter) {
    lambda=makeLambdaSteps(centerLambda,stepLambda)
    sigma=makeSigmaSteps(centerSigma,stepSigma)
    # Grid is lambda X sigma, low to high for both
    grid=fit1(X,lambda,sigma,kernelFunc,folds,processMean,scl,Y,display) 
    if (!any(!is.infinite(grid))) {
      if(display) print(paste("Step: ",i,".Only infinities found. Zooming out.",sep=""))
      stepLambda=stepLambda*2
      stepSigma=stepSigma*2
      next
    }
    
    best=which(grid==max(grid),arr.ind=T)
    scr=grid[best[1],best[2]]
    gridCpy=grid
    gridCpy[best[1],best[2]]=-Inf
    secondBest=which(gridCpy==max(gridCpy),arr.ind=T)[1,] # If multiple equal
    secondScr=gridCpy[secondBest[1],secondBest[2]]
    if (scr-secondScr < sensitivity && stepLambda<.1 && stepSigma<.1) {
      message("Search converged.")
      if (display) print(paste("Step ",i,"/",maxIter,": ",scr))
      found=c(lambda[best[1]],sigma[best[2]])
      break
    }
    else if (i==maxIter) {
      warning("Covergence failed. Exiting after max iterations.")
      found=c(lambda[best[1]],sigma[best[2]])
      break
    }
    else if (display) {
      print(paste("Step ",i,"/",maxIter,": ",scr))
      print(paste("Sigma:",sigma[best[2]]))
      print(paste("Lambda:",lambda[best[1]]))
      print(paste("Difference:",scr-secondScr))
      if (detailed) {
        print("Sigmas:")
        print(sigma)
        print("Lambdas:")
        print(lambda)
        print("Grid:")
        print(grid)
      }
    }
    if (best[1]==1) {
      centerLambda=centerLambda-2*stepLambda
      stepLambda=stepLambda*1.5
    }
    if (best[1]==2) {
      centerLambda=centerLambda-stepLambda
      stepLambda=stepLambda*.75
    }
    if (best[1]==3) {
      stepLambda=stepLambda*.75
    }
    if (best[1]==4) {
      centerLambda=centerLambda+stepLambda
      stepLambda=stepLambda*.9
    }
    if (best[1]==5) {
      centerLambda=centerLambda+2*stepLambda
      stepLambda=stepLambda*1.1
    }
    if (best[2]==1) {
      centerSigma=centerSigma-2*stepSigma
      stepSigma=stepSigma*1.1
    }
    if (best[2]==2) {
      centerSigma=centerSigma-stepSigma
      stepSigma=stepSigma*.9
    }
    if (best[2]==3) {
      stepSigma=stepSigma*.9
    }
    if (best[2]==4) {
      centerSigma=centerSigma+stepSigma
      stepSigma=stepSigma*.9
    }
    if (best[2]==5) {
      centerSigma=centerSigma+2*stepSigma
      stepSigma=stepSigma*1.1
    }    
  }
  message(paste("Sigma:",found[2]))
  message(paste("Lambda:",found[1]))
  model=gp(X,kernelFunc(found[2]),found[1],processMean,scl,Y)
}
#' Find the log probability of the data using ABO CV
#' 
#' Find the log probability of the data using ABO CV
#' @param X Input data, or input and target variables, with Y the last column
#' @param k The kernel function
#' @param lambda The lambda parameter
#' @param processMean The process mean of the models
#' @param scl Whether the data should be scaled and centered
#' @param Y The target variables, if not included in X. If NULL it is assumed they are included in X.
#' @return The sum of the log probability of each data point using ABO CV
#' @keywords internal
abo.gp=function(
  X,
  k=makeLaplacianKernel(1),
  lambda=0,
  processMean=0,
  scl=T,
  Y=NULL,
  display=F
) {
  if (is.null(Y)) {
    Y=X[,ncol(X)]
    X=X[,-ncol(X),drop=F]
  }
  p=sapply(1:nrow(X), function(i) {
    model=gp(X[-i,],k,lambda,processMean,scl,Y[-i])
    p=predict(model,X[i,])
    dnorm(Y[i],p[1,1],p[1,2])
  })
#  if (display) print(p)
  sum(log(p))
}

#' Find the probability of points
#' 
#' Find the probability of points
#' @param model The gp model
#' @param X The input data
#' @param Y The target data
#' @param logprob Whether the probabilities should be log
#' @return A scalar giving the probability or log probability
#' @keywords internal
probability=function (model,X,Y,logprob=T) {
  p=predict(model,X)
  probs=sapply(1:length(Y),function(i) dnorm(Y[i],p[i,1],p[i,2]))  
  if (logprob)
    sum(log(probs))  
  else
    prod(probs)
}
#' Perform K-folds fitting on GP
#' 
#' Perform K-folds fitting on GP
#' @param X Input data, or input and target variables, with Y the last column
#' @param folds The number of folds to use
#' @param k The kernel function
#' @param lambda The lambda parameter
#' @param processMean The process mean of the models
#' @param scl Whether the data should be scaled and centered
#' @param Y The target variables, if not included in X. If NULL it is assumed they are included in X.
#' @return The sum of the log probability of each data point using K-folds CV
#' @keywords internal
kfolds.gp = function (
  X,
  folds=10,
  k=makeLaplacianKernel(1),
  lambda=0,
  processMean=0,
  scl=T,
  Y=NULL 
) {
  if (is.null(Y)) {
    Y=X[,ncol(X)]
    X=X[,-ncol(X),drop=F]
  }
  sum(sapply(1:folds,function(i){
    holdout=(((i-1)/folds)*nrow(X)):((i/folds)*nrow(X))
    model=gp(X[-holdout,,drop=F],k,lambda,processMean,scl,Y[-holdout])
    probability(model,X[holdout,,drop=F],Y[holdout])
  }))
}