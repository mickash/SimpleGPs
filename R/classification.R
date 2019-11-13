#' Produce a logistic GP classifier
#' 
#' Produce a logistic GP classifier.
#' Uses a Monte Carlo (Metropolis Hastings) Engine.
#' @param X The X data as a matrix, or able to be case to a matrix
#' @param Y The Y data as a vector. If NULL, the last column of X is used as Y.
#' @param samples The number of samples to use
#' @param burn The burn period
#' @param method The MCMC method to use ("M" is Metropolis; "MG" is Metrolopolis in Gibbs)
#' @param denom Denominator for the standard Gaussian random walk in the Metropolis steps
#' @param fmean The mean of the process for the latent variable. If NA the mean of the 
#' estimate points is used.
#' @param progress Whether progress should be printed
#' @return A gpc object
#' @export
gpc = function (
  X,
  Y=NULL,
  samples=1000,
  burn=samples/10,
  method="MG",
  kernel="Laplacian",
  kernelParams=NULL,
  fDenom=ifelse(method=="MG",10,length(Y)+1),
  sDenom=10,
  fmean=0, # -5 good for low prob, -10 for tiny prob mean
  progress=F,
  createKernel=makeLaplacianKernel
  ) {
  # Deal with X and Y
  if (is.data.frame(X)) X=as.matrix(X)
  if (is.null(dim(X))) X=matrix(X,ncol=1)
  if (is.null(Y)) {
    Y=X[,ncol(X)]
    X=X[,-ncol(X),drop=FALSE]
  }

  
  # Make latent variable params
  params=list()
  params$f=rep(0,length(Y))
  params$unknownIndices=which(is.na(Y))
  params$fDenom=fDenom
  params$mult=1
  
  # Make kernel packet
  params$kPacket=makeKernelPacket(kernel,kernelParams,X)
  
  # Sample from:
  #   p(Y|f)p(f|X,B)p(B)
  #   B = lambda,sigma: log-normal, uniform?
  #   f|X,B = g(X,B)
  #   Y|f = log(f)
  #   
  # MH Step: Adjust lambda, sigma, f values
  
  for (i in 1:burn) {
    params=evaluate(params,X,Y,method)
    if (i %% 10 == 0 && progress) cat(paste("\nBurn:",i,"of",burn))    
  }
  f=matrix(rep(NA,length(Y)*samples),nrow=samples)
  sigma=rep(NA,samples)
  mult=rep(NA,samples)
  sampled=list()
  for (i in 1:samples) {
    params=evaluate(params,X,Y,method)
    f[i,]=params$f
    sigma[i]=params$kPacket$params$sigma
    mult[i]=params$mult
    if (i %% 10 == 0 && progress) cat(paste("\nSample:",i,"of",samples))    
  }

#   accepted=sapply(sampled,function(s)s$first)
#   acceptance=c(length(which(accepted)),length(which(!accepted)))
#   f=sapply(1:length(Y),function(i) sapply(sampled,function(s)s$f[i]))
#   sigma=sapply(sampled,function(s) s$kPacket$sigma)
#  gpModel=gp(X,params$kPacket$create(mean(sigma)),processMean=fmean,Y=colMeans(f))
  acceptanceRateF=sum(sapply(2:nrow(f),
      function (row) sapply(1:ncol(f),
      function(col)as.numeric(f[row,col]!=f[row-1,col]))))/length(f) 
  acceptanceRateS=sum(sapply(2:length(sigma),
      function(i)as.numeric(sigma[i]!=sigma[i-1])))/length(sigma)

  out=list (
    sampled=sampled,
    unknownIndices=params$unknownIndices,
    accF=acceptanceRateF,
    accS=acceptanceRateS,
    f=f,
    sigma=sigma,
    mult=mult,
#    gpModel=gpModel,
    X=X,
    Y=Y
    )
  class(out)="gpc"
  return (out)
}

candidate = function (params,denom) {
  params$sigma=params$sigma+rnorm(1)/denom
  params$f=params$f+mvtnorm::rmvnorm(1,rep(0,length(params$f)))/denom
  params$Y=sample(c(0,1),length(params$Y))
  params$first=T
  return (params)
}
logProb = function (params,X,Y) {
  log(dnorm(params$sigma))+
    log(emptyProbability(X,makeLaplacianKernel(params$sigma),Y=params$f))+
    probY_log(params$f,Y)
}

evaluate=function(params,X,Y,method) {
  if (method=="M") evaluate1(params,X,Y)
  else if (method=="MG") evaluate3(params,X,Y)
  else stop("Invalid method.")
}
evaluate1 = function (params,X,Y,denom) {
  params2=candidate(params,denom)
  pOld_log=logProb(params,X,Y)
  pNew_log=logProb(params2,X,Y)
  update=log(runif(1))<pNew_log-pOld_log
  if (update) return (params2)
  params$first=F
  return (params)
}
updateCheck=function(params,params2,X,Y) 
  log(runif(1))<logProb(params2,X,Y)-logProb(params,X,Y)  

logProbF = function (params,X,Y,i) {
  p=1/(1+exp(-f[i]))
  if (Y[i]==1) p=log(p)
  else p=log(1-p)
  p+log(emptyProbability(X,makeLaplacianKernel(params$sigma),Y=params$f))
}

probY_log = function (f,Y) {
  p=1/(1+exp(-f))  
  ones=which(Y==1)
  zeros=which(Y==0)
  sum(log(p[ones]))+sum(log(1-p[zeros]))
}



evaluate3 = function (
  params ,
  X ,
  Y
  ) {
  # Kernel
  params$kPacket$params=params$kPacket$update(params$kPacket,X,params$f)
  
  # Unknown Fs. f_i is conditional on X_i,sigma,f_-i 
  logProbf=mvtnorm::dmvnorm(params$f,rep(0,length(params$f)),params$kPacket$params$K,log=T)
  newF = params$f
  for ( i in 1:length(params$f) ) {
    # Peturb
    newF[i] = params$f[i] + rnorm(1)/params$fDenom

    # Alter if appropriate
    newPf_log=mvtnorm::dmvnorm(newF,rep(0,length(newF)),params$kPacket$params$K,log=T)
    p= newPf_log - logProbf 
    if (!is.na(Y[i])) {
      oldPy=1/(1+exp(-(params$mult*params$f[i])))  
      newPy=1/(1+exp(-(params$mult*newF[i])))  
      p=p+log(ifelse(Y[i]==1,newPy,1-newPy)) - log(ifelse(Y[i]==1,oldPy,1-oldPy))            
    }
    
    if (log(runif(1))<p) {
      params$f[i]=newF[i]
      logProbf=newPf_log
    } 
    else {
      newF[i]=params$f[i]
    }
  }

  newMult=params$mult+rnorm(1)
  if (newMult>1 && newMult<100) {
    pMultOld=sum(sapply(1:length(params$f),
                        function(i)ifelse(is.na(Y[i]),0,log(1/(1+exp(-(params$mult*params$f[i])))))))  
    pMultNew=sum(sapply(1:length(params$f),
                        function(i)ifelse(is.na(Y[i]),0,log(1/(1+exp(-(newMult*params$f[i])))))))  
    if (log(runif(1))<pMultNew-pMultOld)
      params$mult=newMult    
  }
  
  return (params)
}

