#' Make Exponential Kernel
#' 
#' Constructor function for the exponential kernel
#' @param sigma The sigma parameter
#' @return The specified kernel function
#' @export
makeExponentialKernel = function (sigma) function(a,b) exp(-(sqrt(sum((a-b)^2)))/(2*sigma^2))
#' Make Laplacian Kernel
#' 
#' Constructor function for the Laplacian kernel
#' @param sigma The sigma parameter
#' @return The specified kernel function
#' @export
makeLaplacianKernel = function (sigma) function(a,b) exp(-(sqrt(sum((a-b)^2)))/sigma)

#' Make Kernel Package for GPC
#' 
#' Make Kernel Package for GPC
#' @param type Kernel type: "Exponential", "Laplacian"
#' @return Kernel package for gpc function
#' @export 
makeKernelPacket = function (type,params,X) {
  if (type=="Exponential") {
  }
  else if (type=="Laplacian") {
    sigma=1
    sDenom=1
    if (!is.null(params$sigma)) sigma=params$sigma 
    if (!is.null(params$sDenom)) sDenom=params$sDenom
    out=list(
      update=updateLaplacianKernel,
      create=createLapacianKernelFromPacket,
      params=list(sigma=sigma),
      hyperparams=list(sDenom=sDenom)
    )
    out$params$K=createKernelMatrix(X,out$create(out))
    return(out)    
  } 
}

#' Perform a Metropolis step on the Laplacian kernel
#' 
#' Perform a Metropolis step on the Laplacian kernel
#' @param kPacket
#' @export
updateLaplacianKernel = function (kPacket,X,Y) {
  newSigma = kPacket$params$sigma + rnorm(1)/kPacket$hyperparams$sDenom
  newk = makeLaplacianKernel(newSigma)
  newK = createKernelMatrix(X,newk)
  pNew=mvtnorm::dmvnorm(Y,rep(0,length(Y)),newK,log=T) + dlnorm(newSigma,log=T)
  pOld=mvtnorm::dmvnorm(Y,rep(0,length(Y)),kPacket$params$K,log=T) + dlnorm(kPacket$params$sigma,log=T)
  if (log(runif(1))<pNew-pOld) {
    kPacket$params$sigma=newSigma
    kPacket$params$K=newK
  }
  return(kPacket$params)
}
createLapacianKernelFromPacket = function (kPacket) 
  makeLaplacianKernel(kPacket$params$sigma)

initializeLaplacianKernel = function (kPacket,X) {
  k=makeLaplacianKernel(kPacket$params$sigma)
  kPacket$params$K=createKernelMatrix(X,k)
  return (kPacket$params)
}
