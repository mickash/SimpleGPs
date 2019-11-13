quicker=function(model,x) {
  # Make sure x is a matrix
  if (is.null(nrow(x))) x=matrix(x,nrow=1)
  a=(nrow(model$X)+nrow(x))^2
  b=nrow(x)*(nrow(model$X)+1)^2
  a<b
}

#' Predict
#' 
#' Predict
#' @param model The gp object
#' @param x The X data, as a matrix. Vectors will be interpreted as single row matrices.
#' @param full Should the full covariance matrix be used.
#' @return A matrix, with the first column giving predictions (means) and the second
#' column giving variance.
#' @export
predict.gp=function(model,x,full=quicker(model,x)){
  # Make sure x is a matrix
  if (is.null(nrow(x))) x=matrix(x,nrow=1)
  # Scale Data, and make sure it is a matrix
  sX=NULL
  if (model$scale)
    sX=sapply(1:ncol(x),function(i)(x[,i]-attr(model$X,"scaled:center")[i])/attr(model$X,"scaled:scale")[i])
  else
    sX=x
  if (nrow(x)==1) sX=matrix(sX,nrow=1)
  
  if (full) {
    # Create covariance matrix components
    k_xX=createKernelMatrix2(sX,model$X,model$k)
    k_xx=createKernelMatrix(sX,model$k)+diag(model$lambda,nrow(sX))
    # K-Matrix & Dimensions:
    # KXX  KXx    n*n   n*m 
    # KxX  Kxx    m*n   m*m
    #      m*m  m*n     n*n       n*m 
    cond_K=diag(k_xx-k_xX%*%model$K_%*%t(k_xX))
    cond_M=model$M+k_xX%*%model$K_%*%(model$Y-model$M)  
    # Correct for any negative variances due to numerical instability
    cond_K[which(cond_K<0)]=0
    # Unscale output
    if (model$scale)
      cbind(cond_M*attr(model$Y,"scaled:scale")+attr(model$Y,"scaled:center"),
            cond_K*attr(model$Y,"scaled:scale"))
    else
      cbind(cond_M,cond_K)    
  }
  else {
    out=sapply(1:nrow(sX),function(i) {
      # Create covariance matrix components
      k_xX=createKernelMatrix2(sX[i,,drop=F],model$X,model$k)
      k_xx=createKernelMatrix(sX[i,,drop=F],model$k)+model$lambda
      # K-Matrix & Dimensions:
      # KXX  KXx    n*n   n*m 
      # KxX  Kxx    m*n   m*m
      #      m*m  m*n     n*n       n*m 
      cond_K=max(0,k_xx-k_xX%*%model$K_%*%t(k_xX))
      cond_M=model$M+k_xX%*%model$K_%*%(model$Y-model$M)  
      c(cond_M,cond_K)
    })
    out[1,]=out[1,]*attr(model$Y,"scaled:scale")+attr(model$Y,"scaled:center")
    out[2,]=out[2,]*attr(model$Y,"scaled:scale")
    t(out)
  }  
}


#' Sample
#' 
#' Sample
#' @param model The gp object
#' @param x The X data, as a matrix. Vectors will be interpreted as single row matrices.
#' This gives the points to sample Y values at.
#' @return A matrix, giving Y values.
#' @export
#' @export
sample=function(model,X) {
  # Make sure X is a matrix
  if (is.null(nrow(X))) 
    X=matrix(X,nrow=1)
  # Scale Data, and make sure it is a matrix
  sX=NULL
  if(model$scale)
    sX=sapply(1:ncol(X),function(i)(X[,i]-attr(model$X,"scaled:center")[i])/attr(model$X,"scaled:scale")[i])
  else
    sX=X
  if (nrow(X)==1) 
    sX=matrix(sX,nrow=1)
  
  # Create covariance matrix components
  k_xX=createKernelMatrix2(sX,model$X,model$k)
  k_xx=createKernelMatrix(sX,model$k)+diag(model$lambda,nrow(X))
  cond_K=k_xx-k_xX%*%model$K_%*%t(k_xX)
  cond_M=model$M+k_xX%*%model$K_%*%(model$Y-model$M) 
  if (model$scale)
    mvtnorm::rmvnorm(1,cond_M,cond_K,"chol")[1,]*attr(model$Y,"scaled:scale")+attr(model$Y,"scaled:center")
  else
    mvtnorm::rmvnorm(1,cond_M,cond_K,"chol")[1,]
  # The following could be faster in certain circumstances...
  #   Create output vector
  #   out=rep(NA,nrow(X))
  #   for (i in 1:nrow(sX)) {
  #     # Create covariance matrix components
  #     k_xX=createKernelMatrix2(sX[i,,drop=F],model$X,model$k)
  #     k_xx=createKernelMatrix(sX[i,,drop=F],model$k)+model$lambda
  #     # K-Matrix & Dimensions:
  #     # KXX  KXx    n*n   n*m 
  #     # KxX  Kxx    m*n   m*m
  #     #      m*m  m*n     n*n       n*m 
  #     cond_K=k_xx-k_xX%*%model$K_%*%t(k_xX)
  #     cond_M=model$M+k_xX%*%model$K_%*%(model$Y-model$M)  
  #     out[i]=rnorm(1,cond_M,cond_K)*scaleY+centerY
  #     model$X=rbind(model$X,sX[i,])
  #     model$Y=rbind(model$Y,out[i])
  #     model$K=cbind(model$K,t(k_xX))
  #     model$K=rbind(model$K,cbind(k_xX,k_xx))
  #     model$K_=solve(model$K)
  #   }
  #   return(out)
}

#' Prediction function for classification model
#' 
#' S3 override prediction function for classification model
#' @param model gpc model
#' @param X Input data
#' @param ... Other arguments, which are all ignored.
#' @return A three column matrix consisting of the probabilitys predicted 
#' @export
predict.gpc = function (model,X,...) {
  f=predict(model$gpModel,X)
  f_=1/(1+exp(-f[,1]))
  cbind(f_,as.numeric(f_>.5),f)
}