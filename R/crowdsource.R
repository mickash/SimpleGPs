csProj=function(dim,mean,sigma,lambda) {
  dens=makeDensity(dim,mean,sigma)
  dens=crowdsource(dens)
  dens=crowdModel(dens,lambda)
  return(dens)
}
gen0grid=function(lim){
  x1=c()
  x2=c()
  y=c()
  for (i in 0:lim) {
    for (j in 0:lim) {
      if (i==0 || i==lim || j==0 || j==lim) {
        x1=c(x1,i)
        x2=c(x2,j)
        y=c(y,0)
      }
    }
  }  
  data.frame(x1=x1,x2=x2,y=y)
}
genNAgrid=function(lim){
  x1=c()
  x2=c()
  y=c()
  for (i in 0:lim) {
    for (j in 0:lim) {
      if (i!=0 && i!=lim && j!=0 && j!=lim) {
        x1=c(x1,i)
        x2=c(x2,j)
        y=c(y,NA)
      }
    }
  }  
  data.frame(x1=x1,x2=x2,y=y)
}

makeDensity=function(lim=20,mean=3,sigma=5,neg=F) {
  border=gen0grid(lim)
  df=genNAgrid(lim)
  model=gp(border,k=makeLaplacianKernel(sigma),scl=F,processMean=mean)
  df$y=sample(model,df[-3])[1,]
  if (!neg)
    df$y[which(df$y<0)]=0
  df=rbind(border,df)
  df=df[order(df$x1,df$x2),]
  rgl::plot3d(df)
  rgl::surface3d(0:lim,0:lim,matrix(df$y,nrow=21,byrow=T),col="red",alpha=.6)
  list(model=model,df=df,border=border,lim=lim,mean=base::mean(df$y),sigma=sigma,neg=neg)
}

crowdsource=function(dens,people=20,goes=5,sd=.5,odd=.2) {
  oddities=runif(people)<odd
  x1=base::sample(0:dens$lim,people*goes,T)
  x2=base::sample(0:dens$lim,people*goes,T)
  y=rep(NA,people*goes)
  for (i in 1:people){
    for (j in 1:goes) {
      n=(i-1)*goes+j
      row=which(dens$df$x1==x1[n]&dens$df$x2==x2[n])
      real=dens$df[row,3]
      y[n]=rnorm(1,real,sd)
      if (oddities[i])
        y[n]=y[n]/2
    }
  }
  if (!dens$neg)
    y[which(y<0)]=0
  dens$cs=data.frame(x1,x2,y)
  return (dens)
}
crowdModel=function(dens,lambda) {
  cmodel=gp(
    rbind(dens$cs,dens$border),
    k=makeLaplacianKernel(sigma),
    scl=F,
    processMean=dens$mean,
    lambda=c(rep(lambda,nrow(dens$cs)),rep(0,nrow(dens$border)))
    )
  rgl::plot3d(dens$df)
  rgl::surface3d(0:dens$lim,0:dens$lim,matrix(dens$df$y,nrow=dens$lim+1,byrow=T),col="red",alpha=.5)
  rgl::points3d(dens$cs,col="blue")
  y=predict(cmodel,dens$df[1:2])[,1]
  if (!dens$neg)
    y[which(y<0)]=0
  rgl::surface3d(0:dens$lim,0:dens$lim,matrix(y,nrow=dens$lim+1,byrow=T),col="blue",alpha=.5)
  dens$cmodel=cmodel
  dens$est=y
  return(dens)
}

reestimatePoints=function(dens,pnts,sd=.3,num=3) {
  est=apply(pnts,1,function(row){
    r=which(dens$df$x1==row[1]&dens$df$x2==row[2])
    real=dens$df[r,3]
    rnorm(num,real,sd)/num   
  })
  dens$expert=data.frame(x1=pnts[,1],x2=pnts[,2],est=est)
  return(dens)
}

# Strategy 1: Reestimate 5 points at random
# Strategy 2: Reestimate 5 most unlikely points given all
# Strategy 3: Reestimate 5 most unlikely points given others
# Strategy 4: Reestimate 5 points with largest KL-divergence from others
# Strategy 5: Reestimate 5 most anomolous points

# KL divergence of density
# mean absolute distance of grid points
# root mean squared distance of grid points