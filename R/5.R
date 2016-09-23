three.RC.1<-function(X,Y,kapa){

  n=length(X[2,])

  if(length(kappa)==1){
    kappa[2]=kappa[1]
  }

  A=X[seq(along=Y)[Y==1],] # selecting class +1
  B=X[seq(along=Y)[Y==-1],] #selecting class -1
  
  ## %%% Medidas estadisticas
  mu=rbind(colMeans(A),colMeans(B))
  Ylab=c(1,-1)
  
  sigma=list(cov(A)+1e-7*diag(x=1,nrow = n , ncol = n),
             cov(B)+1e-7*diag(x=1,nrow = n , ncol = n))
  
  Mchol=list(t(chol(sigma[[1]])),t(chol(sigma[[2]])))

    ##  linear coeficient
  bb=-cbind(c(-1,numeric(n+1))) 
  
  ##  Linear constraint
  Dt=rbind(diag(x=1,nrow=n+1,ncol=n+2),cbind(0, -diag(x=1,nrow=n,ncol=n+1)))
  f=cbind(c(0,rep(x=1, 2*n)))
  
  ## %% Building the 2nd and 3rd constraint
  At2=matrix(0,nrow=n+1,ncol = n+2)
  c2=matrix(0,nrow=n+1,ncol = 1)
  At3=matrix(0,nrow=n+1,ncol = n+2)
  c3=matrix(0,nrow=n+1,ncol = 1)
  
  At2[1,1]=-1
  At2[1,1:n+1]=mu[1,]
  At2[1,n+2]=1
  At2[1:n+1,1:n+1]=kappa[1]*t(Mchol[[1]])
  
  At3[1,2]=-1
  At3[1,1:n+1]=-mu[2,]
  At3[1,n+1]=-1
  At3[1:n+1,1:n+1]=kappa[2]*t(Mchol[[2]])

  ## Solve the SOC-problem with SCS
  At=-rbind(Dt,At2,At3) #
  ct=rbind(f,c2,c3) ## f,c1,c2
  K.q=c(n+1,n+1) #dimension of cone
  cone=list( q = K.q , l=(2*n+1))

  scs=scs(At,ct,-bb,cone)

  w=cbind(scs$x[1:n+1])
  b=scs$x[n+2]
  return(list(w=w,b=b))


}