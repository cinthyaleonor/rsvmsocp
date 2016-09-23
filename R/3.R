two.RR.1<-function(X,Y,kappa, C) # Socp_ R12 svmV1_scs  ##  what's r12????
{
 
  if(length(kappa)==1){
    kappa[2]=kappa[1]
  }

 ## calc 
  A=X[seq(along=Y)[Y==1],] # selecting class +1
  B=X[seq(along=Y)[Y==-1],] #selecting class -1

  ## %%% Medidas estadisticas
  mu=rbind(colMeans(A),colMeans(B))
  n=dim(mu)[2]
  Ylab=c(1,-1)

  sigma=list(cov(A)+1e-6*diag(x=1,nrow = n , ncol = n),
             cov(B)+1e-6*diag(x=1,nrow = n , ncol = n))

  Mchol=list(t(chol(sigma[[1]])),t(chol(sigma[[2]])))
  #### 
  ## %% linear coeficient

  bb=-cbind(c(-1,-C,numeric(n+1)))

  ## %% Linear constraint
  Dt=rbind(c(1, numeric(n+2)),c(0, 1,numeric(n+1)))
  #f=cbind(c(0,0))
  #K.l=2
  
   ##%% Building the 2nd constraint
  At2=matrix(0, nrow=n+1, ncol=n+3) 
  At3=matrix(0, nrow=n+1, ncol=n+3)
  c2=cbind(numeric(n+1))
  c3=cbind(numeric(n+1))

  At2[1,1]=-1
  At2[1,1:n+2]=mu[1,]
  At2[1,n+3]=1
  At2[1:n+1,1:n+2]=kappa[1]*t(Mchol[[1]])

  At3[1,2]=-1
  At3[1,1:n+2]=-mu[2,]
  At3[1,n+3]=-1
  At3[1:n+1,1:n+2]=kappa[2]*t(Mchol[[2]])

  At=-rbind(Dt,At2,At3) #
  ct=rbind(cbind(c(0,0)),c2,c3) ## f,c1,c2
  K.q=c(n+1, n+1) #dimension of cone
  cone <- list( q = K.q , l=2)

  scs=scs(At,ct,-bb,cone)

  w=cbind(scs$x[1:n+2])
  b=scs$x[n+3]
  return(list(w=w,b=b))

}
