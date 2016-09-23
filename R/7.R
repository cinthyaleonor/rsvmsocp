four.SC.1<-function(X,Y,kappa,C){

  n=length(X[2,]) # number of dimentions/attributes of X
  if(length(kappa)==1){
    kappa[2]=kappa[1]
  }
  
  A=X[seq(along=Y)[Y==1],] # selecting class +1
  B=X[seq(along=Y)[Y==-1],] #selecting class -1
  
  ## Statistics Measures
  mu=rbind(colMeans(A),colMeans(B))
  Ylab=c(1,-1)
  
  sigma=list(cov(A)+1e-7*diag(x=1,nrow = n , ncol = n),
             cov(B)+1e-7*diag(x=1,nrow = n , ncol = n))
  
  Mchol=list(t(chol(sigma[[1]])),t(chol(sigma[[2]])))
  
  #  %% linear coeficient
  bb=-cbind(c(C,1,numeric(n+1)))

  ## % Building the 1st constraint
  At1=cbind(0,diag(x=1,nrow = n+1, ncol=n+2))
  c1=cbind(numeric(n+1))
  # % Building the 2nd constraint
  At2=matrix(0,nrow=(n+1)*2,ncol=n+3)
  c2=cbind(numeric(2*(n+1)))
  At2[1,1]= 1
  
  K.q<-length(At1[,1]) #dimension of cone
  
  
  for(i in 1:2){
    f=(n+1)*i-n #starting row
    At2[f,]= c(1,0,Ylab[i]*mu[i,],Ylab[i])
    At2[(f+1):((n+1)*i),3:(n+2)]=kappa[i]*t(Mchol[[i]])
    K.q[(i+1)]=n+1
    c2[f]=-1
  }
  
  Dt=c(1, numeric(n+2))
  f=0
  K.l=1
  
  At=-rbind(Dt,At1,At2)
  ct=rbind(f,c1,c2)
  cone <- list( l=K.l, q = K.q  )
  
  scs <- scs(At, ct, bb , cone)
  
  w=cbind(scs$x[3:(n+2)])
  b=scs$x[(n+3)]

  return(list(w=w,b=b))
  
}
