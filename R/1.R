one.HC.1 <-function(X,Y,kappa) # Socp_H Chol - svmV1_scs
{


  n=length(X[2,])

  if(length(kappa)==1){
    kappa[2]=kappa[1] #
  }

  A=X[seq(along=Y)[Y==1],] # selecting class +1
  B=X[seq(along=Y)[Y==-1],] #selecting class -1

  ## %%% Medidas estadisticas
  mu=rbind(colMeans(A),colMeans(B))
  Ylab=c(1,-1)

  sigma=list(cov(A)+1e-7*diag(x=1,nrow = n , ncol = n),
             cov(B)+1e-7*diag(x=1,nrow = n , ncol = n))

  Mchol=list(t(chol(sigma[[1]])),t(chol(sigma[[2]])))

  ## linear coeficient
  bb=-cbind(c(1,numeric(n+1))) #check symbol

  ## %% Building the 1st constraint
  At1=diag(x=1,nrow = n+1, ncol=n+2)
  c1=numeric(n+1)

  ## %% Building the 2nd and 3rd constraint
  At2=matrix(0,nrow=2*(n+1),ncol = (n+2))
  c2=numeric(2*(n+1)) #%matrix(0,nrow=2*(n+1),ncol = 1)

  
  for(i in 1:2){
    f=(i-1)*(n+1)+1 #starting row
    At2[f,2:(n+2)]=c(Ylab[i]*mu[i,], Ylab[i])
    At2[f+1:n,2:(n+1)]=kappa[i]*t(Mchol[[i]])
    c2[f]=-1
  }

   At=-rbind(At1,At2) ## sparce
   ct=cbind(c(c1,c2))

  ## %% Solve the SOC-problem with SCS

  K.q=c(dim(At1)[1],n+1,n+1) #dimension of cone

  cone <- list( q = K.q)
  scs <- scs(At, ct, -bb , cone)

    w=cbind(scs$x[2:(n+1)])
    b=scs$x[n+2]


return(list(w=w,b=b))
}
