 five.Tw.1<-function(X,Y,kappa,theta1,theta2){
  
    n=length(X[2,]) # number of dimentions/attributes of X
    if(length(kappa)==1){
      kappa[2]=kappa[1]
    }
    
      ###calcmedv
    A=X[seq(along=Y)[Y==1],] # selecting class +1
    B=X[seq(along=Y)[Y==-1],] #selecting class -1
    m1=nrow(A)    # number of +1 class
    m2=nrow(B)  
      ## Statistics Measures
    mu=rbind(colMeans(A),colMeans(B))
  
    sigma=list(cov(A)+1e-7*diag(x=1,nrow = n , ncol = n),
               cov(B)+1e-7*diag(x=1,nrow = n , ncol = n))
    
    Mchol=list(t(chol(sigma[[1]])),t(chol(sigma[[2]]))) ##Sr but seem to be the same
    
    #######
    e1=cbind(rep(1,m1))
    e2=cbind(rep(1,m2))
    H= cbind(A,e1)
    G= cbind(B,e2)
    
    HH=t(H)%*%(H)
    HH = HH + theta1*diag(dim(HH)[2]) #regularization
    R1=chol(HH)
    rm(HH)
    
    GG=t(G)%*%G
    GG=GG + theta2*diag(dim(GG)[2])#%regularization
    R2=chol(GG)
    rm(GG)
    
    bb=-cbind(c(1,1,numeric(n+1)))
     ## First Problem
    ##  Building the 1st constraint
  
    At1=rbind(c(1,numeric(n+1)),cbind(0,R1))
    c1=numeric(n+2)
  
  ## Building the 2nd constraint
    At2=matrix(0,nrow=n+1,ncol=n+2)
    At2[1,]=c( 0,- mu[2,], -1)
    At2[1:n+1,1:n+1]=kappa[2]*t(Mchol[[2]])
    c2=c(-1,numeric(n+1))
  
    At=-rbind(At1,At2)
    ct=rbind(c1,c2)
    K.q=c(n+2,n+1)
   
  ## Solve the SOC-problem with SCS
    cone <- list( q = K.q)
    
    scs <- scs(At, ct, -bb , cone)
    w1=cbind(scs$x[1:n+1])
    b1=scs$x[n+2]
    
    rm(At,At1, At2)
  
    ## Second Problem
    ## Building the 1st constraint
    At1=rbind(c(1,numeric(n+1)),cbind(0,R2))
 
    ## Building the 2nd constraint
    At2=matrix(0,nrow=n+1,ncol=n+2)
    At2[1,]=c( 0, mu[1,], 1)
    At2[1:n+1,1:n+1]=kappa[1]*t(Mchol[[1]])

  
    At=-rbind(At1,At2)
    #ct=rbind(c1,c2)
   
    ## Solve the SOC-problem with SCS
    cone <- list( q = K.q )
       
    scs <- scs(At, ct, -bb , cone)
    rm(At,At1, At2, c1, c2)
  
    w2=cbind(scs$x[1:n+1])
    b2=scs$x[n+2]
  
  
  return(list(w1=w1,b1=b1,w2=w2,b2=b2))
 
}