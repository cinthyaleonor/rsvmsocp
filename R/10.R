  five.Tw.2<-function(X,Y,kappa, theta1,theta2) {
    
    n=length(X[2,])
    ### calc_med
    A=X[seq(along=Y)[Y==1],] # selecting class +1
    B=X[seq(along=Y)[Y==-1],] #selecting class -1
    m1=nrow(A)    # number of +1 class
    m2=nrow(B)  
    ## Statistics Measures
    mu=rbind(colMeans(A),colMeans(B))
    
    Mchol=list((t(A)-mu[1,])*as.vector(array(1,m1))/sqrt(m1),(t(B)-mu[2,])*as.vector(array(1,m2))/sqrt(m2))
     
    #end calcmed
    #######
    e1=cbind(rep(1,m1))
    e2=cbind(rep(1,m2))
    H= cbind(A,e1)
    G= cbind(B,e2)
    
    HH=t(H)%*%H
    HH = HH + theta1*diag(dim(HH)[2]) #regularization
    R1=chol(HH)
    rm(HH)
    
    GG=t(G)%*%G
    GG=GG + theta2*diag(dim(GG)[2])#%regularization
    R2=chol(GG)
    rm(GG)
    
      ## linear coeficient
    bb=-cbind(c(1,1,numeric(n+1)))
    
    ## First Problem
    ##  Building the 1st constraint
    At1=rbind(c(1,numeric(n+1)),cbind(0,R1))
    c1=numeric(n+2)
    
    ## Building the 2nd constraint
    At2=matrix(0,nrow=m2+1,ncol=n+2)
    At2[1,]=c( 0,- mu[2,], -1)
    At2[2:(m2+1),2:(n+1)]=kappa[2]*t(Mchol[[2]])
    
    c2=c(-1,numeric(m2+1))
    
    At=-rbind(At1,At2)
    ct=cbind(c(c1,c2))
    K.q=c(n+2,m2+1)
    
    ## Solve the SOC-problem with SCS
    cone <- list( q = K.q)
    
    scs <- scs(At, ct, -bb , cone)
    w1=cbind(scs$x[2:(n+1)])
    b1=scs$x[(n+2)]
    
    rm(At,At1, At2,c2,K.q)
    
    ## Second Problem
    ## Building the 1st constraint
    At1=rbind(c(1,numeric(n+1)),cbind(0,R2))
    
    ## Building the 2nd constraint
    At2=matrix(0,nrow=m1+1,ncol=n+2)
    At2[1,]=c( 0, mu[1,], 1)
    At2[2:(m1+1),2:(n+1)]=kappa[1]*t(Mchol[[1]])
    c2=c(-1,numeric(m1+1))
    
    
    At=-rbind(At1,At2)
    ct=cbind(c(c1,c2))
    K.q=c(n+2,m1+1)
   
    ## Solve the SOC-problem with SCS
    cone <- list( q = K.q )
    
    scs <- scs(At, ct, -bb , cone)
    rm(At,At1, At2, c1, c2)
    
    w2=cbind(scs$x[2:(n+1)])
    b2=scs$x[(n+2)]
    
    
    return(list(w1=w1,b1=b1,w2=w2,b2=b2))
    }