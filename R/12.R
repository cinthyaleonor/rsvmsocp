six.NH.2<-function(X,Y,kappa,theta1){
  
  n=length(X[2,])
  
  A=X[seq(along=Y)[Y==1],] # selecting class +1
  B=X[seq(along=Y)[Y==-1],] #selecting class -1
  m1=nrow(A)    # number of +1 class
  m2=nrow(B)  
  ## Statistics Measures
  mu=rbind(colMeans(A),colMeans(B))
  
  Mchol=list( (t(A)-mu[1,])*as.vector(array(1,m1))/sqrt(m1),(t(B)-mu[2,])*as.vector(array(1,m2))/sqrt(m2) )
  #Mchol_2=(t(B)-mu[2,])*as.vector(array(1,m2))/sqrt(m2)
  
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
  GG=GG + theta2*diag(dim(GG)[2])# regularization
  R2=chol(GG)
  rm(GG)
  
  ## linear coeficient
  bb=-cbind(c(1,1,numeric(2*n+2)))
  
 
  ## Building the 1st constraint
  At1=rbind(c(1,numeric(2*n+3)),cbind(matrix(0,nrow=n+1,ncol=2),R1,matrix(0,nrow=n+1,ncol=n+1)))
  c1=numeric(n+2)
  
  ### Building the 2nd constraint
  At2=rbind(c(0,1,numeric(2*n+2)),cbind(matrix(0,nrow=n+1,ncol=n+3),R2))
  c2=numeric(n+2)
  ## Building the 3rd constraint
  At3=matrix(0, nrow=m1+1, ncol=2*n+4)
  c3=c(-1, numeric(m1))
 
  At3[1,]=c(0,0,mu[1,],1,-mu[1,],-1) #first row
  At3[-1,]=cbind(0,0,kappa[1]*t(Mchol[[1]]),0,-kappa[1]*t(Mchol[[1]]),0) 
  
  
  ## Building the 4rt constraint
  At4=matrix(0, nrow= m2+1,ncol=2*n+4)
  c4=c(-1,numeric(m2+1)) 
  
  At4[1,3:(2*n+4)]=c(-mu[2,],-1,mu[2,],1) #first Row
  At4[-1,]=cbind(0,0,kappa[2]*t(Mchol[[2]]),0,-kappa[2]*t(Mchol[[2]]),0) 
 
 
  ### 
  K.q=c(n+2, n+2, m1+1, m2+1) 
  ct=cbind(c(c1,c2,c3,c4))
  
  At=-rbind(At1,At2,At3,At4)
  ## Solve the SOC-problem with SCS
  cone <- list( q = K.q )
  
  scs <- scs(At, ct, -bb , cone)
  
  rm(At,At1, At2, c1, c2)
  
  w1=scs$x[3:(n+2)]
  b1=scs$x[n+3]
  w2=scs$x[(n+4):(2*n+3)]
  b2=scs$x[(2*n+4)]
  
  return(list(w1=w1,b1=b1,w2=w2,b2=b2))
  
}