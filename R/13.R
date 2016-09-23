SOCP_Ker_Hsvm_scs<-function(X,Y,kappa,P){


  ### calc_med
  A=X[seq(along=Y)[Y==1],] # selecting class +1
  B=X[seq(along=Y)[Y==-1],] #selecting class -1
  m_i=c(nrow(A), nrow(B))
  D= rbind(A,B)
  ## Statistics Measures
  mu=rbind(colMeans(A),colMeans(B))
  Kg=type_kern(D,P)
  m=sum(m_i)

  g=matrix(nrow=m,ncol=2)
  
  for(r in 1:2){ 
    tr=c((ifelse((r-1)>0,sum(m_i[1:(r-1)])+1,0)):sum(m_i[1:r]))
  
      for(j in 1:2){
        tj= c(ifelse(j-1>0,sum(m_i[1:(j-1)])+1,0):sum(m_i[1:j]))
        g[tj,r]=Kg[tj,tr]%*%rep(1,m_i[r])/m_i[r]  ## [g1',g2']
      }
    }

  
  Kt=list((diag(x=1,m_i[1]) - rep(1,m_i[1]) / m_i[1]) %*% t(Kg[,1:m_i[1]]),
          (diag(x=1,m_i[2]) - rep(1,m_i[2]) / m_i[2]) %*% t(Kg[,(m_i[1]+1):m]))
  
  Mchol=list(Kt[[1]]/sqrt(m_i[1]), Kt[[2]]/sqrt(m_i[1]))
  
  
    
  MK=Kg+(1e-7)*diag(dim(Kg)[2])
  R_chol=chol(MK)
  rm(Kt, MK)

  ## %% linear coeficient
  bb=-cbind(c(1,numeric(m+1)))
  ## %% Building the 1st constraint
  At1=rbind(c(1,numeric(m+1)), cbind(matrix(0,nrow=m,ncol=1),R_chol,matrix(0,nrow=m,ncol=1)))
    c1= numeric(m+1)
  

  ## %% Building the 2nd constraint
  At2=matrix(0, nrow=(m_i[1]+1),ncol=m+2) 
  c2=c(-1, numeric(m_i[1]))

  At2[1,]=c(1,t(g[,1]),1) ## check!!! 
  At2[-1,2:(m+1)]=kappa[1]*Mchol[[1]]

  ## %% Building the 3rd constraint
  At3=matrix(0, nrow=m_i[2]+1,ncol=m+2)
  c3=c(-1,numeric(m_i[2]))
  At3[1,]=c(-1,-t(g[,2]),-1)
  At3[-1,2:(m+1)]=kappa[2]*Mchol[[2]]

  At=-rbind(At1,At2, At3)
  ct=cbind(c(c1,c2,c3))
  K.q=c(m+1,m_i[1]+1, m_i[2]+1)
  rm(At1, At2, At3, c1, c2, c3)
## Solve the SOC-problem with SCS
  cone <- list( q = K.q )

  scs <- scs(At, ct, -bb , cone)
  
  s=cbind(scs$x[2:(m+1)])
  b=scs$x[(m+2)]

return(list(w=s,b=b))
  }