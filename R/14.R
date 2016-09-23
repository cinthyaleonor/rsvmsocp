# % Solver that solving the problem (Xi-SOC-SVM con Kernel )
# 
# %   min t+C*xi
# %   s.t.
# %         |Rv|<=t
# %         (v'*h_1+b)>=1+kapa_1|D^T_1v|,  clase +1
#            %         -(v'*h_2+b)>=1+kapa_2|D^T_2v|, clase -1
# %         Xi>=0      
# % Datos salida: 
#   %       s = L^(-1)v
#   %       b - bias
#   %       Kv - funcion evaluado en el Testing
  
SOCP_Ker_Ssvm_scs<-function(X,Y,kappa,P,C){
  
  A=X[seq(along=Y)[Y==1],] # selecting class +1
  B=X[seq(along=Y)[Y==-1],] #selecting class -1
  m_i=c(nrow(A), nrow(B))
  D=rbind(A,B)
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
  
  ##linear coeficient
  bb=-cbind(C,1,numeric(m+1))
  
  ## Building the 1st constraint
  At1=rbind(c(0,1,numeric(m+1)), cbind(matrix(0,nrow=m,ncol=2),R_chol,matrix(0,nrow=m,ncol=1)))
  c1= numeric(m+1)
  
  At2=matrix(0, nrow=(m_i[1]+1),ncol=m+3)
  c2=c(-1, numeric(m_i[1]))
  
  At2[1,]=c(1,0,t(g[,1]),1)
  At2[-1,3:(m+2)]=kappa[1]*Mchol[[1]]
  
  ## %% Building the 3rd constraint
  At3=matrix(0, nrow=m_i[2]+1,ncol=m+3)
  c3=c(-1,numeric(m_i[2]))
  At3[1,]=c(1,0,-t(g[,2]),-1)
  At3[-1,3:(m+2)]=kappa[2]*Mchol[[2]]
  
  Dt=c(1,numeric(m+2))
  f=0
  K.l=1
  
  At=-rbind(Dt,At1,At2,At3)
  ct=cbind(c(f,c1,c2,c3))
  K.q=c(m+1,m_i[1]+1, m_i[2]+1)
  rm(At1, At2, At3, c1, c2, c3,Dt,f)
  cone <- list( q = K.q , l=K.l )
  
  scs <- scs(At, ct, -bb , cone)
  
  s=cbind(scs$x[3:(m+2)])
  b=scs$x[(m+3)]
  
  return(list(w=s,b=b))
}