two.RR.2<-function(X,Y,kappa, C) # Socp_ R12 svmV1_scs  ##  what's r12????
{
  
  if(length(kappa)==1){
    kappa[2]=kappa[1]
  }
  ##########################################################
  cm=calc_med(X,Y)
  ###########################################################
  n=dim(cm$mu)[2] # number of dimentions/attributes of X
  ## %% linear coeficient
  
  bb=-cbind(c(-1,-C,numeric(n+1)))
  
  ## %% Linear constraint
  Dt=rbind(c(1, numeric(n+2)),c(0, 1,numeric(n+1)))
  
  ##%% Building the 2nd constraint
  At2=matrix(0, nrow=cm$npos+1, ncol=n+3) 
  At3=matrix(0, nrow=cm$nneg+1, ncol=n+3)
  c2=numeric(cm$npos+1)
  c3=numeric(cm$nneg+1)
  
  At2[1,]=c(-1,0,cm$mu[1,], 1)
  At2[2:(cm$npos+1),3:(n+2)]=kappa[1]*t(cm$mchol1)
  
  At3[1,]=c(0,-1,cm$mu[2,], -1)
  At3[2:(cm$nneg+1),3:(n+2)]=kappa[2]*t(cm$mchol2)
  
  At=-rbind(Dt,At2,At3) #
  ct=cbind(c(0,0,c2,c3)) ## f,c1,c2
  K.q=c(cm$npos+1, cm$nneg+1) #dimension of cone
  cone <- list( q = K.q , l=2)
  
  scs=scs(At,ct,-bb,cone)
  
  w=cbind(scs$x[3:(n+2)])
  b=scs$x[(n+3)]
  return(list(w=w,b=b))
  
}
