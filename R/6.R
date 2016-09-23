three.RV.2<-function(X,Y,kappa){

  if(length(kappa)==1){
    kappa[2]=kappa[1]
  }
  ##########################################################
  ## This is cal_med  maybe is better to use that function
  cm=calc_med(X,Y)
  
  ## end calc_med
  ###########################################################
  n=dim(cm$mu)[2] # number of dimentions/attributes of X
  
    ##  linear coeficient
  bb=-cbind(c(-1,numeric(n+1))) 
  
  ##  Linear constraint
  Dt=rbind(diag(x=1,nrow=n+1,ncol=n+2),cbind(0, -diag(x=1,nrow=n,ncol=n+1)))
  
           
  f=cbind(c(0,rep(x=1,2*n)))
  
  ## %% Building the 2nd and 3rd constraint
  ##  Building the 2nd and 3rd constraint
  At2=matrix(0, nrow=cm$npos+1, ncol=n+2)
  At2[1,1]=-1
  At2[1,1:n+1]=cm$mu[1,]
  At2[1,n+2]=1
  At2[1:cm$npos+1,1:n+1]=kappa[1]*t(cm$mchol1)
  
  At3=matrix(0, nrow=cm$nneg+1, ncol=n+2)
  At3[1,n+2]=-1
  At3[1,1:n+1]=-cm$mu[2,]
  At3[1,n+2]=-1
  At3[1:cm$nneg+1,1:n+1]=kappa[2]*t(cm$mchol2)
  
  c2=cbind(numeric(cm$npos+1))
  c3=cbind(numeric(cm$nneg+1))

  ## Solve the SOC-problem with SCS
  At=-rbind(Dt,At2,At3) #
  ct=rbind(f,c2,c3) ## f,c1,c2
  K.q=c(cm$npos+1,cm$nneg+1) #dimension of cone
  cone=list( q = K.q , l=(2*n+1))

  scs=scs(At,ct,-bb,cone)

  w=cbind(scs$x[1:n+1])
  b=scs$x[n+2]
  return(list(w=w,b=b))


}