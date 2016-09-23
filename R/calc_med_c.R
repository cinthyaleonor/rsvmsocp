
calc_med<-function(X,Y){


  A=X[seq(along=Y)[Y==1],] # selecting class +1
  B=X[seq(along=Y)[Y==-1],] #selecting class -1
  #DataTrain=list("Positive"=A,"Negative"=B)
  npos=nrow(A)    # number of +1 class
  nneg=nrow(B)   # number of -1 class

  #Statistics measures
  mu= rbind(colMeans(A),colMeans(B))

  Mchol_1=(t(A)-mu[1,])*as.vector(array(1,npos))/sqrt(npos)
  Mchol_2=(t(B)-mu[2,])*as.vector(array(1,nneg))/sqrt(nneg)


  return(list(mu=mu,mchol1=Mchol_1,mchol2=Mchol_2,npos=npos,nneg=nneg))
}




