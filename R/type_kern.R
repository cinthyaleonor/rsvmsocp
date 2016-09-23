type_kern<-function(X,P,Xt){
    
  if(missing(Xt)){Xt=X}
    
 
  if(P[[1]][1]=="linear"){ K=X%*%t(Xt)}
    else if(P[[1]][1]=="rbf"){
      ps = X%*%t(Xt) 
      nps=dim(ps)[1]
      pps=dim(ps)[2]
      
      normx = rowSums(X^2)## 
      normxsup = rowSums(Xt^2) ##
      ps = -2*ps + t(repmat(normx,pps,1)) + repmat(t(normxsup),nps,1)
      K = exp(-ps/(2*P$sigma^2))
    }
  
   return(K)

}