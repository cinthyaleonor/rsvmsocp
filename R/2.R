one.HV.2<-function(X,Y,kappa) ## Socp_Hard svmV2_scs
{

  if(length(kappa)==1){
      kappa[2]=kappa[1]
  }
##########################################################
 cm=calc_med(X,Y)

###########################################################
n=dim(cm$mu)[2] # number of dimentions/attributes of X

## %% linear coeficient
bb=-cbind(c(1,numeric(n+1)))

##  Building the 1st constraint
At1=diag(x=1,nrow = n+1, ncol=n+2)
c1=numeric(n+1)

##  Building the 2nd and 3rd constraint
At2=matrix(0, nrow=cm$npos+1, ncol=n+2)
At2[1,2:(n+1)]=cm$mu[1,]
At2[1,(n+2)]=1
At2[2:(cm$npos+1),2:(n+1)]=kappa[1]*t(cm$mchol1)

At3=matrix(0, nrow=cm$nneg+1, ncol=n+2)
At3[1,2:(n+1)]=-cm$mu[2,]
At3[1,(n+2)]=-1
At3[2:(cm$nneg+1),2:(n+1)]=kappa[2]*t(cm$mchol2)  

c2=c(-1,numeric(cm$npos+1))
c3=c(-1,numeric(cm$nneg+1))

At=-rbind(At1,At2,At3) 
ct=cbind(c(c1,c2,c3))

###%% Solve the SOC-problem with SCS
K.q=c(dim(At1)[1],cm$npos +1,cm$nneg+1) #dimension of cone
cone <- list( q = K.q)

scs <- scs(At, ct, -bb , cone)

w=cbind(scs$x[2:(n+1)])
b=scs$x[(n+2)]


 return(list(w=w,b=b))
}

