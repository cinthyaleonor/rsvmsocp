four.SC.2<-function(X,Y,kappa,C){

    if(length(kappa)==1){
      kappa[2]=kappa[1]
    }
    
    ##########################################################
    cm=calc_med(X,Y)
    ###########################################################
    n=dim(cm$mu)[2] # number of dimentions/attributes of X
    
    # linear coeficient
    bb=-cbind(c(C,1,numeric(n+1)))
    
    ## Building the 1st constraint
    At1=cbind(0,diag(x=1,ncol=n+2,nrow = n+1)) 
    c1=numeric(n+1)
    
    ## Building the 2nd and 3rd constraints
    At2=matrix(0, nrow=cm$npos+1, ncol=n+3)
    At3=matrix(0, nrow=cm$nneg+1, ncol=n+3)
    c2=c(-1,numeric(cm$npos-1))
    c3=c(-1,numeric(cm$nneg-1))
    
    
    At2[1,]=c(1, 0,  cm$mu[1,], 1)
    At2[2:(cm$npos+1),3:(n+2)]=kappa[1]*t(cm$mchol1)
    
    At3[1,]=c(1, 0,  -cm$mu[2,], -1)
    At3[2:(cm$nneg+1),3:(n+2)]=kappa[2]*t(cm$mchol2) 
    
    K.q=c(dim(At1)[1],cm$npos +1,cm$nneg+1) #dimension of cone
    
    ### Linear constraint
    Dt=c(1, numeric(n+2))
    f=0
    K.l=1
    
    At=-rbind(Dt,At1,At2,At3)
    ct=cbind(c(f,c1,c2,c3))
    cone <- list( q = K.q , l=K.l)
    
    scs <- scs(At, ct, -bb , cone)
    w=cbind(scs$x[3:(n+2)])
    b=scs$x[(n+3)]
    
    
    return(list(w=w,b=b))


}