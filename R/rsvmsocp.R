# Robust Support Vector Machines based on Second-Order Cone
#

# This librery solves prediction problems using supervised model
# Input must have an X matrix with attributes and a Y vector with
# the real class of data set.

#' @param data : X matrix
#' @param data : Y vector \in \{1,-1\}^n
#' @param  nu :
#' @param  c : solver chosed to solve the optimization problem
#' @param  Method: selection of sigma method estimation where 1 is Cholesky Factorization and 2 Variance Estimation D*D
#' @return Prediction
#' @export
#' @keywords SVM, SOCP, Optimization, Predictive Model
#' @references
#' @seealso \code{\link{sensitivity}}, \code{\link{specificity}}, \code{\link{accuracy}},  \code{\link{auc}}.
#' @return A numeric value between zero and one denoting the area under the curve

setClass(Class="RSVMData",
         representation(nu='numeric', cost='numeric',kappa='numeric', X='matrix', Y='matrix'))


rsvmsocp<-function(X,Y,nu,type=1,cost=0,theta1,theta2)
{  #rsocp(rsvmdata,apptype,sigmatype,solver)
  #rsvmdata=rsvmdata
  #X=rvsmdata@X
  #Y=rvsmdata@Y
  #kappa=rvsmdata@kappa
  if(missing(type)){ type=1}
  if(missing(cost)){ cost=0}
  if(missing(theta1)){ theta1=2^1}
  if(missing(theta2)){ theta2=theta1}


  data<-new('RSVMData', nu=nu, cost=2^cost,kappa=sqrt(nu/(1-nu)), X=X, Y=Y)


  if (!requireNamespace("scs", quietly = TRUE)) {
    stop("SCS needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix needed for this function to work. Please install it.",
         call. = FALSE)
  }

    if(type==9||type==10)
    {
      if(type==9){ wb=five.Tw.1(data@X,data@Y,data@kappa, theta1,theta2)}
      else if(type==10){ wb= five.Tw.2(data@X,data@Y,data@kappa, theta1,theta2)}
      
      mt=dim(data@X)[1]
      w11=sqrt(t(wb$w1)%*%wb$w1)
      w22=sqrt(t(wb$w2)%*%wb$w2)
      y1=data@X%*%wb$w1+wb$b1*cbind(rep(1,mt))
      y2=data@X%*%wb$w2+wb$b2*cbind(rep(1,mt))
      Predict = sign(abs(y2/w22[1])-abs(y1/w11[1]))
      

    }else{
         if(type==1){ wb=one.HC.1(data@X,data@Y,data@kappa) } 
    else if(type==2){ wb=one.HV.2(data@X,data@Y,data@kappa) }
    else if(type==3){ wb=two.RR.1(data@X,data@Y,data@kappa, data@cost)}
    else if(type==4){ wb=two.RR.2(data@X,data@Y,data@kappa, data@cost)}
    else if(type==5){ wb=three.RC.1(data@X,data@Y,data@kappa)}
    else if(type==6){ wb=three.RV.2(data@X,data@Y,data@kappa) } 
    else if(type==7){ wb= four.SC.1(data@X,data@Y,data@kappa,data@cost)}
    else if(type==8){ wb=four.SC.2(data@X,data@Y,data@kappa,data@cost) }
  
    Predict=sign(data@X%*%wb$w+wb$b) # 
    }
    
  # Label prediction class and measure model capacity
  
  acc=medi_auc_accu(Predict,Y)


  return(list(wb=wb, acc=acc, Yp=Predict))


}

