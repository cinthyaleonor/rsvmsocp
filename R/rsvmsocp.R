# Robust Support Vector Machines based on Second-Order Cone
#

# This librery solves prediction problems using supervised model
# Input must have an X matrix with attributes and a Y vector with
# the real class of data set.

#' @param data : X matrix
#' @param data : Y vector \in \{1,-1\}^n
#' @param  nu :
#' @param  c : solver chosed to solve the optimization problem
#' @param  Method = type: selection of sigma method estimation where 1 is Cholesky Factorization and 2 Variance Estimation D*D
#' @return Prediction
#' @export
#' @keywords SVM, SOCP, Optimization, Predictive Model
#' @references
#' @seealso \code{\link{sensitivity}}, \code{\link{specificity}}, \code{\link{accuracy}},  \code{\link{auc}}.
#' @return A numeric value between zero and one denoting the area under the curve


#  1) Hard-margin: Eq. (1)
#  2) R1-R2-margin
#  3) R-margin
#  4) Soft-margin: Eq.(2) 
#  5) Tw
#  6) NH 
#  7) L1-socp


setClass(Class="RSVMData",
         representation(nu='numeric', cost='numeric',kappa='numeric', X='matrix', Y='matrix'))


rsvmsocp<-function(X,Y,nu,cost=0,theta1,theta2, P, type, opt)
{  
  if(missing(type)){ type=1}
  if(missing(cost)){ cost=0}
  if(missing(theta1)){ theta1=2^1}
  if(missing(theta2)){ theta2=theta1}
  if(missing(P)){  P=list(kernel='rbf', sigma=2^2)}
  if(missing(opt)){ opt=1}
 
  data<-new('RSVMData', nu=nu, cost=2^cost,kappa=sqrt(nu/(1-nu)), X=X, Y=Y)

  
          if(opt==1){
              switch(type,
                { wb = one.HC.1(data@X, data@Y, data@kappa) },  #1
                { wb = two.RR.1(data@X, data@Y, data@kappa, data@cost) }, #2 
                { wb = three.RC.1(data@X, data@Y, data@kappa) }, #3 
                { wb = four.SC.1(data@X, data@Y, data@kappa, data@cost)}, #4
                { wb = five.Tw.1(data@X, data@Y, data@kappa, theta1, theta2)}, #5
                { wb = six.NH.1(data@X, data@Y, data@kappa, theta1)}, #6
                { wb = SOCP_Ker_Hsvm_scs(data@X, data@Y, data@kappa, P)}, #7
                { wb = SOCP_Ker_Ssvm_scs(data@X, data@Y, data@kappa, P, C)}, #8
                stop("Invalid Model"))
           } 
           else{
                switch(type,
                {wb = one.HV.2(data@X, data@Y, data@kappa)},
                {wb = two.RR.2(data@X, data@Y, data@kappa, data@cost) },
                {wb = three.RV.2(data@X, data@Y, data@kappa) },
                {wb = four.SC.2(data@X, data@Y, data@kappa, data@cost) },
                {wb= five.Tw.2(data@X, data@Y, data@kappa, theta1, theta2)},
                {wb= six.NH.2(data@X, data@Y, data@kappa, theta1)},
                stop("Invalid Model"))
            } 
  
      if( any(type %in% c(1:4)) == TRUE)
      { Predict=sign(data@X%*%wb$w+wb$b) } 
      else if(type %in% c(7,8)){
        Predict=sign(wb$w+wb$b)
      
      }
      else{
      
        mt=dim(data@X)[1]
        w11=sqrt(t(wb$w1)%*%wb$w1)
        w22=sqrt(t(wb$w2)%*%wb$w2)
        y1=data@X%*%wb$w1+wb$b1*cbind(rep(1,mt))
        y2=data@X%*%wb$w2+wb$b2*cbind(rep(1,mt))
      
        Predict = sign(abs(y2/w22[1])-abs(y1/w11[1]))}
  
  # Label prediction class and measure model capacity
  
  acc=medi_auc_accu(Predict,Y)
  return(list(wb=wb, acc=acc, Yp=Predict))

}


