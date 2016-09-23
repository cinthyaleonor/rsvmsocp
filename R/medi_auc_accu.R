# Function for compute: AUC, Accuracy, sensitivity, specificity 
# Output: 
#      Yt         - real label
#        Prediction - predictive  label

medi_auc_accu= function(Predict,Y){

# Calculate accuracy identifying  where de predicted class is equal to real class
accu <- (length(subset(Predict, Predict==Y))/length(Predict)) 
tPos=length(Predict[Y==1  & Predict==Y]) #True Positive
tNeg=length(Predict[Y==-1 & Predict==Y]) #True Negative
fPos=length(Predict[Y==-1 & Predict!=Y]) #False Positive
fNeg=length(Predict[Y==1  & Predict!=Y]) #False Negative
sens=tPos/(tPos+fNeg) # sensitivity 
spec=tNeg/(fPos+tNeg) # specificity
AUC=(sens+spec)/2

return(list(AUC=AUC,Accuracy=accu,Sensitivity=sens,Specificity=spec))
}
