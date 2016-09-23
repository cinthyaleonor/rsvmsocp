# Function for compute: AUC, Accuracy, sensitivity, specificity 
# Output: 
#      Yt         - real label
#        Prediction - predictive  label

medi_auc_accu= function(Predict,Y){

# Calculate accuracy identifying  where de predicted class is equal to real class
accu <- (1-length(subset(Predict, Predict==Y))/length(Predict)) 
tPos=length(Predict[Y==1&Predict==1]) #True Positive
tNeg=length(Predict[Y==-1 & Predict==-1]) #True Negative
fPos=length(Predict[Y==-1 & Predict==1]) #False Positive
fNeg=length(Predict[Y==1 & Predict==-1]) #False Negative
sens=tPos/(tPos+fNeg) # sensitivity 
spec=tNeg/(fPos+tNeg) # specificity
AUC=(sens+spec)/2

return(list(AUC=AUC,Accuracy=accu,Sensitivity=sens,Specificity=spec))
}
