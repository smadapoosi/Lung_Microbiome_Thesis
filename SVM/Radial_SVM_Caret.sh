#Runs RF on stool data set and generates table of OTU importances, selection probabilities, and models

Rscript ~/Radial_SVM-RFE_Caret.R > ~/T_Flavan3ols_Caret_SVM_Output_Model.CAARS_ONLY_Normalized.txt #output model info to text file (also generates table of ROC sensitivity and specificity in a separate doc and plots ROC curves)

#auc
grep "Max AUC" ~/T_Flavan3ols_Caret_SVM_Output_Model.CAARS_ONLY_Normalized.txt | awk '{print $4}' > ~/AUC.txt

#cost
grep "cost" ~/T_Flavan3ols_Caret_SVM_Output_Model.CAARS_ONLY_Normalized.txt | awk '{print $6}' > ~/Cost.txt

#sigma
grep "Tuning" ~/T_Flavan3ols_Caret_SVM_Output_Model.CAARS_ONLY_Normalized.txt | awk '{print $11}' > ~/Sigma.txt

#optimal number of predictors from RFE
grep "Optimal Number of Predictors" ~/T_Flavan3ols_Caret_SVM_Output_Model.CAARS_ONLY_Normalized.txt | awk '{print $6}' > ~/Optim.txt

#number of support vectors
grep "Vectors" ~/T_Flavan3ols_Caret_SVM_Output_Model.CAARS_ONLY_Normalized.txt | awk '{print $6}' >  ~/Num_Support_Vectors.txt

#otus
grep "Otu" ~/T_Flavan3ols_Caret_SVM_Output_Model.CAARS_ONLY_Normalized.txt > ~/Otus.txt

#training error
grep "Training error" ~/T_Flavan3ols_Caret_SVM_Output_Model.CAARS_ONLY_Normalized.txt | awk '{print $4}' > ~/Training_Errors.txt

#type 1 error rate on training
grep "Training Type 1 Error" ~/T_Flavan3ols_Caret_SVM_Output_Model.CAARS_ONLY_Normalized.txt | awk '{print $6}' > ~/Training_Type_1_Error.txt

#type 2 error rate on training
grep "Training Type 2 Error" ~/T_Flavan3ols_Caret_SVM_Output_Model.CAARS_ONLY_Normalized.txt | awk '{print $6}' > ~/Training_Type_2_Error.txt

#type 1 error rate on testing set 
grep "Testing Type 1 Error" ~/T_Flavan3ols_Caret_SVM_Output_Model.CAARS_ONLY_Normalized.txt | awk '{print $6}' > ~/Testing_Type_1_Error.txt

#type 2 error rate on testing set 
grep "Testing Type 2 Error" ~/T_Flavan3ols_Caret_SVM_Output_Model.CAARS_ONLY_Normalized.txt | awk '{print $6}' > ~/Testing_Type_2_Error.txt

Rscript Radial_SVM_Caret_Packager.R

rm ~/Cost.txt
rm ~/AUC.txt
rm ~/Sigma.txt
rm ~/Num_Support_Vectors.txt
rm ~/Otus.txt
rm ~/Training_Errors.txt
rm ~/Testing_Type_1_Error.txt
rm ~/Testing_Type_2_Error.txt
rm ~/Training_Type_1_Error.txt
rm ~/Optim.txt
rm ~/Training_Type_2_Error.txt

mv ~/*.txt ~/Desktop/CAARS_Files/16S_SVM_V1/Output/TMP/