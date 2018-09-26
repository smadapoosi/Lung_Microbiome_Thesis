#Runs RF on stool data set and generates table of OTU importances, selection probabilities, and models

Rscript ~/Regression_RF_Caret.R > ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model.CAARS_ONLY_Normalized.txt #output model info to text file (also generates table of ROC sensitivity and specificity in a separate doc and plots ROC curves)

#Features (OTU) Table - includes otus, prob.select, importance, number of models included in, and percent of models included in
awk '{print $2}' ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model.CAARS_ONLY_Normalized.txt | grep "Otu" > ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_Otus.CAARS_ONLY_Normalized.txt #outputs otus to a file
awk '{print $1}' ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_Otus.CAARS_ONLY_Normalized.txt | sort | uniq -c > ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_OTU_Counts.CAARS_ONLY_Normalized.txt #list of OTUs and how many models they come up in 

paste ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_Otus.CAARS_ONLY_Normalized.txt ~/IL_5_Continuous_Caret_Continuous_Output_Model_Importance.CAARS_ONLY_Normalized.txt > ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_OTU_Table.CAARS_ONLY_Normalized.txt #generate aggregate table for packaging by later Rscript

#MSE Average
grep "Mean of squared residuals:" ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model.CAARS_ONLY_Normalized.txt | awk '{print $5}' > ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_Stats.CAARS_ONLY_Normalized.txt
sed 's/.$//' ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_Stats.CAARS_ONLY_Normalized.txt > ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_MSE.CAARS_ONLY_Normalized.txt #extracts OOB error rates for all models and outputs to a file

#% of Variance Explained by Model
grep "% Var explained:"  ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model.CAARS_ONLY_Normalized.txt | awk '{print $4}' > ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_Var.CAARS_ONLY_Normalized.txt

#generates files for type 1 and type 2 error data
awk '{print $4}' ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model.CAARS_ONLY_Normalized.txt | grep -Eo "[0-9]+\.[0-9]+" > ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_Stats.CAARS_ONLY_Normalized.txt

#Ntree average - used as control to make sure all 100 have 500 ntrees

grep "Number of trees:" ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model.CAARS_ONLY_Normalized.txt | awk '{print $4}' > ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_Ntree.CAARS_ONLY_Normalized.txt

#Mtry average
grep "No. of variables tried at each split:" ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model.CAARS_ONLY_Normalized.txt | awk '{print $8}' > ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_Mtry.CAARS_ONLY_Normalized.txt #outputs mtrys from each model into text file

#Kopt (optimal number of predictors)
grep "predictors" ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model.CAARS_ONLY_Normalized.txt | awk '{print $1}' > ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_Kopt.CAARS_ONLY_Normalized.txt

#Packager of all of individual files into two main output files
Rscript ~/Regression_RF_Caret_Packager.R #package and make table of OTUs, importance, and probability of selection and new files including the average OOB error rate and AUC

#remove temp files
rm ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_Otus.CAARS_ONLY_Normalized.txt
rm ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_Prob_Select.CAARS_ONLY_Normalized.txt
rm ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_Importance.CAARS_ONLY_Normalized.txt
rm ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_Stats.CAARS_ONLY_Normalized.txt
rm ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_Stats2.CAARS_ONLY_Normalized.txt 
rm ~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_Stats1.CAARS_ONLY_Normalized.txt

mv ~/*.CAARS_ONLY_Normalized.txt ~/TMP

