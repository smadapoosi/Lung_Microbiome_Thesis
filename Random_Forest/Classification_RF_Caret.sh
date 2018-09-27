#Runs RF on stool data set and generates table of OTU importances, selection probabilities, and models

Rscript ~/Classification_RF_Caret.R > ~/HA_Healthy_Asthma_Caret_RF_Output_Model.CAARS_ONLY_Normalized.txt #output model info to text file (also generates table of ROC sensitivity and specificity in a separate doc and plots ROC curves)

#Features (OTU) Table - includes otus, prob.select, importance, number of models included in, and percent of models included in
awk '{print $2}' ~/HA_Healthy_Asthma_Caret_RF_Output_Model.CAARS_ONLY_Normalized.txt | grep "Otu" > ~/HA_Healthy_Asthma_Caret_RF_Output_Model_Otus.CAARS_ONLY_Normalized.txt #outputs otus to a file
awk '{print $1}' ~/HA_Healthy_Asthma_Caret_RF_Output_Model_Otus.CAARS_ONLY_Normalized.txt | sort | uniq -c > ~/HA_Healthy_Asthma_Caret_RF_Output_Model_OTU_Counts.CAARS_ONLY_Normalized.txt #list of OTUs and how many models they come up in 

paste ~/HA_Healthy_Asthma_Caret_RF_Output_Model_Otus.CAARS_ONLY_Normalized.txt ~/HA_Healthy_Asthma_Caret_RF_Output_Model_Importance.CAARS_ONLY_Normalized.txt > ~/HA_Healthy_Asthma_Caret_RF_Output_Model_OTU_Table.CAARS_ONLY_Normalized.txt #generate aggregate table for packaging by later Rscript

#Out of box error average
grep "OOB" ~/HA_Healthy_Asthma_Caret_RF_Output_Model.CAARS_ONLY_Normalized.txt | awk '{print $6}' > ~/HA_Healthy_Asthma_Caret_RF_Output_Model_Stats.CAARS_ONLY_Normalized.txt
grep "%" ~/HA_Healthy_Asthma_Caret_RF_Output_Model_Stats.CAARS_ONLY_Normalized.txt | sed 's/.$//' > ~/HA_Healthy_Asthma_Caret_RF_Output_Model_OOB.CAARS_ONLY_Normalized.txt #extracts OOB error rates for all models and outputs to a file

awk '{print $2}' ~/HA_Healthy_Asthma_Caret_RF_Output_Model.CAARS_ONLY_Normalized.txt | awk '/Accuracy/ {for(i=1; i<=5; i++) {getline; print}}' | grep '0.' > ~/HA_Healthy_Asthma_Caret_RF_Output_Model_AUC_CV.CAARS_ONLY_Normalized.txt

#generates files for type 1 and type 2 error data
awk '{print $4}' ~/HA_Healthy_Asthma_Caret_RF_Output_Model.CAARS_ONLY_Normalized.txt | grep -Eo "[0-9]+\.[0-9]+" > ~/HA_Healthy_Asthma_Caret_RF_Output_Model_Stats.CAARS_ONLY_Normalized.txt

#Type 2 error
sed '1d; n; d' ~/HA_Healthy_Asthma_Caret_RF_Output_Model_Stats.CAARS_ONLY_Normalized.txt > ~/HA_Healthy_Asthma_Caret_RF_Output_Model_T_Type_2_Error.CAARS_ONLY_Normalized.txt

#Type 1 error
sed 'n; d' ~/HA_Healthy_Asthma_Caret_RF_Output_Model_Stats.CAARS_ONLY_Normalized.txt > ~/HA_Healthy_Asthma_Caret_RF_Output_Model_T_Type_1_Error.CAARS_ONLY_Normalized.txt #grabs #selects type 1 errors from errors file

#Testing type 1 error
awk '{print $6}' ~/HA_Healthy_Asthma_Caret_RF_Output_Model_Stats1.CAARS_ONLY_Normalized.txt | sed '/^$/d' > ~/HA_Healthy_Asthma_Caret_RF_Output_Model_S_Type_1_Error.CAARS_ONLY_Normalized.txt

#Testing type 2 error
awk '{print $6}' ~/HA_Healthy_Asthma_Caret_RF_Output_Model_Stats2.CAARS_ONLY_Normalized.txt | sed '/^$/d' > ~/HA_Healthy_Asthma_Caret_RF_Output_Model_S_Type_2_Error.CAARS_ONLY_Normalized.txt

#Ntree average - used as control to make sure all 100 have 500 ntrees

grep "Number of trees:" ~/HA_Healthy_Asthma_Caret_RF_Output_Model.CAARS_ONLY_Normalized.txt | awk '{print $4}' > ~/HA_Healthy_Asthma_Caret_RF_Output_Model_Ntree.CAARS_ONLY_Normalized.txt

#Mtry average
grep "No. of variables tried at each split:" ~/HA_Healthy_Asthma_Caret_RF_Output_Model.CAARS_ONLY_Normalized.txt | awk '{print $8}' > ~/HA_Healthy_Asthma_Caret_RF_Output_Model_Mtry.CAARS_ONLY_Normalized.txt #outputs mtrys from each model into text file

#Kopt (optimal number of predictors)
grep "predictors" ~/HA_Healthy_Asthma_Caret_RF_Output_Model.CAARS_ONLY_Normalized.txt | awk '{print $1}' > ~/HA_Healthy_Asthma_Caret_RF_Output_Model_Kopt.CAARS_ONLY_Normalized.txt

#Packager of all of individual files into two main output files
Rscript ~/Classification_RF_Caret_Packager.R #package and make table of OTUs, importance, and probability of selection and new files including the average OOB error rate and AUC

#remove temp files
rm ~/HA_Healthy_Asthma_Caret_RF_Output_Model_Otus.CAARS_ONLY_Normalized.txt
rm ~/HA_Healthy_Asthma_Caret_RF_Output_Model_Prob_Select.CAARS_ONLY_Normalized.txt
rm ~/HA_Healthy_Asthma_Caret_RF_Output_Model_Importance.CAARS_ONLY_Normalized.txt
rm ~/HA_Healthy_Asthma_Caret_RF_Output_Model_Stats.CAARS_ONLY_Normalized.txt
rm ~/HA_Healthy_Asthma_Caret_RF_Output_Model_Stats2.CAARS_ONLY_Normalized.txt 
rm ~/HA_Healthy_Asthma_Caret_RF_Output_Model_Stats1.CAARS_ONLY_Normalized.txt

mv ~/*.CAARS_ONLY_Normalized.txt ~/TMP

