#Extracts AUCs from all 100 models and places them in independent files by phenotype
#!/bin/bash
cd  ~/SVM_RFE_Model_Output
for file in *
do 
	grep "Max AUC" $file | awk '{print $4}' > ~/$file.AUCs.txt
done
