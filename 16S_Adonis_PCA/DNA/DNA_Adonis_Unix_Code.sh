#output DNA Adonis Script to Text File
Rscript DNA_Adonis.R > ~/DNA_Adonis_Output.txt
#select out adonis p-values
grep "new_data" ~/DNA_Adonis_Output.txt > ~/DNA_Adonis_X+Formulas.txt #print rows from adonis output containing "new_data" to a new file (includes formula and output)
grep -v "adonis" ~/DNA_Adonis_X+Formulas.txt > ~/DNA_Adonis_X.txt #remove lines containing formulas
awk '{print $8}' ~/DNA_Adonis_X.txt > ~/DNA_Adonis_P_Vals.txt #print just p-values to a new file
#bind environmental variables to their respective p-values and output to csv file
Rscript DNA_Adonis_Variable_P_Value_Matcher.R 
#sort csv file in reverse numerical order to show which variables have p < 0.05
sort -k3 -n -t, ~/DNA_Adonis_Matched_PVals.csv > ~/DNA_Adonis_Matched_PVals-sorted.csv #final file with p-vals from adonis test of all variables, sorted in ascending order of p-value