#remove row names from list of variables
awk '{print $2}' ~/DNA_Adonis_LT_Variables.txt > ~/DNA_Adonis_LT_Variables-norownames.txt
tail -n +2 ~/Variable_Full_Names.txt > ~/Variable_Full_Names-nofirstrow.txt 
paste -s -d' \n' ~/Variable_Full_Names-nofirstrow.txt > ~/Variable_Full_Names-nofirstrow-concatenated.txt #outputs rows of trait and their description from original file