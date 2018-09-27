variables <- read.table("~/DNA_Adonis_Variables.txt")
pvals <- read.table("~/DNA_Adonis_P_Vals.txt")
variable_pval_matrix <- cbind(variables, pvals)
write.csv(variable_pval_matrix, file="~/DNA_Adonis_Matched_PVals.csv")