library(vegan)
library(dplyr)
env <- read.csv("~/sputDNA.lowsquam.env.V1.csv", header=T) #rows are study IDs and columns are clinical data
otus <- read.csv("~/sputDNA.lowsquam.otu.V1.csv", header=T) #rows are study ids and columns are microbiota abundances
env_variables <- colnames(env)
env_variables <- env_variables[c(4, 8:96)] #generates list of potential environmental variables
env_variables <- env_variables[-c(16)]
otu_matrix <- otus[, -c(1)] #remove first column (non-numerical) from otu matrix
otu.hel <-decostand(otu_matrix,"hellinger") #hellinger transformation of matrix
env_subset <- env[, c(env_variables)] #remove first 7 columns from environmental matrix so column numbers align
#env_subset <- env_subset[, -16] #remove sputum squamous 80 (already accounted for in row selection)
env_variables <- colnames(env_subset) #reset elements to account for sputum squamous trait removal
write.table(env_variables, file="~/DNA_Adonis_Variables.txt", quote=F)
#output adonis tables for all variables
set.seed(20) #ensures reproducibility of results
for (var in colnames(env_subset)) {
	print(var) #identifies which trait is X
	index <- match((var), env_variables) #returns index of variable to substitute in the adonis formula
	print(index)
	not_na_rows <- which(!is.na(env_subset[, index])) #find which patients have NAs for the given index
	otu.hel.noNA <- otu.hel[c(not_na_rows), ] #subset otu table to rows for patients without NAs for the specified environmental variable
	new_data <- env_subset[c(not_na_rows), ] #remove NAs from env table if they exist
	print(adonis(otu.hel.noNA~new_data[, index], method="euclidean"))
}



