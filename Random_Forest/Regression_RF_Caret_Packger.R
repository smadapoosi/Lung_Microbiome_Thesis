library(plyr)
library(dplyr)

#Model features file
options(scipen = 999) #convert scientific notation to decimal
model_table <- read.table("~/IL_5_Continuous_Caret_Continuous_Output_Model_Importance.CAARS_ONLY_Normalized.txt", stringsAsFactors = FALSE) #read in OTU and package prob select and importance data
otu_counts <- read.table("~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_OTU_Counts.CAARS_ONLY_Normalized.txt") #reads in table with otus and how many models they were included in
otu_counts <- otu_counts[, c(2,1)]
colnames(model_table) <- c("OTU", "MDA", "MDA_Magnitude") #colnames for otu table with importance metrics
model_table <- model_table[-1, ] #remove title row
model_table$MDA <- as.numeric(model_table$MDA) #convert to numeric for aggregation
model_table$MDA_Magnitude <- as.numeric(model_table$MDA_Magnitude) #same as above
colnames(otu_counts) <- c("OTU", "Percent_Inclusion") #includes percent of models that had a given OTU
model_table <- ddply(model_table, "OTU", numcolwise(mean)) #collapses OTUs to unique ones and aggregates their parameters
model_table <- left_join(model_table, otu_counts, by="OTU") #joins counts of OTUs among models to table of otus and their probability of selection
model_table <- na.omit(model_table) #remove rows with NAs, which means there are importance metrics (since it calculates for everything) but they're not featured in any model
write.csv(model_table, file="~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_OTU_Table_Packaged.CAARS_ONLY_Normalized.csv", row.names=F, quote=F)

#Model statistics file
mse <- read.table("~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_MSE.CAARS_ONLY_Normalized.txt") #read in OOB error tables
mean_mse <- signif(mean(mse$V1), digits=10) #calculates mean out of box error rate
mean_mse <- cbind("Mean MSE", mean_mse)
 
kopt <- read.table("~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_Kopt.CAARS_ONLY_Normalized.txt") #read in OOB error tables
mean_kopt <- signif(mean(kopt$V1), digits=10) #calculates mean out of box error rate
mean_kopt <- cbind("Total Number of Predictors", mean_kopt)

perc_var <- read.table("~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_Var.CAARS_ONLY_Normalized.txt")
mean_perc_var<- mean(perc_var$V1) #calculates mean cross-validated AUC over 100 models
mean_perc_var <- cbind("Mean % Variance Explained by Model", mean_perc_var)

ntree <- read.table("~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_Ntree.CAARS_ONLY_Normalized.txt")
mean_ntree <- signif(mean(ntree$V1), digits=10) #calculates mean ntree over 100 models
mean_ntree <- cbind("Mean Ntree", mean_ntree)

mtry <- read.table("~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_Mtry.CAARS_ONLY_Normalized.txt")
mean_mtry <- signif(mean(mtry$V1), digits=10) #calculates mean mtry over 100 models
mean_mtry <- cbind("Mean Mtry", mean_mtry)

mean_stats <- rbind(mean_mse, mean_perc_var, mean_ntree, mean_mtry, mean_kopt) #generate df with all of the summary stats
colnames(mean_stats) <- c("Parameter", "Mean Value")

write.csv(mean_stats, file="~/IL_5_Continuous_Caret_RF_Continuous_Output_Model_Stats_Packaged.CAARS_ONLY_Normalized.csv", row.names=F, quote=F)