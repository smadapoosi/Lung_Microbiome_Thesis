library(plyr)
library(dplyr)

#Model features file
options(scipen = 999) #convert scientific notation to decimal
model_table <- read.table("~/HA_Healthy_Asthma_Caret_RF_Output_Model_Importance.CAARS_ONLY_Normalized.txt", stringsAsFactors = FALSE) #read in OTU and package prob select and importance data
otu_counts <- read.table("~/HA_Healthy_Asthma_Caret_RF_Output_Model_OTU_Counts.CAARS_ONLY_Normalized.txt") #reads in table with otus and how many models they were included in
otu_counts <- otu_counts[, c(2,1)]
colnames(model_table) <- c("OTU", "MDA", "MDA_Magnitude") #colnames for otu table with importance metrics
model_table <- model_table[-1, ] #remove title row
model_table$MDA <- as.numeric(model_table$MDA) #convert to numeric for aggregation
model_table$MDA_Magnitude <- as.numeric(model_table$MDA_Magnitude) #same as above
colnames(otu_counts) <- c("OTU", "Percent_Inclusion") #includes percent of models that had a given OTU
model_table <- ddply(model_table, "OTU", numcolwise(mean)) #collapses OTUs to unique ones and aggregates their parameters
model_table <- left_join(model_table, otu_counts, by="OTU") #joins counts of OTUs among models to table of otus and their probability of selection
model_table <- na.omit(model_table) #remove rows with NAs, which means there are importance metrics (since it calculates for everything) but they're not featured in any model
write.csv(model_table, file="~/HA_Healthy_Asthma_Caret_RF_Output_Model_OTU_Table_Packaged.CAARS_ONLY_Normalized.csv", row.names=F, quote=F)

#Model statistics file
oob <- read.table("~/HA_Healthy_Asthma_Caret_RF_Output_Model_OOB.CAARS_ONLY_Normalized.txt") #read in OOB error tables
mean_oob <- signif(mean(oob$V1), digits=10) #calculates mean out of box error rate
mean_oob <- cbind("Mean OOB error rate", mean_oob)
 
kopt <- read.table("~/HA_Healthy_Asthma_Caret_RF_Output_Model_Kopt.CAARS_ONLY_Normalized.txt") #read in OOB error tables
mean_kopt <- signif(mean(kopt$V1), digits=10) #calculates mean out of box error rate
mean_kopt <- cbind("Mean Optimal Number of Predictors", mean_kopt)

auc_cv <- read.table("~/HA_Healthy_Asthma_Caret_RF_Output_Model_AUC_CV.CAARS_ONLY_Normalized.txt")
mean_auc_cv <- mean(auc_cv$V1) #calculates mean cross-validated AUC over 100 models
mean_auc_cv <- cbind("Mean AUC CV", mean_auc_cv)

ntree <- read.table("~/HA_Healthy_Asthma_Caret_RF_Output_Model_Ntree.CAARS_ONLY_Normalized.txt")
mean_ntree <- signif(mean(ntree$V1), digits=10) #calculates mean ntree over 100 models
mean_ntree <- cbind("Mean Ntree", mean_ntree)

mtry <- read.table("~/HA_Healthy_Asthma_Caret_RF_Output_Model_Mtry.CAARS_ONLY_Normalized.txt")
mean_mtry <- signif(mean(mtry$V1), digits=10) #calculates mean mtry over 100 models
mean_mtry <- cbind("Mean Mtry", mean_mtry)

type_1_error <- read.table("~/HA_Healthy_Asthma_Caret_RF_Output_Model_T_Type_1_Error.CAARS_ONLY_Normalized.txt")
mean_type_1_error <- signif(mean(type_1_error$V1), digits=10) #calculates mean type 1 error over 100 models
mean_type_1_error <- cbind("Mean CV Type I Error", mean_type_1_error)

type_2_error <- read.table("~/HA_Healthy_Asthma_Caret_RF_Output_Model_T_Type_2_Error.CAARS_ONLY_Normalized.txt")
mean_type_2_error <- signif(mean(type_2_error$V1), digits=10) #calculates mean type 1 error over 100 models
mean_type_2_error <- cbind("Mean CV Type II Error", mean_type_2_error)

type_1_error_testing <- read.table("~/HA_Healthy_Asthma_Caret_RF_Output_Model_S_Type_1_Error.CAARS_ONLY_Normalized.txt")
#mean_type_1_error_testing <- signif(mean(type_1_error_testing$V1), digits=10) #calculates mean type 1 error over 100 models
mean_type_1_error_testing <- t(c("Mean Testing Type I Error", type_1_error_testing))

type_2_error_testing <- read.table("~/HA_Healthy_Asthma_Caret_RF_Output_Model_S_Type_2_Error.CAARS_ONLY_Normalized.txt")
#mean_type_2_error_testing <- signif(mean(type_2_error_testing$V1), digits=10) #calculates mean type 1 error over 100 models
mean_type_2_error_testing <- t(c("Mean Testing Type II Error", type_2_error_testing))

mean_stats <- rbind(mean_oob, mean_auc_cv, mean_ntree, mean_mtry, mean_kopt, mean_type_1_error, mean_type_2_error, mean_type_1_error_testing, mean_type_2_error_testing) #generate df with all of the summary stats
colnames(mean_stats) <- c("Parameter", "Mean Value")

write.csv(mean_stats, file="~/HA_Healthy_Asthma_Caret_RF_Output_Model_Stats_Packaged.CAARS_ONLY_Normalized.csv", row.names=F, quote=F)