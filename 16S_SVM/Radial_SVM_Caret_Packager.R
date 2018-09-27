library(plyr)
library(dplyr)

#Model features file
otu_counts <- read.table('~/Otus.txt')
colnames(otu_counts) <- c("Otu", "Area_Under_ROC_Curve")
otu_aggregate <- otu_counts %>% select(Otu) %>% group_by(Otu) %>% summarize(percent_inclusion= n()) #count number of times each OTU appears in the optimized model
model_table <- ddply(otu_counts, "Otu", numcolwise(mean)) #sum all occurrances' importances by mean for each OTU
model_table <- left_join(model_table, otu_aggregate)
write.csv(model_table, file="~/T_Flavan3ols_Caret_SVM_Output_Model_OTU_Table_Packaged.CAARS_ONLY_Normalized.csv", row.names=F, quote=F)

#Model statistics file
auc <- read.table('~/AUC.txt')
mean_auc <- mean(auc$V1)
mean_auc <- cbind("Mean AUC CV", mean_auc)

n_predictors <- read.table('~/Optim.txt')
mean_n_predictors <- mean(n_predictors$V1)
mean_n_predictors <- cbind("Mean Optimal Number of Predictors", mean_n_predictors)

cost <- read.table('~/Cost.txt') #output model info to text file (also generates table of ROC sensi
mean_cost <- mean(cost$V1)
mean_cost <- cbind("Mean Cost", mean_cost)

sigma <- read.table('~/Sigma.txt') #output model info to text file (also generates table of ROC sensi
mean_sigma <- mean(sigma$V1)
mean_sigma <- cbind("Mean Sigma", mean_sigma)

num_support_vectors <- read.table('~/Num_Support_Vectors.txt') #output model info to text file (also generates table of ROC sensi
mean_num_support_vectors <- mean(num_support_vectors$V1)
mean_num_support_vectors <- cbind("Mean Number of Support Vectors", mean_num_support_vectors)

training_errors <- read.table('~/Training_Errors.txt') #output model info to text file (also generates table of ROC sensi
mean_training_errors <- mean(training_errors$V1)
mean_training_errors <- cbind("Mean Training Error", mean_training_errors)

training_type_1_error <- read.table('~/Training_Type_1_Error.txt')
training_type_1_error <- training_type_1_error[which(training_type_1_error$V1 != "<NA>"), ]
mean_training_type_1_error <- mean(as.numeric(paste(training_type_1_error))) #convert to numerics
mean_training_type_1_error <- cbind("Mean Training Type 1 Error", mean_training_type_1_error)

training_type_2_error <- read.table('~/Training_Type_2_Error.txt')
training_type_2_error <- training_type_2_error[which(training_type_2_error$V1 != "<NA>"), ]
mean_training_type_2_error <- mean(as.numeric(paste(training_type_2_error))) #convert to numerics
mean_training_type_2_error <- cbind("Mean Training Type 2 Error", mean_training_type_2_error)

testing_type_1_error <- read.table('~/Testing_Type_1_Error.txt')
testing_type_1_error <- testing_type_1_error[which(testing_type_1_error$V1 != "<NA>"), ]
mean_testing_type_1_error <- mean(as.numeric(paste(testing_type_1_error))) #convert to numerics
mean_testing_type_1_error <- cbind("Mean Testing Type 1 Error", mean_testing_type_1_error)

testing_type_2_error <- read.table('~/Testing_Type_2_Error.txt')
testing_type_2_error <- testing_type_2_error[which(testing_type_2_error$V1 != "<NA>"), ]
mean_testing_type_2_error <- mean(as.numeric(paste(testing_type_2_error))) #convert to numerics
mean_testing_type_2_error <- cbind("Mean Testing Type 2 Error", mean_testing_type_2_error)

mean_stats <- rbind(mean_auc, mean_cost, mean_sigma, mean_n_predictors, mean_num_support_vectors, mean_training_errors, mean_training_type_1_error, mean_training_type_2_error, mean_testing_type_1_error, mean_testing_type_2_error) #generate df with all of the summary stats
colnames(mean_stats) <- c("Parameter", "Value")

write.csv(mean_stats, file="~/T_Flavan3ols_Caret_SVM_Output_Model_Stats_Packaged.CAARS_ONLY_Normalized.csv", row.names=F, quote=F)