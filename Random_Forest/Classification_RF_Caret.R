require(caret)
require(ggplot2)
require(randomForest)
library(pROC)
library(vegan)
library(plyr)
library(dplyr)
library(tidyverse)

#prepare data set for rf model generation
otu_matrix <- read.csv("~/Desktop/CAARS_Files/Random_Forest_V1/Data/DNA/HA_Healthy_Asthma_Caret_RF.dna.otu.caarsonly.normalized.csv")
otu_matrix_otus <- otu_matrix[, 3:ncol(otu_matrix)] #assign only otus to this matrix
sample_category <- otu_matrix[, 2] #assign environmental variable to this matrix
sample_category <- factor(sample_category)
data <- data.frame(sample_category, otu_matrix_otus) #generates otu matrix with environmental variable category as the first column
data <- data[, -nearZeroVar(data)] #remove otus that have a variance near 0, which slow algorithm time and are not usually predictive

roc_summary_train <- NULL #store roc scores on training set
roc_summary_test <- NULL #store roc scores on testing set
type_1_error_list <- NULL #store list of type 1 errors on testing set
type_2_error_list <- NULL

for (i in 1:100) {
	iteration_number <- i #reference iteration number for files

	#Generate training set
	training_prop <- 0.8 #proportion of total data set to be used for training
	training_samples_count <- training_prop*(nrow(data)) #number of samples to be used in training data set
	prop_cases <- length(which(data$sample_category == 1))/nrow(data) #finds proportion of patients that are cases
	training_cases_count <- prop_cases*training_samples_count #number of cases to be included in training set
	cases <- which(data$sample_category == 1) #identifies cases in original data set
	cases_sample <- sample(cases, training_cases_count) #samples cases to make up identical distribution of cases in the training set as in the original data set

	prop_controls <- length(which(data$sample_category == 0))/nrow(data) #finds proportion of patients that are controls
	training_controls_count <- prop_controls*training_samples_count #number of controls to be included in training set
	controls <- which(data$sample_category == 0) #identifies controls in original data set
	controls_sample <- sample(controls, training_controls_count) #samples controls to make up identical distribution of cases in the training set as in the original data set
	data_cases <- data[cases_sample, ] #subset bmi data for cases to make up half of training data set
	data_controls <- data[controls_sample, ] #subset environmental variable data for controls to make up half of training data set
	data_training <- rbind(data_cases, data_controls) #create training data set 

	set.seed(20) #reproducibility
	#Run random forest model on training set along with 10-fold cross validation
	trcontrol <- trainControl(method="repeatedcv", #cross-validation
							  number=10, #10-fold cross validation
							  repeats = 3, #reduces overfitting on training set
							  savePredictions = TRUE) #saves sensitivity and specificity measures for ROC curve
	#mtry <- sqrt(ncol(data_training))
	model <-train(sample_category~.,data=data_training, 
		method="rf",
		trControl = trcontrol,
		prox=TRUE,
		#tuneLength = 15, #select 15 random numbers of variables to try at each split (mtry)
		saveDetails = TRUE,
		allowParallel=FALSE,
		importance = TRUE)

	#generate ROC data for model on training set
	probabilities <- predict(model$finalModel, type='prob')[, 2] #generate fitted values for roc curve
	roc <- roc(data_training$sample_category~probabilities) #generate roc curve
	roc_summary_train <- rbind(roc_summary_train, data.frame(sensitivity = roc$sensitivities, specificity = roc$specificities))

	#generate testing set for model to test on
	unsampled_cases <- setdiff(cases, cases_sample) #generate list of cases for testing data set not in the training sample
	unsampled_controls <- setdiff(controls, controls_sample) #generate list of controls for testing data set not in the training sample
	testing_data_cases <- data[unsampled_cases, ] #subsets bmi data set by cases
	testing_data_controls <- data[unsampled_controls, ] #subsets bmi data set by controls
	data_testing <- rbind(testing_data_cases, testing_data_controls) #combines cases and controls to make a testing data set
	test <- predict(model, data_testing, type="prob", norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE) #test random forest model on remaining 3/4 of data set
	test_results <- as.factor(colnames(test)[max.col(test,ties.method="first")]) #convert default output to one for confusion matrix

	confusion_matrix <- table(data_testing$sample_category, test_results) #confusion matrix
	if (length(confusion_matrix) == 4) {
		type_1_error <- signif(as.numeric(confusion_matrix["0","1"])/(as.numeric(confusion_matrix["0","0"])+as.numeric(confusion_matrix["0","1"])), digits=10) #calculates type 1 error on testing set
		type_2_error <- signif(as.numeric(confusion_matrix["1","0"])/(as.numeric(confusion_matrix["1","1"]) + as.numeric(confusion_matrix["1","0"])), digits=10) #calculates type 2 error on testing set
	} else {
		type_1_error <- 0
		type_2_error <- 1
	}
	type_1_error_list <- c(type_1_error_list, type_1_error)
	type_2_error_list <- c(type_2_error_list, type_2_error)

	#generate ROC curve on testing data set (out of box)
	probs <- test[, 2] #positives among testing set
	roc <- roc(data_testing$sample_category~probs) #generate roc curve
	roc_summary_test<- rbind(roc_summary_test, data.frame(sensitivity = roc$sensitivities, specificity = roc$specificities))

	#print out model statistics
	print(model$finalModel) #outputs final model statistics (type of RF, ntree, mtry, and confusion matrix)
	print(model) #prints accuracy measures and info on each resampling
	importance <- varImp(model$finalModel, scale = FALSE, type = 1) #outputs predictors and their respective importance measure (MDA)
	abs_importance <- abs(importance) #total change in accuracy (increase or decrease)
	options(scipen = 999) #convert scientific notation to decimal
	abs_importance <- data.frame(abs_importance)
	abs_importance$"OTU" <- rownames(abs_importance)
	importance <- data.frame(importance)
	importance$"OTU" <- rownames(importance)
	importance <- left_join(importance, abs_importance, by = "OTU")
	importance <- select(importance, OTU, Overall.x, Overall.y) #make table of otus, MDA raw, and magnitude of MDA
	colnames(importance) <- c("OTU", "MDA", "MDA_Magnitude")
	write.table(importance, file="~/HA_Healthy_Asthma_Caret_RF_Output_Model_Importance.CAARS_ONLY_Normalized.txt", quote=F, row.names=F)
	sig_importances <- arrange(importance, desc(MDA)) #sort OTUs by descending magnitude to subset important ones
	sig_importances <- head(sig_importances, 25) #select out top 25 OTUs (most important) - TODO CHANGE TO 30 WITH RERUN********
	print(sig_importances)

	#generate AUC curve for model
	file_path <- paste("~/HA_Healthy_Asthma_Caret_RF_Plot_IT", iteration_number, "_CAARS_ONLY_Normalized.pdf", sep="")
	pdf(file_path)
	plot <- plot(model, main = "Random Forest Model Variable Selection for HA_Healthy_Asthma_Median - Caret Control")
	print(plot)
	dev.off()
}
#using ggplot
#generate cumulative ROC curve (all specificities and sensitivities across n = 100 cv models) -> TRAINING SET
colnames(roc_summary_train) <- c("sensitivity", "specificity")
roc_summary_train <- as.data.frame(roc_summary_train)
write.table(roc_summary_train, file="~/HA_Healthy_Asthma_Caret_RF_Output_Model_ROC_Training_CAARS_ONLY_Normalized.txt", quote=F, row.names=F)
roc_summary_train$specificity <- parse_number(roc_summary_train$specificity) #convert factor to number
roc_summary_train <- roc_summary_train %>% group_by(sensitivity) %>% summarize(specificity = max(specificity)) #summarize by max specificity per sensitivity
perfect_sensitivity <- c(1, 0) #add in data point for perfect sensitivity
perfect_specificity <- c(0, 1) #add in data point for perfect specificity
roc_summary_train <- rbind(roc_summary_train, perfect_sensitivity, perfect_specificity)
pdf("~/HA_Healthy_Asthma_Caret_RF_ROC_Curve_Training_CAARS_ONLY_Normalized.pdf")
ggplot(data = roc_summary_train) + 
	geom_step(mapping = aes(x = (1-specificity), y = parse_number(roc_summary_train$sensitivity)), color = "blue", direction = "vh") + #generate roc curve
	geom_abline(intercept = 0, slope = 1) + #line of equal specificity and sensitivity
	labs(x = "False Positive Rate (1 - Specificity)", y = "True Positive Rate (Sensitivity)") + 
	scale_y_continuous(lim=c(0, 1))
	#scale_x_reverse(lim=c(1, 0)) #flip x axis
dev.off()

#ROC curve on testing data set
colnames(roc_summary_test) <- c("sensitivity", "specificity")
roc_summary_test <- as.data.frame(roc_summary_test)
write.table(roc_summary_test, file="~/HA_Healthy_Asthma_Caret_RF_Output_Model_ROC_Testing_CAARS_ONLY_Normalized.txt", quote=F, row.names=F)
roc_summary_test$specificity <- parse_number(roc_summary_test$specificity) #convert factor to number
roc_summary_test <- roc_summary_test %>% group_by(sensitivity) %>% summarize(specificity = max(specificity)) #summarize by max specificity per sensitivity
perfect_sensitivity <- c(1, 0) #add in data point for perfect sensitivity
perfect_specificity <- c(0, 1) #add in data point for perfect specificity
roc_summary_test <- rbind(roc_summary_test, perfect_sensitivity, perfect_specificity)
pdf("~/HA_Healthy_Asthma_Caret_RF_ROC_Curve_Testing_CAARS_ONLY_Normalized.pdf")
ggplot(data = roc_summary_test) + 
	geom_step(mapping = aes(x = (1-specificity), y = parse_number(roc_summary_test$sensitivity)), color = "blue", direction = "vh") + #generate roc curve
	geom_abline(intercept = 0, slope = 1) + #line of equal specificity and sensitivity
	labs(x = "False Positive Rate (1 - Specificity)", y = "True Positive Rate (Sensitivity)") + 
	scale_y_continuous(lim=c(0, 1))
	#scale_x_reverse(lim=c(1, 0)) #flip x axis
dev.off()

#generate mean type 1 and 2 error table
x <- cbind("Testing Type 1 Error Mean", mean(type_1_error_list))
y <- cbind("Testing Type 2 Error Mean", mean(type_2_error_list))
write.table(x, file="~/HA_Healthy_Asthma_Caret_RF_Output_Model_Stats1.CAARS_ONLY_Normalized.txt", quote=F, row.names=F)
write.table(y, file="~/HA_Healthy_Asthma_Caret_RF_Output_Model_Stats2.CAARS_ONLY_Normalized.txt", quote=F, row.names=F)
