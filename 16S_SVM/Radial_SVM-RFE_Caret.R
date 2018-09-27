require(caret)
require(ggplot2)
library(pROC)
library(vegan)
library(plyr)
library(tidyverse)
library(stringr)

#prepare data set for rbf svm model generation
otu_matrix <- read.csv("~/Desktop/CAARS_Files/Random_Forest_V1/Data/DNA/PF_SEAFD_LOW_Caret_RF.dna.otu.caarsonly.normalized.csv")
otu_matrix_otus <- otu_matrix[, 3:ncol(otu_matrix)] #assign only otus to this matrix
sample_category <- otu_matrix[, 2] #assign environmental variable to this matrix
sample_category <- factor(sample_category)
data <- data.frame(sample_category, otu_matrix_otus) #generates otu matrix with environmental variable category as the first column

#Subsets features by mean relative abundance - call in command line with "M" and the #of features to truncate at
mean_rel_abundance_selector <- function(df, n) {
	otu_mean_relative_abundances <- colMeans(df[2:ncol(df)])
	otu_mean_relative_abundances <- sort(otu_mean_relative_abundances, decreasing = TRUE)
	otu_mean_relative_abundances_top_100 <- names(otu_mean_relative_abundances[1:n]) #vector of top 100 otus by mean relative abundance
	df <- df[, which(colnames(df) %in% otu_mean_relative_abundances_top_100)] #filter data
	return(df)
}

#Subsets features with near zero variance - call in command line with "NZV 0"
near_zero_var_selector <- function(df) {
	df <- df[, -nearZeroVar(df)] #remove otus that have a variance near 0, which slow algorithm time and are not usually predictive
	return(df)
}

#Selected feature selection scheme - passed in from command-line
args <- commandArgs(trailingOnly = TRUE) 
if (args[1] == "M") {
	data <- mean_rel_abundance_selector(data, args[2]) #select top 100 OTUs by mean relative abundance, subset data, 
	data <- cbind(sample_category, data) #bind to category
}
if (args[1] == "NZV") {
	data <- near_zero_var_selector(data) #select otus that vary significantly across samples by nearZeroVar
}

#Feature F-Scores
f_score_calculator <- function(data) { #sourced from Chen and Lin 2011 (same implementation as libsvm)
	f_score_record <- NULL
	for(i in 2:ncol(data)) { #iterate through all variables
		df <- data[, c(1, i)] #select variable of interest
		positives <- df[which(df$sample_category == 1), ]
		negatives <- df[which(df$sample_category == 0), ]
		mean_pos <- mean(positives[, 2]) #mean rel abundance for all positive patients
		mean_neg <- mean(negatives[, 2]) 
		sd_pos <- sd(positives[, 2])
		sd_neg <- sd(negatives[, 2])
		f_score <- abs(((mean_pos-mean(df[, 2]))^2 + (mean_neg-mean(df[, 2]))^2)/(sd_pos^2 + sd_neg^2))
		f_score_record <- c(f_score_record, f_score)
	}
	f_score_record <- data.frame(colnames(data)[2:ncol(data)], f_score_record)
	write.csv(f_score_record, file = "~/Desktop/CAARS_Files/16S_SVM_V1/Output/PF_SEAFD_LOW/OTU_NZV_F-Scores.csv", quote = F, row.names = F)
}
f_score_calculator(data)

roc_summary_test <- NULL
type_1_error_list <- data.frame(NULL)
type_2_error_list <- data.frame(NULL)

for (i in 1:100) {
	iteration_number <- i #reference iteration number for files

	#************#
	#Training Set#
	#************#

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

	levels(data_training$sample_category) <- make.names(c("Non_Obese", "Obese")) #manually set factor levels of independent variable

	#Training protocol for RBF model on Training Data Set with 10-fold cross-validation
	trcontrol <- trainControl(method="repeatedcv", #cross-validation
							  number=10, #10-fold cross validation
							  repeats = 3, #reduces overfitting on training set
							  classProbs = TRUE,
							  summaryFunction = twoClassSummary,
							  savePredictions = TRUE) #saves sensitivity and specificity measures for ROC curve
	#Train model on all predictors
	model <-train(sample_category~.,data=data_training, #Run model on all predictors
						method="svmRadial",
						trControl = trcontrol,
						prox=TRUE,
						tuneLength = 15, #select 15 random numbers of variables to try at each split (cost)
						saveDetails = TRUE,
						allowParallel=FALSE,
						metric = 'ROC', #select model based on highest AUC
						importance = TRUE)
	importance <- varImp(model, scale = FALSE)$importance[1] #predictors and their respective importance measure (area under ROC curve)
	
	#Recursive Feature Elimination
	feature_subset_sizes <- c(100, 50, 25, 20, 15, 10, 5, 1) #list of subsets of features to test

	models <- vector("list", length(feature_subset_sizes) + 1) #initialize list to store models based on number of feature subsets
	models[[1]] <- model #initialize with current (naive) model
	max_model_metric <- max(model$results$ROC) #store best model's max AUC - initialize with naive model
	optimal_predictors <- ncol(data_training) - 1 #initiate to naive model
	j <- 1 #counter for number of times RFE has been performed (initialize to 1 to count for naive model)
	
	for (i in feature_subset_sizes) { #iterate through list of feature subset sizes
		j <- j + 1
		#Select features
		prior_importances <- data.frame(importance, OTU = rownames(importance)) #convert rownames into a column
		prior_importances <- prior_importances[order(prior_importances$Non_Obese, decreasing = TRUE)[1:i], ] #select top i predictors by ROC score
		feature_subset <- data_training[, c(1, which(colnames(data_training) %in% prior_importances$OTU))] #select top n features
		
		rfe_model <-train(sample_category~.,data=feature_subset, 
						method="svmRadial",
						trControl = trcontrol,
						prox=TRUE,
						tuneLength = 15, #select 15 random numbers of variables to try at each split (cost)
						saveDetails = TRUE,
						allowParallel=FALSE,
						metric = 'ROC', #select model based on highest AUC
						importance = TRUE)
		
		importance <- varImp(model, scale = FALSE)$importance[1] #update feature importances
		models[[j]] <- rfe_model #store model in list for later access
		if (max(rfe_model$results$ROC) > max_model_metric) {
			best_model_element <- j #store list element corresponding to current best model
			max_model_metric <- max(rfe_model$results$ROC) #replace best model AUC with current model max AUC
			optimal_predictors <- ncol(feature_subset) - 1 #replace with current model feature count
		}
	}
	if (max_model_metric == max(model$results$ROC)) { #if none of the subsets yielded a better model than naive
		optimal_model <- models[[1]] 
	}
	else {
		optimal_model <- models[[best_model_element]] #retrieve the best model from RFE
	}
	#Type Errors from Training Set
	training_type_1_error <- 1 - optimal_model$results[which(optimal_model$results$ROC == max(optimal_model$results$ROC)), ]$Spec
	if (training_type_1_error < 0) {
		training_type_1_error <- 0 #adjust for rounding in R
	}
	training_type_2_error <- 1 - optimal_model$results[which(optimal_model$results$ROC == max(optimal_model$results$ROC)), ]$Sens
	if (training_type_2_error < 0) {
		training_type_2_error <- 0 #adjust for rounding in R
	}

	#***********#
	#Testing Set#
	#***********#

	unsampled_cases <- setdiff(cases, cases_sample) #generate list of cases for testing data set not in the training sample
	unsampled_controls <- setdiff(controls, controls_sample) #generate list of controls for testing data set not in the training sample
	testing_data_cases <- data[unsampled_cases, ] #subsets bmi data set by cases
	testing_data_controls <- data[unsampled_controls, ] #subsets bmi data set by controls
	data_testing <- rbind(testing_data_cases, testing_data_controls) #combines cases and controls to make a testing data set
	levels(data_testing$sample_category) <- c("Non_Obese", "Obese")

	#Test Optimal Model
	test <- predict(optimal_model, data_testing, type="prob", norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE) #test random forest optimal_model on remaining 3/4 of data
	test_2 <- predict(optimal_model, data_testing) #just output predictions not probabilities
	confusion_matrix <- t(confusionMatrix(test_2, data_testing$sample_category)$table)
	print(confusion_matrix)
	print(c("Training Type 1 Error", training_type_1_error))
	print(c("Training Type 2 Error", training_type_2_error))
	testing_type_1_error <- confusion_matrix[[3]]/(confusion_matrix[[1]] + confusion_matrix[[3]])
	testing_type_1_error <- c("Testing Type 1 Error", testing_type_1_error) #label in output for selection by shell script
	type_1_error_list <- rbind(type_1_error_list, testing_type_1_error) #append to list of type 1 errors
	testing_type_2_error <-  confusion_matrix[[2]]/(confusion_matrix[[2]] + confusion_matrix[[4]])
	testing_type_2_error <- c("Testing Type 2 Error", testing_type_2_error) #label in output for selection by shell script
	type_2_error_list <- rbind(type_2_error_list, testing_type_2_error)

	#print out optimal_model statistics
	print(optimal_model$finalModel) #outputs final optimal model statistics (type of RF, ntree, mtry, and confusion matrix)
	print(c("Max AUC", max(optimal_model$results$ROC))) #print maximum AUC (used to select best-fit optimal_model)
	print(c("Optimal Number of Predictors", ncol(optimal_model$trainingData) - 1)) #print optimal number of predictors based on RFE
	print(optimal_model) #prints accuracy measures and info on each resampling
	importance <- varImp(optimal_model, scale = FALSE) #outputs predictors and their respective importance measure (area under ROC curve)
	print(importance$importance[1]) #print all variables from RFE updated model and their importance for class discrimination

	#Testing set ROC data
	probs <- test[, 2] #positives among testing set
	roc <- roc(data_testing$sample_category~probs) #generate roc curve
	roc_summary_test<- rbind(roc_summary_test, data.frame(sensitivity = roc$sensitivities, specificity = roc$specificities))
	
	#generate AUC curve for optimal model
	file_path <- paste("~/Desktop/CAARS_Files/16S_SVM_V1/Output/AUC_Plots/PF_SEAFD_LOW_Caret_SVM_Plot_IT", iteration_number, "_CAARS_ONLY_Normalized.pdf", sep="")
	pdf(file_path)
	plot <- plot(optimal_model, main = "Radial RBF SVM Model Variable Selection for PF_SEAFD_LOW")
	print(plot)
	dev.off()
}

#Cumulative ROC Curve
colnames(roc_summary_test) <- c("sensitivity", "specificity")
roc_summary_test <- as.data.frame(roc_summary_test)
write.table(roc_summary_test, file="~/Desktop/CAARS_Files/16S_SVM_V1/Output/PF_SEAFD_LOW_Caret_SVM_ROC_Testing_CAARS_ONLY_Normalized.txt", quote=F, row.names=F)
roc_summary_test$specificity <- parse_number(roc_summary_test$specificity) #convert factor to number
roc_summary_test <- roc_summary_test %>% group_by(sensitivity) %>% summarize(specificity = max(specificity)) #summarize by max specificity per sensitivity
perfect_sensitivity <- c(1, 0) #add in data point for perfect sensitivity
perfect_specificity <- c(0, 1) #add in data point for perfect specificity
roc_summary_test <- rbind(roc_summary_test, perfect_sensitivity, perfect_specificity)
pdf("~/Desktop/CAARS_Files/16S_SVM_V1/Output/AUC_Plots/PF_SEAFD_LOW_Caret_SVM_ROC_Curve_Testing_CAARS_ONLY_Normalized.pdf")
ggplot(data = roc_summary_test) + 
	geom_step(mapping = aes(x = (1-specificity), y = parse_number(roc_summary_test$sensitivity)), color = "blue", direction = "vh") + #generate roc curve
	geom_abline(intercept = 0, slope = 1) + #line of equal specificity and sensitivity
	labs(x = "False Positive Rate (1 - Specificity)", y = "True Positive Rate (Sensitivity)") + 
	scale_y_continuous(lim=c(0, 1))
	#scale_x_reverse(lim=c(1, 0)) #flip x axis
dev.off()
#Cumulative Type Errors
print(type_1_error_list)
print(type_2_error_list)