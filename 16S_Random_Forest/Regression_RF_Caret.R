require(caret)
require(ggplot2)
require(randomForest)
library(pROC)
library(vegan)
library(plyr)
library(dplyr)

#prepare data set for rf model generation
otu_matrix <- read.csv("IL_5_Continuous_Caret_RF.dna.otu.caarsonly.normalized.csv")
otu_matrix_otus <- otu_matrix[, 3:ncol(otu_matrix)] #assign only otus to this matrix
response <- otu_matrix[, 2] #assign environmental variable to this matrix
#sample_category <- factor(sample_category)
data <- data.frame(response, otu_matrix_otus) #generates otu matrix with environmental variable category as the first column
data <- data[, -nearZeroVar(data)] #remove otus that have a variance near 0, which slow algorithm time and are not usually predictive
data <- data %>% arrange(., response) #rearrange in ascending order of response to promote proper rf grouping

roc_summary <- NULL
model_summary <- NULL
type_1_error_list <- NULL
type_2_error_list <- NULL
predicted <- c()

data_names <- rownames(data)
for (i in 1:100) {
	iteration_number <- i #reference iteration number for files
	training_prop <- 0.8 #proportion of total data set to be used for training
	training_samples_count <- training_prop*(nrow(data)) #number of samples to be used in training data set
	training_sample <- sample(data_names, training_samples_count) #select training cases randomly
	data_training <- data[training_sample, ] 

	set.seed(20) #reproducibility
	#Run random forest model on training set along with 10-fold cross validation
	trcontrol <- trainControl(method="repeatedcv", #cross-validation
							  number=10, #10-fold cross validation
							  repeats = 3, #reduces overfitting on training set
							  savePredictions = TRUE) #saves sensitivity and specificity measures for ROC curve
	#mtry <- sqrt(ncol(data_training))
	model <-train(response~.,data=data_training, 
		method="rf",
		trControl = trcontrol,
		prox=TRUE,
		tuneLength = 15, #select 15 random numbers of variables to try at each split (mtry)
		saveDetails = TRUE,
		allowParallel=FALSE,
		importance = TRUE)

	#generate testing set for model to test on
	unsampled <- setdiff(data_names, training_sample)
	data_testing <- data[unsampled, ] #selects remaining data for testing
	test <- predict(model, data_testing, predict.all=FALSE, proximity=TRUE, nodes=FALSE, se.fit = TRUE) #test random forest model on remaining 3/4 of data set
	
	test_results <- cbind(data_testing[, "response"], test)
	print(test_results)
	predicted <- rbind(predicted, test_results)

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
	write.table(importance, file="IL_5_Continuous_Caret_Continuous_Output_Model_Importance.CAARS_ONLY_Normalized.txt", quote=F, row.names=F)
	sig_importances <- arrange(importance, desc(MDA)) #sort OTUs by descending magnitude to subset important ones
	sig_importances <- head(sig_importances, 25) #select out top 25 OTUs (most important)
	print(sig_importances) #print importances for top 25 OTUs

	#generate AUC curve for model
	file_path <- paste("AUC_Plots/IL_5_Continuous_Continuous_Caret_Continuous_Plot_IT", iteration_number, "_CAARS_ONLY_Normalized.pdf", sep="")
	pdf(file_path)
	plot <- plot(model, main = "Random Forest Model Variable Selection for IL_5_Continuous Continuous")
	print(plot)
	dev.off()
}
#output mean residuals plot
colnames(predicted) <- c("Actual", "Predicted")
predicted <- data.frame(predicted)
aggregated <- ddply(predicted, 'Actual', summarize, Mean_Predicted = mean(Predicted)) #aggregate mean of predicted IL_5_Continuous_Continuous for each level of actual IL_5_Continuous_Continuous
write.table(aggregated, file = "IL_5_Continuous_Continuous_Caret_Prediction_vs._Actual.txt")
#write.table(predicted, file = "IL_5_Continuous_Continuous_Caret_Prediction_vs._Actual.txt")
pdf("IL_5_Continuous_Continuous_Caret_Prediction_vs._Actual_Plot.pdf")
ggplot(data = aggregated) +
geom_point(mapping = aes(x = Actual, y = Mean_Predicted), color = "blue") +
geom_smooth(mapping = aes(x = Actual, y = Mean_Predicted), method = "lm") +
geom_abline(intercept = 0, slope = 1) + #line of equal specificity and sensitivity
scale_x_continuous(limits = c(min(aggregated$Actual), max(aggregated$Actual))) + 
scale_y_continuous(limits = c(min(aggregated$Actual), max(aggregated$Actual))) + 
ggtitle(paste("Plot of Predicted Residuals over 100 iterations - ", "IL_5_Continuous_Continuous", sep = " ")) + 
labs(x = paste("Actual", "IL_5_Continuous_Continuous", sep = " "), y = paste("Mean Predicted", "IL_5_Continuous_Continuous", "over 100 RF Models", sep = " "))
dev.off()

