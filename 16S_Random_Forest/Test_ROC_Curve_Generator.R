require(caret)
require(ggplot2)
require(randomForest)
library(pROC)
library(vegan)
library(plyr)
library(dplyr)
library(tidyverse)
library(Bolstad2) #AUC calculations

roc_summary_test <- read.table('~/Gender_Caret_RF_Output_Model_ROC_CAARS_ONLY_Normalized.txt')
colnames(roc_summary_test) <- c("sensitivity", "specificity")
roc_summary_test$specificity <- parse_number(roc_summary_test$specificity) #convert factor to number
roc_summary_test <- roc_summary_test[-1, ] #remove column labels
roc_summary_test <- roc_summary_test %>% arrange(sensitivity) %>% group_by(sensitivity) %>% summarize(specificity = mean(specificity)) #summarize by max specificity per sensitivity
perfect_sensitivity <- c(1, 0) #add in data point for perfect sensitivity
perfect_specificity <- c(0, 1) #add in data point for perfect specificity
roc_summary_test <- rbind(roc_summary_test, perfect_sensitivity, perfect_specificity)

#AUC calculations
x_axis <- 1 - roc_summary_test$specificity
y_axis <- as.numeric(levels(roc_summary_test$sensitivity))[roc_summary_test$sensitivity]
auc <- sintegral(x_axis, y_axis)$int #integrate under ROC curve using Simpson method to get AUC

pdf("~/Gender_Caret_RF_ROC_Curve_Testing_CAARS_ONLY_Normalized.pdf")
ggplot(data = roc_summary_test) + 
	geom_step(mapping = aes(x = (1-specificity), y = parse_number(roc_summary_test$sensitivity)), color = "blue", direction = "vh") + #generate roc curve
	geom_abline(intercept = 0, slope = 1) + #line of equal specificity and sensitivity
	labs(x = "False Positive Rate (1 - Specificity)", y = "True Positive Rate (Sensitivity)") + 
	ggtitle(paste("AUC = ", auc, sep = '')) + 
	scale_y_continuous(lim=c(0, 1)) +
	theme(plot.title = element_text(hjust = 0.5, size = 12), #modify plot title
				legend.position = "none") #remove legend
dev.off()

