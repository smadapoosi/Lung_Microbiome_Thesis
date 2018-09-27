library(ggpubr)
library(ggplot2)
library(car)
library(dunn.test) #Dunn post-hoc test
source('games_howell.R') # Games-Howell post hoc test https://gist.github.com/aschleg/ea7942efc6108aedfa9ec98aeb6c2096#file-games_howell-r

rf <- read.table("~/T_FLAVAN3OLS_Caret_RF_Output_Model_AUC_CV.CAARS_ONLY_Normalized.txt")
rf <- cbind(rf, rep("RF", n = nrow(rf)))
colnames(rf) <- c("AUC", "Model")
svm <- read.table("~/T_FLAVAN3OLS_Caret_SVM_Output_Model.CAARS_ONLY_Normalized.txt.AUCs.txt")
svm <- cbind(svm, rep("SVM", n = nrow(svm))) #add label to these
colnames(svm) <- c("AUC", "Model")
svm_rfe <- read.table("~/T_FLAVAN3OLS_Caret_SVM_RFE_Output_Model.CAARS_ONLY_Normalized.txt.AUCs.txt")
svm_rfe <- cbind(svm_rfe, rep("SVM-RFE", n = nrow(svm_rfe))) #add label to these
colnames(svm_rfe) <- c("AUC", "Model")
auc_table <- rbind(rf, svm, svm_rfe) #aggregate data table with column 1 being method and column 2 being AUCs

asterisk_selection <- function(n) {
	label <- NULL
	if(n > 0.05) {
		label <- "ns" #relative abundance difference not significant
	}
	else if (n <= 0.05 && n > 0.01) {
		label <- "*"
	}
	else if (n <= 0.01 && n > 0.001) {
		label <- "**"
	}
	else if (n <= 0.001 && n > 0.0001) {
		label <- "***"
	}
	else if (n <= 0.0001) {
		label <- "****"
	}
	return(label)
}

test_name <- "ANOVA"
stat.test <- aov(AUC ~ Model, data = auc_table) #perform anova
stat.test.p.value <- summary(aov(AUC ~ Model, data = auc_table))[[1]][["Pr(>F)"]][[1]]
p_val <- signif(p.adjust(stat.test.p.value, method = "BH", n = 8), digits = 3) #adjust p-value based on all features in model

asterisk_naive <- asterisk_selection(stat.test.p.value)
asterisk <- asterisk_selection(p_val) #determine which significance indicator to be selected

post.hoc.test <- TukeyHSD(stat.test) #post hoc Tukey HSD test
post.hoc.test.p <- data.frame("Comparison" = c("SVM-RF","SVM-RFE-RF", "SVM-RFE-SVM"),
							  "p.value" = c(post.hoc.test$Model[[10]], post.hoc.test$Model[[11]], post.hoc.test$Model[[12]]))

annot_1 <- post.hoc.test.p[[2]][[1]]
annot_2 <- post.hoc.test.p[[2]][[2]]
annot_3 <- post.hoc.test.p[[2]][[3]]
asterisk_annot_1 <- asterisk_selection(annot_1)
asterisk_annot_2 <- asterisk_selection(annot_2)
asterisk_annot_3 <- asterisk_selection(annot_3)

#Check ANOVA assumptions and redo using other test if failed
if (leveneTest(AUC ~ Model, auc_table)$'Pr(>F)'[[1]] < 0.05) { #if it fails the homoscedacity assumption
	if (shapiro.test(x = residuals(stat.test))$p.value < 0.05) { #if it fails the normality of residuals assumption as well
		stat.test <- kruskal.test(AUC ~ Model, auc_table) #use Kruskal-Wallis nonparametric test instead
		stat.test.p.value <- stat.test$p.value
		test_name <- "Kruskal-Wallis" #update test name for plot
		p_val <- signif(p.adjust(stat.test$p.value, method = "BH", n = 8), digits = 3) #adjust p-value based on all features in model
		asterisk_naive <- asterisk_selection(stat.test$p.value)
		asterisk <- asterisk_selection(p_val) #determine which significance indicator to be selected
		post.hoc.test <- dunn.test(auc_table$AUC, auc_table$Model)$P #run post-hoc Dunn test
		post.hoc.test.p <- data.frame("Comparison" = c("SVM-RF","SVM-RFE-RF", "SVM-RFE-SVM"),
							  "p.value" = post.hoc.test)
		annot_1 <- post.hoc.test.p[[2]][[1]]
		annot_2 <- post.hoc.test.p[[2]][[2]]
		annot_3 <- post.hoc.test.p[[2]][[3]]
		asterisk_annot_1 <- asterisk_selection(annot_1)
		asterisk_annot_2 <- asterisk_selection(annot_2)
		asterisk_annot_3 <- asterisk_selection(annot_3)
	}
	else {
		stat.test <- oneway.test(AUC ~ Model, auc_table) #Welch One-Way Test with No Homoscedacity Assumption
		test_name <- "Welch" 
		stat.test.p.value <- stat.test$p.value
		p_val <- signif(p.adjust(stat.test$p.value, method = "BH", n = 8), digits = 3) #adjust p-value based on all features in model
		asterisk_naive <- asterisk_selection(stat.test$p.value)
		asterisk <- asterisk_selection(p_val) #determine which significance indicator to be selected
		post.hoc.test <- games.howell(auc_table$Model, auc_table$AUC)$p #run post-hoc Games-Howell test
		post.hoc.test.p <- data.frame("Comparison" = c("SVM-RF","SVM-RFE-RF", "SVM-RFE-SVM"),
							  "P-Value" = post.hoc.test)
	}
}

print(post.hoc.test.p)

#boxplot significance positions
df1 <- data.frame(a = c(1, 1, 2, 2), b = c(1.1, 1.25, 1.25, 1.1)) #location of brackets
df2 <- data.frame(a = c(1, 1:3,3), b = c(1.25, 1.30, 1.30, 1.30, 1.25)) #location of brackets (middle comparison)
df3 <- data.frame(a = c(2, 2, 3, 3), b = c(1.1, 1.25, 1.25, 1.1)) #location of brackets

#Generate boxplot
pdf(paste("~/T_FLAVAN3OLS_AUCs_Boxplot.pdf", sep=""))
p <- ggplot(data = auc_table) +
	geom_boxplot(mapping = aes(x = Model, y = AUC, fill = Model)) +
 	geom_line(data = df1, aes(x = a, y = b)) + annotate("text", x = 1.5, y = 1.25, label = asterisk_annot_1, size = 8, color = "red") +
    geom_line(data = df2, aes(x = a, y = b)) + annotate("text", x = 2, y = 1.30, label = asterisk_annot_2, size = 8, color = "red") +
    geom_line(data = df3, aes(x = a, y = b)) + annotate("text", x = 2.5, y = 1.25, label = asterisk_annot_3, size = 8, color = "red") + 
	annotate("text", x = 0.25, y = 1.4, label = paste(test_name, " P-Value = ", signif(stat.test.p.value, digits = 3), " (", asterisk_naive, ")", sep = ''), hjust = 0, vjust = 1, size = 4, color = 'red') + 
	annotate("text", x = 0.25, y = 1.39, label = paste("BH-Corrected P-Value = ", p_val, " (", asterisk, ")", sep = ''), hjust = 0, vjust = 2.5, size = 4, color = 'red') + 
	labs(x = "Classification Method", y = "Cross-Validated Training AUC") +
	ggtitle("Comparison of Daily Total Flavan-3-ols (Above/Below Median) Model AUCs by Classification Method") + #plot title
	scale_fill_manual(values = c("blue", "yellow", "orange")) +
	theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 8), #modify plot title
			axis.title = element_text(face="bold", size = 10), #set axes to bold
			legend.position = "none") #remove legend
plot(p)