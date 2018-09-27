library(ggplot2)
x <- read.csv("~/T_FLAVAN3OLS_Caret_RF_Output_Model_OTU_Table_Packaged.CAARS_ONLY_Normalized.csv")
small_prob <- which(x[, "MDA"] < 0) #remove otus that make the model worse
x <- x[-small_prob, ]
otu <- x$"OTU"
otu <- gsub("\\..*","",otu) #remove everything after "OTU"
importance <- x$"MDA" #importance metric for feature ranking
formatted_name <- "Daily Total Flavan-3-ols (Above/Below Median)"
plot <- ggplot(x, aes(x=reorder(otu, -importance), y=importance)) + 
geom_point() + 
	geom_smooth(method = "lm") + 
	#manually draw rectangle around selected otus
	theme(text = element_text(size=7.5)) + 
	theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
	labs(title = paste(formatted_name, "Caret Selected OTU MDA Importance Stripchart - Relative Abundance", sep = " "),
		x = "OTU", y="Mean MDA") + 
	theme(plot.title = element_text(face="bold", hjust = 0.5, size = 7))
	#geom_rect(mapping = aes(xmin = "Otu161", xmax = "Otu104", ymin = 0.004, ymax = 0.019), color = 'blue', alpha = 0)
pdf("~/T_FLAVAN3OLS_Caret_RF_OTU_Stripchart.pdf")
plot(plot)
dev.off()