library(ggplot2)
x <- read.csv("~/OTU_NZV_F-Scores.csv")
small_otus <- which(x$"f_score_record" < 0.1)
x <- x[-small_otus, ]
otu <- x$"colnames.data..2.ncol.data.."
otu <- gsub("\\..*","",otu) #remove everything after "OTU"
importance <- x$"f_score_record" 
formatted_name <- "ICS Use Among Asthmatics (Yes/No)"
plot <- ggplot(x, aes(x=reorder(otu, -importance), y=importance)) + 
geom_point() + 
	geom_smooth(method = "lm") + 
	#manually draw rectangle around selected otus
	theme(text = element_text(size=10)) + 
	theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
	labs(title = paste(formatted_name, "F-Score Importance Stripchart - Relative Abundance", sep = " "),
		x = "OTU", y="Mean F-Score") + 
	theme(plot.title = element_text(face="bold", hjust = 0.5, size = 10))
	#geom_rect(mapping = aes(xmin = "Otu161", xmax = "Otu104", ymin = 0.004, ymax = 0.019), color = 'blue', alpha = 0)
pdf("~/Asthma_ICS_Yes_No_F_Score_NZV_OTU_Stripchart.pdf")
plot(plot)
dev.off()