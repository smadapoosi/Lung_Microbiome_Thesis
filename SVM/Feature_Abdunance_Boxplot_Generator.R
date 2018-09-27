library(ggplot2)
library(stringr)
library(gridExtra) #for arranging grids
#library(ggsignif)
library(ggpubr)

env <- read.csv("~/Asthma_ICS_Yes_No_Caret_RF.dna.otu.caarsonly.normalized.csv", header=T) #env table with first-pass adonis significant variables log-transformed

variable <- "Asthma_ICS_Yes_No" #selected variable to view levels
formatted_variable_name <- "ICS Use Among Asthmatics (Yes/No)" #variable name for display on plot

otu <- c("Otu0008.Firmicutes.Bacilli.100..Lactobacillales.Streptococcaceae.Streptococcus")

myplot <- function(i) { #plots individual boxplots using ggplot
	df <- data.frame(factor(env[, 2]), env[, i]) #create new data frame with just variable and otu
	colnames(df) <- c("variable", "otu_relative_abundance")
	levels(df$"variable") <- c("No", "Yes") #levels of variable
	stat.test <- compare_means(otu_relative_abundance ~ variable, data = df, method = "t.test", p.adjust.method = "BH")
	p_val <- signif(p.adjust(stat.test$p.adj, method = "BH", n = 112), digits = 3) #adjust p-value based on all features in model
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
	asterisk_naive <- asterisk_selection(stat.test$p.adj)
	asterisk <- asterisk_selection(p_val) #determine which significance indicator to be selected
	#Generate boxplot
	pdf(paste("~/Asthma_ICS_Yes_No/", variable, "_", i, ".pdf", sep=""))
	p <- ggplot(data = df) +
		geom_boxplot(mapping = aes(x = variable, y = otu_relative_abundance, fill = variable)) +
		annotate("text", -Inf, Inf, label = paste("P-Value = ", signif(stat.test$p.adj, digits = 3), " (", asterisk_naive, ")", sep = ''), hjust = 0, vjust = 1.5, size = 5, color = 'red') + 
		annotate("text", -Inf, Inf, label = paste("BH-Corrected P-Value = ", p_val, " (", asterisk, ")", sep = ''), hjust = 0, vjust = 3, size = 5, color = 'red') + 
		labs(x = formatted_variable_name, y = "Relative Abundance") +
		ggtitle(str_wrap(paste(i, " Relative Abundance by ", formatted_variable_name, sep=""), width = 40)) + #plot title
		scale_fill_manual(values = c("blue", "yellow")) +
		theme(plot.title = element_text(hjust = 0.5, size = 15), #modify plot title
				axis.title = element_text(face="bold", size = 15), #set axes to bold
				axis.text.x = element_text(size = 10),
				axis.text.y = element_text(size = 10),
				legend.position = "none") #remove legend
	plot(p) #plot individual graph
	return(p)
}
p <- lapply(otu, myplot) #call function on all otus of interest
pdf(paste("~/", variable, "_RFE_All_Hits.pdf", sep=""), width = 10, height = 10)
do.call(grid.arrange, c(p, ncol = 2)) #arrange plots in a grid
dev.off()