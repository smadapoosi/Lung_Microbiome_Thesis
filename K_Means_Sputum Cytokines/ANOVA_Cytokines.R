library(tidyverse)
library(stringr)
library(magrittr)

#Gather list of patients in each K-means group
group_1_ids <- read.csv('Group_1.csv') %>% select('X') %>% head(3) %$% X %>% str_sub(., 1, 3) #group 1
group_2_ids <- read.csv('Group_2.csv') %>% select('X') %>% head(6) %$% X %>% str_sub(., 1, 3) #group 2
group_3_ids <- read.csv('Group_3.csv') %>% select('X') %>% head(12) %$% X %>% str_sub(., 1, 3) #group 3

otus <- cytokines

rownames(otus) <- str_sub(rownames(otus), 1, 3) #remove V1 tag from rownames

group_1_otus <- otus[which(rownames(otus) %in% group_1_ids), ]
group_1_otus <- group_1_otus %>% mutate(group = 1)
group_2_otus <- otus[which(rownames(otus) %in% group_2_ids), ]
group_2_otus <- group_2_otus %>% mutate(group = 2)
group_3_otus <- otus[which(rownames(otus) %in% group_3_ids), ]
group_3_otus <- group_3_otus %>% mutate(group = 3)

group_updated <- rbind(group_1_otus, group_2_otus, group_3_otus)

#CYTOKINE DIFFERENCES BETWEEN GROUPS
cyt_anova <- function(data) { #pass in cytokines table with k-means groups as last column, output which species differ by group using ANOVA
		species_names <- colnames(data)[1:ncol(data) - 1] #list all bacterial species in otu table
		output_list <- list(NULL)
		p_val_list <- c(NULL)
		j <- 0
		for (i in species_names) {
			model <- aov(log10(data[[i]])~as.factor(data[, ncol(data)])) #run anova for selected species comparing all groups
			means <- model.tables(model, "means")
			p_val_list <- c(p_val_list, summary(model)[[1]][["Pr(>F)"]][1])
			if (summary(model)[[1]][["Pr(>F)"]][1] #bonferroni correction given number of tests based on number of OTUs
				 < 0.1) { #ANOVA test using bonferroni correction if if otu rel abundance is significantly different among groups
				j <- j + 1 #counter for list elements
				output_list[[i]] <- c(summary(model)[[1]], means) #append to list
			}
		}
	df <- tibble(species_names, p.adjust(p_val_list, "holm", length(p_val_list)), #holm p-values
									p.adjust(p_val_list, "BH", length(p_val_list)), #benjamini-hochberg p-values
									p.adjust(p_val_list, "bonferroni", length(p_val_list))) #bonferroni p-values
	colnames(df) <- c("Species", "Holm P", "BH P", "Bonferroni P")
	write.csv(df, file = "Cyt_P_Vals.csv", quote = F, row.names = F)
	print(output_list) #bonferroni p-values
}
cyt_anova(group_updated)