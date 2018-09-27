library(tidyverse)
library(stringr)
library(magrittr)
otus <- read.csv("~/sputDNA.lowsquam.WGS.species.caarsonly.top_100.csv")
rownames(otus) <- otus[, 1]
otus <- otus[, -1]

#Gather list of patients in each K-means group
group_1_ids <- read.csv('Group_1.csv') %>% select('X') %>% head(3) %$% X %>% str_sub(., 1, 3) #group 1
group_2_ids <- read.csv('Group_2.csv') %>% select('X') %>% head(6) %$% X %>% str_sub(., 1, 3) #group 2
group_3_ids <- read.csv('Group_3.csv') %>% select('X') %>% head(12) %$% X %>% str_sub(., 1, 3) #group 3

bacteria_anova <- function(otus) { #outputs anova results

	group_1_otus <- otus[which(rownames(otus) %in% group_1_ids), ]
	group_2_otus <- otus[which(rownames(otus) %in% group_2_ids), ]
	group_3_otus <- otus[which(rownames(otus) %in% group_3_ids), ]

	#************************************#
	#Looking at broader taxa (genus-level)
	#************************************#

	colnames(group_1_otus) <- gsub("\\..*","", colnames(group_1_otus)) #remove everything related to species
	g1_aggr <- as.data.frame(do.call(cbind, by(t(group_1_otus),INDICES=names(group_1_otus),FUN=colSums)))
	rownames(g1_aggr) <- group_1_ids

	#Looking at broader taxa (genus-level)
	colnames(group_2_otus) <- gsub("\\..*","", colnames(group_2_otus)) #remove everything related to species
	g2_aggr <- as.data.frame(do.call(cbind, by(t(group_2_otus),INDICES=names(group_2_otus),FUN=colSums)))
	rownames(g2_aggr) <- group_2_ids

	#Looking at broader taxa (genus-level)
	colnames(group_3_otus) <- gsub("\\..*","", colnames(group_3_otus)) #remove everything related to species
	g3_aggr <- as.data.frame(do.call(cbind, by(t(group_3_otus),INDICES=names(group_3_otus),FUN=colSums)))
	rownames(g3_aggr) <- group_3_ids

	g1 <- g1_aggr[,order(colMeans(g1_aggr))] #sort in ascending order by mean relative abundance at genus level
	g2 <- g2_aggr[,order(colMeans(g2_aggr))] #sort in ascending order by mean relative abundance at genus level
	g3 <- g3_aggr[,order(colMeans(g3_aggr))] #sort in ascending order by mean relative abundance at genus level
	#print(g1)
	#print(g2)
	#print(g3)


	group_updated <- rbind(g1_aggr, g2_aggr, g3_aggr) #add in data frames
	group_updated <- group_updated %>% mutate(group = ifelse(rownames(group_updated) %in% group_1_ids, 1, 
											   			ifelse(rownames(group_updated) %in% group_2_ids, 2, 3)))

	species_anova <- function(data) { #pass in otu table with k-means groups as last column, output which species differ by group using ANOVA
		species_names <- colnames(data)[1:ncol(data) - 1] #list all bacterial species in otu table
		output_list <- list(NULL)
		p_val_list <- c(NULL)
		j <- 0
		for (i in species_names) {
			model <- aov(data[[i]]~as.factor(data[, ncol(data)])) #run anova for selected species comparing all groups
			means <- model.tables(model, "means")
			p_val_list <- c(p_val_list, summary(model)[[1]][["Pr(>F)"]][1])
			print(summary(model)[[1]])
			#if (p.adjust(summary(model)[[1]][["Pr(>F)"]][1], "bonferroni", length(species_names)) #bonferroni correction given number of tests based on number of OTUs
			if (summary(model)[[1]][["Pr(>F)"]][1]
				 < 0.05) { #ANOVA test using bonferroni correction if if otu rel abundance is significantly different among groups
				j <- j + 1 #counter for list elements
				output_list[[i]] <- c(summary(model)[[1]], means) #append to list
			}
		}
		df <- tibble(species_names, p.adjust(p_val_list, "holm", length(p_val_list)), #holm p-values
									p.adjust(p_val_list, "BH", length(p_val_list)), #benjamini-hochberg p-values
									p.adjust(p_val_list, "bonferroni", length(p_val_list))) #bonferroni p-values
		colnames(df) <- c("Species", "Holm P", "BH P", "Bonferroni P")
		write.csv(df, file = "Genus_P_Vals.csv", quote = F, row.names = F)
		print(output_list) #bonferroni p-values
	}
	species_anova(group_updated) #output anova results of species
}
bacteria_anova(otus)


