
#ENVIRONMENTAL VARIABLE DIFFERENCES
env <- read.csv('~/sputDNA.lowsquam.WGS.env.caarsonly.csv')
rownames(env) <- env[, 1]
env <- env[, -1]

rownames(env) <- str_sub(rownames(env), 1, 3) #remove V1 tag from rownames

group_1_data <- env[which(rownames(env) %in% group_1_ids), ]
group_1_data <- group_1_data %>% mutate(group = 1)
group_2_data <- env[which(rownames(env) %in% group_2_ids), ]
group_2_data <- group_2_data %>% mutate(group = 2)
group_3_data <- env[which(rownames(env) %in% group_3_ids), ]
group_3_data <- group_3_data %>% mutate(group = 3)
#group_4_data <- env[which(rownames(env) %in% group_4_ids), ]
#group_4_data <- group_4_data %>% mutate(group = 4)

group_updated <- rbind(group_1_data, group_2_data, group_3_data)

#Features that are selected by K-Means PCA
lists <- c("Age_asthmadiag",
			  "Methpc20_log10",
			  "SUP_EPA",
			  "SUP_OM_6",
			  "SUP_OM_3",
			  "SUP_OLEIC",
			  "SUP_ALA",
			  "SUP_VITE",
			  "SUP_VITD",
			  "SUP_VITA",
			  "A_DRINKS",
			  "DT_ALCO",
			  "Sputum..eos",
			  "absol_eos..K.uL.",
			  "SNQscore_avg")

log10_list <- c("Age", 
				"ige_total", 
				"eos_percent",
				"Sputum..neutr",
				"DT_KCAL",
				"DT_PROT",
				"DT_CARB",
				"DT_TFAT",
				"DT_SUG_T",
				"DT_FIBE",
				"DT_MOIS",
				"DT_SFAT",
				"DT_MFAT",
				"DT_PFAT",
				"DT_CHOL",
				"DT_CAFFN",
				"QUERCE",
				"T_ANTHOCYADNS",
				"T_FLAVAN3OLS",
				"T_FLAVANONES",
				"T_FLAVONES",
				"T_FLAVONOLS",
				"T_ISOFLAVONES",
				"T_FLAVONOIDS",
				"GI",
				"F_TOTAL",
				"V_TOTAL",
				"G_TOTAL",
				"PF_TOTAL",
				"PF_MPS_TOTAL",
				"PF_MEAT",
				"PF_SEAFD_HI",
				"PF_SEAFD_LOW",
				"D_TOTAL",
				"ADD_SUGARS",
				"KCAL_EXPENDITURE_ALL",
				"METMINS",
				"TRP_G",
				"BETN_C",
				"DT_FIBER_INSOL",
				"DT_FIBER_SOL",
				"DT_PROT_ANIMAL",
				"DT_PROT_VEGETABLE",
				"bmi_v1",
				"fev1_v1",
				"fev1_percent_v1",
				"Shannon")

#Functions to run ANOVA
env_anova <- function(data, list) {
	#list <- list
	p_val_list <- c(NULL)
	for (i in list) { #iterate through intended list of values
		model <- aov(data[[i]]~as.factor(data[, ncol(data)])) #run anova for selected species comparing all groups
		means <- model.tables(model, "means")
		print(i)
		print(means) #print in-group means for recording
		print(summary(model)[[1]][["Pr(>F)"]][1])
		p_val_list <- c(p_val_list, summary(model)[[1]][["Pr(>F)"]][1]) #append to list of p-values
	}
	df <- tibble(lists, p.adjust(p_val_list, "holm", length(p_val_list)), #holm p-values
							p.adjust(p_val_list, "BH", length(p_val_list)), #benjamini-hochberg p-values
							p.adjust(p_val_list, "bonferroni", length(p_val_list))) #bonferroni p-values
	colnames(df) <- c("Species", "Holm P", "BH P", "Bonferroni P")
	write.csv(df, file = "Env_P_Vals.csv", quote = F, row.names = F)
}
env_anova_log10 <- function(data, list) { #pass in data frame of environmental variables and list of variables to run ANOVA on
	#list <- list
	p_val_list <- c(NULL)
	for (i in list) { #iterate through intended list of values
		model <- aov(log10(data[[i]])~as.factor(data[, ncol(data)])) #run anova for selected species comparing all groups
		model_1 <- aov(data[[i]]~as.factor(data[, ncol(data)]))
		means <- model.tables(model_1, "means")
		print(i)
		print(means) #print in-group means for recording
		print(summary(model)[[1]][["Pr(>F)"]][1])
		p_val_list <- c(p_val_list, summary(model)[[1]][["Pr(>F)"]][1]) #append to list of p-values
	}
	df <- tibble(log10_list, p.adjust(p_val_list, "holm", length(p_val_list)), #holm p-values
							p.adjust(p_val_list, "BH", length(p_val_list)), #benjamini-hochberg p-values
							p.adjust(p_val_list, "bonferroni", length(p_val_list))) #bonferroni p-values
	colnames(df) <- c("Species", "Holm P", "BH P", "Bonferroni P")
	write.csv(df, file = "~/Env_Log10_P_Vals.csv", quote = F, row.names = F)
}
#call env_anova on selected variables
env_anova_log10(group_updated, log10_list)
env_anova(group_updated, lists)