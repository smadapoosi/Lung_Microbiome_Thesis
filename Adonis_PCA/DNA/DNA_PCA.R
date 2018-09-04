library(vegan)
library(dplyr)
env <- read.csv("~/sputDNA.lowsquam.env.V1.caarsonly.csv", header=T)
otus <- read.csv("~/sputDNA.lowsquam.otu.V1.caarsonly.csv", header=T)
otus <- otus[, -1] #remove sample id column from otu file
rownames(env) <- env[, 1]
env <- env[, -1]

sig_traits <- read.csv("~/DNA_Adonis_Matched_PVals-sorted.csv") #gather list of traits and their p-values
sig_traits <- subset(sig_traits, sig_traits$V1 <= 0.05) #subset to just the significant traits
colnames(sig_traits) <- c("Variable", "Adonis P-Value")
sig_variables <- sig_traits$"Variable"

true_names <- read.csv("~/RNA_Adonis_LT_Variables_Fullnames.csv", header=T) 

#for loop to convert numerical traits to categorical with respect to the median
for(sample in sig_variables) {
	if (is.numeric(env[1, sample]) == TRUE) { #determines if trait is continuous and needs to be converted to categorical
		sample_true_name <- true_names[true_names$x == sample, "full_name"]
		sample_units <- true_names[true_names$x == sample, "units"]
		adonis_p_value <- sig_traits[sig_traits$"Variable" == sample, "Adonis P-Value"]
		titleplot <- sample
		mypath <- file.path(paste("~/DNA_V1_Numerical_Hist_Post_LT_", titleplot, ".pdf", sep = ""))
		pdf(file=mypath) #sets output path for the plot
		if (sample == "BMI.category"){
			hist(env[, sample], main=paste("Histogram of ", sample, " after Log Transformation", sep=""), xlab=paste(sample, " (", sample_true_name, ")", sep="")) #produces histogram of BMI (no units involved)
		}
		else {
			hist(env[, sample], main=paste("Histogram of ", sample, " after Log Transformation", sep=""), xlab=paste(sample, " (", sample_true_name, " (", sample_units, ")", ")", sep="")) #produces histogram of all significant numeric variables for reference
		}
		dev.off()
		medians <- median(na.omit(as.numeric(env[, sample]))) #identify the mean of the given column
		keep_patients <- which(!is.na(env[, sample])) #identify which cases should be reatained for the PCA (not NA)
		env_copy <- env[keep_patients, ]
		otus <- otus[keep_patients, ] #subset OTU table to just the patients to be retained
		otu.hel <- decostand(otus,"hellinger") #hellinger transformation of matrix

		sample.pca <- rda(otu.hel) #generate pca plot and space to color
		pc1 <- summary(eigenvals(sample.pca))[[1]][[2]]*100 #convert proportion of variance explained by pc1 to a percentage
		pc2 <- summary(eigenvals(sample.pca))[[1]][[5]]*100 #convert proportion of variance explained by pc2 to a percentage

		index <- match((sample), colnames(env_copy)) #grabs the index of the column in the original environmental variable data set
		median <- env[, index] %>% median()
		greater <- which(env_copy[, index] > median)
		lower <- which(env_copy[, index] <= median)
		i <- rownames(env_copy)
		if (sample != "BMI.category") {
			env_copy[greater, index] <- 2 #replace with 2 if above median
			env_copy[lower, index] <- 1 #replace with 1 if below median
		}
		titleplot <- sample
		colors <- as.factor(env_copy[, index]) #convert categories of the elements of the column to colors
		mypath <- file.path(paste(paste("~/DNA_V1_PCA_", titleplot, ".pdf", sep = "")))
		pdf(file=mypath) #sets output path for the plot
		#ordiplot(sample.pca, type="n", font=2, font.lab=2, main=paste("DNA_V1 PCA by ", sample, " (", sample_true_name, ")", sep = ""), xlab = paste("PC1 (", pc1, "% of variance explained)", sep = ""), 
			#ylab = paste("PC2 (", pc2, "% of variance explained)", sep = ""))
		ordiplot(sample.pca, type="n", font=2, font.lab=2, main="Sputum 16S DNA PCA by BMI", xlab = paste("PC1 (", pc1, "% of variance explained)", sep = ""), 
			ylab = paste("PC2 (", pc2, "% of variance explained)", sep = ""))
		if (sample == "BMI.category") {
			for (id in rownames(env_copy)) {
				if (env_copy[id, "BMI.category"] == 1) {
					env_copy[id, "BMI.category"] = c("Normal") 
				}
				else if (env_copy[id, "BMI.category"] == 2) {
					env_copy[id, "BMI.category"] = c("Overweight") 
				}
				else if (env_copy[id, "BMI.category"] == 3) {
					env_copy[id, "BMI.category"] = c("Obese")
				}
			}
			#beta dispersion calculation
			print(sample)
			dist_matrix <- vegdist(otu.hel, method = "euclidean")
			beta_dispersion <- betadisper(dist_matrix, group = as.factor(env_copy[, "BMI.category"]))
			print(permutest(beta_dispersion))

			colors <- as.factor(env_copy[, "BMI.category"])
			points(sample.pca, col=as.factor(colors), pch=21, bg=as.factor(colors)) #fills in points with colors based on category
			ordispider(sample.pca, groups = env_copy[, "BMI.category"], show.groups = c("Normal", "Overweight", "Obese"), col = c("black", "red", "green"))
			for (i in seq(1:4)) {
				normal_counts <- length(which(env_copy$"BMI.category" == "Normal"))
				overweight_counts <- length(which(env_copy$"BMI.category" == "Overweight"))
				obese_counts <- length(which(env_copy$"BMI.category" == "Obese"))
				legend("topright", legend=paste(c("Normal (n = ", "Obese (n = ", "Overweight (n = "), c(normal_counts, obese_counts, overweight_counts), c(")", ")", ")"), sep=""), col=c("black", "red", "green"), fill=c("black", "red", "green"), border=c("black", "red", "green"), bty="n")
				legend("bottomright", legend=paste("Adonis P-Value = ", adonis_p_value))
			}
		}
		else {
			#beta dispersion calculation
			print(sample)
			dist_matrix <- vegdist(otu.hel, method = "euclidean")
			beta_dispersion <- betadisper(dist_matrix, group = as.factor(env_copy[, sample]))
			print(permutest(beta_dispersion))

			points(sample.pca, col=as.factor(colors), pch=21, bg=as.factor(colors)) #fills in points with colors based on category
			for (i in seq(1:2)) {
				ordispider(sample.pca, groups = as.factor(env_copy[, sample]), show.groups = i, col = i)
			}
			legend_median <- round(median, digits=2) #truncate at hundreth's place for better visualization
			for (i in seq(1:4)) {
				above_median <- length(which(env_copy[, sample] == 2))
				below_median <- length(which(env_copy[, sample] == 1))
				legend("topright", legend=c(paste("=< ", legend_median, " ", sample_units, " ", "(n = ", below_median, ")", sep=""), 
				paste("> ", legend_median, " ", sample_units, " ", "(n = ", above_median, ")", sep="")), col=c("black", "red"), fill=c("black", "red"), border=c("black", "red"), bty="n")
				legend("bottomright", legend=paste("Adonis P-Value = ", adonis_p_value))
			}
		}
		dev.off()
	}
}

sig_traits <- read.csv("~/DNA_Adonis_Matched_PVals-sorted.csv") #gather list of traits and their p-values
colnames(sig_traits) <- c("Variable", "Adonis P-Value")

#for variables with > 2 factors & NAs - have to do manually :( 

#1 - Study
otus <- read.csv("~/sputDNA.lowsquam.otu.V1.csv", header=T)
otus <- otus[, -1] #remove sample id column from otu file
otu.hel <- decostand(otus,"hellinger") #hellinger transformation of matrix
study.pca <- rda(otu.hel)
pc1 <- summary(eigenvals(study.pca))[[1]][[2]]*100 #convert proportion of variance explained by pc1 to a percentage
pc2 <- summary(eigenvals(study.pca))[[1]][[5]]*100 #convert proportion of variance explained by pc2 to a percentage

#beta dispersion
if (nlevels(env[, "Study"]) > 1) { #if looking at multiple studies and not just CAARS
	print("Study")
	dist_matrix <- vegdist(otu.hel, method = "euclidean")
	beta_dispersion <- betadisper(dist_matrix, group = env[, "Study"])
	print(permutest(beta_dispersion))

	colors <- as.factor(env[, "Study"])
	abba <- length(which(env[, "Study"] == "ABBA"))
	caars <- length(which(env[, "Study"] == "CAARS"))
	mampa <- length(which(env[, "Study"] == "MAMPA"))
	op <- par(cex = 0.5)
	pdf("~/DNA_PCA_Study.pdf")
	ordiplot(study.pca, type="n", font=2, font.lab=2, 
		main="Sputum 16S DNA PCA by Study",
		xlab = paste("PC1 (", pc1, "% of variance explained)", sep = ""), 
		ylab = paste("PC2 (", pc2, "% of variance explained)", sep = ""))
	points(study.pca, col=as.factor(colors), pch=21, bg=as.factor(colors)) #fills in points with colors based on category
	ordispider(study.pca, env[, "Study"], show.groups = c("ABBA", "CAARS", "MAMPA"), col = c("black", "red", "green"))
	legend("topright", legend=paste(c("ABBA", "CAARS", "MAMPA"), c(" (n= ", " (n= ", " (n= "), c(abba, caars, mampa), c(")", ")", ")"), sep = "") , col=c("black", "red", "green"), fill=c("black", "red", "green"), border=c("black", "red", "green"), bty="n", text.font = 2)
	legend("bottomright", legend=paste("Adonis P-Value = ", sig_traits[sig_traits$"Variable" == "Study", "Adonis P-Value"]), text.font = 2)
	dev.off()
}

#2 - ICSdose_HiMedLow
ics_dose.pca <- rda(otu.hel)
pc1 <- summary(eigenvals(ics_dose.pca))[[1]][[2]]*100 #convert proportion of variance explained by pc1 to a percentage
pc2 <- summary(eigenvals(ics_dose.pca))[[1]][[5]]*100 #convert proportion of variance explained by pc2 to a percentage

#beta dispersion
print("ICSdose_HiMedLow")
dist_matrix <- vegdist(otu.hel, method = "euclidean")
beta_dispersion <- betadisper(dist_matrix, group = env[, "ICSdose_HiMedLow"])
print(permutest(beta_dispersion))

colors <- as.factor(env[, "ICSdose_HiMedLow"])
low_count <- length(which(env[, "ICSdose_HiMedLow"] == "Low"))
med_count <- length(which(env[, "ICSdose_HiMedLow"] == "Med"))
none_count <- length(which(env[, "ICSdose_HiMedLow"] == "None"))
op <- par(cex = 0.5)
pdf("~/DNA_PCA_ICSdose_HiMedLow.pdf")
ordiplot(ics_dose.pca, type="n", font=2, font.lab=2, 
	main="Sputum 16S DNA PCA by ICS Dose Category",
	xlab = paste("PC1 (", pc1, "% of variance explained)", sep = ""), 
	ylab = paste("PC2 (", pc2, "% of variance explained)", sep = ""))
points(ics_dose.pca, col=as.factor(colors), pch=21, bg=as.factor(colors)) #fills in points with colors based on category
ordispider(ics_dose.pca, env[, "ICSdose_HiMedLow"], show.groups = c("Low", "Med", "None", "NA"), col = c("black", "red", "green"))
legend("topright", legend=paste(c("Low", "Med", "None"), c(" (n= ", " (n= ", " (n= "), c(low_count, med_count, none_count), c(")", ")", ")"), sep = "") , col=c("black", "red", "green"), fill=c("black", "red", "green"), border=c("black", "red", "green"), bty="n", text.font = 2)
legend("bottomright", legend=paste("Adonis P-Value = ", sig_traits[sig_traits$"Variable" == "ICSdose_HiMedLow", "Adonis P-Value"]), text.font = 2)
dev.off()

#3 - For Asthmatics, ICS Yes/No
ics_dose.pca <- rda(otu.hel)
pc1 <- summary(eigenvals(ics_dose.pca))[[1]][[2]]*100 #convert proportion of variance explained by pc1 to a percentage
pc2 <- summary(eigenvals(ics_dose.pca))[[1]][[5]]*100 #convert proportion of variance explained by pc2 to a percentage

#beta dispersion
print("ICS_use_yes_no")
dist_matrix <- vegdist(otu.hel, method = "euclidean")
beta_dispersion <- betadisper(dist_matrix, group = env[, "ICS_use_yes_no"])
print(permutest(beta_dispersion))

colors <- as.factor(env[, "ICS_use_yes_no"])
yes_count <- length(which(env[, "ICS_use_yes_no"] == "Yes"))
no_count <- length(which(env[, "ICS_use_yes_no"] == "No"))
op <- par(cex = 0.5)
pdf("~/DNA_PCA_ICS_Yes_No.pdf")
ordiplot(ics_dose.pca, type="n", font=2, font.lab=2, 
	main="Sputum 16S DNA PCA by ICS use among Asthmatics",
	xlab = paste("PC1 (", pc1, "% of variance explained)", sep = ""), 
	ylab = paste("PC2 (", pc2, "% of variance explained)", sep = ""))
points(ics_dose.pca, col=as.factor(colors), pch=21, bg=as.factor(colors)) #fills in points with colors based on category
ordispider(ics_dose.pca, env[, "ICS_use_yes_no"], show.groups = c("Yes", "No", "NA"), col = c("black", "red"))
legend("topright", legend=paste(c("Yes", "No"), c(" (n= ", " (n= "), c(yes_count, no_count), c(")", ")"), sep = "") , col=c("red", "black"), fill=c("red", "black"), border=c("red", "black"), bty="n", text.font = 2)
legend("bottomright", legend=paste("Adonis P-Value = ", sig_traits[sig_traits$"Variable" == "ICS_use_yes_no", "Adonis P-Value"]), text.font = 2)
dev.off()

#4 - ICStype
otus <- read.csv("~/sputDNA.lowsquam.otu.V1.csv", header=T)

remove_patients <- which(is.na(env[, "ICStype"])) #identify which cases should be reatained for the PCA (not NA)
none <- which(env[, "ICStype"] == "None")
remove_patients <- c(remove_patients, none)
names <- rownames(env)
env_copy <- env %>% filter(!names %in% rownames(env[remove_patients, ]))
rownames(otus) <- otus[, 1]
names <- rownames(otus)
otus <- otus %>% filter(!names %in% rownames(otus[remove_patients, ])) #subset OTU table to just the patients to be retained
otus <- otus[, -1]
otu.hel <- decostand(otus,"hellinger") #hellinger transformation of matrix
ics_type.pca <- rda(otu.hel)
pc1 <- summary(eigenvals(ics_type.pca))[[1]][[2]]*100 #convert proportion of variance explained by pc1 to a percentage
pc2 <- summary(eigenvals(ics_type.pca))[[1]][[5]]*100 #convert proportion of variance explained by pc2 to a percentage

#beta dispersion
print("ICStype")
dist_matrix <- vegdist(otu.hel, method = "euclidean")
beta_dispersion <- betadisper(dist_matrix, group = env[, "ICStype"])
print(permutest(beta_dispersion))

colors <- as.factor(env_copy[, "ICStype"])
bf_count <- length(which(env_copy[,  "ICStype"] == "budesone and fluticasone"))
budesonide_count <- length(which(env_copy[,  "ICStype"] == "budesonide"))
fluticasone_count <- length(which(env_copy[,  "ICStype"] == "fluticasone"))
mb_count <- length(which(env_copy[,  "ICStype"] == "budesone and fluticasone"))
pdf("~/DNA_PCA_ICStype.pdf")
ordiplot(ics_type.pca, type="n", font=2, font.lab=2, 
	main="Sputum 16S DNA PCA by ICS Type",
	xlab = paste("PC1 (", pc1, "% of variance explained)", sep = ""), 
	ylab = paste("PC2 (", pc2, "% of variance explained)", sep = ""))
points(ics_type.pca, col=as.factor(colors), pch=21, bg=as.factor(colors)) #fills in points with colors based on category
ordispider(ics_type.pca, env_copy[, "ICStype"], show.groups = c("budesone and fluticasone", "budesonide", "fluticasone", "mometasone and budesonide", "NA"), col = c("black", "red", "green", "blue"))
legend("bottomleft", legend=paste(c("budesone and fluticasone", "budesonide", "fluticasone", "mometasone and budesonide"), c(" (n=", " (n=", " (n=", " (n="), c(bf_count, budesonide_count, fluticasone_count, mb_count), c(")", ")", ")", ")"), sep="") , col = c("black", "red", "green", "blue"), fill=c("black", "red", "green", "blue"), border=c("black", "red", "green", "blue"), bty="n", text.font = 2)
legend("bottomright", legend=paste("Adonis P-Value = ", sig_traits[sig_traits$"Variable" == "ICStype", "Adonis P-Value"]), text.font = 2)
dev.off()

#5 - Other Respiratory Diseases
otus <- read.csv("~/sputDNA.lowsquam.otu.V1.csv", header=T)
otus <- otus[, -1]
otu.hel <- decostand(otus,"hellinger") #hellinger transformation of matrix
resp_dis.pca <- rda(otu.hel)
pc1 <- summary(eigenvals(resp_dis.pca))[[1]][[2]]*100 #convert proportion of variance explained by pc1 to a percentage
pc2 <- summary(eigenvals(resp_dis.pca))[[1]][[5]]*100 #convert proportion of variance explained by pc2 to a percentage

#beta dispersion
print("Other_respir_dis")
dist_matrix <- vegdist(otu.hel, method = "euclidean")
beta_dispersion <- betadisper(dist_matrix, group = env[, "other_respir_dis"])
print(permutest(beta_dispersion))

colors <- as.factor(env[, "other_respir_dis"])
ar_count <- length(which(env[,  "other_respir_dis"] == "allergic rhinitis"))
sa_count <- length(which(env[,  "other_respir_dis"] == "sinusitis, alpha1-AT without lung emphysema"))
sp_count <- length(which(env[,  "other_respir_dis"] == "sinus polyps(Samter's triad)"))
sn_count <- length(which(env[,  "other_respir_dis"] == "sinusitis, nasal polyps"))
pdf("~/DNA_PCA_other_respir_dis.pdf")
ordiplot(resp_dis.pca, type="n", font=2, font.lab=2, 
	main="Sputum 16S DNA PCA by Other Respiratory Diseases",
	xlab = paste("PC1 (", pc1, "% of variance explained)", sep = ""), 
	ylab = paste("PC2 (", pc2, "% of variance explained)", sep = ""))
points(resp_dis.pca, col=as.factor(colors), pch=21, bg=as.factor(colors)) #fills in points with colors based on category
ordispider(resp_dis.pca, env[, "other_respir_dis"], show.groups = c("allergic rhinitis", "sinusitis, alpha1-AT without lung emphysema", "sinus polyps(Samter's triad)", "sinusitis, nasal polyps", "NA"), col = c("black", "red", "green", "blue"))
legend("bottomleft", legend=paste(c("allergic rhinitis", "sinusitis, alpha1-AT without lung emphysema", "sinus polyps(Samter's triad)", "sinusitis, nasal polyps"), 
	c(" (n=", " (n=", " (n=", " (n="), 
	c(ar_count, sa_count, sp_count, sn_count), 
	c(") ", ") ", ") ", ") "), sep=""), 
	col = c("black", "red", "green", "blue"), 
	fill=c("black", "red", "green", "blue"), 
	border=c("black", "red", "green", "blue"), 
	bty="n", text.font = 2)
legend("bottomright", legend=paste("Adonis P-Value = ", sig_traits[sig_traits$"Variable" == "other_respir_dis", "Adonis P-Value"]), text.font = 2)
dev.off()
