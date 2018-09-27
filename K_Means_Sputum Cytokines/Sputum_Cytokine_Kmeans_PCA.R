library(tidyverse)
library(vegan)
otus <- read.csv("~/sputDNA.lowsquam.WGS.species.caarsonly.top_100.csv")
rownames(otus) <- otus[, 1]
otus <- otus[, -1]

otu.hel <- decostand(otus,"hellinger") #hellinger transformation of matrix
otus.pca <- rda(otu.hel) #generate naive pca plot and space to color
pc1.otu <- summary(eigenvals(otus.pca))[[1]][[2]]*100 #convert proportion of variance explained by pc1 to a percentage
pc2.otu <- summary(eigenvals(otus.pca))[[1]][[5]]*100 #convert proportion of variance explained by pc2 to a percentage
scores <- scores(otus.pca)$sites

#Elbow method for optimal k
k.max <- 20 #max clusters = number of patient samples
wss <- sapply(1:k.max, function(k) {
              	kmeans(scores, k, nstart=50,iter.max = 20)$tot.withinss})
pdf("~/OTUs_Kmeans_Scree_Plot.pdf")
df <- data.frame(k = 1:k.max, wss)
ggplot (data = df) +
geom_point(mapping = aes(x = k, y = wss)) +
geom_line(mapping = aes(x = k, y = wss)) +
ggtitle("Scree Plot of OTU PCA K-means Clustering") +
labs(x = "Number of clusters K", y = "Total within-clusters sum of squares") +
scale_x_continuous(breaks = c(1:20))
dev.off()

#k means clustering of otu PCA
km.otu <- kmeans(scores, centers=3, nstart=50) #select 3 based on scree plot

#plot PCA and color by k means clustering
mypath <- file.path("~/WGS_K_Means_PCA/OTUs_K_Means_PCA.pdf")
pdf(file=mypath) #sets output path for the plot
ordiplot(otus.pca, type = "n", main = "PCA of WGS Sputum Microbial Communities Clustered by K-Means")
points(otus.pca, col=as.factor(km.otu$cluster), pch=21, bg=as.factor(km.otu$cluster)) #fills in points with colors based on category
ordispider(otus.pca, groups = km.otu$cluster, show.groups = km.otu$cluster)
dev.off()

cytokines <- read.csv("~/sputDNA.lowsquam.cytokine.caarsonly.cytokines_subsetted.csv")
rownames(cytokines) <- cytokines[, 1]
cytokines <- cytokines[, -1]

cyt.hel <- decostand(cytokines,"hellinger") #hellinger transformation of matrix
cyt.pca <- rda(cyt.hel) #generate naive pca plot and space to color - only includes numerically defined cytokine measurements
pc1.cyt <- summary(eigenvals(cyt.pca))[[1]][[2]]*100 #convert proportion of variance explained by pc1 to a percentage
pc2.cyt <- summary(eigenvals(cyt.pca))[[1]][[5]]*100 #convert proportion of variance explained by pc2 to a percentage
scores <- scores(cyt.pca)$sites

#Elbow method for optimal k (i.e. scree plot)
k.max <- 20 #max clusters = number of patient samples
wss <- sapply(1:k.max, function(k) {
              	kmeans(scores, k, nstart=50,iter.max = 20)$tot.withinss})
pdf("~/WGS_K_Means_PCA/Cytokines_Kmeans_Scree_Plot.pdf")
df <- data.frame(k = 1:k.max, wss)
ggplot (data = df) +
geom_point(mapping = aes(x = k, y = wss)) +
geom_line(mapping = aes(x = k, y = wss)) +
ggtitle("Scree Plot of Sputum Cytokine PCA K-means Clustering") +
labs(x = "Number of clusters K", y = "Total within-clusters sum of squares") +
scale_x_continuous(breaks = c(1:20))
dev.off()

#k means clustering of cytokine PCA
km.cyt <- kmeans(scores, centers=3, nstart=50) #select 5 based on scree plot

#plot cytokine PCA and color by k means clustering
mypath <- file.path("~/Cytokines_K_Means_PCA.pdf")
pdf(file=mypath) #sets output path for the plot
ordiplot(cyt.pca, type = "n", main = "PCA of WGS Sputum Cytokine Samples Clustered by K-Means")
points(cyt.pca, col=as.factor(km.cyt$cluster), pch=21, bg=as.factor(km.cyt$cluster)) #fills in points with colors based on category
ordispider(cyt.pca, groups = km.cyt$cluster, show.groups = km.cyt$cluster)
dev.off()

#plot cytokine PCA and color by k means clustering of OTUs to show relationship between OTUs and cytokine measurements
mypath <- file.path("~/Cytokines_PCA_Colored_By_OTU_K_Means_PCA.pdf")
pdf(file=mypath) #sets output path for the plot
ordiplot(cyt.pca, type = "n", cex.main = 0.75, font = 2, font.lab = 2,
	main = "PCA of Sputum Cytokine Measurements, Colored by WGS Microbial Community K Means Cluster", 
	xlab = paste("PC1 (", pc1.cyt, "% of variance explained)", sep = ""), 
	ylab = paste("PC2 (", pc2.cyt, "% of variance explained)", sep = ""))
points(cyt.pca, col=as.factor(km.otu$cluster), pch=21, bg=as.factor(km.otu$cluster)) #fills in points with colors based on category
ordispider(cyt.pca, groups = km.otu$cluster, show.groups = km.otu$cluster, col = c("Black", "Red", "Green", "Blue"))
dev.off()

#plot backwards (cytokine k means used to describe PCA variations in microbial communities)
mypath <- file.path("~/OTUs_PCA_Colored_By_Cytokines_K_Means_PCA.pdf")
pdf(file=mypath) #sets output path for the plot
ordiplot(otus.pca, type = "n", cex.main = 0.75, font = 2, font.lab = 2,
	main = "PCA of WGS Microbial Communities, Colored by Sputum Cytokine Measurements K Means Cluster", 
	xlab = paste("PC1 (", pc1.otu, "% of variance explained)", sep = ""), 
	ylab = paste("PC2 (", pc2.otu, "% of variance explained)", sep = ""))
points(otus.pca, col=as.factor(km.cyt$cluster), pch=21, bg=as.factor(km.cyt$cluster)) #fills in points with colors based on category
ordispider(otus.pca, groups = km.cyt$cluster, show.groups = km.cyt$cluster, col = c("Black", "Red", "Green", "Blue"))
dev.off()

