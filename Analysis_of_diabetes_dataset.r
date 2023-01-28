################### INFO ######################
#    DATA ANALYSIS AND STATISTICAL LEARNING   #
#     DATA SCIENCE FOR MANAGEMENT - LM 91     #

##### TITLE: ANALYSIS OF DIABETES DATASET #####

################## LIBRARIES ##################
library(e1071)
library(cluster)
library(factoextra)
library(stats)
library(clustertend)
library(NbClust)
library(fpc)
library(clValid)
library(gridExtra)
library(grid)
library(gamlss)
set.seed(123)


################ READ THE DATA ################
library(mclust)
attach(diabetes)
str(diabetes)
head(diabetes)


############ UNIVARIATE ANALYSIS ##############
#### Class ####
table(diabetes$class)
round(table(diabetes$class)/length(diabetes$class)*100, digits=2)
pie(table(diabetes$class)/length(diabetes$class)*100, col = c("green","blue", "red"))


#### Glucose ####
summary(diabetes$glucose)
sd(diabetes$glucose)
var(diabetes$glucose)
(quantile(diabetes$glucose,0.75)) - (quantile(diabetes$glucose,0.25))
hist(diabetes$glucose, freq=FALSE, main = "Histogram of Glucose", xlab = "mg/100ml·h")
lines(density(diabetes$glucose), col="red")
skewness(diabetes$glucose)
kurtosis(diabetes$glucose, type = 1)
# Fitting Glucose
par(mfrow=c(2,2)) # Set the plot area 
fit.gamma <- histDist(diabetes$glucose, family=GA, nbins=21, main="Gamma distribution")
fit.t = histDist(diabetes$glucose, family=TF, nbins=21, main="t-Student distribution")
fit.sep2 = histDist(diabetes$glucose, family=SEP2, nbins=21, main="Skew Power Exponential type 2 distribution")
fit.st5 = histDist(diabetes$glucose, family=ST5, nbins=21, main="Skew t type 5")
# AIC and BIC values
data.frame(row.names=c("Gamma", "t-Student", "SEP2", "ST5"),
           AIC=c(AIC(fit.gamma), AIC(fit.t), AIC(fit.sep2),AIC(fit.st5)),
           BIC=c(fit.gamma$sbc, fit.t$sbc, fit.sep2$sbc, fit.st5$sbc))
par(mfrow=c(1,1)) # Reset the plot area

# Boxplot of Glucose
boxplot(diabetes$glucose, main = "Boxplot of Glucose", ylab = "mg/100ml·h")
boxplot.stats(diabetes$glucose)$out


#### Insulin ####
summary(diabetes$insulin)
sd(diabetes$insulin)
var(diabetes$insulin)
(quantile(diabetes$insulin,0.75)) - (quantile(diabetes$insulin,0.25))
hist(diabetes$insulin, freq=FALSE, main = "Histogram of Insulin", xlab = "μU/ml·h")
lines(density(diabetes$insulin), col="red")
skewness(diabetes$insulin)
kurtosis(diabetes$insulin, type = 1)
# Fitting Insulin
par(mfrow=c(2,2)) # Set the plot area 
fit.gamma <- histDist(diabetes$insulin, family=GA, nbins=21, main="Gamma distribution")
fit.t = histDist(diabetes$insulin, family=TF, nbins=21, main="t-Student distribution")
fit.sep2 = histDist(diabetes$insulin, family=SEP2, nbins=21, main="Skew Power Exponential type 2 distribution")
fit.st5 = histDist(diabetes$insulin, family=ST5, nbins=21, main="Skew t type 5")
# AIC and BIC values
data.frame(row.names=c("Gamma", "t-Student", "SEP2", "ST5"),
           AIC=c(AIC(fit.gamma), AIC(fit.t), AIC(fit.sep2),AIC(fit.st5)),
           BIC=c(fit.gamma$sbc, fit.t$sbc, fit.sep2$sbc, fit.st5$sbc))
par(mfrow=c(1,1)) # Reset the plot area

# Boxplot of Insulin
boxplot(diabetes$insulin, main = "Boxplot of Insulin", ylab = "μU/ml·h")
boxplot.stats(diabetes$insulin)$out


#### Sspg ####
summary(diabetes$sspg)
sd(diabetes$sspg)
var(diabetes$sspg)
(quantile(diabetes$sspg,0.75)) - (quantile(diabetes$sspg,0.25))
hist(diabetes$sspg, freq=FALSE, main = "Histogram of Sspg", xlab = "mg/100ml")
lines(density(diabetes$sspg), col="red")
skewness(diabetes$sspg)
kurtosis(diabetes$sspg, type = 1)
# Fitting Sspg
par(mfrow=c(2,2)) # Set the plot area 
fit.gamma <- histDist(diabetes$sspg, family=GA, nbins=21, main="Gamma distribution")
fit.t = histDist(diabetes$sspg, family=TF, nbins=21, main="t-Student distribution")
fit.sep2 = histDist(diabetes$sspg, family=SEP2, nbins=21, main="Skew Power Exponential type 2 distribution")
fit.st5 = histDist(diabetes$sspg, family=ST5, nbins=21, main="Skew t type 5")
# AIC and BIC values
data.frame(row.names=c("Gamma", "t-Student", "SEP2", "ST5"),
           AIC=c(AIC(fit.gamma), AIC(fit.t), AIC(fit.sep2),AIC(fit.st5)),
           BIC=c(fit.gamma$sbc, fit.t$sbc, fit.sep2$sbc, fit.st5$sbc))
par(mfrow=c(1,1)) # Reset the plot area

# Boxplot of Sspg
boxplot(diabetes$sspg, main = "Boxplot of Sspg", ylab = "mg/100ml")
boxplot.stats(diabetes$sspg)$out


############ UNIVARIATE ANALYSIS FOR EACH LEVEL OF CLASS ############

# Group by class
normal_data = diabetes[diabetes$class == "Normal", ]
chemical_data = diabetes[diabetes$class == "Chemical", ]
overt_data = diabetes[diabetes$class == "Overt", ]

# Normal - Glucose
summary(normal_data$glucose)
sd(normal_data$glucose)
var(normal_data$glucose)
(quantile(normal_data$glucose,0.75)) - (quantile(normal_data$glucose,0.25))
hist(normal_data$glucose, freq=FALSE,  main = "Histogram of Normal Glucose", xlab = "mg/100ml·h")
lines(density(normal_data$glucose), col="red")
skewness(normal_data$glucose)
kurtosis(normal_data$glucose, type = 1)
boxplot(normal_data$glucose, main = "Boxplot Normal of Glucose", ylab = "mg/100ml·h")

# Normal - Insulin
summary(normal_data$insulin)
sd(normal_data$insulin)
var(normal_data$insulin)
(quantile(normal_data$insulin,0.75)) - (quantile(normal_data$insulin,0.25))
hist(normal_data$insulin, freq=FALSE, main = "Histogram of Normal Insulin", xlab = "μU/ml·h")
lines(density(normal_data$insulin), col="red")
skewness(normal_data$insulin)
kurtosis(normal_data$insulin, type = 1)
boxplot(normal_data$insulin, main = "Boxplot of Normal Insulin", ylab = "μU/ml·h")

# Normal - SSPG
summary(normal_data$sspg)
sd(normal_data$sspg)
var(normal_data$sspg)
(quantile(normal_data$sspg,0.75)) - (quantile(normal_data$sspg,0.25))
hist(normal_data$sspg, freq=FALSE, main = "Histogram of Normal Sspg", xlab = "mg/100ml")
lines(density(normal_data$sspg), col="red")
skewness(normal_data$sspg)
kurtosis(normal_data$sspg, type = 1)
boxplot(normal_data$sspg, main = "Boxplot of Normal Sspg", ylab = "mg/100ml")
boxplot.stats(normal_data$sspg)$out

# Chemical - Glucose
summary(chemical_data$glucose)
sd(chemical_data$glucose)
var(chemical_data$glucose)
(quantile(chemical_data$glucose,0.75)) - (quantile(chemical_data$glucose,0.25))
hist(chemical_data$glucose, freq=FALSE, main = "Histogram of Chemical Glucose", xlab = "mg/100ml·h")
lines(density(chemical_data$glucose), col="red")
skewness(chemical_data$glucose)
kurtosis(chemical_data$glucose, type = 1)
boxplot(chemical_data$glucose, main = "Boxplot Chemical of Glucose", ylab = "mg/100ml·h")

# Chemical - Insulin
sd(chemical_data$insulin)
var(chemical_data$insulin)
(quantile(chemical_data$insulin,0.75)) - (quantile(chemical_data$insulin,0.25))
hist(chemical_data$insulin, freq=FALSE, main = "Histogram of Chemical Insulin", xlab = "μU/ml·h")
lines(density(chemical_data$insulin), col="red")
skewness(chemical_data$insulin)
kurtosis(chemical_data$insulin, type = 1)
boxplot(chemical_data$insulin, main = "Boxplot of Chemical Insulin", ylab = "μU/ml·h")
boxplot.stats(chemical_data$insulin)$out

# Chemical - SSPG
sd(chemical_data$sspg)
var(chemical_data$sspg)
(quantile(chemical_data$sspg,0.75)) - (quantile(chemical_data$sspg,0.25))
hist(chemical_data$sspg, freq=FALSE, main = "Histogram of Chemical Sspg", xlab = "mg/100ml")
lines(density(chemical_data$sspg), col="red")
skewness(chemical_data$sspg)
kurtosis(chemical_data$sspg, type = 1)
boxplot(chemical_data$sspg, main = "Boxplot of Chemical Sspg", ylab = "mg/100ml")
boxplot.stats(chemical_data$sspg)$out

# Overt - Glucose
summary(overt_data$glucose)
sd(overt_data$glucose)
var(overt_data$glucose)
(quantile(overt_data$glucose,0.75)) - (quantile(overt_data$glucose,0.25))
hist(overt_data$glucose, freq=FALSE, main = "Histogram of Overt Glucose", xlab = "mg/100ml·h")
lines(density(overt_data$glucose), col="red")
skewness(overt_data$glucose)
kurtosis(overt_data$glucose, type = 1)
boxplot(overt_data$glucose, main = "Boxplot Overt of Glucose", ylab = "mg/100ml·h")

# Overt - Insulin
sd(overt_data$insulin)
var(overt_data$insulin)
(quantile(overt_data$insulin,0.75)) - (quantile(overt_data$insulin,0.25))
hist(overt_data$insulin, freq=FALSE, main = "Histogram of Overt Insulin", xlab = "μU/ml·h")
lines(density(overt_data$insulin), col="red")
skewness(overt_data$insulin)
kurtosis(overt_data$insulin, type = 1)
boxplot(overt_data$insulin, main = "Boxplot of Overt Insulin", ylab = "μU/ml·h")

# Overt - SSPG
sd(overt_data$sspg)
var(overt_data$sspg)
(quantile(overt_data$sspg,0.75)) - (quantile(overt_data$sspg,0.25))
hist(overt_data$sspg, freq=FALSE, main = "Histogram of Overt Sspg", xlab = "mg/100ml")
lines(density(overt_data$sspg), col="red")
skewness(overt_data$sspg)
kurtosis(overt_data$sspg, type = 1)
boxplot(overt_data$sspg, main = "Boxplot of Overt Sspg", ylab = "mg/100ml")
boxplot.stats(overt_data$sspg)$out


######### PRINCIPAL COMPONENT ANALYSIS ########

# Correlation first
################ CORRELATION ##################
# Simple correlation matrix
(cor(diabetes[2:4]))

# By default, prcomp() centers the variables to have mean zero.
# By using the option scale = TRUE, we scale the variables to have standard deviation one.
pr.out = prcomp(diabetes[2:4], scale = TRUE)
pr.out

# Single components of prcomp
# Standard deviations of principal components
pr.out$sdev

# Means of original variables
pr.out$center

# Standard deviation of original variables
pr.out$scale

# Rotation are the principal component loadings (with sign flip) 
pr.out$rotation = -pr.out$rotation
pr.out$rotation

# Out are the scores (with sign flip)
pr.out$x = -pr.out$x
head(pr.out$x)

# Visualize PCA
fviz_pca_ind(pr.out, title = "PCA - Diabetes Data", legend = "top",
             habillage = diabetes$class, palette = c("green","blue","red"),
             geom = "point", ggtheme = theme_classic())

################## BIPLOT ####################
biplot(pr.out, scale=0, cex = 0.8)

############## KAISER'S RULE #################
pr.var = pr.out$sdev^2 # variances of the PCs 
pr.var

################ SCREE PLOT ##################
pve = pr.var/sum(pr.var)
pve
par(mfrow=c(1,2)) #Set the plotting area into a 1*2 array
plot(pve, main = "Scree plot", xlab="Principal Component", ylab="Proportion of Variance Explained", ylim=c(0,1), type='b')
plot(cumsum(pve), main ="Cumulative Scree plot", xlab="Principal Component", ylab="Cumulative Proportion of Variance Explained", ylim=c(0,1),type='b')
par(mfrow=c(1,1)) #Reset plot


############### CLUSTER ANALYSIS ##############
########## CLUSTER TENDENCY ############
random_df = apply(diabetes[2:4], 2, function(x){runif(length(x), min(x), (max(x)))})
random_df = as.data.frame(random_df)
scaled.random_df = scale(random_df)

# Graphical representation of the data
# Standardized df data
scaled.df = scale(diabetes[2:4])
pairs(scaled.df, pch = 21)

# Standardized uniform random df data
pairs(scaled.random_df, pch = 21)

# Plot the standardized df data 
fviz_pca_ind(pr.out, title = "PCA - Diabetes Data", legend = "bottom",
             habillage = diabetes$class, palette = c("green","blue","red"),
             geom = "point", ggtheme = theme_classic())

# Plot the random df data
fviz_pca_ind(prcomp(scaled.random_df), title = "PCA - Random data",
             geom = "point", ggtheme = theme_classic())

# Hopkins statistic
hopkins(scaled.df, n = nrow(scaled.df)-1)

# Compute Hopkins statistic for a random dataset
hopkins(scaled.random_df, n = nrow(random_df)-1)

# Computing the dissimilarity matrix with Euclidean distance
dist.eucl = dist(scaled.df, method = "euclidean")

# Subset the first 6 columns and rows and round the values
round(as.matrix(dist.eucl)[1:6, 1:6], 2)

# (VAT algorithm)
## Visualizing distance matrix
fviz_dist(dist.eucl, show_labels = FALSE)

# Visualizing distance matrix of random df
fviz_dist(dist(scaled.random_df, method = "euclidean"), show_labels = FALSE)


########## OPTIMAL NUMBER OF CLUSTERS ############
#### Hierarchical Clustering ####
# Elbow method
fviz_nbclust(scaled.df, hcut, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)+
  labs(subtitle = "Hierarchical: Elbow method")

# Silhouette method
fviz_nbclust(scaled.df, hcut, method = "silhouette") +
  labs(subtitle = "Hierarchical: Silhouette method")

# Gap statistic
fviz_nbclust(scaled.df, hcut, method = "gap_stat", nboot = 500) +
  labs(subtitle = "Hierarchical: Gap statistic method")

# NbClust function
grid.arrange(arrangeGrob(fviz_nbclust(NbClust(scaled.df, distance="euclidean", method="single")),
                         top=textGrob("Single linkage method with euclidean distance",gp=gpar(fontsize=17))),
             arrangeGrob(fviz_nbclust(NbClust(scaled.df, distance="euclidean", method="average")),
                         top=textGrob("Average linkage method with euclidean distance",gp=gpar(fontsize=17))),
             arrangeGrob(fviz_nbclust(NbClust(scaled.df, distance="euclidean", method="complete")),
                         top=textGrob("Complete linkage method with euclidean distance",gp=gpar(fontsize=17))),
             arrangeGrob(fviz_nbclust(NbClust(scaled.df, distance="euclidean", method="ward.D")),
                         top=textGrob("Ward.D linkage method with euclidean distance",gp=gpar(fontsize=17))),
             arrangeGrob(fviz_nbclust(NbClust(scaled.df, distance="euclidean", method="ward.D2")),
                         top=textGrob("Ward.D2 linkage method with euclidean distance",gp=gpar(fontsize=17))),
             arrangeGrob(fviz_nbclust(NbClust(scaled.df, distance="euclidean", method="centroid")),
                         top=textGrob("Centroid linkage method with euclidean distance",gp=gpar(fontsize=17))), ncol = 2, nrow = 3)

#### K-means ####
# Elbow method
fviz_nbclust(scaled.df, kmeans, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)+
  labs(subtitle = "K-Means: Elbow method")

# Silhouette method
fviz_nbclust(scaled.df, kmeans, method = "silhouette") +
  labs(subtitle = "K-Means: Silhouette method")

# Gap statistic
fviz_nbclust(scaled.df, kmeans, nstart = 25, method = "gap_stat", nboot = 500) +
  labs(subtitle = "K-Means: Gap statistic method")

# NbClust function
fviz_nbclust(NbClust(scaled.df, distance = "euclidean", min.nc = 2, max.nc = 10, method = "kmeans"))

#### K-Medoids ####
# Elbow method
fviz_nbclust(scaled.df, cluster::pam, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)+
  labs(subtitle = "K-Medoids: Elbow method")

# Silhouette method
fviz_nbclust(scaled.df, cluster::pam, method = "silhouette") +
  labs(subtitle = "K-Medoids: Silhouette method")

# Gap statistic
fviz_nbclust(scaled.df, cluster::pam, method = "gap_stat", nboot = 500) +
  labs(subtitle = "K-Medoids: Gap statistic method")


########## HIERARCHICAL CLUSTERING ############
# Single linkage method
hc.single = hclust(d = dist.eucl, method = "single")
fviz_dend(hc.single, cex = 0.5, main="Single linkage method")
cor(dist.eucl, cophenetic(hc.single))

# Average linkage method
hc.average = hclust(d = dist.eucl, method = "average")
fviz_dend(hc.average, cex = 0.5, main="Average linkage method")
cor(dist.eucl, cophenetic(hc.average))

# Complete linkage method
hc.complete = hclust(d = dist.eucl, method = "complete")
fviz_dend(hc.complete, cex = 0.5, main="Complete linkage method")
cor(dist.eucl, cophenetic(hc.complete))

# Ward.D linkage method
hc.wardd = hclust(d = dist.eucl, method = "ward.D")
fviz_dend(hc.wardd, cex = 0.5, main="Ward.D linkage method")
cor(dist.eucl, cophenetic(hc.wardd))

# Ward.D2 linkage method
hc.wardd2 = hclust(d = dist.eucl, method = "ward.D2")
fviz_dend(hc.wardd2, cex = 0.5, main="Ward.D2 linkage method")
cor(dist.eucl, cophenetic(hc.wardd2))

# Centroid linkage method
hc.centroid = hclust(d = dist.eucl, method = "centroid")
fviz_dend(hc.centroid, cex = 0.5, main="Centroid linkage method")
cor(dist.eucl, cophenetic(hc.centroid))

# Cut the dendrogam (average)
grp.average = cutree(hc.average, k = 3)
fviz_dend(hc.average, k = 3, cex = 0.5, k_colors = c("green", "blue", "red"),
          color_labels_by_k = TRUE, rect = TRUE, main="Average linkage method")

# Cut the dendrogam (complete)
grp.complete = cutree(hc.complete, k = 3)
fviz_dend(hc.complete, k = 3, cex = 0.5, k_colors = c("green", "blue", "red"),
          color_labels_by_k = TRUE, rect = TRUE, main="Complete linkage method")

# Cut the dendrogam (wardd2)
grp.wardd2 = cutree(hc.wardd2, k = 3)
fviz_dend(hc.wardd2, k = 3, cex = 0.5, k_colors = c("green", "blue", "red"),
          color_labels_by_k = TRUE, rect = TRUE, main="Ward.D2 linkage method") 


############# SCATTERPLOT MATRIX ##############
# Scatterplot matrix on the original variables using the Class Variable
pairs(scaled.df, main = "Glucose, Insulin and Sspg for Classes", pch = 21, bg = c("blue", "green", "red")[unclass(diabetes$class)])

# Scatterplot matrix on the original variables with Average
pairs(scaled.df, main = "Glucose, Insulin and Sspg with Average", pch = 19, col=c("green", "blue", "red")[grp.average])

# Scatterplot matrix on the original variables with Complete
pairs(scaled.df, main = "Glucose, Insulin and Sspg with Complete", pch = 19, col=c("green", "blue", "red")[grp.complete])

# Scatterplot matrix on the original variables with Ward.D2
pairs(scaled.df, main = "Glucose, Insulin and Sspg with Ward.D2", pch = 19, col=c("green", "blue", "red")[grp.wardd2])

# Cluster plot of the first 2 PCs Average
fviz_cluster(list(data = scaled.df, cluster = grp.average),
             palette = c("green", "blue", "red"), ellipse.type = "convex",
             repel = TRUE, show.clust.cent = FALSE, ggtheme = theme_classic(),
             main = "Cluster plot the first 2 PCs with Average method")

# Cluster plot of the first 2 PCs Complete
fviz_cluster(list(data = scaled.df, cluster = grp.complete),
             palette = c("green", "blue", "red"), ellipse.type = "convex",
             repel = TRUE, show.clust.cent = FALSE, ggtheme = theme_classic(),
             main = "Cluster plot the first 2 PCs with Complete method")

# Cluster plot of the first 2 PCs Ward.D2
fviz_cluster(list(data = scaled.df, cluster = grp.wardd2),
             palette = c("green", "blue", "red"), ellipse.type = "convex",
             repel = TRUE, show.clust.cent = FALSE, ggtheme = theme_classic(),
             main = "Cluster plot the first 2 PCs with Ward.D2 method")


########## PARTITIONING CLUSTERING ############
################## K-MEANS ####################
# Computing K-means clustering
km.res = eclust(scaled.df, "kmeans", k = 3, nstart = 50, graph = FALSE)

# Print the results
print(km.res)

# Visualizing K-means clusters on the original space
cl <- km.res$cluster
pairs(scaled.df, main = "Cluster Plot K-Means", pch = 19, col=c("green", "blue", "red")[cl])

# Visualize K-means clusters in the PC space
fviz_cluster(km.res, geom = "point", ellipse.type = "norm",
             palette = c("green", "blue", "red"), ggtheme = theme_classic(),
             main = "Cluster plot on the first 2 PCs with K-Means")


################# K-MEDOIDS ###################
# Computing PAM clustering
pam.eucl = pam(scaled.df, 3, metric = "euclidean")
print(pam.eucl)

# Visualizing K-medoids clusters on the original space
cl = pam.eucl$clustering
pairs(scaled.df, main = "Cluster Plot K-Medoids", pch = 19, col=c("green", "blue","red")[cl])

# Visualizing PAM clusters on the space of first 2 PCs
fviz_cluster(pam.eucl, geom = "point", ellipse.type = "norm",
             palette = c("green","blue","red"), ggtheme = theme_classic(),
             main = "Cluster plot on the first 2 PCs with K-Medoids")


############## CLUSTER VALIDATION ################
############ HIERARCHICAL CLUSTERING #############
# Internal validation for hc
hc.res = eclust(scaled.df, "hclust", k = 3, hc_metric = "euclidean", hc_method = "ward.D2", graph = FALSE)

# Silhouette plot for hc
fviz_silhouette(hc.res, palette = c("green","blue","red"), ggtheme = theme_classic())

# Silhouette information for hc
hc.res$silinfo$clus.avg.widths
hc.res$silinfo$avg.width

# Units with negative silhouette widths
neg_sil_index = which(hc.res$silinfo$widths[, "sil_width"] < 0)
hc.res$silinfo$widths[neg_sil_index, , drop = FALSE] 

# Statistics for hc
hc_stats = cluster.stats(dist(scaled.df), hc.res$cluster)

# Dunn index hc
hc_stats$dunn

# External clustering validation
# Confusion matrix hc
table(diabetes$class, hc.res$cluster)

# Agreement between class and HC clusters
class = as.numeric(diabetes$class)
hc.clust_stats = cluster.stats(d = dist(scaled.df), class, hc.res$cluster)

# Corrected Rand index hc
hc.clust_stats$corrected.rand 

# Meila's VI index hc
hc.clust_stats$vi


############ K-MEANS CLUSTERING #############
# Internal validation k-means
km.res = eclust(scaled.df, "kmeans", k = 3, nstart = 50, graph = FALSE)

# Silhouette plot K-means
fviz_silhouette(km.res, palette = c("green","blue","red"), ggtheme = theme_classic())

# Silhouette information
head(km.res$silinfo$widths[, 1:3])
km.res$silinfo$clus.avg.widths
km.res$silinfo$avg.width

# Units with negative silhouette widths
neg_sil_index = which(km.res$silinfo$widths[, "sil_width"] < 0)
km.res$silinfo$widths[neg_sil_index, , drop = FALSE]   

# Statistics for K-means clustering
km_stats = cluster.stats(dist(scaled.df), km.res$cluster)

# Dunn index k-means
km_stats$dunn   

# External clustering validation
# Confusion matrix k-means
table(diabetes$class, km.res$cluster)

# Compute cluster stats
km_stats.clust = cluster.stats(d = dist(scaled.df), class, km.res$cluster)

# Corrected Rand index k-means
km_stats.clust$corrected.rand  

# Meila's VI index k-means
km_stats.clust$vi


############ K-MEDOIDS CLUSTERING #############
# Internal validation k-medoids
pam.res = eclust(scaled.df, "pam", k = 3, graph = FALSE)

# Silhouette plot k-medoids
fviz_silhouette(pam.res, palette = c("green","blue","red"), ggtheme = theme_classic())

# Silhouette information
head(pam.res$silinfo$widths[, 1:3])
pam.res$silinfo$clus.avg.widths
pam.res$silinfo$avg.width

# Units with negative silhouette widths
neg_sil_index = which(pam.res$silinfo$widths[, "sil_width"] < 0)
pam.res$silinfo$widths[neg_sil_index, , drop = FALSE]   

# Statistics for K-medoids clustering
pam_stats = cluster.stats(dist(scaled.df), pam.res$cluster)

# Dunn index k-medoids
pam_stats$dunn   

# External clustering validation
# Confusion matrix k-medoids
table(diabetes$class, pam.res$cluster)

# Compute cluster stats
pam_stats.clust = cluster.stats(d = dist(scaled.df), class, pam.res$cluster)

# Corrected Rand index k-medoids
pam_stats.clust$corrected.rand  

# Meila's VI index k-medoids
pam_stats.clust$vi


########### BEST CLUSTERING ALGORITHMS ############
# Internal measures
intern = clValid(scaled.df, nClust = 2:10, clMethods = c("hierarchical","kmeans","pam"), validation = "internal")
summary(intern)

# Stability measures
stab = clValid(scaled.df, nClust = 2:10, clMethods = c("hierarchical","kmeans","pam"), validation = "stability")
summary(stab)


############# MODEL BASED CLUSTERING ##############
# External information
class = factor(diabetes$class)

# Try if there is relation between Class and underlying clustering
pairs(scaled.df, gap=0, pch = 16, col = as.numeric(class))

# Fit the parsimonious models
mod = Mclust(scaled.df, G = 1:9, modelNames = NULL)

# Visualize the top three BIC models
summary(mod$BIC)

# Plot of the BIC values for all the fitted models
fviz_mclust(mod, "BIC", palette = "jco")

# Visualize the main output
summary(mod)

# Useful quantities
head(round(mod$z, 6))    # Probality to belong to a given cluster

# pairs plot with classification
pairs(scaled.df, gap=0, pch = 16, col = mod$classification)

# Confusion matrix (External Cluster Validation)
table(class, mod$classification)

# Adjusted (or Corrected) Rand Index
adjustedRandIndex(class, mod$classification)

# Visulaing the obtained results in the PC space
# Classification: plot showing the clustering
fviz_mclust(mod, "classification", geom = "point", pointsize = 1.5, palette = c("green","blue","red"))

# Classification uncertainty
# larger symbols indicate the more uncertain observations
fviz_mclust(mod, "uncertainty", palette = c("green","blue","red"))

# Plot the data using only two variables of interest, let say here "glucose" and "insulin"
fviz_mclust(mod, "uncertainty", palette = c("green","blue","red"), choose.vars = c("glucose", "insulin"))

# Plot the data using only two variables of interest, let say here "glucose" and "sspg"
fviz_mclust(mod, "uncertainty", palette = c("green","blue","red"), choose.vars = c("glucose", "sspg"))

# Plot the data using only two variables of interest, let say here "sspg" and "insulin"
fviz_mclust(mod, "uncertainty", palette = c("green","blue","red"), choose.vars = c("sspg", "insulin"))