################################################################################
# Worked example using FlowSOM for clustering and Rtsne for visualization of a 
# high-dimensional mass cytometry (CyTOF) data set
#
# References:
#
# FlowSOM:
# - van Gassen et al. (2015), "FlowSOM: Using self-organizing maps for 
# visualization and interpretation of cytometry data", Cytometry Part A, 
# http://www.ncbi.nlm.nih.gov/pubmed/25573116
# - Bioconductor package: http://www.ncbi.nlm.nih.gov/pubmed/25573116
#
# Rtsne:
# - CRAN package: https://cran.r-project.org/web/packages/Rtsne/index.html
# - additional details: https://github.com/jkrijthe/Rtsne
#
# Data set:
# - We use "sample 01" from the following paper:
# - Samusik et al. (2016), "Automated mapping of phenotype space with 
# single-cell data", Nature Methods, http://www.ncbi.nlm.nih.gov/pubmed/27183440
# - Data file (Samusik_01_notransform.fcs) is also available for download from 
# https://flowrepository.org/id/FR-FCM-ZZPH
#
# Lukas Weber, September 2016
################################################################################


# install packages (if not already installed)

# from Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("flowCore")
biocLite("FlowSOM")

# from CRAN
install.packages("Rtsne")
install.packages("ggplot2")


# load packages

library(flowCore)
library(FlowSOM)
library(Rtsne)
library(ggplot2)




#############################
### LOAD AND PREPARE DATA ###
#############################

# load data using flowCore package (download from link above or data/ folder in
# GitHub repository for this example)

filename <- "data/Samusik_01_notransform.fcs"

data <- flowCore::exprs(flowCore::read.FCS(filename, transformation = FALSE, truncate_max_range = FALSE))

head(data)
dim(data)


# select protein marker columns to use for clustering
marker_cols <- 9:47


# apply arcsinh transformation

# (with standard scale factor of 5 for CyTOF data; alternatively 150 for flow 
# cytometry data; see Bendall et al. 2011, Science, Supplementary Figure S2)

asinh_scale <- 5
data[, marker_cols] <- asinh(data[, marker_cols] / asinh_scale)

summary(data)


# create flowFrame object (required input format for FlowSOM)
data_FlowSOM <- flowCore::flowFrame(data)




###################
### RUN FLOWSOM ###
###################

# set seed for reproducibility
set.seed(123)

# run FlowSOM (initial steps prior to meta-clustering)
out <- FlowSOM::ReadInput(data_FlowSOM, transform = FALSE, scale = FALSE)
out <- FlowSOM::BuildSOM(out, colsToUse = marker_cols)
out <- FlowSOM::BuildMST(out)

# optional visualization
FlowSOM::PlotStars(out)

# extract cluster labels (pre meta-clustering) from output object
labels_pre <- out$map$mapping[, 1]

# specify final number of clusters for meta-clustering (can also be selected 
# automatically, but this often does not perform well)
k <- 40

# run meta-clustering
set.seed(123)
out <- FlowSOM::metaClustering_consensus(out$map$codes, k = k)

# extract cluster labels from output object
labels <- out[labels_pre]

# summary of cluster sizes and number of clusters
table(labels)
length(table(labels))




####################
### SAVE RESULTS ###
####################

# save cluster labels
res <- data.frame(cluster = labels)

write.table(res, file = "results/cluster_labels_FlowSOM.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")




#################
### RUN RTSNE ###
#################

# subsampling (required due to runtime)
n_sub <- 10000

set.seed(123)
ix <- sample(1:length(labels), n_sub)

# prepare data for Rtsne (matrix format required)
data_Rtsne <- data[ix, marker_cols]
data_Rtsne <- as.matrix(data_Rtsne)

head(data_Rtsne)
dim(data_Rtsne)

# remove any near-duplicate rows (required by Rtsne)
dups <- duplicated(data_Rtsne)
data_Rtsne <- data_Rtsne[!dups, ]

dim(data_Rtsne)


# run Rtsne (Barnes-Hut-SNE algorithm; runtime ~2 min)

# note initial PCA step is not required, since we do not have too many
# dimensions (i.e. not thousands, which may be the case in other domains)

set.seed(123)
out_Rtsne <- Rtsne(data_Rtsne, pca = FALSE, verbose = TRUE)




###################
### CREATE PLOT ###
###################

# load cluster labels (if not still loaded)
file_labels <- "results/cluster_labels_FlowSOM.txt"
data_labels <- read.table(file_labels, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
labels <- data_labels[, "cluster"]

# select points used by Rtsne
labels_plot <- labels[ix][!dups]
length(labels_plot)  ## should be same as number of rows in data_Rtsne

# prepare Rtsne output data for plot
data_plot <- as.data.frame(out_Rtsne$Y)
colnames(data_plot) <- c("tSNE_1", "tSNE_2")
head(data_plot)
dim(data_plot)  ## should match length of labels_plot (otherwise labels will not match up correctly)

data_plot[, "cluster"] <- as.factor(labels_plot)

head(data_plot)


# plot 2D t-SNE projection

ggplot(data_plot, aes(x = tSNE_1, y = tSNE_2, color = cluster)) + 
  geom_point(size = 0.2) + 
  coord_fixed(ratio = 1) + 
  ggtitle("t-SNE projection with FlowSOM clustering") + 
  theme_bw()

ggsave("plots/FlowSOM_Rtsne_plot.pdf", height = 6, width = 7)



