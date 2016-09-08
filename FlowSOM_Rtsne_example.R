#########################################################################################
# Worked example showing how to use FlowSOM for clustering and Rtsne for visualization of
# a mass cytometry (CyTOF) data set
#
# References:
#
# FlowSOM:
# - van Gassen et al. (2015), "FlowSOM: Using self-organizing maps for visualization and 
# interpretation of cytometry data", Cytometry Part A, 87/7, 636-645, 
# http://www.ncbi.nlm.nih.gov/pubmed/25573116
# - Bioconductor package: http://www.ncbi.nlm.nih.gov/pubmed/25573116
#
# Rtsne:
# - available from CRAN: https://cran.r-project.org/web/packages/Rtsne/index.html
# - additional details: https://github.com/jkrijthe/Rtsne
#
# Data set:
# - Healthy human bone marrow data set "Marrow1" from Figure 1b in the following paper:
# - Amir et al. (2013), "viSNE enables visualization of high dimensional single-cell data
# and reveals phenotypic heterogeneity of leukemia", Nature Biotechnology, 
# http://www.ncbi.nlm.nih.gov/pubmed/23685480
# - Link to raw data: http://www.c2b2.columbia.edu/danapeerlab/html/viSNE-data.html 
# (filename "visne_marrow1.fcs")
#
# Lukas Weber, September 2016
#########################################################################################


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

# load data using flowCore package (download from link above or data/ folder in this repository)

filename <- "data/visne_marrow1.fcs"

data_flowFrame <- flowCore::read.FCS(filename, transformation = FALSE, truncate_max_range = FALSE)
data_flowFrame

data <- flowCore::exprs(data_flowFrame)
params <- flowCore::parameters(data_flowFrame)
desc <- flowCore::description(data_flowFrame)

head(data)
dim(data)


# select following protein expression columns to use in FlowSOM calculations:
# CD11b, CD123, CD19, CD20, CD3, CD33, CD34, CD38, CD4, CD45, CD45RA, CD8, CD90
# (see Amir et al. 2013, Supplementary Tables 1 and 2)

cols_to_use <- c(11, 23, 10, 16, 7, 22, 14, 28, 12, 6, 8, 13, 30)
unname(colnames(data))[cols_to_use]  # check


# apply arcsinh transformation

cols_proteins <- 4:36

asinh_scale <- 5
data[, cols_proteins] <- asinh(data[, cols_proteins] / asinh_scale)

summary(data)


# create flowFrame object (required as input format for FlowSOM)

data_FlowSOM <- flowCore::flowFrame(exprs = data, parameters = params, description = desc)




###################
### RUN FLOWSOM ###
###################

# run initial steps in FlowSOM workflow (prior to meta-clustering)

# set seed for reproducibility
set.seed(123)

# run FlowSOM
out <- FlowSOM::ReadInput(data_FlowSOM, transform = FALSE, scale = FALSE)
out <- FlowSOM::BuildSOM(out, colsToUse = cols_to_use)
out <- FlowSOM::BuildMST(out)

# optional visualization
FlowSOM::PlotStars(out)

# extract cluster labels (pre meta-clustering) from output object
labels_pre <- out$map$mapping[, 1]


# run meta-clustering (final step in FlowSOM workflow)

# set final number of clusters (can also be set automatically, but this does not perform 
# well for many data sets)
k <- 20

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

# subsample points for plot (required for runtime)

n_sub <- 10000

set.seed(123)
ix <- sample(1:length(labels), n_sub)


# prepare data for Rtsne (matrix format required)

data_Rtsne <- data[ix, cols_to_use]
data_Rtsne <- as.matrix(data_Rtsne)

class(data_Rtsne)
head(data_Rtsne)
dim(data_Rtsne)


# remove near-duplicate rows (required by Rtsne)

dups <- duplicated(data_Rtsne)
data_Rtsne <- data_Rtsne[!dups, ]

dim(data_Rtsne)


# run Rtsne (Barnes-Hut-SNE algorithm; runtime ~2 min)
# without PCA step (see Amir et al. 2013, Online Methods, "viSNE analysis")

set.seed(123)
out_Rtsne <- Rtsne(data_Rtsne, pca = FALSE, verbose = TRUE)




###################
### CREATE PLOT ###
###################

# load cluster labels from FlowSOM (if not still loaded)
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
  geom_point(size = 0.75, alpha = 0.5) + 
  coord_fixed(ratio = 1) + 
  ggtitle("t-SNE projection with FlowSOM clustering") + 
  theme_bw()

ggsave("plots/FlowSOM_Rtsne_plot.pdf", height = 6, width = 7)



