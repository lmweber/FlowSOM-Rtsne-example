#########################################################################################
# Short example showing how to run FlowSOM on a mass cytometry data set
#
# References:
#
# FlowSOM paper:
# Van Gassen et al. (2015), "FlowSOM: Using self-organizing maps for visualization and 
# interpretation of cytometry data", Cytometry Part A, 87/7, 636-645, 
# http://www.ncbi.nlm.nih.gov/pubmed/25573116
#
# FlowSOM R package on Bioconductor:
# http://bioconductor.org/packages/release/bioc/html/FlowSOM.html
#
# Data set:
# Healthy human bone marrow data set "Marrow1" from Figure 1b in the following paper:
# Amir et al. (2013), "viSNE enables visualization of high dimensional single-cell data 
# and reveals phenotypic heterogeneity of leukemia", Nature Biotechnology, 
# http://www.ncbi.nlm.nih.gov/pubmed/23685480
#
# The data set can be downloaded from the Dana Pe'er lab page at: 
# http://www.c2b2.columbia.edu/danapeerlab/html/viSNE-data.html
# (filename "visne_marrow1.fcs")
#
# Lukas M. Weber, December 2015
#########################################################################################


# install packages from Bioconductor (skip if already installed)

source("https://bioconductor.org/biocLite.R")
biocLite("flowCore")
biocLite("FlowSOM")


# install packages from CRAN (skip if already installed)

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

# load data using flowCore package (data source: see above)

fn <- "data/visne_marrow1.fcs"

data_ff <- flowCore::read.FCS(fn, transformation = FALSE)
data_ff

data <- flowCore::exprs(data_ff)
params <- flowCore::parameters(data_ff)
desc <- flowCore::description(data_ff)

head(data)
dim(data)


# select 13 protein expression columns to use in FlowSOM calculations
# CD11b, CD123, CD19, CD20, CD3, CD33, CD34, CD38, CD4, CD45, CD45RA, CD8, CD90
# (see Amir et al. 2013, Supplementary Tables 1 and 2)

cols_to_use <- c(11, 23, 10, 16, 7, 22, 14, 28, 12, 6, 8, 13, 30)
unname(colnames(data))[cols_to_use]  # check


# arcsinh transformation

cols_proteins <- 4:36

asinh_scale <- 5
data[, cols_proteins] <- asinh(data[, cols_proteins] / asinh_scale)


# normalization

# Note that I have not included a normalization step, since each dimension is on roughly
# the same scale, making normalization unnecessary. However an optional normalization
# step could be included here (percentile normalization, z-score normalization, etc).

summary(data)


# create flowFrame object

# Note: FlowSOM requires input data to be in the form of a flowFrame object, which can
# also be created automatically with "flowCore::read.FCS". However this makes it more
# complicated to apply the arcsinh transform (see flowCore vignette).

data <- flowCore::flowFrame(exprs = data, parameters = params, description = desc)



###################
### RUN FLOWSOM ###
###################

# run FlowSOM (with set.seed for reproducibility)

set.seed(123)

out_fSOM <- FlowSOM::ReadInput(data, transform = FALSE, scale = FALSE)
out_fSOM <- FlowSOM::BuildSOM(out_fSOM, colsToUse = cols_to_use)
out_fSOM <- FlowSOM::BuildMST(out_fSOM)


# plot

FlowSOM::PlotStars(out_fSOM)


# extract cluster labels from output object

str(out_fSOM$map)
head(out_fSOM$map$mapping)
dim(out_fSOM$map$mapping)

labels <- out_fSOM$map$mapping[, 1]


# cluster sizes and number of clusters (before metaclustering)

table(labels)
length(table(labels))


# run metaclustering step

# set number of metaclusters (can also be set automatically)
k <- 15

set.seed(123)
out_meta <- FlowSOM::metaClustering_consensus(out_fSOM$map$codes, k = k)


# extract metacluster labels from output object

out_meta
labels_meta <- out_meta[labels]


# cluster sizes and number of clusters (with metaclustering)

table(labels_meta)
length(table(labels_meta))



####################
### SAVE RESULTS ###
####################

# save cluster labels

res <- data.frame(cluster = labels_meta)

write.table(res, file = "results/cluster_labels_FlowSOM.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")



#################################
### VISUALIZE WITH T-SNE PLOT ###
#################################

# subsample points for plot

n <- length(labels_meta)
n_sub <- 10000
set.seed(123)
ix <- sample(1:n, n_sub)

# prepare data for Rtsne

data_rtsne <- flowCore::exprs(data)
data_rtsne <- data_rtsne[ix, cols_to_use]
data_rtsne <- as.matrix(data_rtsne)  # Rtsne requires matrix format

dim(data_rtsne)

# remove near-duplicate rows (required by Rtsne)

dups <- duplicated(data_rtsne)
data_rtsne <- data_rtsne[!dups, ]

dim(data_rtsne)


# run Rtsne (Barnes-Hut-SNE algorithm)
# without PCA step (see Amir et al. 2013, Online Methods, "viSNE analysis")

set.seed(123)
out_rtsne <- Rtsne(data_rtsne, pca = FALSE, verbose = TRUE)

# plot 2D t-SNE projection

data_plot <- as.data.frame(out_rtsne$Y)
colnames(data_plot) <- c("tSNE_1", "tSNE_2")

labels_plot <- labels_meta[ix][!dups]
length(labels_plot)

data_plot$cluster <- as.factor(labels_plot)

head(data_plot)

ggplot(data_plot, aes(x = tSNE_1, y = tSNE_2, color = cluster)) + 
  geom_point(size = 1.25) + 
  coord_fixed(ratio = 1) + 
  ggtitle("tSNE projection: FlowSOM clustering") + 
  theme_bw()

ggsave("plots/tSNE_plot.pdf", height = 7, width = 8)


