rm(list=ls())

# Importing data.
load('hw2_data.RData')

# Work on a time series of a single person.
data <-  asd_sel$caltech_0051472

# Each row is a time series and each column is a region of the brain
head(data)
library(corrplot)
data_cor <- cor(data)
corrplot(data_cor, method='color', tl.pos='n')
corrplot(data_cor, method='color', tl.pos='n', order="hclust")
