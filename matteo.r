rm(list=ls())

# Importing data.
load('data/hw2_data.RData')

# Work on a time series of a single person.
data <-  asd_sel$caltech_0051472

# Each row is a time series and each column is a region of the brain
head(data)
library(corrplot)
data_cor <- cor(data)
corrplot(data_cor, method='color', tl.pos='n')
corrplot(data_cor, method='color', tl.pos='n', order="hclust")





















rm(list=ls())


# PACKAGES ----------------------------------------------------------------

library(jpeg)
library(ppcor)
library(igraph)
library(htmlwidgets)
library(animation)



# LOAD DATA ---------------------------------------------------------------

load('data/hw2_data.RData')



# DATA CLEANING -----------------------------------------------------------

scale_datasets_list <- function(ls){
  scaled_list <- list()
  names <- names(ls)
  for(i in 1:length(ls)){
    
    #### Scaling by column.
    patient <- data.frame(apply(ls[[i]], 2, scale))
    colnames(patient) <- colnames(ls[[i]])
    scaled_list[[names[i]]] <-patient
  }
  return(scaled_list)
}

td_sel_scale <- scale_datasets_list(td_sel)
asd_sel_scale <- scale_datasets_list(asd_sel)



# POOL DATASET ------------------------------------------------------------

#### Mean by cell given list of datasets
cells_value_array <- function(ls, i, j){
  cells_value <- c()
  for(patient in ls){
    cells_value <- c(cells_value, patient[i,j])
  }
  return(cells_value)
}


#### Choose metric =c('mean', 'median', 'sd')
summary_dataset <- function(ls, metric='mean'){
  if(metric=='mean'){fun <- mean}
  if(metric=='median'){fun <- median}
  if(metric=='sd'){fun <- sd}
  
  n <- nrow(ls[[1]])
  m <- ncol(ls[[1]])
  mean_data <- matrix(rep(NA,n*m), n, m)
  for(i in 1:n){
    for(j in 1:m){
      mean_data[i,j] <- fun(cells_value_array(ls, i, j))
    }
  }
  data_frame <-data.frame(mean_data)
  colnames(data_frame) <- colnames(ls[[1]])
  return(data_frame)
}


TD <- summary_dataset(td_sel_scale, metric='mean')
ASD <- summary_dataset(asd_sel_scale, metric='mean')




# FUNCTIONS ---------------------------------------------------------------

lower_or_upper <- function(data, bound, cor_type='normal', bonferroni=TRUE){
  
  #### Setting Parameters
  n <- dim(data)[1]
  D <- dim(data)[2]
  alpha <- .05
  m <- choose(D, 2)   # Binomial coefficient
  
  #### Bonferroni Correction
  if(bonferroni == TRUE){ alpha <- alpha / m }
  
  #### Use "Correlation" or "Partial Correlation"
  if(cor_type == 'normal'){
    g <- 0
    corr_matrix <- cor(data) }
  if(cor_type == 'partial'){
    g <- D-2
    corr_matrix <- pcor(data) }
  
  #### Computing Fisher Z-Transform
  Z_j_k_td <- (1/2)*log((1+corr_matrix)/(1-corr_matrix))
  
  #### Confidence intervals for theta
  se <- sqrt(1/( n - g - 3))
  Log_lower <- Z_j_k_td - qnorm(1 - (alpha/2)) * se
  Log_upper <- Z_j_k_td + qnorm(1 - (alpha/2)) * se
  
  #### Confidence intervals for rho
  Lower_bound <- (exp(2*Log_lower) - 1 ) /   ((exp(2*Log_lower) + 1))
  Upper_bound <- (exp(2*Log_upper) - 1 ) /   ((exp(2*Log_upper) + 1)) 
  
  #### Remove NA (on diagonal)
  Lower_bound[is.na(Lower_bound)] = 1
  Upper_bound[is.na(Upper_bound)] = 1
  
  if(bound == 'L'){ return(Lower_bound) }
  if(bound == 'U'){ return(Upper_bound) }
}


adj_matrix_func <- function(list, t, cor_type='normal', bonferroni=TRUE){
  
  #### Find lower and upper for each person
  lowers <- lapply(list, lower_or_upper, 'L', cor_type, bonferroni)
  uppers <- lapply(list, lower_or_upper, 'U', cor_type, bonferroni)
  
  #### Compute mean lower matrix
  Y <- do.call(cbind, lowers)
  Y <- array(Y, dim=c(dim(lowers[[1]]), length(lowers)))
  LOWER <- apply(Y, c(1, 2), mean, na.rm = TRUE)
  
  #### Compute mean upper matrix 
  Y <- do.call(cbind, uppers)
  Y <- array(Y, dim=c(dim(uppers[[1]]), length(uppers)))
  UPPER <- apply(Y, c(1, 2), mean, na.rm = TRUE)
  
  #### Adjacency Matrix
  adj_matrix <- as.matrix(UPPER < -t | LOWER > t)
  
  #### Rename rows and columns as ROI zones
  colnames(adj_matrix) <- colnames(list[[1]])
  rownames(adj_matrix) <- colnames(list[[1]])
  
  return(adj_matrix)
}


cor_matrix_function <- function(list){
  cor_matrices <- lapply(list, cor)
  Y <- do.call(cbind, cor_matrices)
  Y <- array(Y, dim=c(dim(cor_matrices[[1]]), length(cor_matrices)))
  mean_cor_matrix <- apply(Y, c(1, 2), mean, na.rm = TRUE)
  return(mean_cor_matrix)
}



# MATRICES ----------------------------------------------------------------

t <- 0.7

td_adj_normal_cor <- adj_matrix_func(td_sel_scale, t, bonferroni=T)
asd_adj_normal_cor <- adj_matrix_func(asd_sel_scale, t, bonferroni=T)

td_adj_td_partial_cor <- adj_matrix_func(td_sel_scale, t, 'partial')
asd_adj_partial_cor <- adj_matrix_func(asd_sel_scale, t, 'partial')

td_cor_matrix <- cor_matrix_function(td_sel_scale)
asd_cor_matrix <- cor_matrix_function(asd_sel_scale)



# GRAPHS ------------------------------------------------------------------

plot_graphs <- function(adj_mat_1, adj_mat_2,dimensions=2){
  
  #### Check if there are edges
  if(sum(adj_mat_1) == 116 & sum(adj_mat_2)==116){
    return('There are no edges. Try with lower value of t')
  }
  
  #### Create Graphs
  g1 <- graph.adjacency(adj_mat_1, mode = "undirected", diag = FALSE )
  g2 <- graph.adjacency(adj_mat_2, mode = "undirected", diag = FALSE )
  
  
  #### Set nodes colors based on ROIs
  colors <- c('aquamarine', 'chartreuse', 'darkorchid1', 'gold', 'lightpink1', 'orangered', 'salmon', 'seashell2', 'skyblue2')
  
  for(i in 1:length(V(g1)$name)){
    V(g1)$color[i] <- colors[strtoi(substr(V(g1)$name[i],1,1))]
    V(g2)$color[i] <- colors[strtoi(substr(V(g2)$name[i],1,1))] }
  
  
  #### Change edges colors based on if they are present in both graphs or in only one
  for(i in 1:length(E(g1))){
    edge <- as_ids(E(g1)[i])
    if(edge %in% as_ids(E(g2))){E(g1)[i]$color = 'black'}
    else{E(g1)[i]$color = 'blue'} }
  
  
  for(i in 1:length(E(g2))){
    edge <- as_ids(E(g2)[i])
    if(edge %in% as_ids(E(g1))){E(g2)[i]$color = 'black'}
    else{E(g2)[i]$color = 'red'} }
  
  
  #### Import ROI coordinates
  load('data/coordinates.RData')
  coord <- aal116coordinates
  
  if(dimensions == 2){
    #### 2D Plot
    par(mfrow=c(1,2))
    layout <- matrix(c(coord$x.mni, coord$y.mni), 116,2)
    
    my_image <- readJPEG("images/brain.jpg")
    main='3D brain'
    for(graph in list(g1,g2)){
      
      if(identical_graphs(graph, g1)){main = paste0("Brain's ROI correlation\n of TD patients (t=", t, ')')}
      if(identical_graphs(graph, g2)){main = paste0("Brain's ROI correlation\n of ASD patients (t=", t, ')')}
      
      plot(0,0, type='n', xlim=c(-1.2, 1.2), ylim=c(-1.2, 1.2), axes=F, main=main, xlab='', ylab='')
      
      rasterImage(my_image, xleft=-1.2, xright=1.2, ybottom=-1.2, ytop=1.3)
      
      plot(graph, vertex.size=10, vertex.label.cex=.5, vertex.color=V(graph)$color, vertex.shape='circle',
           edge.width=4, edge.color=E(graph)$color, vertex.label.col = 'black',
           layout=layout, add=T) }}
  
  
  if(dimensions == 3){
    #### 3D Plot
    layout <- matrix(c(coord$x.mni, coord$y.mni, coord$z.mni), 116,3)
    
    for(graph in list(g1, g2)){
      rglplot(graph,
              vertex.size=7, vertex.label.cex=.5, vertex.color=V(graph)$color,
              edge.width=4, edge.color=E(graph)$color,
              layout=layout) }}
}


plot_graphs(td_adj_normal_cor, asd_adj_normal_cor,dimensions=2)





















# GIF ---------------------------------------------------------------------
# -------------------------------------------------------------------------


rm(list=ls())


# PACKAGES ----------------------------------------------------------------

library(jpeg)
library(ppcor)
library(igraph)
library(htmlwidgets)



# LOAD DATA ---------------------------------------------------------------

load('data/hw2_data.RData')



# DATA CLEANING -----------------------------------------------------------

scale_datasets_list <- function(ls){
  scaled_list <- list()
  names <- names(ls)
  for(i in 1:length(ls)){
    
    #### Scaling by column.
    patient <- data.frame(apply(ls[[i]], 2, scale))
    colnames(patient) <- colnames(ls[[i]])
    scaled_list[[names[i]]] <-patient
  }
  return(scaled_list)
}

td_sel_scale <- scale_datasets_list(td_sel)
asd_sel_scale <- scale_datasets_list(asd_sel)



# POOL DATASET ------------------------------------------------------------

#### Mean by cell given list of datasets
cells_value_array <- function(ls, i, j){
  cells_value <- c()
  for(patient in ls){
    cells_value <- c(cells_value, patient[i,j])
  }
  return(cells_value)
}


#### Choose metric =c('mean', 'median', 'sd')
summary_dataset <- function(ls, metric='mean'){
  if(metric=='mean'){fun <- mean}
  if(metric=='median'){fun <- median}
  if(metric=='sd'){fun <- sd}
  
  n <- nrow(ls[[1]])
  m <- ncol(ls[[1]])
  mean_data <- matrix(rep(NA,n*m), n, m)
  for(i in 1:n){
    for(j in 1:m){
      mean_data[i,j] <- fun(cells_value_array(ls, i, j))
    }
  }
  data_frame <-data.frame(mean_data)
  colnames(data_frame) <- colnames(ls[[1]])
  return(data_frame)
}


TD <- summary_dataset(td_sel_scale, metric='mean')
ASD <- summary_dataset(asd_sel_scale, metric='mean')






lower_or_upper <- function(data, bound, cor_type='normal', bonferroni=TRUE){
  
  #### Setting Parameters
  n <- dim(data)[1]
  D <- dim(data)[2]
  alpha <- .05
  m <- choose(D, 2)   # Binomial coefficient
  
  #### Bonferroni Correction
  if(bonferroni == TRUE){ alpha <- alpha / m }
  
  #### Use "Correlation" or "Partial Correlation"
  if(cor_type == 'normal'){
    g <- 0
    corr_matrix <- cor(data) }
  if(cor_type == 'partial'){
    g <- D-2
    corr_matrix <- pcor(data) }
  
  #### Computing Fisher Z-Transform
  Z_j_k_td <- (1/2)*log((1+corr_matrix)/(1-corr_matrix))
  
  #### Confidence intervals for theta
  se <- sqrt(1/( n - g - 3))
  Log_lower <- Z_j_k_td - qnorm(1 - (alpha/2)) * se
  Log_upper <- Z_j_k_td + qnorm(1 - (alpha/2)) * se
  
  #### Confidence intervals for rho
  Lower_bound <- (exp(2*Log_lower) - 1 ) /   ((exp(2*Log_lower) + 1))
  Upper_bound <- (exp(2*Log_upper) - 1 ) /   ((exp(2*Log_upper) + 1)) 
  
  #### Remove NA (on diagonal)
  Lower_bound[is.na(Lower_bound)] = 1
  Upper_bound[is.na(Upper_bound)] = 1
  
  if(bound == 'L'){ return(Lower_bound) }
  if(bound == 'U'){ return(Upper_bound) }
}


adj_matrix_func <- function(mat , t, bonf=T){
  L <-  lower_or_upper(mat , "L", bonferroni=bonf)
  U <-  lower_or_upper(mat , "U", bonferroni= bonf)
  adj <- as.matrix(L > t | U < -t)
  return(adj)
}




plot_graphs <- function(mat_1, mat_2, t, dimensions=2, bonf=TRUE){
  
  adj_mat_1 <- adj_matrix_func(mat_1, t = t, bonf=bonf)
  adj_mat_2 <- adj_matrix_func(mat_2, t = t, bonf=bonf)
  
  #Check if there are edges
  if(sum(adj_mat_1) == 116 & sum(adj_mat_2)==116){
    return('There are no edges. Try with lower value of t')
  }
  
  #### Create Graphs
  g1 <- graph.adjacency(adj_mat_1, mode = "undirected", diag = FALSE )
  g2 <- graph.adjacency(adj_mat_2, mode = "undirected", diag = FALSE )
  
  
  #### Set nodes colors based on ROIs
  colors <- c('0','#77DD77', '#C8A2C8', '#FFFF66', '#CB3234', '#E1E1E2', '#ff8c69', '#7fffd4', '#ABCDEF')
  
  for(i in 1:length(V(g1)$name)){
    V(g1)$color[i] <- colors[strtoi(substr(V(g1)$name[i],1,1))]
    V(g2)$color[i] <- colors[strtoi(substr(V(g2)$name[i],1,1))] }
  
  
  #### Change edges colors based on if they are present in both graphs or in only one
  for(i in 1:length(E(g1))){
    edge <- as_ids(E(g1)[i])
    if(edge %in% as_ids(E(g2))){E(g1)[i]$color = 'black'}
    else{E(g1)[i]$color = 'blue'} }
  
  
  for(i in 1:length(E(g2))){
    edge <- as_ids(E(g2)[i])
    if(edge %in% as_ids(E(g1))){E(g2)[i]$color = 'black'}
    else{E(g2)[i]$color = 'red'} }
  
  
  #### Import ROI coordinates
  load('data/coordinates.RData')
  coord <- aal116coordinates
  
  if(dimensions == 2){
    #### 2D Plot
    par(mfrow=c(1,2))
    layout <- matrix(c(coord$x.mni, coord$y.mni), 116,2)
    
    my_image <- readJPEG("images/brain.jpg")
    
    for(graph in list(g1,g2)){
      
      if(identical_graphs(graph, g1)){main = paste0("Brain's ROI correlation\n of TD patients (t=", t, ')')}
      if(identical_graphs(graph, g2)){main = paste0("Brain's ROI correlation\n of ASD patients (t=", t, ')')}
      
      plot(0,0, type='n', xlim=c(-1.2, 1.2), ylim=c(-1.2, 1.2), axes=F, main=main, xlab='', ylab='')
      
      rasterImage(my_image, xleft=-1.2, xright=1.2, ybottom=-1.2, ytop=1.3)
      
      plot(graph, vertex.size=10, vertex.label.cex=.5, vertex.color=V(graph)$color, vertex.shape='circle',
           edge.width=4, edge.color=E(graph)$color, vertex.label.col = 'black',
           layout=layout, add=T) }}
  
  
  if(dimensions == 3){
    #### 3D Plot
    layout <- matrix(c(coord$x.mni, coord$y.mni, coord$z.mni), 116,3)
    
    for(graph in list(g1, g2)){
      rglplot(graph,
              vertex.size=7, vertex.label.cex=.5, vertex.color=V(graph)$color,
              edge.width=4, edge.color=E(graph)$color,
              layout=layout, main=main) }}
}


plot_graphs(TD, ASD ,t = .3 , dimensions=2, bonf=F)






# GIFs --------------------------------------------------------------------




library(magick)
library(dplyr)
library(tidyr)





#### WITH BONF ADJ
grid <- seq(0, 0.6, by=.05)

for(t in grid){
  path <- paste0("images/gif_normal_cor_bonf/t=0_",t*100,".jpg")
  if(t == .05){path <-  paste0("images/gif_normal_cor_bonf/t=0_0",t*100,".jpg")}
  jpeg(file=path, width = 900, height = 560)
  plot_graphs(TD, ASD ,t = t , dimensions=2)
  dev.off()
}



imgs <- list.files("images/gif_normal_cor_bonf/", full.names = TRUE)
img_list <- lapply(imgs, image_read)


## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 1)

## view animated image
img_animated

image_write(image = img_animated,
            path = "images/brain_graph_normal_corr_bonf.gif")








### WITHOUT BONF ADJ

grid <- seq(0, 0.6, by=.05)

for(t in grid){
  path <- paste0("images/gif_normal_cor_no_bonf/t=0_",t*100,".jpg")
  if(t == .05){path <-  paste0("images/gif_normal_cor_no_bonf/t=0_0",t*100,".jpg")}
  jpeg(file=path, width = 900, height = 560)
  plot_graphs(TD, ASD ,t = t , dimensions=2, bonf=F)
  dev.off()
}



imgs <- list.files("images/gif_normal_cor_no_bonf/", full.names = TRUE)
img_list <- lapply(imgs, image_read)


## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 1)

## view animated image
img_animated

image_write(image = img_animated,
            path = "images/brain_graph_normal_corr_no_bonf.gif")














