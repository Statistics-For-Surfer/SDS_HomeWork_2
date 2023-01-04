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
library(magick)
library(dplyr)
library(tidyr)


par(col='white')
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






lower_or_upper <- function(data, data2 = NULL , bound, cor_type='normal', bonferroni=TRUE, delta = F){
  
  
  #### Setting Parameters
  n <- dim(data)[1]
  D <- dim(data)[2]
  alpha <- .05
  m <- choose(D, 2)   # Binomial coefficient
  
  #### Bonferroni Correction
  if(bonferroni == TRUE){ alpha <- alpha / m }
  if(isFALSE(delta)){
    #### Use "Correlation" or "Partial Correlation"
    if(cor_type == 'normal'){
      g <- 0
      corr_matrix <- cor(data) }
    if(cor_type == 'partial'){
      g <- D-2
      corr_matrix <- pcor(data)$estimate
    }}
  
  if(delta){
    #### Use "Correlation" or "Partial Correlation"
    if(cor_type == 'normal'){
      g <- 0
      corr_matrix1 <- cor(data)
      corr_matrix2 <- cor(data2)
      corr_matrix <- corr_matrix1 - corr_matrix2}
    if(cor_type == 'partial'){
      g <- D-2
      corr_matrix1 <- pcor(data)$estimate
      corr_matrix2 <- pcor(data2)$estimate
      corr_matrix <- corr_matrix1 - corr_matrix2}}
  
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
  
  if(bound == 'L'){ 
    Log_lower <- Z_j_k_td - qnorm(1 - (alpha/2)) * se
    
    #### Confidence intervals for rho
    Lower_bound <- (exp(2*Log_lower) - 1 ) /   ((exp(2*Log_lower) + 1))
    
    
    #### Remove NA (on diagonal)
    Lower_bound[is.na(Lower_bound)] = 1
    
    return(Lower_bound) }
  if(bound == 'U'){ 
    
    Log_upper <- Z_j_k_td + qnorm(1 - (alpha/2)) * se
    
    #### Confidence intervals for rho
    Upper_bound <- (exp(2*Log_upper) - 1 ) /   ((exp(2*Log_upper) + 1)) 
    
    #### Remove NA (on diagonal)
    Upper_bound[is.na(Upper_bound)] = 1
    return(Upper_bound) }
}


adj_matrix_func <- function(mat , t, bonf=T, cor_type='normal', delta = FALSE , data2 = NULL){
  L <-  lower_or_upper(mat , "L", bonferroni= bonf, cor_type=cor_type , delta = delta, data2 = data2)
  U <-  lower_or_upper(mat , "U", bonferroni= bonf, cor_type=cor_type , delta = delta, data2 = data2)
  adj <- as.matrix(L > t | U < -t)
  return(adj)
}


plot_graphs <- function(mat_1, mat_2, t, dimensions=2, bonf=TRUE, cor_type='normal', delta=F){
  
  if(delta){
    adj_mat_1 <- adj_matrix_func(mat_1, t = t, bonf=bonf, cor_type = cor_type, delta = delta, data2 = mat_2)
    
    #Check if there are edges
    # if(sum(adj_mat_1) == 116 & sum(adj_mat_2)==116){
    #   return('There are no edges. Try with lower value of t')
    # }
    
    colnames(adj_mat_1) <- colnames(mat_1)
    
    #### Create Graphs
    g1 <- graph.adjacency(adj_mat_1, mode = "undirected", diag = FALSE )
    
    #### Set nodes colors based on ROIs
    colors <- c('0','#77DD77', '#C8A2C8', '#FFFF66', '#CB3234', '#E1E1E2', '#ff8c69', '#7fffd4', '#ABCDEF')
    
    for(i in 1:length(V(g1)$name)){
      V(g1)$color[i] <- colors[strtoi(substr(V(g1)$name[i],1,1))] }
    
    
    load('data/coordinates.RData')
    coord <- aal116coordinates
    
    if(dimensions == 2){
      #### 2D Plot
      par(mfrow=c(1,1))
      layout <- matrix(c(coord$x.mni, coord$y.mni), 116,2)
      
      my_image <- readJPEG("images/brain.jpg")
      
      main = paste0("Brain's ROI correlation\n of DELTA (t=", t, ')')
        
      plot(0,0, type='n', xlim=c(-1.2, 1.2), ylim=c(-1.2, 1.2), axes=F, main=main, xlab='', ylab='')
      
      rasterImage(my_image, xleft=-1.2, xright=1.2, ybottom=-1.2, ytop=1.3)
      
      plot(g1, vertex.size=10, vertex.label.cex=.5, vertex.color=V(g1)$color, vertex.shape='circle',
           edge.width=4, edge.color='black', vertex.label.col = 'black',
           layout=layout, add=T) }
    
    
    if(dimensions == 3){
      #### 3D Plot
      layout <- matrix(c(coord$x.mni, coord$y.mni, coord$z.mni), 116,3)
      
      
      rglplot(g1,
              vertex.size=7, vertex.label.cex=.5, vertex.color=V(g1)$color,
              edge.width=4, edge.color=E(g1)$color,
              layout=layout, main=main) }
    
  }
  else{
    
    adj_mat_1 <- adj_matrix_func(mat_1, t = t, bonf=bonf, cor_type = cor_type)
    adj_mat_2 <- adj_matrix_func(mat_2, t = t, bonf=bonf, cor_type=cor_type)
    
    #Check if there are edges
    # if(sum(adj_mat_1) == 116 & sum(adj_mat_2)==116){
    #   return('There are no edges. Try with lower value of t')
    # }
    
    colnames(adj_mat_1) <- colnames(mat_1)
    colnames(adj_mat_2) <- colnames(mat_2)
    
    #### Create Graphs
    g1 <- graph.adjacency(adj_mat_1, mode = "undirected", diag = FALSE )
    g2 <- graph.adjacency(adj_mat_2, mode = "undirected", diag = FALSE )
  
    #### Set nodes colors based on ROIs
    colors <- c('0','#77DD77', '#C8A2C8', '#FFFF66', '#CB3234', '#E1E1E2', '#ff8c69', '#7fffd4', '#ABCDEF')
    
    for(i in 1:length(V(g1)$name)){
      V(g1)$color[i] <- colors[strtoi(substr(V(g1)$name[i],1,1))]
      V(g2)$color[i] <- colors[strtoi(substr(V(g2)$name[i],1,1))] }
    
    #### Change edges colors based on if they are present in both graphs or in only one
    
    if(sum(adj_mat_1) != 116){
    for(i in 1:length(E(g1))){
      edge <- as_ids(E(g1)[i])
      if(edge %in% as_ids(E(g2))){E(g1)[i]$color = 'black'}
      else{E(g1)[i]$color = 'blue'} }}
    
    if(sum(adj_mat_2)!=116){
    for(i in 1:length(E(g2))){
      edge <- as_ids(E(g2)[i])
      if(edge %in% as_ids(E(g1))){E(g2)[i]$color = 'black'}
      else{E(g2)[i]$color = 'red'} }}
    
    load('data/coordinates.RData')
    coord <- aal116coordinates
    
    if(dimensions == 2){
      #### 2D Plot
      par(mfrow=c(1,2))
      layout <- matrix(c(coord$x.mni, coord$y.mni), 116,2)
      
      my_image <- readJPEG("images/brain.jpg")
      
      for(graph in list(g1,g2)){
        
        if(identical_graphs(graph, g1)){main = paste0("Brain's ROI correlation\n TD patients (t=", t, ')')}
        if(identical_graphs(graph, g2)){main = paste0("Brain's ROI correlation\n ASD patients (t=", t, ')')}
        
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
                layout=layout, main=main) }}}
  
}

plot_graphs(TD, ASD, 0, delta=F, bonf=T, cor_type = 'partial')






TD <- summary_dataset(td_sel_scale, metric='median')
ASD <- summary_dataset(asd_sel_scale, metric='median')



grid_t <- seq(0, 0.3, by=.05)

for(t in grid_t){
  path <- paste0("images/handmade/", t*100, ".jpeg")
  if(t > 0 & t < .1){path <-  paste0("images/handmade/0", t*100, ".jpeg")}
  jpeg(file=path, width = 450, height = 560)
  plot_graphs(TD, ASD, t, delta=T, bonf=T, cor_type = 'partial')
  dev.off()
  }

imgs <- list.files("images/handmade/", full.names = TRUE)
img_list <- lapply(imgs, image_read)
img_joined <- image_join(img_list)
img_animated <- image_animate(img_joined, fps = 1)
unlink(paste0("images/handmade/", '*'))

# image_browse(img_animated)


image_write(image = img_animated,
            path = paste0("images/delta_partial_bonf.gif")) 













