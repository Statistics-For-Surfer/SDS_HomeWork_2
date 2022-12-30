rm(list=ls())

# Importing data.
load('data/hw2_data.RData')
# PACKAGES ----------------------------------------------------------------

library(jpeg)
library(ppcor)
library(igraph)

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



adj_matrix_func <- function(mat , t){
  L <-  lower_or_upper(mat , "L")
  U <-  lower_or_upper(mat , "U")
  adj <- as.matrix(L > t | U < -t)
  return(adj)
}








# MATRICES ----------------------------------------------------------------



td_adj_normal_cor <- adj_matrix_func(td_sel_scale, t, bonferroni=T)
asd_adj_normal_cor <- adj_matrix_func(asd_sel_scale, t, bonferroni=T)

td_adj_td_partial_cor <- adj_matrix_func(td_sel_scale, t, 'partial')
asd_adj_partial_cor <- adj_matrix_func(asd_sel_scale, t, 'partial')




# GRAPHS ------------------------------------------------------------------

plot_graphs <- function(mat_1 , mat_2 ,t, dimensions=2){
  
  adj_mat_1 <- adj_matrix_func(mat_1, t = t)
  adj_mat_2 <- adj_matrix_func(mat_2, t = t)
  
  #Check if there are edges
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




plot_graphs(TD, ASD ,t = .6 , dimensions=2)

