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






rm(list=ls())


# PACKAGES ----------------------------------------------------------------

library(jpeg)
library(ppcor)
library(igraph)
library(JointNets)



# DATA --------------------------------------------------------------------

load('hw2_data.RData')



# FUNCTIONS ---------------------------------------------------------------

lower_or_upper <- function(data, bound, cor_type='normal'){
  
  #### Setting Parameters
  n <- dim(data)[1]
  D <- dim(data)[2]
  alpha <- .05
  m <- choose(D, 2)   # Binomial coefficient
  
  #### Bonferroni Correction
  bon_alpha <- alpha / m  
  
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
  Log_lower <- Z_j_k_td - qnorm(1 - (bon_alpha/2)) * se
  Log_upper <- Z_j_k_td + qnorm(1 - (bon_alpha/2)) * se
  
  #### Confidence intervals for rho
  Lower_bound <- (exp(2*Log_lower) - 1 ) /   ((exp(2*Log_lower) + 1))
  Upper_bound <- (exp(2*Log_upper) - 1 ) /   ((exp(2*Log_upper) + 1)) 
  
  #### Remove NA (on diagonal)
  Lower_bound[is.na(Lower_bound)] = 1
  Upper_bound[is.na(Upper_bound)] = 1
  
  if(bound == 'L'){ return(Lower_bound) }
  if(bound == 'U'){ return(Upper_bound) }
}


adj_matrix_func <- function(list, t, cor_type='normal'){
  
  #### Find lower and upper for each person
  lowers <- lapply(list, lower_or_upper, 'L', cor_type)
  uppers <- lapply(list, lower_or_upper, 'U', cor_type)
  
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

t <- 0.4

td_adj_normal_cor <- adj_matrix_func(td_sel, t)
td_adj_td_partial_cor <- adj_matrix_func(td_sel, t, 'partial')

asd_adj_normal_cor <- adj_matrix_func(asd_sel, t)
asd_adj_partial_cor <- adj_matrix_func(asd_sel, t, 'partial')

td_cor_matrix <- cor_matrix_function(td_sel)
asd_cor_matrix <- cor_matrix_function(asd_sel)



# GRAPHS ------------------------------------------------------------------

plot_graphs <- function(adj_mat_1, adj_mat_2,dimensions=2){
    
  
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
  data(aal116coordinates)
  coord <- aal116coordinates
  
  if(dimensions == 2){
  #### 2D Plot
  par(mfrow=c(1,2))
  layout <- matrix(c(coord$x.mni, coord$y.mni), 116,2)
  
  my_image <- readJPEG("brain.jpg")
  
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


plot_graphs(td_adj_normal_cor, asd_adj_normal_cor,dimensions=2)
