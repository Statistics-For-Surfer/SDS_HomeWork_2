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
    corr_matrix <- cor(data)
  }
  if(cor_type == 'partial'){
    g <- D-2
    v <- diag(var(data))
    corr_matrix <- - cov(data) / sqrt(v %*% t(v))
    }

  
  
  #### Computing Fisher Z-Transform
  Z_j_k_td <- (1/2)*log((1+corr_matrix)/(1-corr_matrix))
  
  #### Confidence intervals for theta
  se <- sqrt(1/( n - g - 3))
  Log_lower <- Z_j_k_td - qnorm(1 - (bon_alpha/2)) * se
  Log_upper <- Z_j_k_td + qnorm(1 - (bon_alpha/2)) * se

  #### Confidence intervals for rho
  Lower_bound <- (exp(2*Log_lower) - 1 ) /   ((exp(2*Log_lower) + 1))
  Upper_bound <- (exp(2*Log_upper) - 1 ) /   ((exp(2*Log_upper) + 1)) 
  
  #### Remove na (on diagonal)
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

t <- 0.03

td_adj_normal_cor <- adj_matrix_func(td_sel, t)
td_adj_td_partial_cor <- adj_matrix_func(td_sel, t, 'partial')

asd_adj_normal_cor <- adj_matrix_func(asd_sel, t)
asd_adj_partial_cor <- adj_matrix_func(asd_sel, t, 'partial')

td_cor_matrix <- cor_matrix_function(td_sel)
asd_cor_matrix <- cor_matrix_function(asd_sel)



# GRAPHS ------------------------------------------------------------------

#library(igraph)

g <- graph.adjacency(td_adj_normal_cor, mode = "undirected", diag = FALSE )
set.seed(69)
plot.igraph(g, vertex.size = 1, edge.size = 20)


g <- graph.adjacency(td_adj_td_partial_cor, mode = "undirected", diag = FALSE )
set.seed(69)
plot.igraph(g, vertex.size = 1, edge.size = 20)


