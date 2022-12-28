#' ---
#' title:  "Statistical Methods for Data Science I"
#' author: 'Barba Paolo , Candi Matteo'
#' ---

# *********************************************************************** #
#                ** Statistical Methods for Data Science I **             #
#                      ** Homework 2**                 #
# *********************************************************************** #

# Load data ---------------------------------------------------------------
rm(list=ls())
load("~/Documents/GitHub/SDS_HomeWork_2/hw2_data.RData")
load("hw2_data.RData")


# Pool-Together -------------------------------------------------------------------------------------------------------------------------------------

#[TODO]

#### Normal



# per ogni area bisogna fa l'utente medio



#### autistich (like Matteo)


hist(asd_sel$caltech_0051472$`2001`)
hist(asd_sel$trinity_0050234$`2001`)



# Association graph , Confidence intervals-----------------------------------------------------------------------------------------------------------------------------

# Function

L_and_U <- function(data, type){
  n <- dim(data)[1]
  D <- dim(data)[2]
  g <-  D - 2
  alpha <- .05
  m <- choose(D, 2)   # Binomial coefficient
  
  bon_alpha <- alpha / m  # Apply Bonferroni 
  
  correlation_matrix <- cor(data)
  
  Z_j_k_td <- (1/2)*log((1+correlation_matrix)/(1-correlation_matrix))
  
  #### Confidence intervals for theta
  se <- sqrt(1/( n - 3))
  
  Log_lower <- Z_j_k_td - qnorm(1 - (bon_alpha/2)) * se
  Log_upper <- Z_j_k_td + qnorm(1 - (bon_alpha/2)) * se
  
  #### Confidence intervals for rho
  Lower_bound <- (exp(2*Log_lower) - 1 ) /   ((exp(2*Log_lower) + 1))
  Upper_bound <- (exp(2*Log_upper) - 1 ) /   ((exp(2*Log_upper) + 1)) 
  
  Lower_bound[is.na(Lower_bound)] = 1
  Upper_bound[is.na(Upper_bound)] = 1
  
  if(type == 'L'){
    return(Lower_bound)
  }
  if(type == 'U'){
    return(Upper_bound)
  }
}


adj_matrix_func <- function(list, t){
  lowers <- lapply(list, L_and_U, 'L')
  uppers <- lapply(list, L_and_U, 'U')
  
  Y <- do.call(cbind, lowers)
  Y <- array(Y, dim=c(dim(lowers[[1]]), length(lowers)))
  LOWER <- apply(Y, c(1, 2), mean, na.rm = TRUE)
  
  Y <- do.call(cbind, uppers)
  Y <- array(Y, dim=c(dim(uppers[[1]]), length(uppers)))
  UPPER <- apply(Y, c(1, 2), mean, na.rm = TRUE)
  
  adj_matrix <- as.matrix(UPPER < -t | LOWER > t)
  
  colnames(adj_matrix) <- colnames(list[[1]])
  rownames(adj_matrix) <- colnames(list[[1]])
  
  return(adj_matrix)
}


t <- 0.3
adj_matrix_td <- adj_matrix_func(td_sel, t)
adj_matrix_asd <- adj_matrix_func(asd_sel, t)



#library(igraph)
g <- graph.adjacency(adj_matrix_td, mode = "undirected", diag = FALSE )
set.seed(69)
plot.igraph(g, vertex.size = 1, edge.size = 20)


g <- graph.adjacency(adj_matrix_asd, mode = "undirected", diag = FALSE )
set.seed(69)
plot.igraph(g, vertex.size = 1, edge.size = 20)





