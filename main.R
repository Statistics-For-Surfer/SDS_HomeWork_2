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


# Pool-Together -------------------------------------------------------------------------------------------------------------------------------------

#[TODO]

#### Normal



# per ogni area bisogna fa l'utente medio



#### autistich (like Matteo)


hist(asd_sel$caltech_0051472$`2001`)
hist(asd_sel$trinity_0050234$`2001`)



# Association graph , Confidence intervals-----------------------------------------------------------------------------------------------------------------------------

# only use one subjects for group since pool together isn't done






td <- td_sel$caltech_0051487
asd <- asd_sel$caltech_0051472

dim(td)
n <- dim(td)[1]
D <- dim(td)[2]
alpha <- .05
m <-  choose(D, 2)   # Binomial coefficient
m

bon_alpha <- alpha / m 


correlation_matrix_td <- cor(td)
correlation_matrix_asd <- cor(asd)
Z_j_k_td <- (1/2)*log((1+correlation_matrix_td)/(1-correlation_matrix_td))


g <-  D - 2

#### Confidence intervals for theta

Log_lower <- Z_j_k_td - qnorm(1 - (bon_alpha/2)) * sqrt(1/( n - g - 3))
Log_upper <- Z_j_k_td + qnorm(1 - (bon_alpha/2)) * sqrt(1/( n - g - 3))




#### Confince intervals for rho

Lower_bound <- (exp(2*Log_lower) - 1 ) /   ((exp(2*Log_lower) + 1))
Upper_bound <- (exp(2*Log_upper) - 1 ) /   ((exp(2*Log_upper) + 1)) 

Lower_bound[1,2]
Upper_bound[1,2]


# Build the adjcency matrix ---------------------------------------------------------------------------------------------------------------------

# Idea: Set  a threshold t , if abs(t) is inside the interval ---> put an edge between the two ROI.
t <- .3



t <- .004


adj_matrix <- as.matrix(Upper_bound < -t | Lower_bound > t)


adj_matrix[is.na(adj_matrix)] = 1





g <- graph.adjacency(adj_matrix, mode = "undirected", 
                     diag = FALSE , )
set.seed(69)
plot(g, vertex.size = 5 , vertex.size2 = 0.2, edge.size = 5)



#install.packages("igraph")
library(igraph)
edges <- graph_from_adjacency_matrix(
  adj_matrix,
  mode = "undirected",
  weighted = NULL,
  diag = TRUE,
  add.colnames = NULL,
  add.rownames = NA
)
h <- graph.empty() + names(asd)
h <- h + edges
plot(h)


# analysis together -----------------------------------------------------------------------------------------------------------------------------
rm(list = ls())


matrices <-  lapply(td_sel, cor)
Y <- do.call(cbind, matrices)


Y <- array(Y, dim=c(dim(matrices[[1]]), length(matrices)))


mean_rho <- apply(Y, c(1, 2), mean, na.rm = TRUE)









