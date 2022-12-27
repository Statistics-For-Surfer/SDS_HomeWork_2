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



# Association graph , Confidence intervals-----------------------------------------------------------------------------------------------------------------------------

# only use one subjects for group since pool together isn't done


td <- td_sel$caltech_0051487
asd <- asd_sel$caltech_0051472

dim(td)
n <- dim(td)[1]
D <- dim(td)[2]





alpha <- .05

m <-  choose(D, 2)
bonferroni_alpha <- alpha / m 

Z_j_k <- (1/2)*log((1+cor(td))/(1-cor(td)))
g <-  D - 2

Log_lower <- Z_j_k - qnorm(1 - alpha/2) * sqrt(1/( n - g - 3))
Log_upper <- Z_j_k + qnorm(1 - alpha/2) * sqrt(1/( n - g - 3))

Lower_bound <- (exp(2*Log_lower) - 1 ) /   ((exp(2*Log_lower) + 1))
Upper_bound <- (exp(2*Log_upper) - 1 ) /   ((exp(2*Log_upper) + 1)) 





# Build the adjcency matrix ---------------------------------------------------------------------------------------------------------------------

# Idea: Set  a threshold t , if abs(t) is inside the interval ---> put an edge between the two ROI.
t <- .3


edges <- function( l, u,t){
  if(is.nan(l)){
    return(0)
  }
  t1 <- -t
  t2 <- t
  l <- Lower_bound[10,8]
  u <- Upper_bound[10,8]

  
  if (t2 > l ){
    return(0)
  }
  if (t1 < u){
    return(0)
  }
  return(1)
}

edges(Lower_bound[10,8],Upper_bound[10,8] , .01 )


adj <- matrix(NA , D, D)
t <- .01
for ( i in 1:D){
  for( j in 1:D){
    result <- edges(Lower_bound[i,j] , Upper_bound[i,j], t)
    adj[i,j] <- result
    
    
    
  }   
  }
adj

