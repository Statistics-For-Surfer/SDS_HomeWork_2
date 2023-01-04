# PACKAGES ----------------------------------------------------------------

library(jpeg)
library(ppcor)
library(igraph)
library(magick)
library(dplyr)
library(tidyr)



# DATA --------------------------------------------------------------------

### Create a function to scale by column the datasets
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


### Pool the data together creating 2 dataframes: TD and ASD
cells_value_array <- function(ls, i, j){
  cells_value <- c()
  for(patient in ls){
    cells_value <- c(cells_value, patient[i,j])
  }
  return(cells_value)
}

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



# DATA DESCRIPTION ----------------------------------------------------------------

### Histograms of correlations coefficient
corr_distr_fun <- function(data_1 , data_2 , metric){
  pTD  <-  summary_dataset(data_1 , metric = metric)
  pASD <-  summary_dataset( data_2, metric = metric)
  mat  <-  cor(pTD)
  mat2 <-  cor(pASD)
  upper_tri_1 <- upper.tri(mat)
  mat_upper_tri_1 <- mat[upper_tri_1]
  upper_tri_2 <- upper.tri(mat2)
  mat_upper_tri_2 <- mat[upper_tri_2]
  par(mfrow = c(1,2))
  hist( mat_upper_tri_1,
        main = "Histograms of correlation \n coefficients for TD subject",
        xlab = expression(hat(rho)) , 
        border = "white" , 
        col = "#E52B50"  , freq = FALSE )
  box()
  
  hist(mat_upper_tri_2 ,
       main = "Histograms of correlation\n coefficients for ASD subject",
       xlab = expression(hat(rho)) , 
       border = "white" , 
       col = "#E52B50", freq = FALSE )
  box()
  par(mfrow = c(1,1))
}

### Quantile distribution function
quantile_distr_func <- function(data1 , data2 , metric){
  pTD  <-  summary_dataset(data1 , metric = metric)
  pASD <-  summary_dataset(data2, metric = metric)
  mat  <-  cor(pTD)
  mat2 <-  cor(pASD)
  upper_tri_1 <- upper.tri(mat)
  mat_upper_tri_1 <- mat[upper_tri_1]
  upper_tri_2 <- upper.tri(mat2)
  mat_upper_tri_2 <- mat[upper_tri_2]
  q1   <-  sort(abs(mat_upper_tri_1))
  q2   <-  sort(abs(mat_upper_tri_2))
  par(mfrow = c(1,2))
  plot(q1 , main = "Quantile distribution of correlation \n  coefficient for TD subject", ylab = "Quantile" , 
       col = "#ffff66" , cex = .3)
  grid()
  plot(q2 , main = "Quantile distribution of correlation \n  coefficient for ASD subject", ylab = "Quantile",
       col = "#ffff66" , cex = .3)
  grid()
  par(mfrow = c(1,1))
}




# ANALYSIS ----------------------------------------------------------------

### Confidence intervals
lower_or_upper <- function(data, data2 = NULL , bound, cor_type='normal', bonferroni=TRUE, delta = F ){
  
  #### Setting Parameters
  n <- dim(data)[1]
  D <- dim(data)[2]
  alpha <- .05
  m <- choose(D, 2)   # Binomial coefficient
  
  #### Bonferroni Correction
  if(bonferroni == TRUE){ alpha <- alpha / m }
  if (delta == F){
    #### Use "Correlation" or "Partial Correlation"
    if(cor_type == 'normal'){
      g <- 0
      corr_matrix <- cor(data) }
    if(cor_type == 'partial'){
      g <- D-2
      corr_matrix <- pcor(data)$estimate
    }}
  if (delta == T){
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

### Adjacency matrix
adj_matrix_func <- function(mat , t, bonf=T, cor_type='normal', delta = NULL , data2 = NULL){
  L <-  lower_or_upper(mat , "L", bonferroni= bonf, cor_type=cor_type , delta = delta, data2 = data2 )
  U <-  lower_or_upper(mat , "U", bonferroni= bonf, cor_type=cor_type , delta = delta, data2 = data2)
  adj <- as.matrix(L > t | U < -t)
  return(adj)
}

### Plot function for the graphs
plot_graphs <- function(mat_1, mat_2, t, dimensions=2, bonf=TRUE, cor_type='normal'){
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

### Confidence intervals comparing.
relation_CI <- function(data, cor_type='normal'){
  a_norm <-  lower_or_upper(data , bound ="L" , cor_type = cor_type , bonferroni= F)
  b_norm <-  lower_or_upper(data , bound ="U" , cor_type = cor_type , bonferroni= F)
  a_bonf <-  lower_or_upper(data , bound ="L" , cor_type = cor_type , bonferroni= T)
  b_bonf <- lower_or_upper(data , bound ="U" , cor_type = cor_type , bonferroni= T)
  
  upper_tri_a_norm <- upper.tri(a_norm)
  upper_tri_b_norm <- upper.tri(b_norm)  
  upper_tri_a_bonf <- upper.tri(a_bonf)
  upper_tri_b_bonf <- upper.tri(b_bonf)
  
  a_norm <- a_norm[upper_tri_a_norm]
  b_norm <- b_norm[upper_tri_b_norm]
  a_bonf <- a_bonf[upper_tri_a_bonf]
  b_bonf <- b_bonf[upper_tri_b_bonf]
  
  ls <-  sort(a_norm)
  us <-  sort(b_norm)
  ls_bon <- sort(a_bonf) 
  us_bon <- sort(b_bonf)
  xs <- 1:length(ls)
  
  mat  <-  cor(data)
  upper_tri_1 <- upper.tri(mat)
  mat_upper_tri_1 <- mat[upper_tri_1]
  mat = sort(mat_upper_tri_1)
  
  
  par(mfrow = c(1,1))
  plot(ls, type = "n"  , main = expression(paste("Asyntotic confidence intervals of ", rho , " \n with & without Bonferroni correction")), 
       xlab = "" , ylab = "Lower and Upper Bound" , ylim = c(-1,1))
  segments(x0 = xs , y0 = ls_bon , x1 = xs , y1 = us_bon,  col = "lightblue")
  segments(x0 = xs , y0 = ls , x1 = xs , y1 = us,  col = "gold")
  
  legend("topleft" , c("With Bonferroni correction" , "Without Bonferroni correction" , "Person's coefficient") , col = c("lightblue" , "#FFFF66", "#CB3234" ) , lty = 1 , lwd = 3 , bty = "n" )
  
  points( mat , col = "#CB3234", cex =  .3)
  
  
}
relation_CI_par_pear <- function(data, bonferroni = T){
  a_norm <-  lower_or_upper( data , bound = "L" , cor_type = 'normal' , bonferroni= bonferroni )
  b_norm <-  lower_or_upper(data , bound = "U" , cor_type = 'normal' , bonferroni= bonferroni )
  a_par <-  lower_or_upper(data , bound ="L" , cor_type = 'partial', bonferroni= bonferroni)
  b_par <- lower_or_upper(data , bound ="U" , cor_type =   'partial' , bonferroni= bonferroni)
  
  upper_tri_a_norm <- upper.tri(a_norm)
  upper_tri_b_norm <- upper.tri(b_norm)  
  upper_tri_a_par <- upper.tri(a_par)
  upper_tri_b_par <- upper.tri(b_par)
  
  a_norm <- a_norm[upper_tri_a_norm]
  b_norm <- b_norm[upper_tri_b_norm]
  a_par <- a_par[upper_tri_a_par]
  b_par <- b_par[upper_tri_b_par]
  
  
  ls <-  sort(a_norm)
  us <-  sort(b_norm)
  ls_par <- sort(a_par) 
  us_par <- sort(b_par)
  xs <- 1:length(ls)
  
  mat  <-  cor(data)
  mat2 <-  pcor(data)$estimate
  upper_tri_1 <- upper.tri(mat)
  mat_upper_tri_1 <- mat[upper_tri_1]
  upper_tri_2 <- upper.tri(mat2)
  mat_upper_tri_2 <- mat2[upper_tri_1]
  mat1 <- sort(mat_upper_tri_1)
  mat2 <- sort(mat_upper_tri_2)
  
  
  
  
  main = "Asyntotic confidence intervals for \n Partial  correlation coefficient without \n Bonferroni correction"
  if(bonferroni){main = "Asyntotic confidence intervals for \n Partial  correlation coefficient with \n Bonferroni correction"}
  
  
  
  plot(ls, type = "n" , main =  main , 
       xlab = "" , ylab = "Lower and Upper Bound" , ylim = c(-1,1))
  
  segments(x0 = xs , y0 = ls_par , x1 = xs , y1 = us_par,  col = "#f5bf80")
  segments(x0 = xs , y0 = ls , x1 = xs , y1 = us,  col = "#77DD77")
  legend("topleft" , c("With Partial correlation" , "With Pearson correlation" , "Person's correlation" , "Partial correlation") , col = c("#f5bf80" , "#77DD77","#CB3234" ,"white") , lty = 1 , lwd = 2 , bty = "n" )
  points(mat1, col = "#CB3234", cex =  .3)
  
  points(mat2, col = "white", cex =  .3)
}

### Gif generator function
gif_generator <- function(grid_t, TD, ASD, delta, bonf, cor_type){
  for(t in grid_t){
    path <- paste0("images/handmade/", t*100, ".jpeg")
    if(t > 0 & t < .1){path <-  paste0("images/handmade/0", t*100, ".jpeg")}
    jpeg(file=path, width = 900, height = 560)
    plot_graphs(TD, ASD, t, delta=T, bonf=T, cor_type = 'partial')
    dev.off()
  }
  
  imgs <- list.files("images/handmade/", full.names = TRUE)
  img_list <- lapply(imgs, image_read)
  img_joined <- image_join(img_list)
  img_animated <- image_animate(img_joined, fps = 1)
  unlink(paste0("images/handmade/", '*'))

  image_write(image = img_animated,
              path = paste0("images/gif.gif")) }



