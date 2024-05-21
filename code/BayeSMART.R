
# ***README***
# The following script is used to run the revised iIMPACT method for spatial transcriptomics data on multi-sample
# proposed in the submitted manuscript titled "Integrating Image and Molecular Profiles for Spatial Transcriptomics Analysis"
# ***END***


# load cpp functions
Rcpp::sourceCpp('./BayeSMART.cpp')

# load other packages
library(mvtnorm)
library(scater)
library(scran)
library(MASS)
library(readr)


# Summarize neighborhood information: Extract the indices of all neighbors for each spot
# loc: a n*2 matrix for x and y coordinates of spots, where n is the number of spots
# n_neighbor: number of neighbors for lattice
get.neighbor <- function(loc, n_neighbor, tune = 2){
  P <- matrix(0, nrow = nrow(loc), ncol = n_neighbor)
  loc <- as.matrix(loc)
  if (n_neighbor == 4 | n_neighbor == 8){
    loc <- round(loc)
    aa <- sqrt(2)
  } else if (n_neighbor == 6){
    aa <- sqrt(3)
  } else {aa <- 1.2}

  dist_matrix <- vectorized_pdist(loc, loc)
  min_dist <- min(dist_matrix[dist_matrix > 0])


  if (n_neighbor == 8){
    dist_threshold <- min_dist*(aa - 1)*tune + min_dist*aa
  } else{
    dist_threshold <- min_dist*(aa - 1)*tune + min_dist
  }

  for (i in 1:nrow(loc)){
    vec <- dist_matrix[i,]
    valid_indices <- which(vec > 0 & vec < dist_threshold)
    sorted_valid_indices <- valid_indices[order(vec[valid_indices])]
    final_ind <- rep(0, n_neighbor)

    if(length(sorted_valid_indices) >= n_neighbor){
      smallest_indices <- sorted_valid_indices[1:n_neighbor]
    }else{
      final_ind[seq_along(sorted_valid_indices)] <- sorted_valid_indices
      smallest_indices <- final_ind
    }

    P[i,] <- smallest_indices
    }
  return(P)
}

# Function to concat the spatial information of each sample together
# xys: the list containing all the locations of the spots in each sample
# n_neighbor: number of neighbors of each spot
spatial_concat <- function(xys, n_neighbor, tune = 2){
  G <- matrix(ncol = n_neighbor, nrow = 0)
  G_origin <- matrix(ncol = n_neighbor, nrow = 0)
  each_length <- sapply(xys, nrow)
  
  for (s in 1:length(xys)){
    G1 <- get.neighbor(xys[[s]], n_neighbor, tune = tune)
    G_origin <- rbind(G_origin, G1)
    if (s > 1){
      non_zero_indices <- G1 != 0
      G1[non_zero_indices] <- G1[non_zero_indices] + sum(each_length[1:(s-1)])
    }
    G <- rbind(G, G1)
  }
  return (list(G = G, G_origin = G_origin))
}


# Summarize cell abundance information: Extract the number of cells for different cell types for each spot
# cell_loc: a m*2 matrix for x and y coordinates of cells, where m is the number of cells
# cell_type: a vector with length m to record the cell types for all cells
# spot_loc: a n*2 matrix for x and y coordinates of spots
# lattice: lattice type, 'hexagon' or 'square'
get.cell.abundance <- function(cell_loc, cell_type, spot_loc, lattice = 'hexagon'){
  
  n_cell <- dim(cell_loc)[1]
  n_spot <- dim(spot_loc)[1]
  
  dist_spot <- vectorized_pdist(as.matrix(spot_loc), as.matrix(spot_loc))
  
  min_dist <- min(dist_spot[dist_spot > 0])
  
  if (lattice == 'hexagon'){
    threshold <- min_dist*1.01*sqrt(3)/2
  }
  else if(lattice == 'square'){
    threshold <- min_dist*1.01 / sqrt(2)
  }
  else {stop('Parameter lattice: Please enter a correct value.')}
  
  cell_names <- unique(cell_type)
  V <- matrix(0, ncol = length(cell_names) , nrow = n_spot)
  colnames(V) <- cell_names
  
  # assign cells to spots
  count <- 0
  for (i in 1:n_cell){
    ne_class <- which(cell_names == cell_type[i])
    dist_matrix <- vectorized_pdist(as.matrix(cell_loc[i,]), as.matrix(spot_loc))
    
    list_spot <- which(dist_matrix <= threshold)
    for (index_spot in list_spot){
      V[index_spot, ne_class] <- V[index_spot, ne_class] + 1
    }
    
    if (floor(i*100/n_cell) == count)
    {
      print(paste0(count, '% has been done'))
      count <- count + 10
    }
  }
  return(V)
}


# Function to concat the cell-type information of each sample together from the images
# cell_list: the list containing all the cell information of the spots in each sample obtained from the paired images
# xys: the list containing all the locations of the spots in each sample
# lattice: lattice: lattice type, 'hexagon' or 'square'

image_concat <- function(cell_list, xys, lattice = 'square'){
  V <- NULL
  all_cell <- c()
  for (s in 1:length(cell_list)){
    all_cell <- c(all_cell, cell_list[[s]]$name)
  }
  cell_types <- unique(all_cell)
  
  for (s in 1:length(cell_list)){
    cell_info <- cell_list[[s]]
    
    cell_xy <- cell_info[, 1:2]
    cell_name <- cell_info$name
    
    V1 <- get.cell.abundance(cell_xy, cell_name, xys[[s]], lattice = lattice)
    new_V <- matrix(0, nrow = nrow(V1), ncol = length(cell_types))
    colnames(new_V) <- cell_types
    indices <- match(cell_types, colnames(V1))
    
    for (i in seq_along(indices)) {
      index <- indices[i]
      if (!is.na(index)) {  # If the name is found in the column names of V
        new_V[, i] <- V1[, index]
      }
    }
    if (is.null(V)){
      V <- new_V
    }else{
      V <- rbind(V, new_V)
    }
  }
  return (V)
}



# run BayeSMART
# V: spot-level cell abundance (n*q matrix, q is the number of cell types)
# Y: low-dimensional representation of gene expression by function 'batch_remove'
# G: neighborhood information for each spot
# n_cluster: number of spatial domains
# w: scaling parameter for image profile representing the weight of the image information
# label_switch_refer, an index of column for label switching, usually choose column index of tumor cell type
run.BayeSMART <- function(V, Y, G, n_cluster, f_val = 1, w = 1/20, label_switch_refer = 1, lambda = 0){
  n_spot <- nrow(Y)
  n_cell_type <- ncol(V)
  n_PC <- ncol(Y)
  set.seed(123)
  
  e <- rep(1, n_cluster)
  f <- rep(f_val, n_spot)
  omega_initial <- rep(1/n_cell_type, n_cell_type)
  alpha <- rep(1, n_cell_type)
  tau <- 0.01
  eta <- rep(0, n_PC)
  sigma_initial <- diag(1, ncol = n_PC, nrow = n_PC)
  mu_initial <- rmvnorm(1, eta, sigma_initial/tau)
  Z_initial <- matrix(sample(1:n_cluster, n_spot, replace = T), ncol = 1)
  alpha_gamma <- 0.1
  beta_gamma <- diag(0.1, ncol = n_PC, nrow = n_PC)
  
  # run iIMPACT
  result <- BayeSMART(V, Y, G, n_cluster, e, f, eta, tau,  alpha_gamma, beta_gamma, mu_initial, sigma_initial, alpha, omega_initial, Z_initial,rep(w, n_spot), lambda = lambda)
  
  Z_p <- result[["Z"]]
  Z_p_old <- Z_p
  omega_p <- result[['omega']]
  mean_p <- result[['mu']]
  
  n_iter <- nrow(mean_p)
  K <- n_cluster
  
  omega_p_new <- array(0, dim = c(floor(n_iter/2), K, n_cell_type))
  mean_p_new <- array(0, dim = c(floor(n_iter/2), K, n_PC))
  
  # switch label
  for (i in (n_iter - floor(n_iter/2) + 1):n_iter){
    mean_now <- omega_p[i,seq(label_switch_refer, (K - 1)*n_cell_type + label_switch_refer, n_cell_type)]
    for (j in 1:K){
      Z_p[i, which(Z_p_old[i, ] == order(mean_now)[j])] <- j
      old_index = order(mean_now)[j]
      new_index = j
      omega_p_new[i - (n_iter - floor(n_iter/2)), new_index, ] <- omega_p[i, ((old_index - 1)*n_cell_type + 1):(old_index*n_cell_type)]
      mean_p_new[i - (n_iter - floor(n_iter/2)), new_index, ] <- mean_p[i, ((old_index - 1)*n_PC + 1):(old_index*n_PC)]
    }
  }
  
  print('100% has been done')
  
  return(list(Z = Z_p[(n_iter - floor(n_iter/2) + 1):n_iter, ], omega = omega_p_new, mu = mean_p_new, cell_type = colnames(V), n_cluster = K))
}


# Get spatial domain identification results
get.spatial.domain <- function(result){
  Z_p <- result[['Z']]
  spatial_domain <- apply(Z_p, 2, getmode)
  return(spatial_domain)
}



# Get domain-level cell proportion
get.domain.cell.prop <- function(result){
  omega_p <- result[['omega']]
  domain_cell_proportion <- apply(omega_p, c(2, 3), median)
  colnames(domain_cell_proportion) <- result[['cell_type']]
  
  return(domain_cell_proportion)
}


# Get the list of clustering result of each sample
domain_split <- function(spatial_domain, G_origin, xys){
  domain_list <- NULL
  each_length <- sapply(xys, nrow)
  for (s in 1:length(xys)){
    l = each_length[s] ###### revised
    
    if (s == 1){ ######
      domains <- spatial_domain[1:l]
      G1 <- G_origin[1:l, ]
    }else{
      sum_ = sum(each_length[1: (s-1)])
      domains <- spatial_domain[(sum_ + 1): (sum_ + l)]
      G1 <- G_origin[(sum_ + 1): (sum_ + l), ]
    }
    spatial_domain_refined <- refine.cluster(G1, domains, area_unit = 3)
    domain_list[[s]] <- spatial_domain_refined
  }
  return (domain_list)
}


# function to refine the clustering results
# P: Extracted indices of all neighbors for each spot from get.neighbor function
# cluster: vector of spatial domain result
# area_unit: the threshold of number of spots for small area
refine.cluster <- function(P, cluster, area_unit = 1){
  n_spot <- nrow(P)
  visited <- rep(0, n_spot)
  connected_area <- list()
  k <- 1
  # bfs to find all connected areas
  for (i in 1:n_spot){
    if (visited[i] == 0){
      visited[i] <- 1
      temp <- c(i)
      now_cluster <- cluster[i]
      now_area_list <- c(i)
      while (length(temp) > 0){
        now <- temp[1]
        temp <- temp[-1]
        for (j in P[now, ]){
          if (j != 0){
            if (visited[j] == 0 & cluster[j] == now_cluster){
              visited[j] <- 1
              now_area_list <- c(now_area_list, j)
              temp <- c(temp, j)
            }
          }
        }
      }
      connected_area[[k]] <- now_area_list
      k <- k + 1
    }
  }
  
  n_area <- length(connected_area)
  
  # change the cluster for small areas
  cluster_new <- cluster
  for (i in 1:n_area){
    now_area_list <- connected_area[[i]]
    if (length(now_area_list) <= area_unit){
      # find all neighbors of the current connected area
      neighbor_list <- c()
      for (j in now_area_list){
        neighbor_list <- c(neighbor_list, P[j, P[j, ]!= 0])
      }
      neighbor_list <- setdiff(neighbor_list, now_area_list)
      # cluster of neighbor spots
      neighbor_cluster <- unique(cluster[neighbor_list])
      if (length(neighbor_cluster) == 1){
        cluster_new[now_area_list] <- neighbor_cluster[1]
      }}}
  return(cluster_new)
}



# function to create grids for single-cell technologies and get profiles
# count: the single-cell gene count matrix
# location: the x and y axis of each cell
# x_grid, y_grid: the number of spots along the x_axis and y_axis, respectively
get_profiles <- function(count, location, celltype, x_grid = 30, y_grid = 30){
  x_min <- min(location[,1])
  x_max <- max(location[,1])
  y_min <- min(location[,2])
  y_max <- max(location[,2])
  x_inc <- (x_max - x_min) / x_grid
  y_inc <- (y_max - y_min) / y_grid
  
  # store molecular profile M
  num_of_spots <- x_grid * y_grid
  cells <- unique(celltype)
  M <- matrix(0, nrow = num_of_spots, ncol = ncol(count))
  V <- matrix(0, nrow = num_of_spots, ncol = length(cells))
  spot_xy <- matrix(0, nrow = num_of_spots, ncol = 2)
  cell_assign <- rep(NA, length(celltype))
  
  spot_ind <- 1
  visited <- integer(0)
  for (i in 1:x_grid){
    x_left <- x_min + (i-1)*x_inc
    x_right <- x_left + x_inc
    
    for (j in 1:y_grid){
      y_below <- y_min + (j-1)*y_inc
      y_above <- y_below + y_inc
      
      # save spot location
      spot_xy[spot_ind,] <- c(mean(x_left, x_right), mean(y_below, y_above))
      
      indices_to_check <- setdiff(1:length(celltype), visited)
      
      for (k in indices_to_check){
        x <- location[k, 1]
        y <- location[k, 2]
        if (x >= x_left & x<= x_right & y >= y_below & y <= y_above){
          
          cell_assign[k] <- spot_ind
          
          # update M
          M[spot_ind,] <- M[spot_ind,] + count[k,]
          
          # update V
          for (l in 1:length(cells)){
            if (celltype[k] == cells[l]){
              V[spot_ind, l] <- V[spot_ind, l] + 1
            }
          }
          visited <- c(visited, k)
        }
      }
      
      spot_ind <- spot_ind + 1
    }
  }
  colnames(M) <- colnames(count)
  colnames(V) <- cells
  
  # get empty grids
  empty_spot <- which(rowSums(V) == 0)
  
  # update G
  G <- generate_adjacent_matrix(x_grid, y_grid)
  
  # # remove empty grids
  G <- adjust_matrix(G, empty_spot)
  if (length(empty_spot) > 0){
    M <- M[-empty_spot, ]
    V <- V[-empty_spot, ]
    spot_xy <- spot_xy[-empty_spot, ]
  }
  
  new_assignment <- adjust_cell_assignment(cell_assign)
  
  # normalize M based on No. of cells
  M <- M / rowSums(V)
  
  return(list(M = M, V = V, G = G, empty_spot = empty_spot, spot_xy = spot_xy, cell_assign = new_assignment))
}



# function for preprocessing of the molecular profile, with gene count matrix of all the samples combinded together first
# and library size normalization& log-transformation conducted, then SVGs or HVGs are selected, and PCA is applied and
# HARMONY is used for batch correction across samples
# list_count: a list containing all the count matrix of each sample
# each_list: a vector of the number of spots in each sample
# xys: a list containing all the spot location for each sample
# geneSelect: method used for selecting SVGs or HVGs
# n_gene: number of genes selected as SVGs or HVGs
# pcn: number of principle components kept
batch_remove <- function(list_count, xys=NULL, gene_select = "hvgs", n_gene = 2000, pcn = 3){
  set.seed(42)
  cnt_all <- do.call(cbind, list_count)
  # 1.Library size normalization + log transformation
  cnt_all <- scater::normalizeCounts(cnt_all, log = TRUE)
  each_length <- sapply(list_count, ncol)
  if(gene_select == "sparkx"){
    genes <- lapply(1:length(each_length), function(l){
      capture.output(sparkx.l <- SPARK::sparkx(list_count[[l]], xys[[l]]))
      sparkx.l <- sparkx.l$res_mtest[order(sparkx.l$res_mtest$adjustedPval), ]
      genes.l <- head(rownames(sparkx.l), n = n_gene)
    })
    genes <- unique(unlist(genes))
    cnt_all <- cnt_all[genes, ]
  } else if(gene_select == "hvgs"){
    dec <- scran::modelGeneVar(cnt_all)
    genes <- scran::getTopHVGs(dec, n = n_gene)
    cnt_all <- cnt_all[genes, ]
  }else if(gene_select == "none"){
    cnt_all <- cnt_all
  }
  
  idx_rm <- apply(cnt_all, 1, sum) == 0
  cnt_all <- cnt_all[!idx_rm, ]
  
  # 3.Dimension reduction
  cnt_all <- apply(cnt_all, MARGIN = 1, scale, 
                 center = T, scale = T)
  Q <- prcomp(cnt_all, scale. = F)$rotation
  cnt_all_f <- (cnt_all %*% Q)[, 1:pcn]
  
  # library(irlba)
  # svd_result <- irlba(cnt_all, nv = pcn)
  # Q <- svd_result$v
  # cnt_all_f <- cnt_all %*% Q
  
  # 4.Batch effect correction
  set.seed(42)
  cnt_all_f2 <- harmony::HarmonyMatrix(
    data_mat = cnt_all_f,
    meta_data = as.character(rep(1:length(each_length), each_length)), 
    do_pca = F, verbose = F)
  
  return (cnt_all_f2)
}




######################
# Helper functions
######################

# helper function to calculate distance
vectorized_pdist <- function(A,B){
  an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
  bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))
  m = nrow(A)
  n = nrow(B)
  tmp = matrix(rep(an, n), nrow=m)
  tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
  sqrt( tmp - 2 * tcrossprod(A,B) )
}

getmode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


# helper function for the HER2-positive dataset: process the cell-typing from Hover-Net
new.cell.abundance <- function(V){
  mat <- matrix(NA, nrow = nrow(V), ncol = 3)
  colnames(mat) <- c("tumor", "stroma", "immune")
  mat[,1] <- V[, "neopla"] + V[, "necros"]
  mat[,2] <- V[,"connec"] + V[, "no-neo"]
  mat[,3] <- V[, "inflam"]
  
  return (mat)
}


# helper function to generate spatial profile on single-cell dataset
generate_adjacent_matrix <- function(m1, m2) {
  neighbors <- function(i, j) {
    n <- c()
    for (x in -1:1) {
      for (y in -1:1) {
        if (!(x == 0 && y == 0) && (i + x) >= 1 && (i + x) <= m1 && (j + y) >= 1 && (j + y) <= m2) {
          n <- c(n, (i + x - 1) * m2 + j + y)
        }
      }
    }
    return(n)
  }
  
  G <- list()
  for (i in 1:m1) {
    for (j in 1:m2) {
      G[[length(G) + 1]] <- neighbors(i, j)
    }
  }
  return(G)
}


# helper function to generate spatial profile on single-cell dataset
adjust_matrix <- function(G, nums, ncol = 8) {
  # Remove the rows that correspond to the grids in nums
  if (length(nums) > 0){
    G <- G[-nums]
    
    # Create a mapping for old index to new index
    removed_so_far <- 0
    index_mapping <- list()
    for (i in 1:(length(G) + length(nums))) {
      if (i %in% nums) {
        removed_so_far <- removed_so_far + 1
      } else {
        index_mapping[[as.character(i)]] <- i - removed_so_far
      }
    }
    
    # Adjust the remaining grids based on the new indices
    for (i in 1:length(G)) {
      G[[i]] <- sapply(G[[i]], function(x) ifelse(!(x %in% nums), index_mapping[[as.character(x)]], NA))
      G[[i]] <- na.omit(G[[i]])  # remove NAs (corresponding to removed indices)
    }
  }
  
  # convert G from list to matrix
  matt <- matrix(0, nrow = length(G), ncol = ncol)
  for (i in 1:length(G)){
    v <- G[[i]]
    if (length(v) < ncol){
      v <- c(v, rep(0, ncol - length(v)))
    }
    matt[i,] <- v
  }
  return (matt)
}


# helper function to assign each cell back to the generated grids(spots)
adjust_cell_assignment <- function(cell_assign){
  
  non_empty_spots <- unique(cell_assign)
  spot_mapping <- setNames(seq_along(non_empty_spots), sort(non_empty_spots))
  new_cell_assign <- as.integer(spot_mapping[as.character(cell_assign)])
  return(new_cell_assign)
}