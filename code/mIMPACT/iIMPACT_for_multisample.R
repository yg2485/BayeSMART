
# ***README***
# The following script is used to run the revised iIMPACT method for spatial transcriptomics data on multi-sample
# proposed in the submitted manuscript titled "Integrating Image and Molecular Profiles for Spatial Transcriptomics Analysis"
# ***END***


# load cpp functions
Rcpp::sourceCpp('/Users/yanghongguo/Desktop/Research/iBurst/code/mIMPACT/iIMPACT_for_multisample.cpp')

# load other packages
library(DirichletReg)
library(mvtnorm)
library(LaplacesDemon)
library(SingleCellExperiment)
library(scater)
library(scran)
library(MASS)

# Summarize neighborhood information: Extract the indices of all neighbors for each spot
# loc: a n*2 matrix for x and y coordinates of spots, where n is the number of spots
# n_neighbor: number of neighbors for lattice 
get.neighbor <- function(loc, n_neighbor){
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
    dist_threshold <- min_dist*(aa - 1)*0.9 + min_dist*aa
  } else{
    dist_threshold <- mean_dist*(aa - 1)*0.5 + mean_dist
  }
  
  # print(min_dist)
  # print(dist_threshold)
  
  for (i in 1:nrow(loc)){
    k <- 1
    for (j in 1:nrow(loc)){
      if (dist_matrix[i, j] > 0 & dist_matrix[i, j] < dist_threshold){
        P[i, k] <- j
        k <- k + 1
      }}}
  return(P)
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


# Process the gene expression count matrix: Normalize and get the low-dimensional representation of gene expression profile via PCA
# count: a n*p matrix for gene expression counts, where n is the number of spots and p is the number of genes
# n_PC: number of principle components for the gene expression profile
process.gene.expression <- function(count, n_PC = 3, n_HVG = 2000, method = 'PCA', d){
  n_gene <- ncol(count)
  n_spot <- nrow(count)
  
  if (method == "PCA"){
    if (is.null(colnames(count)) == FALSE){
      rowData <- data.frame(gene = colnames(count))} else {rowData <- data.frame(gene = 1:n_gene)}
    
    colnames(count) <- rowData$gene
    
    sce <- SingleCellExperiment(assays = list(counts = as(t(count), 'dgCMatrix')), rowData = rowData)
    
    sce <- logNormCounts(sce)
    
    n_HVG <- min(n_HVG, n_gene)
    
    set.seed(100)
    dec <- scran::modelGeneVar(sce)
    top <- scran::getTopHVGs(dec, n = n_HVG)
    
    set.seed(101)
    sce <- scater::runPCA(sce, subset_row = top, ncomponents = n_PC)
    
    Y <- sce@int_colData@listData$reducedDims@listData$PCA
  }
  
  else if (method == "t-SNE"){
    # library(FactoMineR)
    # library(Rtsne)
    # res.pca <- PCA(count, ncp = 50, graph = FALSE)
    # YY <- Rtsne(res.pca$ind$coord, dims=3, perplexity=30)
    # Y <- YY$Y
    if (is.null(colnames(count)) == FALSE){
      rowData <- data.frame(gene = colnames(count))} else {rowData <- data.frame(gene = 1:n_gene)}
    
    colnames(count) <- rowData$gene
    
    sce <- SingleCellExperiment(assays = list(counts = as(t(count), 'dgCMatrix')), rowData = rowData)
    
    sce <- logNormCounts(sce)
    
    n_HVG <- min(n_HVG, n_gene)
    
    set.seed(100)
    dec <- scran::modelGeneVar(sce)
    top <- scran::getTopHVGs(dec, n = n_HVG)
    
    set.seed(101)
    sce <- scater::runTSNE(sce, subset_row = top, ncomponents = n_PC)
    
    Y <- reducedDim(sce, "TSNE")
  }
  
  else if (method == 'UMAP'){
    set.seed(d)
    library(uwot)
    Y <- umap(count, n_components = n_PC, n_neighbors = 15, min_dist = 0.1, metric = 'euclidean')
  }
  
  else {stop('Parameter method: Please enter a correct value.')}

  return(Y)
}


# Create square lattice for imaging-based SRT data
# cell_info: a m*3 matrix for x and y coordinates, and the cell type for all cells measured on the tissue, where m is the number of cells.
# size: grid size to create the lattice
create.grid <- function(cell_info, size = NULL){
  
  n_cell <- dim(cell_info)[1]
  if (is.null(size)){
    size <- round(mean(c((max(cell_info[, 'x']) - min(cell_info[, 'x'])) / 15, (max(cell_info[, 'y']) -min(cell_info[, 'y'])) / 15)))
  }
  x_total <- floor(max(cell_info[, 'x']) / size) + 1
  y_total <- floor(max(cell_info[, 'y']) / size) + 1
  
  new_loc <- NULL
  k <- 1
  cell_assignment <- numeric(n_cell)
  
  for (i in 1:x_total){
    for (j in 1:y_total){
      
      cell_list <- c()
      
      for (ii in 1:n_cell){
        if (cell_info[ii, 'x'] >= (i - 1)*size & cell_info[ii, 'x'] <  i * size){
          if (cell_info[ii, 'y'] >= (j - 1)*size  & cell_info[ii, 'y'] < j * size){
            cell_list <- c(cell_list, ii)
          }
        }
      }
      if (length(cell_list) > 0){
        x_mid <- round((2 * i - 1)/2*size)
        y_mid <- round((2 * j - 1)/2*size)
        
        new_loc <- rbind(new_loc, c(x_mid, y_mid))
        cell_assignment[cell_list] <- k
        k <- k + 1
      }
    }
  }
  colnames(new_loc) <- c('x', 'y')
  
  return(list(spot_loc = new_loc, cell_assignment = cell_assignment))
}

# Process the gene expression count matrix and cell information for imaging_based SRT data
# Generate cell abundance and low-dimensional representation of gene expression profile via PCA on spot level
# count: a m*p matrix for gene expression counts, where m is the number of cells and p is the number of genes
# cell_info: a m*3 matrix for x and y coordinates, and the cell type for all cells measured on the tissue, where m is the number of cells.
# cell_assignment: results from create_grid function
# n_PC: number of principle components for the gene expression profile
process.imaging.based.SRT <- function(count, cell_info, cell_assignment, n_PC = 3){
  n_spot <- max(cell_assignment)
  cell_type <- cell_info[ ,3]
  Y_cell <- process.gene.expression(count, n_PC)
  
  cell_names <- unique(cell_type)
  V <- matrix(0, ncol = length(cell_names) , nrow = n_spot)
  colnames(V) <- cell_names
  Y <- matrix(0, ncol = n_PC , nrow = n_spot)
  
  for (i in 1:n_spot){
    for (j in which(cell_assignment == i)){
      V[i, cell_type[j]] <- V[i, cell_type[j]] + 1
    }
    if (length(which(cell_assignment == i)) > 1){
      Y[i, ] <- apply(Y_cell[cell_assignment == i,], 2, mean)
    }
    else {Y[i, ] <- Y_cell[which(cell_assignment == i),]}
  }
  
  return(list(Y = Y, V = V))
}


# run iIMPACT
# V: spot-level cell abundance (n*q matrix, q is the number of cell types)
# Y: low-dimensional representation of gene expression profile via PCA
# G: obtained from get.neighbor function
# n_cluster: number of spatial domains
# w: scaling parameter for image profile
# label_switch_refer, an index of column for label switching, usually choose column index of tumor cell type
run.mIMPACT <- function(V, Y, G, n_cluster, w = 1/20, label_switch_refer = 1){
  n_spot <- nrow(Y)
  n_cell_type <- ncol(V)
  n_PC <- ncol(Y)
  set.seed(123)
  
  e <- rep(1, n_cluster)
  f <- rep(1, n_spot)
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
  result <- mIMPACT(V, Y, G, n_cluster, e, f, eta, tau,  alpha_gamma, beta_gamma, mu_initial, sigma_initial, alpha, omega_initial, Z_initial,rep(w, n_spot))
  
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

# Get spatial domain identification results at single-cell resolution
get.cell.spatial.domain <- function(spatial_domain, cell_assignment){
  n_cell <- length(cell_assignment)
  n_spot <- length(spatial_domain)
  
  spatial_domain_cell <- numeric(n_cell)
  for (i in 1:n_spot){
    spatial_domain_cell[cell_assignment == i] <- spatial_domain[i]
  }
  return(spatial_domain_cell)
}

# Get domain-level cell proportion
get.domain.cell.prop <- function(result){
  omega_p <- result[['omega']]
  domain_cell_proportion <- apply(omega_p, c(2, 3), median)
  colnames(domain_cell_proportion) <- result[['cell_type']]
  
  return(domain_cell_proportion)
}

# Get interactive zone
# prop_cut: threshold to define interactive zone
get.interactive.zone <- function(result, prop_cut = 0.9){
  Z_p <- result[['Z']]
  K <- result[['n_cluster']]
  Z_prop <- matrix(0, ncol = K, nrow = ncol(Z_p))
  for (i in 1:ncol(Z_p)){
    for (k in 1:K){
      Z_prop[i, k] <- mean(Z_p[, i] == k)
    }
  }
  
  interactive_zone <- apply(Z_prop, 1, max) < prop_cut
  
  return(interactive_zone)
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

# Filter out genes with a high proportion of zero counts
# count: a n*p matrix for gene expression counts, where n is the number of spots and p is the number of genes
# min_percentage: gene with non-zero expression proportion <= min_percentage will be filtered out.
#####################
filter.count <- function(count, min_percentage = 0.1){
  gene_num <- ncol(count)
  sample_num <- nrow(count)
  min_total <- 0
  
  if (sum(colSums(count == 0) > (1-min_percentage)*sample_num) > 0){
    count_f <- count[,-(which(colSums(count == 0) > (1-min_percentage)*sample_num))]
  }
  else{
    count_f <- count
  }
  
  return(count_f)
}

# Get the estimated size factor
# count: the n*p matrix for gene expression counts, where n is the number of spots and p is the number of genes
# norm_method: the method to calculate size factor. Choices include 'tss', 'q75', and 'rle'. For details, please refer paper https://onlinelibrary.wiley.com/doi/10.1002/sim.9530
get.size.factor <- function(count, norm_method = 'tss'){
  gene_num <- ncol(count)
  sample_num <- nrow(count)
  
  count_rowsum <- rowSums(count)
  
  if(norm_method == "tss")
  {
    ## TSS(Total Sum Scaling)
    ### scale-factors
    raw_s_factors <- rowSums(count)/mean(rowSums(count))
  }
  else if(norm_method == "q75")
  {
    ## Q75(Upper Quantile normalization)
    ### scale factors
    count_q75 <- apply(count, 1,function(x){quantile(x[x>0],0.75)} )
    count_N <- rowSums(count)/nrow(count)
    raw_s_factors <- count_q75/count_N
  }
  else if(norm_method == "rle")
  {
    ## RLE(Relative Log Expression normalization)
    ### scale_factors
    ### function for calculating the geometric mean
    geo_mean <- function(x){
      exp(sum(log(x[x>0]))/length(x))
    }
    ### function for calculating non-zero median
    non_zero_median <- function(x){
      median(x[x>0])
    }
    ref_sample <- apply(count, 2, geo_mean)
    norm_rle_1 <- sweep(count, 2, ref_sample, FUN = "/")
    raw_s_factors <- apply(as.matrix(norm_rle_1), 1, non_zero_median)
  }
  else{
    stop("Please choose a valid normalization method")
  }
  return(raw_s_factors)
}

# Fit negative binomial regression for domain-specific spatially variable detection
# count: the n*p matrix for gene expression counts, where n is the number of spots and p is the number of genes
# spatial_domain: identified spatial domain from 'get.spatial.domain' or 'refine.cluster' function.
# domain_index: index of the target domain
# si: estimated size factor from 'get.size.factor' function.
detect.domainSVG <- function(count, spatial_domain, domain_index, si){
  n_gene <- ncol(count)
  x <- spatial_domain == domain_index
  if (sum(x) == 0){
    stop('Please enter a valid domain index.')
  }
  results <- data.frame(gene = colnames(count), beta = rep(0, n_gene), p_value = rep(1, n_gene))
  c <- 0
  for (jj in 1:n_gene){
    data_temp <- data.frame(y = count[, jj], si = si, x = x)
    if (sum(data_temp$y == 0) > 2){
      suppressWarnings(m1 <- glm.nb(y ~ x + offset(log(si)), data = data_temp))
      rr <- summary(m1)
      results$p_value[jj] <- rr$coefficients[2, 4]
      results$beta[jj] <- rr$coefficients[2, 1]
      
      if (floor(jj*100/n_gene) == c)
      {
        print(paste0(c, '% has been done'))
        c <- c + 10
      }
    }
  }
  results$adjusted_p_value <- p.adjust(results$p_value, 'BH')
  return(results)
}



# Function to switch colors.
# 'real_region' is the real label for each spot
# 'estimated_region' is the label generated by miIMPACT
switch_color <- function(real_region, estimated_region, K){
  
  estimated_region <- spatial_domain_refined
  color_list <- NULL
  spatial_list <- NULL
  
  color_length <- rep(0, K)
  color_spatial <- rep(0, K)
  
  for (k in 1:K){
    # get the length of each color in real_region
    index1 <- which(real_region == k)
    color_list[[k]] <- index1
    color_length[k] <- length(index1)
    
    # get the length of each color in simulated_region
    index2 <- which(estimated_region == k)
    spatial_list[[k]] <- index2
    color_spatial[k] <- length(index2)
  }
  
  sorted_index_real <- order(color_length)
  sorted_index_spatial <- order(color_spatial)
  
  for (k in 1:K){
    index_to_switch <- spatial_list[[sorted_index_spatial[k]]]
    value_to_assign <- sorted_index_real[k]
    estimated_region[index_to_switch] <- value_to_assign
  }
  
  return (estimated_region)
}





# Helper functions
######################

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