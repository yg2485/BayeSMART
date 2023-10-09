## Function to get the normalization factors
normalize.joint <- function(st_count, bulk_count, norm_method = 'tss'){
  n <- nrow(bulk_count)
  m <- nrow(st_count)
  p <- ncol(bulk_count)
  
  bulk_count_rowsum <- rowSums(bulk_count)
  st_count_rowsum <- rowSums(st_count)
  
  if(norm_method == "tss")
  {
    ## TSS(Total Sum Scaling)
    ### scale-factors
    raw_s_factors <- rowSums(bulk_count)
    raw_v_factors <- rowSums(st_count)
    raw_g_factors <- colSums(bulk_count)
    
    scale_coeff <- exp((-1/(n + m)) * (sum(log(raw_s_factors)) + sum(log(raw_v_factors))))
    scale_coeff_2 <- exp((-1/p) * (sum(log(raw_g_factors))))
    scaled_s_factors <- scale_coeff * raw_s_factors
    scaled_v_factors <- scale_coeff * raw_v_factors
    scaled_g_factors <- scale_coeff_2 * raw_g_factors
  }
  else if(norm_method == "q75")
  {
    ## Q75(Upper Quantile normalization)
    ### scale factors (remove those with 0's)
    raw_s_factors <- apply(bulk_count, 1, function(x){quantile(x[x>0],0.75)} )
    raw_v_factors <- apply(st_count, 1, function(x){quantile(x[x>0],0.75)} )
    raw_g_factors <- apply(bulk_count, 2, function(x){quantile(x[x>0],0.75)} )
    
    scale_coeff <- exp((-1/(n + m)) * (sum(log(raw_s_factors)) + sum(log(raw_v_factors))))
    scale_coeff_2 <- exp((-1/p) * (sum(log(raw_g_factors))))
    scaled_s_factors <- scale_coeff * raw_s_factors
    scaled_v_factors <- scale_coeff * raw_v_factors
    scaled_g_factors <- scale_coeff_2 * raw_g_factors
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
    ref_sample <- apply(bulk_count, 2, geo_mean)
    norm_rle_1 <- sweep(bulk_count, 2, ref_sample, FUN = "/")
    raw_s_factors <- apply(as.matrix(norm_rle_1), 1, non_zero_median)
    
    ref_sample <- apply(st_count, 2, geo_mean)
    norm_rle_1 <- sweep(st_count, 2, ref_sample, FUN = "/")
    raw_v_factors <- apply(as.matrix(norm_rle_1), 1, non_zero_median)
    
    scale_coeff <- exp((-1/(n + m)) * (sum(log(raw_s_factors)) + sum(log(raw_v_factors))))
    scaled_s_factors <- scale_coeff * raw_s_factors
    scaled_v_factors <- scale_coeff * raw_v_factors
  }
  else if(norm_method == "tmm")
  {
    ## TMM(Trimmed Mean Method)
    ### scale_factors
    bulk_count_t <- t(bulk_count)
    raw_s_factors <- calcNormFactors(bulk_count_t, method = "TMM")
    
    st_count_t <- t(st_count)
    raw_v_factors <- calcNormFactors(st_count_t, method = "TMM")
    
    scale_coeff <- exp((-1/(n + m)) * (sum(log(raw_s_factors)) + sum(log(raw_v_factors))))
    scaled_s_factors <- scale_coeff * raw_s_factors
    scaled_v_factors <- scale_coeff * raw_v_factors
  }
  else if(norm_method == "none"){
    scaled_s_factors = rep(1, n)
    scaled_v_factors = rep(1, m)
    scaled_g_factors = rep(1, p)
  }
  else
  {
    stop("Please choose a valid normalization method")
  }
  return(list(s = scaled_s_factors, v = scaled_v_factors, g = scaled_g_factors))
}



##### function to get all the count for each spot'
## use_region specify which label to use
## mode specify the level of region
get_count_region <- function(common_chars, use_region = "real", mode = "cluster"){
  library(readr)
  my_dict <- list("A1", "B1", "C1", "D1", "G1", "H1")
  count_all = matrix(NA, nrow = 0, ncol = length(common_chars))
  region_all <- c()
  
  if (use_region == "real"){
    last_part = "_real_domain.csv"
  } else if (use_region == "estimated"){
    last_part = "_estimated_domain.csv"
  }
  
  if (mode == "spot"){
    for (s in 1:6){
      filename = paste0("../../data/multisample data/processed/", my_dict[[s]], ".RData")
      load(filename)
      
      count <- spot_count[, common_chars]
      count_all <- rbind(count_all, count)
      
      folder_name <- paste0("../../data/each_sample/", my_dict[[s]])
      filename_2 <- paste0(folder_name, "/", my_dict[[s]], last_part)
      
      region <- c(as.matrix((read_csv(filename_2))))
      region_all <- c(region_all, region)
    }
    
    index <- which(region_all >=1 & region_all <=6)
    count_all <- count_all[index,]
    region_all <- region_all[index]
  } else if (mode == "cluster"){
    row_names <- c()
    for (s in 1:6){
      filename = paste0("../../data/multisample data/processed/", my_dict[[s]], ".RData")
      load(filename)
      
      count <- spot_count[, common_chars]
      folder_name <- paste0("../../data/each_sample/", my_dict[[s]])
      filename_2 <- paste0(folder_name, "/", my_dict[[s]], last_part)
      
      region <- c(as.matrix((read_csv(filename_2))))
      
      for (c in 1:6){
        ind = which(region == c)
        if (length(ind) > 0){
          count_c <- count[ind,]
          count_all <- rbind(count_all, colSums(count_c))
          region_all <- c(region_all, c)
          color_name <- c("lb","ye", "or","gr","red", "bl")
          name <- paste0(my_dict[[s]],"_",color_name[c])
          row_names <- c(row_names, name)
          
        }
      }
    }
    rownames(count_all) <- row_names
  }
  
  return(list(count_all = count_all, region_all = region_all))
}


get_syn_bulk <- function(data, label, n_sample = 10, n_spot_min = 150, n_spot_max = 1000, n_region = 6){
  
  bulk_syn <- data.frame(matrix(ncol = 0, nrow = ncol(data)))
  real_prop <- matrix(nrow = n_sample, ncol = n_region)
  rownames(bulk_syn) <- colnames(data)
  n_spot_max = min(n_spot_max, nrow(data))
  
  if (n_spot_min >= n_spot_max){
    stop('n_spot_max must be no less than n_spot_min!')
  }else{
    for (i in 1:n_sample){
      n_spot <- sample(n_spot_min:n_spot_max, 1)
      spot_ind <- sample(1:nrow(data), n_spot, replace = TRUE)
      temp <- data[spot_ind,]
      labels <- label[spot_ind]
      
      all_spot <- colSums(temp)
      name_sample <- paste0("spot", i)
      bulk_syn[[name_sample]] <- c(all_spot)
      
      for (c in 1:6){
        real_prop[i,c] <- sum(labels == c)/n_spot
      }
    }
  }
  return(list(bulk_syn = bulk_syn, real_prop = real_prop))
}
