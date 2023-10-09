library(readr)
library(pROC)
library(caret)

## Binary factorization function:
# If method = 1, set 3 and 5 to 1, the others to 0;
# if method = 2, set 5 to 1, the others to 0
# 3 is the in situ cancer, 5 is invasive cancer
binary_label <- function(arr, method = 1, red_val = 5, orange_val = 3){
  arr <- as.matrix(arr)
  # arr <- real_domain
  if (method == 1){
    ind <- which(arr %in% c(red_val, orange_val))
    arr[ind] <- 1
    arr[-ind] <- 0
  }else if (method == 2){
    ind <- which(arr == red_val)
    arr[ind] <- 1
    arr[-ind] <- 0
  }else{
    stop('Parameter lattice: Please enter a correct value.')
  }
  return(arr)
}


## Get AUC, F1 score, Accuracy of the result
# If method = 1, set 3 and 5 to 1, the others to 0;
# if method = 2, set 5 to 1, the others to 0
get_new_metric <- function(method = 1, red_val = 5, orange_val = 3){
  
  my_dict <- list("A1", "B1", "C1", "D1", "G1", "H1")
  metric_score <- matrix(NA, nrow = 6, ncol = 3)
  colnames(metric_score) <- c('AUC', 'F1', 'Accuracy')
  s <- 1
  for (s in 1:6){
    suppressMessages({
      folder_name <- paste0("../../data/each_sample/", my_dict[[s]])
      real_domain <- as.numeric(as.matrix(read_csv(paste0(folder_name, "/", my_dict[[s]], "_real_domain.csv"))))
      estimated_domain <- as.numeric(as.matrix(read_csv(paste0(folder_name, "/", my_dict[[s]], "_estimated_domain.csv"))))
    })
    
    real_binary <- binary_label(real_domain, method = method)
    estimated_binary <- binary_label(estimated_domain, method = method, red_val = red_val, orange_val = orange_val)
    
    del_ind <- which(real_domain == 7)
    if (length(del_ind) != 0){
      real_binary <- real_binary[-del_ind]
      estimated_binary <- estimated_binary[-del_ind]
    }
    
    # AUC
    roc_obj <- roc(response = real_binary, predictor = estimated_binary)
    auc_value <- auc(roc_obj)
    
    # F1 score
    cm <- confusionMatrix(as.factor(estimated_binary), as.factor(real_binary), mode = 'everything', positive = "1")
    precision <- cm$byClass['Pos Pred Value']
    recall <- cm$byClass['Sensitivity']
    f1_score <- 2 * (precision * recall) / (precision + recall)
    
    # Accuracy
    accuracy <- cm$overall['Accuracy']
    
    metric_score[s, ] <- c(auc_value, f1_score, accuracy)
    
  }
  return (metric_score)
}


#===================================================================#
#===================================================================#
#===================================================================#


######## only calculate cancer (region with label 3 and 5) ########
metric <- get_new_metric(method = 1)



######## only calculate invasive cancer (region with 5) ########
metric <- get_new_metric(method = 2)


