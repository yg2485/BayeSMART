# iBURST: Integrating Bulk RNA-Seq Data and Spatial Transcriptomics

The model consists of two part: 1) Spatial domain identification of spatial transcriptomics(ST) on multi-sample and 2) bulk RNA-seq deconvolution

![iBURST](figure/flowchart.png)

## Spatial Domain Identification for Multi-sample
The code for the first part is in folder code/mIMPACT. Clone the whold repo and run codes in

``` 
main.R
```

For getting the results of six breast cancer patients, source the functions first by
```r
source("../../code/mIMPACT/iIMPACT_for_multisample.R")
```

and get the data ready
```r
my_dict <- list("A1", "B1", "C1", "D1", "G1", "H1")

n_cell_type <- 6 # number of cell types in the image profile
V <- matrix(ncol = n_cell_type, nrow = 0)
Y <- matrix(ncol = 3, nrow = 0)
G <- matrix(ncol = 8, nrow = 0)
G_orgin <- matrix(ncol = 8, nrow = 0)
spot_loc_all <- matrix(ncol = 2, nrow = 0)

each_length <-c()


for (s in 1:6){
  filename = paste0("../../data/multisample data/processed/", my_dict[[s]], ".RData")
  load(filename)
  
  cell_xy <- cell_info[, 1:2]
  cell_name <- cell_info$name
  count <- spot_count[, -1]
  # count <- spot_count[, common_chars]
  
  each_length <- c(each_length, nrow(spot_xy))
  
  V1 <- get.cell.abundance(cell_xy, cell_name, spot_xy, lattice = 'square')
  Y1 <- process.gene.expression(count, n_PC = 3, method = 'PCA',d = 734)
  G1 <- get.neighbor(spot_xy, 8)
  
  G_orgin <- rbind(G_orgin, G1)
  if (s > 1){
    non_zero_indices <- G1 != 0
    G1[non_zero_indices] <- G1[non_zero_indices] + sum(each_length[1:(s-1)])
  }
  V <- rbind(V, V1)
  G <- rbind(G, G1)
  Y <- rbind(Y, Y1)
  
  spot_loc_all = rbind(spot_loc_all, spot_xy)
}
```

Then run the model
```r
# number of spatial domains
K <- 6

# weight of image profile
w <- 1/10

# run the miIMPACT model for spatial domain identification on multi-sample
result <- run.mIMPACT(V, Y, G, n_cluster = K, w)

# get the posterior inference on spatial domains
spatial_domain <- get.spatial.domain(result)
```


Calculate the ARI for each sample
```r
scores = c() # this is for saving the ARI for each sample


for (s in 1:6){
  # get real domain label
  filename_2 = paste0("../../data/multisample data/processed/", my_dict[[s]],"_color.RData") 
  load(filename_2)
  l = each_length[s]
  
  if (s == 1){
    spots <- spot_loc_all[1:l, ]
    domains <- spatial_domain[1:l]
    G1 <- G_orgin[1:l, ]
  }else{
    sum_ = sum(each_length[1: (s-1)])
    spots <- spot_loc_all[(sum_ + 1): (sum_ + l), ]
    domains <- spatial_domain[(sum_ + 1): (sum_ + l)]
    G1 <- G_orgin[(sum_ + 1): (sum_ + l), ]
  }
  
  spatial_domain_refined <- refine.cluster(G1, domains, area_unit = 3)
  
  if (length(undefined) == 0){
    cluster_sim = spatial_domain_refined
    cluster_real = real_region
  }else{
    cluster_sim = spatial_domain_refined[-undefined]
    cluster_real = real_region[-undefined]
  }
  
  ari_score <- adjustedRandIndex(cluster_sim, cluster_real)
  scores <- c(scores, ari_score)
}
```


For visualization, we maunally switch the colors assigned to spots and save the results in data/each_sample. Plot the results by
```r
for (s in 1:6){
  filename = paste0("../../data/multisample data/processed/", my_dict[[s]], ".RData")
  load(filename)
  spots <- spot_xy
  
  filename_2 = paste0("../../data/each_sample/", my_dict[[s]], "/", my_dict[[s]], "_estimated_domain.csv")
  spatial_domain_reorder <- read.csv(filename_2)[[1]]
  
  # spatial_domain_reorder <- switch_color(real_region = real_region, estimated_region = spatial_domain_refined, K = 6)
  df <- data.frame(x = spots[,1], y = spots[,2], domain = as.factor(spatial_domain_reorder))
  x_trans <- min(df$x)
  y_trans <- min(df$y)
  df$x <- df$x - x_trans
  df$y <- df$y - y_trans
  
  directory_path <- "../../fig_result/"
  if (!dir.exists(directory_path)) {
    dir.create(directory_path)
  }
  
  filename = paste0(directory_path, my_dict[[s]],  ".png")
  g <- ggplot(df, aes(x = x, y = y, color = domain)) +
    geom_point(size = 3) + scale_color_manual(values=c('1' = "#7DFFFF", '2' = "#FFFF00", '3' = "darkorange", '4' = '#00FF00', '5' = 'red2' , '6' = '#0000FF'))+
    coord_cartesian(xlim = c(0, 9000), ylim = c(0, 8000))+
    theme(panel.background = element_blank())
  ggsave(filename = filename, plot = g, width = 6, height = 4, dpi = 300)
  
}
```

## Bulk RNA-seq Deconvolution
The code for the second part is in folder code/bulk decov.

```
main_2.R
```
