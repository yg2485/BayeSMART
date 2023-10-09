library(readxl)
library(readr)
library(compositions)
library(ggplot2)


Rcpp::sourceCpp("bulk_decov.cpp")
source("helper_function.R")


my_dict <- list("A1", "B1", "C1", "D1", "G1", "H1")
gene_names <- NULL

for (s in 1:6){
  filename = paste0("../../data/multisample data/processed/", my_dict[[s]], ".RData")
  load(filename)
  count <- spot_count[, -1]
  gene_names[[s]] <- colnames(count)
  
}

common_chars <- Reduce(intersect, gene_names)


### Import Data

# st data
count_region <- get_count_region(common_chars, use_region = "estimated", mode = "cluster")
count_all <- count_region$count_all
region_all <- count_region$region_all

# bulk data
bulk_data <- read_excel("../../data/bulk decov/breast cancer bulk sequencing and clinical data/BRCA.rnaseq__illuminahiseq-Raw-Counts.xlsx")
bulk_data <- bulk_data[-1, ]

bulk_data <- as.data.frame(bulk_data)

common_gene <- intersect(bulk_data$`Sample ID`, common_chars)

bulk_data <- bulk_data[bulk_data$`Sample ID` %in% common_gene, ]
rownames(bulk_data) <- bulk_data[[1]]
bulk_data <- bulk_data[,-1]
sample_id <- colnames(bulk_data)

bulk_data <- t(as.matrix(bulk_data))
bulk_data <- matrix(as.numeric(as.vector(bulk_data)), nrow=nrow(bulk_data))
colnames(bulk_data) <- common_gene
rownames(bulk_data) <- sample_id



marker_gene <- read_tsv("../../data/bulk decov/DEG/domain_markers_estimate.tsv")

# ### get the union of DEG over all samples ###
# marker_list <- marker_gene[marker_gene$p_val_adj < 0.05,]$gene
# genes_union <- intersect(common_gene, marker_list)


### get the intersect of DEG over all samples ###
my_dict <- list("A1", "B1", "C1", "D1", "G2", "H1")
gene_list <- list()
for (i in 1:6){
  gene_list[[i]] <- marker_gene$gene[marker_gene$Sample == my_dict[[i]] & marker_gene$p_val < 0.05]
}
genes <- Reduce(intersect, gene_list)
###


## ================== real data ==========================

common_left <- c()
for (j in 1:length(common_gene)){
  if (common_gene[j] %in% genes){
    common_left <- c(common_left, common_gene[j])
  }
}
gamma <- rep(1, length(common_left))


ref_count <- count_all[,common_left]
bulk_count <- bulk_data[,common_left]
cell_type <- region_all

## ================== simulate data ==========================

common_left <- c()
for (j in 1:length(common_chars)){
  if (common_chars[j] %in% genes){
    common_left <- c(common_left, common_chars[j])
  }
}
gamma <- rep(1, length(common_left))

ref_count <- count_all[,common_left]
bulk_count <- t(as.matrix(bulk_syn))[,common_left]
cell_type <- region_all



# ===============================================================

size.factor = normalize.joint(ref_count = ref_count, bulk_count = bulk_count)
s = size.factor$s
v = size.factor$v
g = size.factor$g


iter = 20000
burn = 10000


# simulated
res = imoscato_new(X=as.matrix(ref_count), 
               Y=as.matrix(bulk_count), 
               cell_type=cell_type, 
               gamma = gamma,
               v=v, 
               s=s, 
               g=g, 
               iter=iter, 
               burn=burn,
               tau_pi = 0.02,
               tau_mu = 0.5,
               d = 42)

# real
res = imoscato_new(X=as.matrix(ref_count), 
                   Y=as.matrix(bulk_count), 
                   cell_type=cell_type, 
                   gamma = gamma,
                   v=v, 
                   s=s, 
                   g=g, 
                   iter=iter, 
                   burn=burn,
                   tau_pi = 0.001,
                   tau_mu = 0.01,
                   d = 42)


Pi_store = lapply(1:iter, function(x) res$Pi_store[ , , x])
proportion = Reduce("+", Pi_store[(burn + 1):iter]) / burn

mu_post <- lapply(1:iter, function(x) res$M_store[ , , x])
mu_post <- Reduce("+", mu_post[(burn + 1):iter]) / burn

# write.csv(proportion, "/Users/yanghongguo/Desktop/Research/iBurst/data/bulk decov/imoscato results/imos_estimated_region.csv", row.names=TRUE)
# write.csv(proportion, "/Users/yanghongguo/Desktop/Research/iBurst/data/bulk decov/imoscato results/imos_real_region.csv", row.names=TRUE)


plot(res$Pi_store[2,5,1:iter])

# real_prop <- matrix(nrow = ncol(bulk_syn), ncol = 6)
# for (s in 1:6){
#   folder_name <- paste0("/Users/yanghongguo/Desktop/each_sample/", my_dict[[s]])
#   # get the truth label (estimated or real)
#   filename_2 <- paste0(folder_name, "/", my_dict[[s]], "_estimated_domain.csv") # defines which as the true label
# 
#   region <- c(as.matrix((read_csv(filename_2))))
#   region <- region[which(region >=1 & region <=6)]
#   for (c in 1:6){
#     real_prop[s,c] <- sum(region == c)/length(region)
#   }
# }

pearson <- c()
for (i in 1:nrow(proportion)){
  pearson[i] <- cor(proportion[i,], real_prop[i,])
}

pearson



#------------------PCA on real bulk with only DEG------------------#
library(readxl)
file_path <- "../../data/bulk decov/breast cancer bulk sequencing and clinical data/"
bulk_clinic <- read_excel(paste0(file_path, "Matched-Clinical-data-with-RNAseq.xlsx"))
bulk_clinic <- as.data.frame(bulk_clinic)


## CLR to map the data from simplex to euclidian space
prop_data <- read.csv("../../data/bulk decov/RCTD results/RCTD_estimated_region.csv",
                      row.names = 1)
library(compositions)
clr_data <- clr(prop_data)
clr_data <- as.matrix(clr_data)


## data cleaning for bulk_clinic
name_list <- c('vital_status', 
               'breast_carcinoma_estrogen_receptor_status',
               'breast_carcinoma_progesterone_receptor_status',
               'ethnicity',
               'gender',
               'histological_type',
               'margin_status',
               'menopause_status',
               'pathologic_stage',
               'person_neoplasm_cancer_status',
               'race')
name_list_2 <- c('vital status', 
                 'breast carcinoma estrogen receptor status',
                 'breast carcinoma progesterone receptor status',
                 'ethnicity',
                 'gender',
                 'histological type',
                 'margin status',
                 'menopause status',
                 'pathologic stage',
                 'person neoplasm cancer status',
                 'race')

bulk_clinic <- bulk_clinic[bulk_clinic$`Sample ID` %in% name_list,]
rownames(bulk_clinic) <- bulk_clinic$`Sample ID`
bulk_clinic <- bulk_clinic[,-1]
bulk_clinic <- t(bulk_clinic)
bulk_clinic <- as.data.frame(bulk_clinic)
bulk_clinic[bulk_clinic == "NA"] <- NA
bulk_clinic$pathologic_stage <- ifelse(bulk_clinic$pathologic_stage %in% c("stage i", "stage ia","stage ib"), "i", 
                                       ifelse(bulk_clinic$pathologic_stage %in% c("stage ii", "stage iia","stage iib"), "ii",
                                              ifelse(bulk_clinic$pathologic_stage %in% c("stage iii", "stage iiia","stage iiib", "stage iiic"), "iii", 
                                                     ifelse(bulk_clinic$pathologic_stage == "iv", "iv", NA))))
bulk_clinic$menopause_status[bulk_clinic$menopause_status == "indeterminate (neither pre or postmenopausal)"] <- "indeterminate"
bulk_clinic$menopause_status[bulk_clinic$menopause_status == "peri (6-12 months since last menstrual period)"] <- "peri"
bulk_clinic$menopause_status[bulk_clinic$menopause_status == "post (prior bilateral ovariectomy or >12 mo since lmp with no prior hysterectomy)"] <- "post"
bulk_clinic$menopause_status[bulk_clinic$menopause_status == "pre (<6 months since lmp and no prior bilateral ovariectomy and not on estrogen replacement)"] <- "pre"

bulk_clinic$histological_type[bulk_clinic$histological_type == "mixed histology (please specify)"] <- "mixed"
bulk_clinic$histological_type[bulk_clinic$histological_type == "other, specify"] <- "other"
bulk_clinic$histological_type[bulk_clinic$histological_type == "infiltrating ductal carcinoma"] <- "infiltrating ductal"
bulk_clinic$histological_type[bulk_clinic$histological_type == "infiltrating lobular carcinoma"] <- "infiltrating lobular"
bulk_clinic$histological_type[bulk_clinic$histological_type == "medullary carcinoma"] <- "medullary"
bulk_clinic$histological_type[bulk_clinic$histological_type == "metaplastic carcinoma"] <- "metaplastic"
bulk_clinic$histological_type[bulk_clinic$histological_type == "mucinous carcinoma"] <- "mucinous"



## CLR to map the data from simplex to euclidian space
proportion <- read.csv("../../data/bulk decov/imoscato results/imos_estimated_region.csv", row.names = 1)

rownames(proportion) <- sample_id

clr_data <- clr(proportion)
clr_data <- as.matrix(clr_data)

## PCA and visualize with ggplot2

# Running PCA
pca_result <- prcomp(clr_data, scale. = TRUE)
# Extracting principal components
pca_data <- as.data.frame(pca_result$x)

# choose clinical variable to visualize
i <- 2
df_A <- bulk_clinic[,name_list[i], drop = FALSE]
df_no_na <- na.omit(df_A)
df_no_indeter <- df_no_na[df_no_na[,1] != "indeterminate", , drop = FALSE]
common_id <- intersect(rownames(df_no_indeter), rownames(pca_data))


# get the final pca and group information
pca_data <- pca_data[common_id,]
pca_data$Group <- df_no_na[common_id,1]


# Plotting the first two principal components
ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 1) +
  ggtitle(name_list_2[i]) +
  xlab("First Principal Component") +
  ylab("Second Principal Component") +
  labs(color = "Receptor Status") +
  theme(panel.background = element_blank(),
        axis.line = element_line(color = "gray"),
        axis.title.x = element_text(size = 12),  # Adjust x-axis label size
        axis.title.y = element_text(size = 12),  # Adjust y-axis label size
        plot.title = element_text(size = 14),
        legend.title = element_text(size = 12),   # Adjust legend title size
        legend.text = element_text(size = 12))



"breast carcinoma estrogen receptor status"
"breast carcinoma progesterone receptor status"







# barplot of first 8 individuals
par(mfrow = c(2, 4), mar = c(1,1,1,1), mai = c(0.3, 0.28, 0.3, 0.1))
groups <- c('adipose tissue', 'immune infiltrate', 'in situ cancer', 'breast glands', 'invasive cancer', 'connective tissue')
colors <- c("cadetblue1", "moccasin", "orange2", "seagreen2", "salmon", "dodgerblue")
for (i in 1:8){
  barplot(unlist(proportion[i,]), 
       col = colors,            # Assign colors to bars
       main = sample_id[i],     # Title of the plot
       xaxt = "n",             # Remove x-axis tick labels
  )
}

# Combine the color legends into a single legend

par(mfrow = c(1, 1))
# par(mar = c(10, 4, 4, 2))
legend("right", legend = groups, fill = colors, title = "Domains", cex = 1.5)



# barplot of 4 positive and 4 negative patients
i <- 2
df_A <- bulk_clinic[,name_list[i], drop = FALSE]
df_no_na <- na.omit(df_A)
df_no_indeter <- df_no_na[df_no_na[,1] != "indeterminate", , drop = FALSE]
common_id <- intersect(rownames(df_no_indeter), rownames(pca_data))
df_left <- df_no_indeter[common_id,, drop = FALSE]

get_positive <- rownames(df_left)[df_left$breast_carcinoma_estrogen_receptor_status == "positive"]
get_negative <- rownames(df_left)[df_left$breast_carcinoma_estrogen_receptor_status == "negative"]

newprop <- matrix(NA, nrow = 8, ncol = 6)

set.seed(731)
newprop[1:4,] <- as.matrix(proportion[sample(get_positive, 4, replace = FALSE),])
newprop[5:8,] <- as.matrix(proportion[sample(get_negative, 4, replace = FALSE),])


par(mfrow = c(2, 4), mar = c(1,1,1,1), mai = c(0.3, 0.28, 0.3, 0.1))
groups <- c('adipose tissue', 'immune infiltrate', 'in situ cancer', 'breast glands', 'invasive cancer', 'connective tissue')
colors <- c("cadetblue1", "moccasin", "orange2", "seagreen2", "salmon", "dodgerblue")
for (i in 1:8){
  barplot(newprop[i,], 
          col = colors,            # Assign colors to bars
          main = sample_id[i],     # Title of the plot
          xaxt = "n",             # Remove x-axis tick labels
  )
}
