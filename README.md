rOCAT
================

-   [1 Installation](#1-installation)
-   [2 Example](#2-example)
    -   [2.1 Working with Single
        Dataset](#21-working-with-single-dataset)
        -   [2.1.1 load package](#211-load-package)
        -   [2.1.2 import data](#212-import-data)
        -   [2.1.3 run OCAT pipeline](#213-run-ocat-pipeline)
        -   [2.1.4 clustering and
            evaluation](#214-clustering-and-evaluation)
    -   [2.2 Working with Multiple
        Datasets](#22-working-with-multiple-datasets)
        -   [2.2.1 import datasets](#221-import-datasets)
        -   [2.2.2 run OCAT pipline](#222-run-ocat-pipline)
        -   [2.2.3 clustering and
            evaluation](#223-clustering-and-evaluation)

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/bowang-lab/rOCAT.svg?branch=main)](https://travis-ci.com/bowang-lab/rOCAT)

<!-- badges: end -->

The goal of rOCAT is to …

# 1 Installation

You can install the development version of rOCAT from
[GitHub](https://github.com/bowang-lab/rOCAT) with:

``` r
# install.packages("devtools")
devtools::install_github("bowang-lab/rOCAT")
```

# 2 Example

This is a basic example which shows you how to solve a common problem:

The dataset Test_5\_Zeisel.mat can be download [here]().

## 2.1 Working with Single Dataset

### 2.1.1 load package

``` r
library(rOCAT)
library(reticulate)
```

### 2.1.2 import data

``` r
# load Test Zeisel data
data <- R.matlab::readMat('../data/Test_5_Zeisel.mat')
# extract the gene feature matrix and make it saved as a sparse matrix
in_X <- as(data$in.X, 'dgCMatrix')
# extract the labels
labels_true <- as.vector(data$true.labs)

# the input data should be a vector of datasets c(datasets1,dataset2,...)
data_list <- c(in_X)
```

### 2.1.3 run OCAT pipeline

``` r
ZW <- run_OCAT(data_list, m_list=list(50), dim=30, 
                p=0.3, log_norm=TRUE, l2_norm=TRUE, if_inference=FALSE, 
                random_seed=42)
#> [1] "Starting Dimension Reduction"
#> [1] "Dimension Reduction Finished"

# first 5x5 dimension of ZW
ZW[1:5,1:5]
#>             [,1]        [,2]        [,3]        [,4]       [,5]
#> [1,] 0.001819295 0.004266508 0.003887158 0.001397503 0.04628960
#> [2,] 0.002862158 0.004538172 0.003428339 0.001378417 0.06061592
#> [3,] 0.002298867 0.004474250 0.004187966 0.001341240 0.05554514
#> [4,] 0.002074808 0.004342816 0.004458587 0.001553739 0.05094255
#> [5,] 0.001530818 0.004259680 0.004892700 0.001822031 0.04234504
```

### 2.1.4 clustering and evaluation

``` r
labels_pred <- evaluate_clusters(np_array(ZW),num_cluster = length(unique(labels_true)))
normalized_mutual_info_score(labels_true, labels_pred)
#> [1] 0.7884375
```

## 2.2 Working with Multiple Datasets

### 2.2.1 import datasets

``` r
# the data is saved in .npz file in this example, so we need scipy from python to help read the data
scipy <- reticulate::import("scipy.sparse")
data_path <- "../data/Pancreas/data/"
files <- c('baron_1','muraro_2', 'seg_3', 'wang_4', 'xin_5')
data_list <- lapply(files,  FUN = function(x) scipy$load_npz(paste0(data_path,x,".npz")))

# import the labels
np <- reticulate::import("numpy") 
label_path <- "../data/Pancreas/label/"
labels <- lapply(files,  FUN = function(x) np$load(paste0(label_path,x,"_label.npy"),allow_pickle = TRUE))
labels_true_combined <- do.call(c, labels)
batch_true_combined <- rep(files,lengths(labels)) #batch labels 
```

### 2.2.2 run OCAT pipline

``` r
m_list <- list(65, 65, 65, 65, 65)
pancreatic_ZW <- run_OCAT(data_list, m_list, dim=60, p=0.3, log_norm=TRUE, l2_norm=TRUE)
#> [1] "Starting Dimension Reduction"
#> [1] "Dimension Reduction Finished"
pancreatic_ZW[0:5,0:5]
#>             [,1]       [,2]        [,3]        [,4]        [,5]
#> [1,] 0.006040815 0.01143108 0.010623352 0.005410251 0.009862246
#> [2,] 0.005507492 0.01095782 0.006144693 0.007255905 0.008626724
#> [3,] 0.006875831 0.01218184 0.005919864 0.014213251 0.008533677
#> [4,] 0.004385151 0.00766135 0.004018047 0.003530303 0.005135434
#> [5,] 0.004746339 0.01214920 0.005628659 0.004717387 0.008688286
```

### 2.2.3 clustering and evaluation

#### 2.2.3.1 cell type clustering

``` r
cell_type_pred <- evaluate_clusters(pancreatic_ZW, num_cluster = length(unique(labels_true_combined)))
normalized_mutual_info_score(labels_true_combined,cell_type_pred)
#> [1] 0.7595895
```

#### 2.2.3.2 batch clustering

``` r
batch_pred <- evaluate_clusters(pancreatic_ZW, num_cluster = length(files))
normalized_mutual_info_score(batch_true_combined,batch_pred)
#> [1] 0.03593664
```
