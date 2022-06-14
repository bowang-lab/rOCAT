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

This is the R package for One Cell At A Time(OCAT), which provides a
fast and memory-efficient framework for analyzing and integrating
large-scale scRNA-seq data. Details of the method can be check
[here](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02659-1).

# 1 Installation

You can install the development version of rOCAT from
[GitHub](https://github.com/bowang-lab/rOCAT) with:

``` r
# install.packages("devtools")
devtools::install_github("bowang-lab/rOCAT")
```

# 2 Example

Here we demonstrate how to sparsely encodes single-cell gene expression
data on both single and multiple datasets

The **single** dataset using 3,005 cells and 4,412 genes in the mouse
somatosensory cortex and hippocampal CA1 region from Zeisel et
al. (2015) and can be download
[here](https://github.com/bowang-lab/OCAT/blob/master/vignettes/Clustering/Test_5_Zeisel.mat).

The **multiple** datasets consist of five scRNA-seq datasets (Baron et
al. 2016, Muraro et al. 2016, Segerstolpe et al. 2016, Wang et al. 2016,
Xin et al. 2016) and can be download
[here](https://data.wanglab.ml/OCAT/Pancreas.zip).

## 2.1 Working with Single Dataset

### 2.1.1 load package

``` r
library(rOCAT)
library(reticulate) # must load this to config the python env

# if no python environment exists in the machine
reticulate::py_config()
# reticulate::install_miniconda() # should prompt to ask for installation

#rOCAT::install_packages() # you could ignore the line if you have 'OCAT', 'scikit-learn', and 'faiss' packages installed in your python env
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
```

### 2.1.4 clustering and evaluation

``` r
labels_pred <- evaluate_clusters(np_array(ZW),num_cluster = length(unique(labels_true)))

# get the NMI score between our predicted cell type and the true cell type 
normalized_mutual_info_score(labels_true, labels_pred)
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

# load the labels from each batch and combine them into one file
labels <- lapply(files,  FUN = function(x) np$load(paste0(label_path,x,"_label.npy"),allow_pickle = TRUE))
labels_true_combined <- do.call(c, labels)

# also save the batch number together
batch_true_combined <- rep(files,lengths(labels)) #batch labels 
```

### 2.2.2 run OCAT pipline

``` r
m_list <- list(65, 65, 65, 65, 65)

# execute the pipeline in these datasets
pancreatic_ZW <- run_OCAT(data_list, m_list, dim=60, p=0.3, log_norm=TRUE, l2_norm=TRUE)
```

### 2.2.3 clustering and evaluation

#### 2.2.3.1 cell type clustering

``` r
cell_type_pred <- evaluate_clusters(pancreatic_ZW, num_cluster = length(unique(labels_true_combined)))

# get the NMI score between our predicted cell type and the true cell type
normalized_mutual_info_score(labels_true_combined,cell_type_pred)
```

#### 2.2.3.2 batch clustering

``` r
batch_pred <- evaluate_clusters(pancreatic_ZW, num_cluster = length(files))

# get the NMI score between our predicted batch and the true batch
normalized_mutual_info_score(batch_true_combined,batch_pred)
```
