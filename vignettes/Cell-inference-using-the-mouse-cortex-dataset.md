In this page, we demonstrated how to run cell type inference with our
rOCAT package

The demonstration dataset using 3,005 cells and 4,412 genes in the mouse
somatosensory cortex and hippocampal CA1 region from Zeisel et
al. (2015) and can be download
[here](https://github.com/bowang-lab/OCAT/blob/master/vignettes/Clustering/Test_5_Zeisel.mat).

# Import Data

    library(rOCAT)
    library(reticulate)

    # load Test Zeisel data
    data <- R.matlab::readMat('../../data/Test_5_Zeisel.mat')

    # extract the gene feature matrix and make it saved as a sparse matrix
    in_X <- as(data$in.X, 'dgCMatrix')

    # extract the labels
    labels_true <- as.vector(data$true.labs)

# split data into reference and inference part

    set.seed(123)
    reference_index <- sort(sample(nrow(in_X), nrow(in_X)*0.9))
         
    in_X_reference <- in_X[reference_index,]
    in_X_inference <- in_X[-reference_index,]

    labels_true_reference <- labels_true[reference_index]
    labels_true_inference <- labels_true[-reference_index]

    # the input data should be a vector of datasets list(datasets1,dataset2,...)
    data_list <- list(in_X_reference)
    data_list_inf  <- list(in_X_inference)

# Run OCAT on the reference dataset

    ZW_db_db_list <-  run_OCAT(data_list, m_list=c(50), dim=30, 
                    p=0.3, log_norm=TRUE, l2_norm=TRUE, if_inference=TRUE, 
                    random_seed=42, labels_true = list(labels_true_reference))
    #> [1] "Starting Dimension Reduction"
    #> [1] "Dimension Reduction Finished"
    #> New m_list based on the true cluster: 5 6 16 13 2 3 3 1 1[1] ""
    #> [1] "Start Sparse Encoding"
    #> [1] "Finished Sparse Encoding ...."
    ZW_db <- ZW_db_db_list[[1]]
    db_list <- ZW_db_db_list[[2]]

# Run OCAT on the inference dataset

    ZW_inf_labels <- run_cell_inference(data_list_inf,labels_db=labels_true_reference,db_list = db_list, true_known = TRUE)
    ZW_inf <- ZW_inf_labels[[1]]
    inf_labels <- ZW_inf_labels[[2]]
