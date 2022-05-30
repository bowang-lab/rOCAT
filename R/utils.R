#' @description  install required python packages
#' @export
load_all <- function(){
  reticulate::py_install('git+https://github.com/bowang-lab/OCAT.git',pip=TRUE)
  reticulate::py_install('faiss')
  reticulate::py_install("scikit-learn")
  }


#' @description  normalize all cells by log10(x+1)
#' @export
normalize_data <- function(data_list){
  data_list <- sapply(data_list, function(x) log10(x+1))
  data_list <- sapply(data_list, function(x) Matrix::Matrix(x,sparse = TRUE))
  return(data_list)
}

#' @description  normalize all cells by L2 normalization
#' @export
l2_normalization <- function(data_list){
  for (i in seq_along(data_list)){
    l2_norm <- sqrt(Matrix::colSums(data_list[[i]]^2))
    l2_norm[l2_norm==0] <- 0.00001
    data_list[[i]] <- sweep(data_list[[i]],2,l2_norm,'/')
    data_list[[i]] <- Matrix::Matrix(data_list[[i]], sparse = TRUE)
  }
  return(data_list)
}

m_estimate <- function(data_list){
  m_list <- c()
  for (i in data_list){
    len <- i@Dim[1]
    if (len<2000)
      m_list <- append(m_list,20)
    else if (len<10000)
      m_list <- append(m_list,40)
    else
      m_list <- append(m_list,60)
  }
  print(paste0('estimated m_list: ', m_list))
  return(m_list)
}

#' @description  normalize cells expression
#' @param data_list   [(a,n)...(z,n)]    -- list of datasets.
#' @param log_norm  flag for log10 normalization.
#' @param l2_norm   flag  for L2 normalization.
#' Out:
#' @return data_list   [(a,dim)...(z,dim)]    -- list of normalized datasets
#'
#' @export
preprocess <- function(data_list, log_norm=log_norm, l2_norm=l2_norm){
  stopifnot("Data list cannot be empty"= length(data_list)>0)
  
  if (log_norm){
    data_list <- normalize_data(data_list)
  }
  if (l2_norm){
    data_list <- l2_normalization(data_list)
  }
  return(data_list)
}


dim_estimate <- function(data_list){
  gene_num <- data_list[[1]]@Dim[1]
  if (gene_num<5000)
    return(50)
  else if (gene_num<10000)
    return(100)
  else
    return(125)
}


#' @description  The function reduces the datasets to dim subspaces
#' @param data_list   [(a,n)...(z,n)]  list of datasets
#' @param dim   desired dimension after dimensionality reduction
#' @param mode  method for finding anchors, default on 'FSM'
#' @param upsample flag for upsample
#' @return data_list   [(a,dim)...(z,dim)] list of datasets
#' @export
apply_dim_reduct <- function(data_list, dim=NULL, mode='FSM', random_seed=42, upsample=FALSE){
  if (is.null(dim))
    dim <- dim_estimate(data_list)

  OCAT <- reticulate::import('OCAT')
  
  print("Starting Dimension Reduction")
  scipy <- reticulate::import("scipy.sparse",convert = FALSE)
  preprocessed_converted_data <- lapply(data_list, function(x) scipy$csr_matrix(reticulate::np_array(x,dtype='float32')))
  reticulate::py_run_string(glue::glue("dim = int({dim});random_seed = {random_seed};upsample={reticulate::r_to_py(upsample)};mode='{mode}'"))
  dim_reducted <- OCAT$apply_dim_reduct(preprocessed_converted_data,dim=py$dim,mode = py$mode, random_seed = py$random_seed,upsample = py$upsample)
  print("Dimension Reduction Finished")
  return(dim_reducted)
}

#' @description  The function generates the corresponding anchors for the datasets
#' @param data_list   [(a,n)...(z,n)] list of datasets
#' @param m_list   lift of num of anchors
#' Out:
#' @return anchor_list [(m,n)...(m,n)] list of anchors
#' @export
find_anchors <- function(data_list,m_list){
  faiss <- reticulate::import('faiss')
  
  anchor_list <- vector(mode = "list", length = length(data_list))
  for (X in seq_along(data_list)){
    dim1 <- dim(data_list[[X]])[[2]]
    temp <- reticulate::np_array(as(data_list[[X]],"matrix"),dtype = 'float32')
    m <- m_list[[X]]
    reticulate::py_run_string(glue::glue("dim={dim1};k={m}",convert=FALSE))
    kmean <- faiss$Kmeans(py$dim,py$k)
    kmean$train(temp)
    anchors <- kmean$centroids
    anchor_list[[X]] <- as(anchors,"matrix")
  }
  return(anchor_list)
}

#'@description  The function generates the sparsed encoding of edges
#' @param data_list   [(a,dim)...(z,dim)]  list of datasets (dim PCs)
#' @param m_list      num of anchors
#' @param s_list      num of anchors to be selected
#' @param p           percentage of NNs to consider
#' @param cn          rounds of optimization
#' @param if_inference     flag for cell inference
#' Out:
#' @return ZW        (a+...+z, m)  OCAT feature matrix
#'
#' @export
sparse_encoding_integration <- function(data_list, m_list, s_list, p=0.3, cn=5, if_inference=TRUE){
  OCAT <- reticulate::import('OCAT')
  anchor_list <- find_anchors(data_list, m_list)
  
  Z_list <- NULL
  for (i in seq_along(data_list)){
    dataset_Z_list <- NULL
    
    dataset <- data_list[[i]]
    dataset <- reticulate::np_array(dataset,dtype='float32')
    
    for (j in seq_along(anchor_list)){
      anchor <- reticulate::np_array(anchor_list[[j]],dtype='float32')
      reticulate::py_run_string(glue::glue("s = {s_list[[j]]};flag=2;cn = {cn}"))
      Z <- OCAT$sparse_encoding$AnchorGraph(dataset$transpose(),anchor$transpose(),py$s,py$flag,py$cn)
      Z <- as(Z,'matrix')
      dataset_Z_list <- cbind(dataset_Z_list,Z)
      dataset_Z_list <- as(dataset_Z_list,'matrix')
    }
    Z_list <- rbind(Z_list,dataset_Z_list)
    dataset_Z_list <- as(Z_list,'matrix') #make sure it keeps the dim with multiple datasets
  }
  Z <- reticulate::np_array(Z_list)
  ZW_Wanchor <- OCAT$sparse_encoding$Z_to_ZW(Z)
  ZW <- ZW_Wanchor[[1]]
  ZW[is.na(ZW)] <- 0
  W_anchor <- ZW_Wanchor[[2]]
  ZW <- OCAT$sparse_encoding$norm(ZW)
  if (if_inference)
    return(list(ZW,anchor_list,s_list,W_anchor))
  else
    return(ZW)
}


#' @description  The function performs OCAT pipeline and return the sparsed encoding oof edges
#' @param object  Seurat Object
#' @param m_list  num of anchors
#' @param s_list  num of anchors to be selected
#' @param dim     default values of dim
#' @param p       percentage of NNs to consider,default = 0.3
#' @param log_norm  num of anchors to be selected
#' @param if_inference  flag for cell inference
#' @param random_seed default = 42
#' Out:
#' @return ZW (a+...+z, m)  OCAT feature matrix
#' @export
run_OCAT <- function(object, group.by.vars, ...) {
  UseMethod("run_OCAT")
}


run_OCAT.default <- function(data_list, m_list=NULL, s_list=NULL, dim=NULL,
                     p=0.3, log_norm=TRUE, l2_norm=TRUE, if_inference=FALSE,
                     random_seed=42){
  if (is.null(m_list))
    m_list <- m_estimate(data_list)
  if (is.null(s_list))
    s_list <- lapply(m_list,"*",p)
    s_list <- lapply(s_list,round)
  if (is.null(dim))
      dim <- dim_estimate(data_list)
  
  data_list <- preprocess(data_list, log_norm=TRUE, l2_norm=TRUE)
  Wm <- NULL
  
  if (if_inference){
    data_list_Wm <- apply_dim_reduct(data_list, dim=dim, mode='FSM', random_seed=random_seed)
    data_list <- data_list_Wm[[1]]
    Wm <- data_list_Wm[[2]]
    sparse_list <- sparse_encoding_integration(data_list, m_list=m_list, s_list=s_list, p=p, cn=5, if_inference=TRUE)
    db_list <- list(sparse_list[[2]],s_list,sparse_list[[4]],Wm)
    return(list(ZW,db_list))
  }
  else {
    data_list <- apply_dim_reduct(data_list, dim=dim, mode='FSM', random_seed=random_seed)[[1]]
    ZW <- sparse_encoding_integration(data_list, m_list=m_list, s_list=s_list, p=p, cn=5, if_inference=FALSE)
    return(ZW)
  }
}



run_OCAT.Seurat <- function(object, reduction.name ='ocat',m_list=NULL, s_list=NULL, dim=NULL,
                             p=0.3, log_norm=TRUE, l2_norm=TRUE, if_inference=FALSE,
                             random_seed=42, assay = "RNA"){
  
  data_list <- c(t(object@assays$RNA@data))
  ZW <- run_OCAT.Default(data_list, m_list, s_list, dim,
                         p, log_norm, l2_norm, if_inference,
                         random_seed) 
  rownames(ZW) <- rownames(data_list[[1]])
  colnames(ZW) <- paste0(reduction.name, "_", seq_len(ncol(ZW)))
  suppressWarnings({
    ocat_data <- Seurat::CreateDimReducObject(embeddings = as(ZW,"matrix"), 
                                                stdev = as.numeric(apply(ZW, 2, stats::sd)), 
                                                assay = assay, key = reduction.name)
  })
  object[[reduction.name]] <- ocat_data
  return(object)
}


#' @description  The function wraps the evaluate_cluster function from OCAT python packagee
#' @param ZW   OCAT feature matrix
#' @param num_clusters  num of clusters
#' Out:
#' @return labels_pred  predicted labels 
#' @export
evaluate_clusters <- function(ZW, num_clusters){
  OCAT <- reticulate::import('OCAT')
  labels_pred <- OCAT$evaluate_clusters(reticulate::np_array(ZW),num_clusters)
  return(labels_pred)
}

#' @description The function wraps the normalized_mutual_info_score from sklearn.metrics.cluster
#' @param labels_true  ground truth labels
#' @param labels_pred  predicted labels
#' Out:
#' @return NMI     NMI score
#' @export
normalized_mutual_info_score <- function(labels_true, labels_pred){
  sklearn <- reticulate::import('sklearn')
  NMI <- sklearn$metrics$cluster$normalized_mutual_info_score(labels_true, labels_pred)
  return(NMI)
}