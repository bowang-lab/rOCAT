#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION

#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[reticulate]{py_install}}
#' @rdname load_all
#' @export 
#' @importFrom reticulate py_install
load_all <- function(){
  reticulate::py_install('git+https://github.com/bowang-lab/OCAT.git',pip=TRUE)
  reticulate::py_install('faiss')
  reticulate::py_install("scikit-learn")
  }


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param data_list PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[Matrix]{Matrix}}
#' @rdname normalize_data
#' @export 
#' @importFrom Matrix Matrix
normalize_data <- function(data_list){
  data_list <- sapply(data_list, function(x) log10(x+1))
  data_list <- sapply(data_list, function(x) Matrix::Matrix(x,sparse = TRUE))
  return(data_list)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param data_list PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[Matrix]{colSums}}, \code{\link[Matrix]{Matrix}}
#' @rdname l2_normalization
#' @export 
#' @importFrom Matrix colSums Matrix
l2_normalization <- function(data_list){
  for (i in seq_along(data_list)){
    l2_norm <- sqrt(Matrix::colSums(data_list[[i]]^2))
    l2_norm[l2_norm==0] <- 0.00001
    data_list[[i]] <- sweep(data_list[[i]],2,l2_norm,'/')
    data_list[[i]] <- Matrix::Matrix(data_list[[i]], sparse = TRUE)
  }
  return(data_list)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param data_list PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname m_estimate
#' @export 
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

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param data_list PARAM_DESCRIPTION
#' @param log_norm PARAM_DESCRIPTION, Default: log_norm
#' @param l2_norm PARAM_DESCRIPTION, Default: l2_norm
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname preprocess
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


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param data_list PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname dim_estimate
#' @export 
dim_estimate <- function(data_list){
  gene_num <- data_list[[1]]@Dim[1]
  if (gene_num<5000)
    return(50)
  else if (gene_num<10000)
    return(100)
  else
    return(125)
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param data_list PARAM_DESCRIPTION
#' @param dim PARAM_DESCRIPTION, Default: NULL
#' @param mode PARAM_DESCRIPTION, Default: 'FSM'
#' @param random_seed PARAM_DESCRIPTION, Default: 42
#' @param upsample PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[reticulate]{import}}, \code{\link[reticulate]{np_array}}, \code{\link[reticulate]{py_run}}
#'  \code{\link[glue]{glue}}
#' @rdname apply_dim_reduct
#' @export 
#' @importFrom reticulate import np_array py_run_string
#' @importFrom glue glue
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

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param data_list PARAM_DESCRIPTION
#' @param m_list PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[reticulate]{import}}, \code{\link[reticulate]{np_array}}, \code{\link[reticulate]{py_run}}
#'  \code{\link[glue]{glue}}
#' @rdname find_anchors
#' @export 
#' @importFrom reticulate import np_array py_run_string
#' @importFrom glue glue
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

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param data_list PARAM_DESCRIPTION
#' @param m_list PARAM_DESCRIPTION
#' @param s_list PARAM_DESCRIPTION
#' @param p PARAM_DESCRIPTION, Default: 0.3
#' @param cn PARAM_DESCRIPTION, Default: 5
#' @param if_inference PARAM_DESCRIPTION, Default: TRUE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[reticulate]{import}}, \code{\link[reticulate]{np_array}}, \code{\link[reticulate]{py_run}}
#'  \code{\link[glue]{glue}}
#' @rdname sparse_encoding_integration
#' @export 
#' @importFrom reticulate import np_array py_run_string
#' @importFrom glue glue
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


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param object PARAM_DESCRIPTION
#' @param group.by.vars PARAM_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname run_OCAT
#' @export 
run_OCAT <- function(object, group.by.vars, ...) {
  UseMethod("run_OCAT")
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param data_list PARAM_DESCRIPTION
#' @param m_list PARAM_DESCRIPTION, Default: NULL
#' @param s_list PARAM_DESCRIPTION, Default: NULL
#' @param dim PARAM_DESCRIPTION, Default: NULL
#' @param p PARAM_DESCRIPTION, Default: 0.3
#' @param log_norm PARAM_DESCRIPTION, Default: TRUE
#' @param l2_norm PARAM_DESCRIPTION, Default: TRUE
#' @param if_inference PARAM_DESCRIPTION, Default: FALSE
#' @param random_seed PARAM_DESCRIPTION, Default: 42
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname run_OCAT.default
#' @export 
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



#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param object PARAM_DESCRIPTION
#' @param reduction.name PARAM_DESCRIPTION, Default: 'ocat'
#' @param m_list PARAM_DESCRIPTION, Default: NULL
#' @param s_list PARAM_DESCRIPTION, Default: NULL
#' @param dim PARAM_DESCRIPTION, Default: NULL
#' @param p PARAM_DESCRIPTION, Default: 0.3
#' @param log_norm PARAM_DESCRIPTION, Default: TRUE
#' @param l2_norm PARAM_DESCRIPTION, Default: TRUE
#' @param if_inference PARAM_DESCRIPTION, Default: FALSE
#' @param random_seed PARAM_DESCRIPTION, Default: 42
#' @param assay PARAM_DESCRIPTION, Default: 'RNA'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[Seurat]{reexports}}
#'  \code{\link[stats]{sd}}
#' @rdname run_OCAT.Seurat
#' @export 
#' @importFrom Seurat CreateDimReducObject
#' @importFrom stats sd
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


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param ZW PARAM_DESCRIPTION
#' @param num_clusters PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[reticulate]{import}}, \code{\link[reticulate]{np_array}}
#' @rdname evaluate_clusters
#' @export 
#' @importFrom reticulate import np_array
evaluate_clusters <- function(ZW, num_clusters){
  OCAT <- reticulate::import('OCAT')
  labels_pred <- OCAT$evaluate_clusters(reticulate::np_array(ZW),num_clusters)
  return(labels_pred)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param labels_true PARAM_DESCRIPTION
#' @param labels_pred PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[reticulate]{import}}
#' @rdname normalized_mutual_info_score
#' @export 
#' @importFrom reticulate import
normalized_mutual_info_score <- function(labels_true, labels_pred){
  sklearn <- reticulate::import('sklearn')
  NMI <- sklearn$metrics$cluster$normalized_mutual_info_score(labels_true, labels_pred)
  return(NMI)
}
