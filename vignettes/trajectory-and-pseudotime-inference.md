In this page, we demonstrated how to run trajectory and pseudotime
inference in [HSMM
dataset](https://github.com/bowang-lab/rOCAT/blob/main/vignettes/HSMM_raw_counts.txt)
with our rOCAT package and HSMM
[labels](https://data.wanglab.ml/OCAT/HSMM.zip)

# load data and apply OCAT pipeline

## load package

    library(rOCAT)
    library(reticulate)
    library(Matrix)

## load data

    # load HSMM data
    data <- read.table("HSMM_raw_counts.txt",row.names = 1)
    data <- Matrix::Matrix(as.matrix(data))
    data <- t(data) # format as cell x Gene

    data_list <- c(data)

    ZW <- run_OCAT(data_list, m_list=list(40), dim=80)
    #> 
    #> Done!
    #> [1] "Starting Dimension Reduction"
    #> [1] "Dimension Reduction Finished"

# Trajectory inference

Load in the annotated labels.

    #load labels
    library(plyr)
    labels <- data.frame(read.table('../../data/HSMM/HSMM_label.txt',row.names = NULL))
    map <- data.frame(cbind(V1=c(1,2,3,4,5), 
              to=c("Fibroblast","Myotubes","Myoblasts","Undiff","Intermediates")))
    labels_mapped <- plyr::join(labels,map,by='V1')$to

OCAT.compute\_lineage() function infers Lineages over clusters with the
OCAT features, predicted/true cluster labels and a user-specified
root\_cluster.

    L <- compute_lineage(ZW,labels_mapped, root_cluster='Myoblasts', name='OE', reverse=0)
    Lineage <- L[[1]]
    root_cluster <- L[[2]]
    cluster_labels <- L[[3]]
    tree <- L[[4]]

# Pseudotime inference

OCAT.compute\_ptime() function infers pseudotime for individual cell
using the OCAT extracted features and the predicted lineage.

    P  <- compute_ptime(ZW, labels_mapped, Lineage, root_cluster)
    Ptime <- P[[1]]
    root_cell_list <- P[[2]]
