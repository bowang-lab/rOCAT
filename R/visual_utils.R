#' @title run_UMAP
#' @description  The function calculate the UMAP embedding
#' @param embeddding  OCAT embedding
#' Out:
#' @return OCAT_umap calculated UMAP embedding
#' @export
run_UMAP <- function(embeddding){
  if (!reticulate::py_module_available("umap"))
    reticulate::py_install("umap-learn")
  umap <- reticulate::import('umap')
  OCAT_umap <- umap$UMAP()$fit_transform(embeddding)
  return(OCAT_umap)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' @title plot_UMAP
#' @description  The function plot the UMAP embedding
#' @param embeddding  UMAP embedding
#' @param labels_true  UMAP embedding
#' @param labels_predicted  UMAP embedding
#' @param legend_labels  UMAP embedding
#' @param title  UMAP embedding
#' @param x  the x and y co-ordinates to be used to position the legend. They can be specified by keyword or in any way which is accepted by xy.coords
#' @param y.intersp  character interspacing factor for vertical (y) line distance
#' @export
plot_UMAP <- function(umap_embedding,labels_true,labels_predicted, legend_labels = NULL, title = '', x='bottomright', y.intersp=1){
  unique_cell_types <- sort(unique(labels_predicted))
  cell_type_color <- gg_color_hue(length(unique_cell_types))

  if (is.null(legend_labels))
    legend_labels <- unique_cell_types

  plot(umap_embedding[,1], umap_embedding[,2], type='n', xlab='UMAP_1', ylab='UMAP_2',
       yaxt  ='n',xaxt='n',
       bty='L', tck=FALSE,fg='grey', mgp=c(1,1,0), main=title)

  for(i in 1:length(unique_cell_types)){
    temp_ind <- which(labels_predicted==unique_cell_types[i])
    points(umap_embedding[temp_ind,1], umap_embedding[temp_ind,2],
           col=cell_type_color[i], pch=20, cex=0.5)
  }
  legend(x,
         legend=legend_labels,
         col=cell_type_color, bty='n', pch=19, y.intersp=y.intersp)
}
