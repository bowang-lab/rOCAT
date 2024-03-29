% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oxy-utils.R, R/utils.R
\name{sparse_encoding_integration}
\alias{sparse_encoding_integration}
\title{FUNCTION_TITLE}
\usage{
sparse_encoding_integration(
  data_list,
  m_list,
  s_list,
  p = 0.3,
  cn = 5,
  if_inference = TRUE
)

sparse_encoding_integration(
  data_list,
  m_list,
  s_list,
  p = 0.3,
  cn = 5,
  if_inference = TRUE
)
}
\arguments{
\item{data_list}{\link{(a,dim)...(z,dim)}  list of datasets (dim PCs)}

\item{m_list}{num of anchors}

\item{s_list}{num of anchors to be selected}

\item{p}{percentage of NNs to consider}

\item{cn}{rounds of optimization}

\item{if_inference}{flag for cell inference
Out:}
}
\value{
OUTPUT_DESCRIPTION

ZW        (a+...+z, m)  OCAT feature matrix
}
\description{
FUNCTION_DESCRIPTION

The function generates the sparsed encoding of edges
}
\details{
DETAILS
}
\examples{
\dontrun{
if(interactive()){
 #EXAMPLE1
 }
}
}
\seealso{
\code{\link[reticulate]{import}}, \code{\link[reticulate]{np_array}}, \code{\link[reticulate]{py_run}}
\code{\link[glue]{glue}}
}
