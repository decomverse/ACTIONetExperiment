% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{import.ace.from.counts}
\alias{import.ace.from.counts}
\title{Constructs an `SingleCellExeriment` object from count matrix, gene names, and sample_metadata}
\usage{
import.ace.from.counts(
  counts_mat,
  feature_metadata = NULL,
  sample_metadata = NULL,
  prefilter = FALSE,
  ...
)
}
\arguments{
\item{counts_mat}{Matrix of counts.mat}

\item{feature_metadata}{Vector or data frame  of annotations for features (rows of the count matrix)}

\item{sample_metadata}{Vector or data frame of annotations for samples (columns of the count matrix)}

\item{prefilter}{Whether to prefilter genes/cells based on the counts.mat}

\item{min.cell.frac.per.gene}{Minimum fraction of cells capturing a gene for it to be retained, if prefilter=TRUE (default=0.005)}

\item{min.genes.per.cell}{Minimum number of required captured genes per cell, if prefilter=TRUE (default=500)}
}
\value{
`SingleCellExeriment` object
}
\description{
Constructs an `SingleCellExeriment` object from count matrix, gene names, and sample_metadata
}
\examples{
ace = import.ace.from.count.matrix(counts.mat.mat, gene_names, prefilter=TRUE)
}
