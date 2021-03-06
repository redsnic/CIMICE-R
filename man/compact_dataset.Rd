% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing.R
\name{compact_dataset}
\alias{compact_dataset}
\title{Compact dataset rows}
\usage{
compact_dataset(mutmatrix)
}
\arguments{
\item{mutmatrix}{input dataset (mutational matrix)}
}
\value{
a list with matrix (the compacted dataset (mutational matrix)), counts 
(frequencies of genotypes) and row_names (comma separated string of sample IDs) fields
}
\description{
Count duplicate rows
and compact the dataset (mutational). The column
'freq' will contain the counts for each row.
}
\examples{
compact_dataset(example_dataset())

}
