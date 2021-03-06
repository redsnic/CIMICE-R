% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing.R
\name{prepare_labels}
\alias{prepare_labels}
\title{Prepare node labels based on genotypes}
\usage{
prepare_labels(samples, genes)
}
\arguments{
\item{samples}{input dataset (mutational matrix) as matrix}

\item{genes}{list of gene names (in the columns' order)}
}
\value{
the computed edge list
}
\description{
Prepare node labels so that each node is labelled with a
comma separated list of the alterated genes
representing its associated genotype.
}
\details{
Note that after this procedure the user is
expected also to run fix_clonal_genotype
to also add the clonal genortype to the
mutational matrix if it is not present.
}
\examples{
require(dplyr) 

# compact
compactedDataset <- compact_dataset(example_dataset())
samples <- compactedDataset$matrix

# save genes' names
genes <- colnames(compactedDataset$matrix)

# keep the information on frequencies for further analysis
freqs <- compactedDataset$counts/sum(compactedDataset$counts)

# prepare node labels listing the mutated genes for each node
labels <- prepare_labels(samples, genes)

}
