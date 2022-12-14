% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getMutations.R
\name{countTableAll}
\alias{countTableAll}
\title{The count table of each SNV, with context length}
\usage{
countTableAll(
  mut_file = NULL,
  vcf_file = NULL,
  ref_genome = NULL,
  context_length = NULL
)
}
\arguments{
\item{mut_file}{is the mutation vector obtained with getMutations.}

\item{vcf_file}{is the loaded VCF file containing the mutations. Needed if no mut_file is provided.}

\item{ref_genome}{is the reference genome. Needed if no mut_file is provided.}

\item{context_length}{is the length of the region displayed. It considers the mutation base and bases upstream and downstream of the mutation. Needed if no mut_file is provided.}
}
\value{
Count table of the absolute frequency of each mutation, considering the context length.
}
\description{
Get a count table of the absolute frequency of each SNV, considering the context length.
}
\examples{
library(VariantAnnotation)
vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
vcf <- readVcf(vcffile, "hg19")[1:20]
vcf <- renameSeqlevels(vcf, 'chr22')
library(BSgenome.Hsapiens.UCSC.hg19)
countTableAll(vcf_file = vcf, ref_genome = Hsapiens, context_length = 5)
mut_file <- getMutations(vcf, Hsapiens, 3)
countTableAll(mut_file)
}
