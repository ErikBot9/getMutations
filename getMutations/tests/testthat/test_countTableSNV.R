test_that('Running countTableSNV with a wrong order of files', {
  library(VariantAnnotation)
  
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  
  vcf <- readVcf(vcffile, "hg19")[1:20]
  
  vcf <- renameSeqlevels(vcf, 'chr22')
  
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  expect_error(countTableSNV(vcf, Hsapiens))
}
)

test_that('Running countTableSNV with useless files addded', {
  library(VariantAnnotation)
  
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  
  vcf <- readVcf(vcffile, "hg19")[1:20]
  
  vcf <- renameSeqlevels(vcf, 'chr22')
  
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  mut_file <- getMutations(vcf, Hsapiens, 3)
  
  expect_output(countTableSNV(mut_file, vcf, Hsapiens))
}
)

test_that('Running countTableSNV with a wrong mutation file', {
  
  mut_file <- 'Not a mutation file'
  
  expect_error(countTableSNV(mut_file))
}
)

test_that('Running countTableSNV without all the parameters', {
  library(VariantAnnotation)
  
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  
  vcf <- readVcf(vcffile, "hg19")[1:20]
  
  vcf <- renameSeqlevels(vcf, 'chr22')
  
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  expect_error(countTableSNV(vcf_file = vcf))
}
)

test_that('Running countTableSNV returns a table', {
  library(VariantAnnotation)
  
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  
  vcf <- readVcf(vcffile, "hg19")[1:20]
  
  vcf <- renameSeqlevels(vcf, 'chr22')
  
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  expect_equal(class(countTableSNV(vcf_file = vcf, ref_genome = Hsapiens)), 'table')
}
)

test_that('Running countTableSNV and checking the names of the table', {
  library(VariantAnnotation)
  
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  
  vcf <- readVcf(vcffile, "hg19")[1:20]
  
  vcf <- renameSeqlevels(vcf, 'chr22')
  
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  table <- countTableSNV(vcf_file = vcf, ref_genome = Hsapiens)
  
  for(i in seq(1, length(table))) {
    expect_match(names(table)[i], '[.>.]')
  }
}
)

test_that('Running countTableSNV gives the correct number of mutations', {
  library(VariantAnnotation)
  
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  
  vcf <- readVcf(vcffile, "hg19")[1:20]
  
  vcf <- renameSeqlevels(vcf, 'chr22')
  
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  expect_equal(sum(countTableSNV(vcf_file = vcf, ref_genome = Hsapiens)), length(vcf))
}
)