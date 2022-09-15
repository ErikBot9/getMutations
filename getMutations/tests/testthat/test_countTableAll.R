test_that('Running countTableAll with a wrong order of files', {
  library(VariantAnnotation)
  
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  
  vcf <- readVcf(vcffile, "hg19")[1:20]
  
  vcf <- renameSeqlevels(vcf, 'chr22')
  
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  expect_error(countTableAll(vcf, Hsapiens, 3))
}
)

test_that('Running countTableAll with useless files addded', {
  library(VariantAnnotation)
  
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  
  vcf <- readVcf(vcffile, "hg19")[1:20]
  
  vcf <- renameSeqlevels(vcf, 'chr22')
  
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  mut_file <- getMutations(vcf, Hsapiens, 3)
  
  expect_output(countTableAll(mut_file, vcf, Hsapiens, 3))
}
)

test_that('Running countTableAll with a wrong mutation file', {
  
  mut_file <- 'Not a mutation file'
    
  expect_error(countTableAll(mut_file))
}
)

test_that('Running countTableAll without all the parameters', {
  library(VariantAnnotation)
  
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  
  vcf <- readVcf(vcffile, "hg19")[1:20]
  
  vcf <- renameSeqlevels(vcf, 'chr22')
  
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  expect_error(countTableAll(vcf_file = vcf, ref_genome = Hsapiens))
}
)

test_that('Running countTableAll returns a table', {
  library(VariantAnnotation)
  
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  
  vcf <- readVcf(vcffile, "hg19")[1:20]
  
  vcf <- renameSeqlevels(vcf, 'chr22')
  
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  expect_equal(class(countTableAll(vcf_file = vcf, ref_genome = Hsapiens, context_length = 3)), 'table')
}
)

test_that('Running countTableAll and checking the names of the table', {
  library(VariantAnnotation)
  
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  
  vcf <- readVcf(vcffile, "hg19")[1:20]
  
  vcf <- renameSeqlevels(vcf, 'chr22')
  
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  table <- countTableAll(vcf_file = vcf, ref_genome = Hsapiens, context_length = 3)
  
  for(i in seq(1, length(table))) {
    expect_match(names(table)[i], '[.>.]')
  }
}
)
 
test_that('Running countTableAll gives the correct number of mutations', {
  library(VariantAnnotation)
  
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  
  vcf <- readVcf(vcffile, "hg19")[1:20]
  
  vcf <- renameSeqlevels(vcf, 'chr22')
  
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  expect_equal(sum(countTableAll(vcf_file = vcf, ref_genome = Hsapiens, context_length = 3)), length(vcf))
}
)