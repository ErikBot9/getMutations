test_that('Running getMutations with a wrong order of files', {
          library(VariantAnnotation)
          
          vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
          
          vcf <- readVcf(vcffile, "hg19")[1:20]
          
          vcf <- renameSeqlevels(vcf, 'chr22')
          
          library(BSgenome.Hsapiens.UCSC.hg19)
          
          expect_error(getMutations(Hsapiens, vcf, 3))
}
)

test_that('Running getMutations with different sequences names between the VCF file and the reference genome', {
  library(VariantAnnotation)
  
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  
  vcf <- readVcf(vcffile, "hg19")[1:20]
  
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  expect_error(getMutations(vcf, Hsapiens, 3))
}
)

test_that('Running getMutations with the wrong VCF file', {
  library(VariantAnnotation)
  
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  expect_error(getMutations(vcffile, Hsapiens, 3))
}
)

test_that('Running getMutations with the wrong reference genome file', {
  library(VariantAnnotation)
  
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  
  vcf <- readVcf(vcffile, "hg19")[1:20]
  
  vcf <- renameSeqlevels(vcf, 'chr22')
  
  expect_error(getMutations(vcf, 'drosophila', 3))
}
)

test_that('Running getMutations with the wrong context length parameter', {
  library(VariantAnnotation)
  
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  
  vcf <- readVcf(vcffile, "hg19")[1:20]
  
  vcf <- renameSeqlevels(vcf, 'chr22')
  
  expect_error(getMutations(vcf, Hsapiens, '3'))
}
)

test_that('Running getMutations with an even context length', {
  library(VariantAnnotation)
  
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  
  vcf <- readVcf(vcffile, "hg19")[1:20]
  
  vcf <- renameSeqlevels(vcf, 'chr22')
  
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  expect_output(getMutations(vcf, Hsapiens, 4))
}
)

test_that('Running getMutations and checking the size of the output', {
  library(VariantAnnotation)
  
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  
  vcf <- readVcf(vcffile, "hg19")[1:20]
  
  vcf <- renameSeqlevels(vcf, 'chr22')
  
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  expect_length(getMutations(vcf, Hsapiens, 3), length(vcf))
}
)

test_that('Running getMutations and checking the structure of the mutations', {
  library(VariantAnnotation)
  
  vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  
  vcf <- readVcf(vcffile, "hg19")[1:20]
  
  vcf <- renameSeqlevels(vcf, 'chr22')
  
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  mut_file <- getMutations(vcf, Hsapiens, 3)
  
  for(i in seq(1, length(vcf))) {
    expect_match(mut_file[i], '[.>.]')
  }
}
)

