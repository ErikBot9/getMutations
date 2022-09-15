#' Structured mutations from a VCF file
#' 
#' Get a structured vector of single nucleotide variants given a VCF file, a reference genome and a context length
#' @param vcf_file is the loaded VCF file containing the mutations
#' @param ref_genome is the reference genome
#' @param context_length is the length of the region displayed. It considers the mutation base and bases upstream and downstream of the mutation.
#' @return Structured vector of single nucleotide variants
#' @examples
#' library(VariantAnnotation)
#' vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
#' vcf <- readVcf(vcffile, "hg19")[1:20]
#' vcf <- renameSeqlevels(vcf, 'chr22')
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' getMutations(vcf, Hsapiens, 5)
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges seqnames
#' @importFrom methods is
#' @importFrom Biostrings reverseComplement
#' @import IRanges
#' @importFrom MatrixGenerics rowRanges
#' @import VariantAnnotation
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @export

getMutations <- function(vcf_file, ref_genome, context_length) {
  if(is(vcf_file, 'BSgenome') & is(ref_genome, 'CollapsedVCF')) {
    stop('You might have inserted the VCF file under \'ref_genome\' and the reference genome under \'vcf_file\'.')
  }
  if(!is(vcf_file, 'CollapsedVCF')) {
    stop('Check if you have correctly loaded the VCF file! The \'vcf_file\' object class should be \'CollapsedVCF\'. Try to use the function \'readVcf\' of the \'VariantAnnotation\' package.')
  }
  if(!is(ref_genome, 'BSgenome')) {
    stop('The reference genome should be of \'BSgenome\' class.')
  }
  if(!is(context_length, 'numeric')){
    stop('The \'context_length\' parameter should be numeric!')
  }
  if(!all(as.vector(seqnames(vcf_file)) %in%  as.vector(ref_genome@user_seqnames))) {
    stop('Check that the sequences names of the VCF file and the reference genome are the same by using the command: all(seqnames(rowRanges(vcf_file)) %in%  ref_genome@user_seqnames). If they are different, you might change them using the function: renameSeqlevels(vcf_file, names)')
  }
  if(context_length %% 2 == 0) {
    print('The context length is even: the extra base will be added at the end of the string, not at the beginning.')
  }
  v <- c()
  gr <- rowRanges(vcf_file)
  ref <- gr$REF
  alt <- unlist(gr$ALT)
  new_gr <- GRanges(seqnames = seqnames(gr), IRanges((start(gr) - (context_length - 1)/2), (end(gr) + (context_length - 1)/2)))
  s <- getSeq(ref_genome, new_gr)
  rev_s <- as.character(reverseComplement(s))
  s <- as.character(s)
  for(i in seq(1,length(gr))) {
    if(as.character(ref[i]) %in% c('A', 'G')) {
      v <- c(v, paste0(substr(rev_s[i], 1, (nchar(rev_s[i])- 1)/2), "[", reverseComplement(ref[i]), ">", reverseComplement(alt[i]), "]", substr(rev_s[i], (nchar(rev_s[i])-((nchar(rev_s[i])- 1)/2)+1), nchar(rev_s[i])), sep = ''))
    }
    else {
      v <- c(v, paste0(substr(s[i], 1, (nchar(s[i])- 1)/2), "[", ref[i], ">", alt[i], "]", substr(s[i], (nchar(s[i])-((nchar(s[i])- 1)/2)+1), nchar(s[i])), sep = ''))
    }
  }
  # to keep only the SNV in the VCF file
  v <- v[nchar(v) == context_length + 4]
  return(v)
}

#' The count table of each SNV, with context length
#' 
#' Get a count table of the absolute frequency of each SNV, considering the context length.
#' @param mut_file is the mutation vector obtained with getMutations.
#' @param vcf_file is the loaded VCF file containing the mutations. Needed if no mut_file is provided.
#' @param ref_genome is the reference genome. Needed if no mut_file is provided.
#' @param context_length is the length of the region displayed. It considers the mutation base and bases upstream and downstream of the mutation. Needed if no mut_file is provided.
#' @return Count table of the absolute frequency of each mutation, considering the context length.
#' @examples
#' library(VariantAnnotation)
#' vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
#' vcf <- readVcf(vcffile, "hg19")[1:20]
#' vcf <- renameSeqlevels(vcf, 'chr22')
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' countTableAll(vcf_file = vcf, ref_genome = Hsapiens, context_length = 5)
#' mut_file <- getMutations(vcf, Hsapiens, 3)
#' countTableAll(mut_file)
#' @importFrom GenomicRanges GRanges
#' @importFrom methods is
#' @importFrom Biostrings reverseComplement
#' @import IRanges
#' @importFrom MatrixGenerics rowRanges
#' @import VariantAnnotation
#' @importFrom stringr str_sub
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @export

countTableAll <- function(mut_file = NULL, vcf_file = NULL, ref_genome = NULL, context_length = NULL) {
  if(!is.null(mut_file)) {
    if(is(mut_file, 'CollapsedVCF')) {
      stop('You have submitted a VCF file instead of a mutation file. If you wanted to use a VCF file, you need to specify \'vcf_file\', \'ref_genome\' and \'context_length\' before the submitted files/values.')
    }
    if(!(length(grep('[*>*]', mut_file) == length(mut_file)))) {
      stop('The mutation file should be a vector obtained via getMutations (or the same structure)!')
    }
    if(!(is.null(vcf_file) & is.null(ref_genome) & is.null(context_length))) {
      print('If you provide the mutation file, you don\'t need to provide the VCF file or the reference genome or the context length. They will be ignored.')
    }
    countTable <- table(mut_file)
  }
  else {
    if(is.null(ref_genome) | is.null(vcf_file) | is.null(context_length)){
      stop('You have either to provide a mutation file or the reference genome, the vcf file AND the context length!')}
    mut_file <- getMutations(vcf_file, ref_genome, context_length)
    countTable <- table(mut_file)
  }
  return(countTable)
}

#' The count table of each SNV, without context length
#' 
#' Get a count table of the absolute frequency of each SNV, not considering the context length.
#' @param mut_file is the mutation vector obtained with getMutations.
#' @param vcf_file is the loaded VCF file containing the mutations. Needed if no mut_file is provided.
#' @param ref_genome is the reference genome. Needed if no mut_file is provided.
#' @return Count table of the absolute frequency of each SNV, not considering the context length.
#' @examples
#' library(VariantAnnotation)
#' vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
#' vcf <- readVcf(vcffile, "hg19")[1:20]
#' vcf <- renameSeqlevels(vcf, 'chr22')
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' countTableSNV(vcf_file = vcf, ref_genome = Hsapiens)
#' mut_file <- getMutations(vcf, Hsapiens, 3)
#' countTableSNV(mut_file)
#' @importFrom GenomicRanges GRanges
#' @importFrom methods is
#' @importFrom Biostrings reverseComplement
#' @import IRanges
#' @importFrom MatrixGenerics rowRanges
#' @import VariantAnnotation
#' @importFrom stringr str_sub
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @export
countTableSNV <- function(mut_file = NULL, vcf_file = NULL, ref_genome = NULL) {
  if(!is.null(mut_file)) {
    if(is(mut_file, 'CollapsedVCF')) {
      stop('You have submitted a VCF file instead of a mutation file. If you wanted to use a VCF file, you need to specify \'vcf_file\', \'ref_genome\' and \'context_length\' before the submitted files/values.')
    }
    if(!(length(grep('[*>*]', mut_file) == length(mut_file)))) {
      stop('The mutation file should be a vector obtained via getMutations (or the same structure)!')
    }
    if(!(is.null(vcf_file) & is.null(ref_genome))) {
      print('If you provide the mutation file, you don\'t need to provide the VCF file or the reference genome. They will be ignored.')
    }
    if(nchar(mut_file[1]) -4 %% 2 == 0) {
      mut_file <- str_sub(mut_file, (((nchar(mut_file[1]) - 4 -2)/2) + 1), ((-(nchar(mut_file[1]) - 4)/2) - 1))}
    else {
      mut_file <- str_sub(mut_file, (((nchar(mut_file[1]) -4 -1)/2) + 1), ((-(nchar(mut_file[1])-4 -1)/2) - 1))}
    countTable <- table(mut_file)
  }
  else{
    if(is.null(vcf_file) | is.null(ref_genome)){
      stop('You have either to provide a mutation file or the reference genome AND the VCF file!')}
    mut_file <- getMutations(vcf_file, ref_genome, 1)
    countTable <- table(mut_file)
  }
  return(countTable)
}
