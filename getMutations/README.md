The getMutations package allows to build a structured representation of single nucleotide variants contained in a VCF file, given a reference genome and a context length. 
It has 3 functions: getMutations, countTableAll and countTableSNV.

# getMutations

It provides a structured representation of the single nucleotide variants in a VCF file, provided a reference genome and a context length, which defines how many bases are shown before and after the SNP. 
The context length is the sum of bases before and after the mutation and the mutation base. If it is even, the number of bases after the mutation will be one more than the number of bases before the mutation.

# countTableAll

It provides a count table of the mutations, considering the context length (how many times a certain mutation appears in your file).

# countTableSNV

It provides a count table of only the SNP mutations, without considering the context length (how many times a certain SNV appears in your file).
