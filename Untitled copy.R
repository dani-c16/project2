library(GenomicRanges)
library(Biostrings)
library(bsseq)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install(c("GenomicFeatures", "Biostrings", "rtracklayer"))
library(Biostrings)
if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE))
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
hg38 <- BSgenome.Hsapiens.UCSC.hg38
BiocManager::install("biomartr")
library(biomartr)
hg38_genome <- genome(db = "refseq", organism = "Homo sapiens", assembly = "GRCh38")

library(rtracklayer)
library(BSgenome)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)

getSeq(BSgenome.Hsapiens.UCSC.hg38, promoters)
promoters <- promoters(hg38_genome, upstream = 2000)

txdb <- TxDb.Hsapiens.UCSC.hg38.knowngene
my_promoters <- GenomicRanges::promoters(x, hg38_genome) #want to pull out the promoter regions mapped to reference genome hg38

#want to define what full coverage is to use in the down-sampling function. Here i wnat to define full coverage as the the result of the CG read numbers for promoters in seq. May start with just the first 15 promoters.
full_coverage <- function()
  ds<- function(full_coverage, fraction = 0.9){
    s <- sample(full_coverage, 0.9)
    return(s)
  }




geno <- seqlengths(BSgenome.Hsapiens.UCSC.hg38) 
wholechromo1 <- GRanges(seqnames = names(geno)
                        , ranges = IRanges(start = 1, end = geno))
wholechromo2 <- GRanges(seqnames = names(geno)
                        , ranges = IRanges(start = 1, end = geno))
wholechromo[1:12]








